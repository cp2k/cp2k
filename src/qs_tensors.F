!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Utility methods to build 3-center integral tensors of various types.
! **************************************************************************************************
MODULE qs_tensors
   USE ai_contraction,                  ONLY: block_add
   USE ai_contraction_sphi,             ONLY: ab_contract,&
                                              libxsmm_abc_contract
   USE ai_overlap,                      ONLY: overlap_ab
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_p_type,&
                                              gto_basis_set_type
   USE block_p_types,                   ONLY: block_p_type
   USE cell_types,                      ONLY: cell_type
   USE cp_array_utils,                  ONLY: cp_2d_r_p_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: dbcsr_filter,&
                                              dbcsr_finalize,&
                                              dbcsr_get_block_p,&
                                              dbcsr_has_symmetry,&
                                              dbcsr_type,&
                                              dbcsr_type_real_8
   USE dbcsr_tensor_api,                ONLY: &
        dbcsr_t_blk_sizes, dbcsr_t_clear, dbcsr_t_copy, dbcsr_t_create, dbcsr_t_default_distvec, &
        dbcsr_t_destroy, dbcsr_t_distribution_destroy, dbcsr_t_distribution_new, &
        dbcsr_t_distribution_type, dbcsr_t_filter, dbcsr_t_get_block, dbcsr_t_get_info, &
        dbcsr_t_get_mapping_info, dbcsr_t_get_num_blocks, dbcsr_t_get_nze_total, &
        dbcsr_t_get_stored_coordinates, dbcsr_t_iterator_blocks_left, dbcsr_t_iterator_next_block, &
        dbcsr_t_iterator_start, dbcsr_t_iterator_stop, dbcsr_t_iterator_type, &
        dbcsr_t_mp_environ_pgrid, dbcsr_t_nblks_total, dbcsr_t_ndims, dbcsr_t_ndims_matrix_column, &
        dbcsr_t_ndims_matrix_row, dbcsr_t_pgrid_type, dbcsr_t_put_block, dbcsr_t_reserve_blocks, &
        dbcsr_t_reserved_block_indices, dbcsr_t_type, dbcsr_t_write_split_info
   USE distribution_1d_types,           ONLY: distribution_1d_type
   USE distribution_2d_types,           ONLY: distribution_2d_type
   USE gamma,                           ONLY: init_md_ftable
   USE hfx_compression_methods,         ONLY: hfx_add_mult_cache_elements,&
                                              hfx_add_single_cache_element,&
                                              hfx_decompress_first_cache,&
                                              hfx_flush_last_cache,&
                                              hfx_get_mult_cache_elements,&
                                              hfx_get_single_cache_element,&
                                              hfx_reset_cache_and_container
   USE hfx_types,                       ONLY: alloc_containers,&
                                              dealloc_containers,&
                                              hfx_cache_type,&
                                              hfx_compression_type,&
                                              hfx_container_type,&
                                              hfx_init_container
   USE input_constants,                 ONLY: do_potential_coulomb,&
                                              do_potential_id,&
                                              do_potential_short,&
                                              do_potential_truncated
   USE input_section_types,             ONLY: section_vals_val_get
   USE kinds,                           ONLY: default_string_length,&
                                              dp,&
                                              int_8
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_type
   USE libint_2c_3c,                    ONLY: cutoff_screen_factor,&
                                              eri_2center,&
                                              eri_3center,&
                                              libint_potential_type
   USE libint_wrapper,                  ONLY: cp_libint_cleanup_2eri,&
                                              cp_libint_cleanup_3eri,&
                                              cp_libint_init_2eri,&
                                              cp_libint_init_3eri,&
                                              cp_libint_set_contrdepth,&
                                              cp_libint_t
   USE molecule_types,                  ONLY: molecule_type
   USE orbital_pointers,                ONLY: ncoset
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_neighbor_list_types,          ONLY: &
        get_iterator_info, get_neighbor_list_set_p, neighbor_list_iterate, &
        neighbor_list_iterator_create, neighbor_list_iterator_p_type, &
        neighbor_list_iterator_release, neighbor_list_set_p_type, nl_sub_iterate, &
        release_neighbor_list_sets
   USE qs_neighbor_lists,               ONLY: atom2d_build,&
                                              atom2d_cleanup,&
                                              build_neighbor_lists,&
                                              local_atoms_type,&
                                              pair_radius_setup
   USE qs_tensors_types,                ONLY: &
        distribution_3d_destroy, distribution_3d_type, neighbor_list_3c_iterator_type, &
        neighbor_list_3c_type, symmetric_ij, symmetric_ijk, symmetric_jk, symmetric_none, &
        symmetrik_ik
   USE t_c_g0,                          ONLY: get_lmax_init,&
                                              init
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_tensors'

   PUBLIC :: build_3c_neighbor_lists, &
             neighbor_list_3c_destroy, neighbor_list_3c_iterate, neighbor_list_3c_iterator_create, &
             neighbor_list_3c_iterator_destroy, get_3c_iterator_info, build_3c_integrals, &
             build_2c_neighbor_lists, build_2c_integrals, cutoff_screen_factor, tensor_change_pgrid, &
             get_tensor_occupancy, compress_tensor, decompress_tensor

   TYPE one_dim_int_array
      INTEGER, DIMENSION(:), ALLOCATABLE    :: array
   END TYPE

   ! cache size for integral compression
   INTEGER, PARAMETER, PRIVATE :: cache_size = 1024

CONTAINS

! **************************************************************************************************
!> \brief Build 2-center neighborlists adapted to different operators
!>        This mainly wraps build_neighbor_lists for consistency with build_3c_neighbor_lists
!> \param ij_list 2c neighbor list for atom pairs i, j
!> \param basis_i basis object for atoms i
!> \param basis_j basis object for atoms j
!> \param potential_parameter ...
!> \param name name of 2c neighbor list
!> \param qs_env ...
!> \param sym_ij Symmetry in i, j (default .TRUE.)
!> \param molecular ...
!> \param dist_2d optionally a custom 2d distribution
!> \param pot_to_rad which radius (1 for i, 2 for j) should be adapted to incorporate potential
! **************************************************************************************************
   SUBROUTINE build_2c_neighbor_lists(ij_list, basis_i, basis_j, potential_parameter, name, qs_env, &
                                      sym_ij, molecular, dist_2d, pot_to_rad)
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: ij_list
      TYPE(gto_basis_set_p_type), DIMENSION(:)           :: basis_i, basis_j
      TYPE(libint_potential_type), INTENT(IN)            :: potential_parameter
      CHARACTER(LEN=*), INTENT(IN)                       :: name
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN), OPTIONAL                      :: sym_ij, molecular
      TYPE(distribution_2d_type), OPTIONAL, POINTER      :: dist_2d
      INTEGER, INTENT(IN), OPTIONAL                      :: pot_to_rad

      CHARACTER(len=*), PARAMETER :: routineN = 'build_2c_neighbor_lists', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: ikind, nkind, pot_to_rad_prv
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: i_present, j_present
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: pair_radius
      REAL(kind=dp)                                      :: subcells
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: i_radius, j_radius
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(distribution_1d_type), POINTER                :: local_particles
      TYPE(distribution_2d_type), POINTER                :: dist_2d_prv
      TYPE(local_atoms_type), ALLOCATABLE, DIMENSION(:)  :: atom2d
      TYPE(molecule_type), DIMENSION(:), POINTER         :: molecule_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      NULLIFY (atomic_kind_set, cell, local_particles, molecule_set, &
               particle_set, dist_2d_prv)

      IF (PRESENT(pot_to_rad)) THEN
         pot_to_rad_prv = pot_to_rad
      ELSE
         pot_to_rad_prv = 1
      ENDIF

      CALL get_qs_env(qs_env, &
                      nkind=nkind, &
                      cell=cell, &
                      particle_set=particle_set, &
                      atomic_kind_set=atomic_kind_set, &
                      local_particles=local_particles, &
                      distribution_2d=dist_2d_prv, &
                      molecule_set=molecule_set)

      CALL section_vals_val_get(qs_env%input, "DFT%SUBCELLS", r_val=subcells)

      ALLOCATE (i_present(nkind), j_present(nkind))
      ALLOCATE (i_radius(nkind), j_radius(nkind))

      i_present = .FALSE.
      j_present = .FALSE.
      i_radius = 0.0_dp
      j_radius = 0.0_dp

      IF (PRESENT(dist_2d)) dist_2d_prv => dist_2d

      !  Set up the radii, depending on the operator type
      IF (potential_parameter%potential_type == do_potential_id) THEN

         !overlap => use the kind radius for both i and j
         DO ikind = 1, nkind
            IF (ASSOCIATED(basis_i(ikind)%gto_basis_set)) THEN
               i_present(ikind) = .TRUE.
               CALL get_gto_basis_set(basis_i(ikind)%gto_basis_set, kind_radius=i_radius(ikind))
            END IF
            IF (ASSOCIATED(basis_j(ikind)%gto_basis_set)) THEN
               j_present(ikind) = .TRUE.
               CALL get_gto_basis_set(basis_j(ikind)%gto_basis_set, kind_radius=j_radius(ikind))
            END IF
         END DO

      ELSE IF (potential_parameter%potential_type == do_potential_coulomb) THEN

         !Coulomb operator, virtually infinite range => set j_radius to arbitrarily large number
         DO ikind = 1, nkind
            IF (ASSOCIATED(basis_i(ikind)%gto_basis_set)) THEN
               i_present(ikind) = .TRUE.
               IF (pot_to_rad_prv == 1) i_radius(ikind) = 1000000.0_dp
            ENDIF
            IF (ASSOCIATED(basis_j(ikind)%gto_basis_set)) THEN
               j_present(ikind) = .TRUE.
               IF (pot_to_rad_prv == 2) j_radius(ikind) = 1000000.0_dp
            END IF
         END DO !ikind

      ELSE IF (potential_parameter%potential_type == do_potential_truncated .OR. &
               potential_parameter%potential_type == do_potential_short) THEN

         !Truncated coulomb/short range: set j_radius to r_cutoff + the kind_radii
         DO ikind = 1, nkind
            IF (ASSOCIATED(basis_i(ikind)%gto_basis_set)) THEN
               i_present(ikind) = .TRUE.
               CALL get_gto_basis_set(basis_i(ikind)%gto_basis_set, kind_radius=i_radius(ikind))
               IF (pot_to_rad_prv == 1) i_radius(ikind) = i_radius(ikind) + cutoff_screen_factor*potential_parameter%cutoff_radius
            END IF
            IF (ASSOCIATED(basis_j(ikind)%gto_basis_set)) THEN
               j_present(ikind) = .TRUE.
               CALL get_gto_basis_set(basis_j(ikind)%gto_basis_set, kind_radius=j_radius(ikind))
               IF (pot_to_rad_prv == 2) j_radius(ikind) = j_radius(ikind) + cutoff_screen_factor*potential_parameter%cutoff_radius
            END IF
         END DO

      ELSE
         CPABORT("Operator not implemented.")
      END IF

      ALLOCATE (pair_radius(nkind, nkind))
      pair_radius = 0.0_dp
      CALL pair_radius_setup(i_present, j_present, i_radius, j_radius, pair_radius)

      ALLOCATE (atom2d(nkind))

      CALL atom2d_build(atom2d, local_particles, dist_2d_prv, atomic_kind_set, &
                        molecule_set, molecule_only=.FALSE., particle_set=particle_set)
      CALL build_neighbor_lists(ij_list, particle_set, atom2d, cell, pair_radius, subcells, &
                                symmetric=sym_ij, molecular=molecular, nlname=TRIM(name))

      CALL atom2d_cleanup(atom2d)

   END SUBROUTINE

! **************************************************************************************************
!> \brief Build a 3-center neighbor list
!> \param ijk_list 3c neighbor list for atom triples i, j, k
!> \param basis_i basis object for atoms i
!> \param basis_j basis object for atoms j
!> \param basis_k basis object for atoms k
!> \param dist_3d 3d distribution object
!> \param potential_parameter ...
!> \param name name of 3c neighbor list
!> \param qs_env ...
!> \param sym_ij Symmetry in i, j (default .FALSE.)
!> \param sym_jk Symmetry in j, k (default .FALSE.)
!> \param sym_ik Symmetry in i, k (default .FALSE.)
!> \param molecular ??? not tested
!> \param op_pos ...
!> \param own_dist ...
! **************************************************************************************************
   SUBROUTINE build_3c_neighbor_lists(ijk_list, basis_i, basis_j, basis_k, &
                                      dist_3d, potential_parameter, name, qs_env, &
                                      sym_ij, sym_jk, sym_ik, molecular, op_pos, own_dist)
      TYPE(neighbor_list_3c_type), INTENT(OUT)           :: ijk_list
      TYPE(gto_basis_set_p_type), DIMENSION(:)           :: basis_i, basis_j, basis_k
      TYPE(distribution_3d_type), INTENT(IN)             :: dist_3d
      TYPE(libint_potential_type), INTENT(IN)            :: potential_parameter
      CHARACTER(LEN=*), INTENT(IN)                       :: name
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN), OPTIONAL                      :: sym_ij, sym_jk, sym_ik, molecular
      INTEGER, INTENT(IN), OPTIONAL                      :: op_pos
      LOGICAL, INTENT(IN), OPTIONAL                      :: own_dist

      CHARACTER(len=*), PARAMETER :: routineN = 'build_3c_neighbor_lists', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, op_pos_prv, sym_level
      TYPE(libint_potential_type)                        :: pot_par_1, pot_par_2

      CALL timeset(routineN, handle)

      IF (PRESENT(op_pos)) THEN
         op_pos_prv = op_pos
      ELSE
         op_pos_prv = 1
      ENDIF

      SELECT CASE (op_pos_prv)
      CASE (1)
         pot_par_1 = potential_parameter
         pot_par_2%potential_type = do_potential_id
      CASE (2)
         pot_par_2 = potential_parameter
         pot_par_1%potential_type = do_potential_id
      END SELECT

      CALL build_2c_neighbor_lists(ijk_list%ij_list, basis_i, basis_j, pot_par_1, TRIM(name)//"_sub_1", &
                                   qs_env, sym_ij=.FALSE., molecular=molecular, &
                                   dist_2d=dist_3d%dist_2d_1, pot_to_rad=1)

      CALL build_2c_neighbor_lists(ijk_list%jk_list, basis_j, basis_k, pot_par_2, TRIM(name)//"_sub_2", &
                                   qs_env, sym_ij=.FALSE., molecular=molecular, &
                                   dist_2d=dist_3d%dist_2d_2, pot_to_rad=2)

      ijk_list%sym = symmetric_none

      sym_level = 0
      IF (PRESENT(sym_ij)) THEN
         IF (sym_ij) THEN
            ijk_list%sym = symmetric_ij
            sym_level = sym_level + 1
         ENDIF
      ENDIF

      IF (PRESENT(sym_jk)) THEN
         IF (sym_jk) THEN
            ijk_list%sym = symmetric_jk
            sym_level = sym_level + 1
         ENDIF
      ENDIF

      IF (PRESENT(sym_ik)) THEN
         IF (sym_ik) THEN
            ijk_list%sym = symmetrik_ik
            sym_level = sym_level + 1
         ENDIF
      ENDIF

      IF (sym_level >= 2) THEN
         ijk_list%sym = symmetric_ijk
      ENDIF

      ijk_list%dist_3d = dist_3d
      IF (PRESENT(own_dist)) THEN
         ijk_list%owns_dist = own_dist
      ELSE
         ijk_list%owns_dist = .FALSE.
      ENDIF

      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Symmetry criterion
!> \param a ...
!> \param b ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION include_symmetric(a, b)
      INTEGER, INTENT(IN)                                :: a, b
      LOGICAL                                            :: include_symmetric

      IF (a > b) THEN
         include_symmetric = (MODULO(a + b, 2) /= 0)
      ELSE
         include_symmetric = (MODULO(a + b, 2) == 0)
      END IF

   END FUNCTION

! **************************************************************************************************
!> \brief Destroy 3c neighborlist
!> \param ijk_list ...
! **************************************************************************************************
   SUBROUTINE neighbor_list_3c_destroy(ijk_list)
      TYPE(neighbor_list_3c_type), INTENT(INOUT)         :: ijk_list

      CALL release_neighbor_list_sets(ijk_list%ij_list)
      CALL release_neighbor_list_sets(ijk_list%jk_list)

      IF (ijk_list%owns_dist) THEN
         CALL distribution_3d_destroy(ijk_list%dist_3d)
      ENDIF

   END SUBROUTINE

! **************************************************************************************************
!> \brief Create a 3-center neighborlist iterator
!> \param iterator ...
!> \param ijk_nl ...
! **************************************************************************************************
   SUBROUTINE neighbor_list_3c_iterator_create(iterator, ijk_nl)
      TYPE(neighbor_list_3c_iterator_type), INTENT(OUT)  :: iterator
      TYPE(neighbor_list_3c_type), INTENT(IN)            :: ijk_nl

      CHARACTER(len=*), PARAMETER :: routineN = 'neighbor_list_3c_iterator_create', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)
      CALL neighbor_list_iterator_create(iterator%iter_ij, ijk_nl%ij_list)
      CALL neighbor_list_iterator_create(iterator%iter_jk, ijk_nl%jk_list, search=.TRUE.)
      iterator%iter_level = 0
      iterator%ijk_nl = ijk_nl

      iterator%bounds_i = 0
      iterator%bounds_j = 0
      iterator%bounds_k = 0

      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief impose atomic upper and lower bounds
!> \param iterator ...
!> \param bounds_i ...
!> \param bounds_j ...
!> \param bounds_k ...
! **************************************************************************************************
   SUBROUTINE nl_3c_iter_set_bounds(iterator, bounds_i, bounds_j, bounds_k)
      TYPE(neighbor_list_3c_iterator_type), &
         INTENT(INOUT)                                   :: iterator
      INTEGER, DIMENSION(2), INTENT(IN), OPTIONAL        :: bounds_i, bounds_j, bounds_k

      IF (PRESENT(bounds_i)) iterator%bounds_i = bounds_i
      IF (PRESENT(bounds_j)) iterator%bounds_j = bounds_j
      IF (PRESENT(bounds_k)) iterator%bounds_k = bounds_k

   END SUBROUTINE

! **************************************************************************************************
!> \brief Destroy 3c-nl iterator
!> \param iterator ...
! **************************************************************************************************
   SUBROUTINE neighbor_list_3c_iterator_destroy(iterator)
      TYPE(neighbor_list_3c_iterator_type), &
         INTENT(INOUT)                                   :: iterator

      CHARACTER(len=*), PARAMETER :: routineN = 'neighbor_list_3c_iterator_destroy', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)
      CALL neighbor_list_iterator_release(iterator%iter_ij)
      CALL neighbor_list_iterator_release(iterator%iter_jk)
      NULLIFY (iterator%iter_ij)
      NULLIFY (iterator%iter_jk)

      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Iterate 3c-nl iterator
!> \param iterator ...
!> \return 0 if successful; 1 if end was reached
! **************************************************************************************************
   RECURSIVE FUNCTION neighbor_list_3c_iterate(iterator) RESULT(iter_stat)
      TYPE(neighbor_list_3c_iterator_type), &
         INTENT(INOUT)                                   :: iterator
      INTEGER                                            :: iter_stat

      CHARACTER(len=*), PARAMETER :: routineN = 'neighbor_list_3c_iterate', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: iatom, iter_level, jatom, jatom_1, &
                                                            jatom_2, katom
      LOGICAL                                            :: skip_this

      iter_level = iterator%iter_level

      IF (iter_level == 0) THEN
         iter_stat = neighbor_list_iterate(iterator%iter_ij)

         IF (iter_stat /= 0) THEN
            RETURN
         ENDIF

         CALL get_iterator_info(iterator%iter_ij, iatom=iatom, jatom=jatom)
         skip_this = .FALSE.
         IF ((iterator%bounds_i(1) > 0 .AND. iatom < iterator%bounds_i(1)) &
             .OR. (iterator%bounds_i(2) > 0 .AND. iatom > iterator%bounds_i(2))) skip_this = .TRUE.
         IF ((iterator%bounds_j(1) > 0 .AND. jatom < iterator%bounds_j(1)) &
             .OR. (iterator%bounds_j(2) > 0 .AND. jatom > iterator%bounds_j(2))) skip_this = .TRUE.

         IF (skip_this) THEN
            iter_stat = neighbor_list_3c_iterate(iterator)
            RETURN
         ENDIF

      ENDIF
      iter_stat = nl_sub_iterate(iterator%iter_jk, iterator%iter_ij)
      IF (iter_stat /= 0) THEN
         iterator%iter_level = 0
         iter_stat = neighbor_list_3c_iterate(iterator)
         RETURN
      ELSE
         iterator%iter_level = 1
      ENDIF

      CPASSERT(iter_stat == 0)
      CPASSERT(iterator%iter_level == 1)
      CALL get_iterator_info(iterator%iter_ij, iatom=iatom, jatom=jatom_1)
      CALL get_iterator_info(iterator%iter_jk, iatom=jatom_2, jatom=katom)

      CPASSERT(jatom_1 == jatom_2)
      jatom = jatom_1

      skip_this = .FALSE.
      IF ((iterator%bounds_k(1) > 0 .AND. katom < iterator%bounds_k(1)) &
          .OR. (iterator%bounds_k(2) > 0 .AND. katom > iterator%bounds_k(2))) skip_this = .TRUE.

      IF (skip_this) THEN
         iter_stat = neighbor_list_3c_iterate(iterator)
         RETURN
      ENDIF

      SELECT CASE (iterator%ijk_nl%sym)
      CASE (symmetric_none)
         skip_this = .FALSE.
      CASE (symmetric_ij)
         skip_this = .NOT. include_symmetric(iatom, jatom)
      CASE (symmetric_jk)
         skip_this = .NOT. include_symmetric(jatom, katom)
      CASE (symmetrik_ik)
         skip_this = .NOT. include_symmetric(iatom, katom)
      CASE (symmetric_ijk)
         skip_this = .NOT. include_symmetric(iatom, jatom) .OR. .NOT. include_symmetric(jatom, katom)
      CASE DEFAULT
         CPABORT("should not happen")
      END SELECT

      IF (skip_this) THEN
         iter_stat = neighbor_list_3c_iterate(iterator)
         RETURN
      ENDIF

   END FUNCTION

! **************************************************************************************************
!> \brief Get info of current iteration
!> \param iterator ...
!> \param ikind ...
!> \param jkind ...
!> \param kkind ...
!> \param nkind ...
!> \param iatom ...
!> \param jatom ...
!> \param katom ...
!> \param rij ...
!> \param rjk ...
!> \param rik ...
!> \param cell_j ...
!> \param cell_k ...
!> \return ...
! **************************************************************************************************
   SUBROUTINE get_3c_iterator_info(iterator, ikind, jkind, kkind, nkind, iatom, jatom, katom, &
                                   rij, rjk, rik, cell_j, cell_k)
      TYPE(neighbor_list_3c_iterator_type), &
         INTENT(INOUT)                                   :: iterator
      INTEGER, INTENT(OUT), OPTIONAL                     :: ikind, jkind, kkind, nkind, iatom, &
                                                            jatom, katom
      REAL(KIND=dp), DIMENSION(3), INTENT(OUT), OPTIONAL :: rij, rjk, rik
      INTEGER, DIMENSION(3), INTENT(OUT), OPTIONAL       :: cell_j, cell_k

      CHARACTER(len=*), PARAMETER :: routineN = 'get_3c_iterator_info', &
         routineP = moduleN//':'//routineN

      INTEGER, DIMENSION(2)                              :: atoms_1, atoms_2, kinds_1, kinds_2
      INTEGER, DIMENSION(3)                              :: cell_1, cell_2
      REAL(KIND=dp), DIMENSION(3)                        :: r_1, r_2

      CPASSERT(iterator%iter_level == 1)

      CALL get_iterator_info(iterator%iter_ij, &
                             ikind=kinds_1(1), jkind=kinds_1(2), nkind=nkind, &
                             iatom=atoms_1(1), jatom=atoms_1(2), r=r_1, &
                             cell=cell_1)

      CALL get_iterator_info(iterator%iter_jk, &
                             ikind=kinds_2(1), jkind=kinds_2(2), &
                             iatom=atoms_2(1), jatom=atoms_2(2), r=r_2, &
                             cell=cell_2)

      IF (PRESENT(ikind)) ikind = kinds_1(1)
      IF (PRESENT(jkind)) jkind = kinds_1(2)
      IF (PRESENT(kkind)) kkind = kinds_2(2)
      IF (PRESENT(iatom)) iatom = atoms_1(1)
      IF (PRESENT(jatom)) jatom = atoms_1(2)
      IF (PRESENT(katom)) katom = atoms_2(2)

      IF (PRESENT(rij)) rij = r_1
      IF (PRESENT(rjk)) rjk = r_2
      IF (PRESENT(rik)) rik = r_1 + r_2

      IF (PRESENT(cell_j)) cell_j = cell_1
      IF (PRESENT(cell_k)) cell_k = cell_1 + cell_2

   END SUBROUTINE

! **************************************************************************************************
!> \brief Allocate blocks of a 3-center tensor based on neighborlist
!> \param t3c empty DBCSR tensor
!>            Should be of shape (1,1) if no kpoints are used and of shape (nimages, nimages)
!>            if k-points are used
!> \param nl_3c 3-center neighborlist
!> \param basis_i ...
!> \param basis_j ...
!> \param basis_k ...
!> \param qs_env ...
!> \param potential_parameter ...
!> \param op_pos ...
!> \param do_kpoints ...
! **************************************************************************************************
   SUBROUTINE alloc_block_3c(t3c, nl_3c, basis_i, basis_j, basis_k, qs_env, potential_parameter, op_pos, do_kpoints)
      TYPE(dbcsr_t_type), DIMENSION(:, :), INTENT(INOUT) :: t3c
      TYPE(neighbor_list_3c_type), INTENT(INOUT)         :: nl_3c
      TYPE(gto_basis_set_p_type), DIMENSION(:)           :: basis_i, basis_j, basis_k
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(libint_potential_type), INTENT(IN)            :: potential_parameter
      INTEGER, INTENT(IN), OPTIONAL                      :: op_pos
      LOGICAL, INTENT(IN), OPTIONAL                      :: do_kpoints

      CHARACTER(LEN=*), PARAMETER :: routineN = 'alloc_block_3c', routineP = moduleN//':'//routineN

      INTEGER :: blk_cnt, handle, i, i_img, iatom, iblk, ikind, iproc, j_img, jatom, jcell, jkind, &
         katom, kcell, kkind, natom, nimg, op_ij, op_jk, op_pos_prv
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: tmp
      INTEGER, DIMENSION(3)                              :: cell_j, cell_k, kp_index_lbounds, &
                                                            kp_index_ubounds
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: do_kpoints_prv, new_block
      REAL(KIND=dp)                                      :: dij, dik, djk, dr_ij, dr_ik, dr_jk, &
                                                            kind_radius_i, kind_radius_j, &
                                                            kind_radius_k
      REAL(KIND=dp), DIMENSION(3)                        :: rij, rik, rjk
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(neighbor_list_3c_iterator_type)               :: nl_3c_iter
      TYPE(one_dim_int_array), ALLOCATABLE, &
         DIMENSION(:, :)                                 :: alloc_i, alloc_j, alloc_k
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)
      NULLIFY (qs_kind_set, atomic_kind_set)

      IF (PRESENT(do_kpoints)) THEN
         do_kpoints_prv = do_kpoints
      ELSE
         do_kpoints_prv = .FALSE.
      ENDIF

      dr_ij = 0.0_dp; dr_jk = 0.0_dp; dr_ik = 0.0_dp

      op_ij = do_potential_id; op_jk = do_potential_id

      IF (PRESENT(op_pos)) THEN
         op_pos_prv = op_pos
      ELSE
         op_pos_prv = 1
      ENDIF

      SELECT CASE (op_pos_prv)
      CASE (1)
         op_ij = potential_parameter%potential_type
      CASE (2)
         op_jk = potential_parameter%potential_type
      END SELECT

      IF (op_ij == do_potential_truncated .OR. op_ij == do_potential_short) THEN
         dr_ij = potential_parameter%cutoff_radius*cutoff_screen_factor
         dr_ik = potential_parameter%cutoff_radius*cutoff_screen_factor
      ELSEIF (op_ij == do_potential_coulomb) THEN
         dr_ij = 1000000.0_dp
         dr_ik = 1000000.0_dp
      ENDIF

      IF (op_jk == do_potential_truncated .OR. op_jk == do_potential_short) THEN
         dr_jk = potential_parameter%cutoff_radius*cutoff_screen_factor
         dr_ik = potential_parameter%cutoff_radius*cutoff_screen_factor
      ELSEIF (op_jk == do_potential_coulomb) THEN
         dr_jk = 1000000.0_dp
         dr_ik = 1000000.0_dp
      ENDIF

      CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, qs_kind_set=qs_kind_set, natom=natom, &
                      dft_control=dft_control, kpoints=kpoints, para_env=para_env)

      IF (do_kpoints_prv) THEN
         nimg = dft_control%nimages
         CALL get_kpoint_info(kpoints, cell_to_index=cell_to_index)
      ELSE
         nimg = 1
      END IF

      ALLOCATE (alloc_i(nimg, nimg))
      ALLOCATE (alloc_j(nimg, nimg))
      ALLOCATE (alloc_k(nimg, nimg))

      IF (do_kpoints_prv) THEN
         kp_index_lbounds = LBOUND(cell_to_index)
         kp_index_ubounds = UBOUND(cell_to_index)
      ENDIF

      CALL neighbor_list_3c_iterator_create(nl_3c_iter, nl_3c)
      DO WHILE (neighbor_list_3c_iterate(nl_3c_iter) == 0)
         CALL get_3c_iterator_info(nl_3c_iter, ikind=ikind, jkind=jkind, kkind=kkind, &
                                   iatom=iatom, jatom=jatom, katom=katom, &
                                   rij=rij, rjk=rjk, rik=rik, cell_j=cell_j, cell_k=cell_k)

         IF (do_kpoints_prv) THEN

            IF (ANY([cell_j(1), cell_j(2), cell_j(3)] < kp_index_lbounds) .OR. &
                ANY([cell_j(1), cell_j(2), cell_j(3)] > kp_index_ubounds)) CYCLE

            jcell = cell_to_index(cell_j(1), cell_j(2), cell_j(3))
            IF (jcell > nimg) CYCLE

            IF (ANY([cell_k(1), cell_k(2), cell_k(3)] < kp_index_lbounds) .OR. &
                ANY([cell_k(1), cell_k(2), cell_k(3)] > kp_index_ubounds)) CYCLE

            kcell = cell_to_index(cell_k(1), cell_k(2), cell_k(3))
            IF (kcell > nimg) CYCLE
         ELSE
            jcell = 1; kcell = 1
         END IF

         djk = NORM2(rjk)
         dij = NORM2(rij)
         dik = NORM2(rik)

         CALL get_gto_basis_set(basis_i(ikind)%gto_basis_set, kind_radius=kind_radius_i)
         CALL get_gto_basis_set(basis_j(jkind)%gto_basis_set, kind_radius=kind_radius_j)
         CALL get_gto_basis_set(basis_k(kkind)%gto_basis_set, kind_radius=kind_radius_k)

         IF (kind_radius_j + kind_radius_i + dr_ij < dij) CYCLE
         IF (kind_radius_j + kind_radius_k + dr_jk < djk) CYCLE
         IF (kind_radius_k + kind_radius_i + dr_ik < dik) CYCLE

         ! tensor is not symmetric therefore need to allocate rows and columns in
         ! correspondence with neighborlist. Note that this only allocates half
         ! of the blocks (since neighborlist is symmetric). After filling the blocks,
         ! tensor will be added to its transposed

         ASSOCIATE (ai=>alloc_i(jcell, kcell))
            ASSOCIATE (aj=>alloc_j(jcell, kcell))
               ASSOCIATE (ak=>alloc_k(jcell, kcell))

                  new_block = .TRUE.
                  IF (ALLOCATED(aj%array)) THEN
                     DO iblk = 1, SIZE(aj%array)
                        IF (ai%array(iblk) == iatom .AND. &
                            aj%array(iblk) == jatom .AND. &
                            ak%array(iblk) == katom) THEN
                           new_block = .FALSE.
                           EXIT
                        ENDIF
                     ENDDO
                  ENDIF
                  IF (.NOT. new_block) CYCLE

                  IF (ALLOCATED(ai%array)) THEN
                     blk_cnt = SIZE(ai%array)
                     ALLOCATE (tmp(blk_cnt))
                     tmp(:) = ai%array(:)
                     DEALLOCATE (ai%array)
                     ALLOCATE (ai%array(blk_cnt + 1))
                     ai%array(1:blk_cnt) = tmp(:)
                     ai%array(blk_cnt + 1) = iatom
                  ELSE
                     ALLOCATE (ai%array(1))
                     ai%array(1) = iatom
                  ENDIF

                  IF (ALLOCATED(aj%array)) THEN
                     tmp(:) = aj%array(:)
                     DEALLOCATE (aj%array)
                     ALLOCATE (aj%array(blk_cnt + 1))
                     aj%array(1:blk_cnt) = tmp(:)
                     aj%array(blk_cnt + 1) = jatom
                  ELSE
                     ALLOCATE (aj%array(1))
                     aj%array(1) = jatom
                  ENDIF

                  IF (ALLOCATED(ak%array)) THEN
                     tmp(:) = ak%array(:)
                     DEALLOCATE (ak%array)
                     ALLOCATE (ak%array(blk_cnt + 1))
                     ak%array(1:blk_cnt) = tmp(:)
                     ak%array(blk_cnt + 1) = katom
                  ELSE
                     ALLOCATE (ak%array(1))
                     ak%array(1) = katom
                  ENDIF

                  IF (ALLOCATED(tmp)) DEALLOCATE (tmp)
               END ASSOCIATE
            END ASSOCIATE
         END ASSOCIATE
      ENDDO

      CALL neighbor_list_3c_iterator_destroy(nl_3c_iter)

      DO i_img = 1, nimg
         DO j_img = 1, nimg
            IF (ALLOCATED(alloc_i(i_img, j_img)%array)) THEN
               DO i = 1, SIZE(alloc_i(i_img, j_img)%array)
                  CALL dbcsr_t_get_stored_coordinates(t3c(i_img, j_img), &
                                                      [alloc_i(i_img, j_img)%array(i), alloc_j(i_img, j_img)%array(i), &
                                                       alloc_k(i_img, j_img)%array(i)], &
                                                      iproc)
                  CPASSERT(iproc .EQ. para_env%mepos)
               ENDDO

               CALL dbcsr_t_reserve_blocks(t3c(i_img, j_img), &
                                           alloc_i(i_img, j_img)%array, &
                                           alloc_j(i_img, j_img)%array, &
                                           alloc_k(i_img, j_img)%array)
            ENDIF
         ENDDO
      ENDDO

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief Build 3-center integral tensor
!> \param t3c empty DBCSR tensor
!>            Should be of shape (1,1) if no kpoints are used and of shape (nimages, nimages)
!>            if k-points are used
!> \param filter_eps Filter threshold for tensor blocks
!> \param qs_env ...
!> \param nl_3c 3-center neighborlist
!> \param basis_i ...
!> \param basis_j ...
!> \param basis_k ...
!> \param potential_parameter ...
!> \param int_eps neglect integrals smaller than int_eps
!> \param op_pos operator position.
!>        1: calculate (i|jk) integrals,
!>        2: calculate (ij|k) integrals
!> \param do_kpoints ...
!> this routine requires that libint has been static initialised somewhere else
!> \param desymmetrize ...
!> \param bounds_i ...
!> \param bounds_j ...
!> \param bounds_k ...
! **************************************************************************************************
   SUBROUTINE build_3c_integrals(t3c, filter_eps, qs_env, &
                                 nl_3c, basis_i, basis_j, basis_k, &
                                 potential_parameter, &
                                 int_eps, &
                                 op_pos, do_kpoints, desymmetrize, &
                                 bounds_i, bounds_j, bounds_k)
      TYPE(dbcsr_t_type), DIMENSION(:, :), INTENT(INOUT) :: t3c
      REAL(KIND=dp), INTENT(IN)                          :: filter_eps
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(neighbor_list_3c_type), INTENT(INOUT)         :: nl_3c
      TYPE(gto_basis_set_p_type), DIMENSION(:)           :: basis_i, basis_j, basis_k
      TYPE(libint_potential_type), INTENT(IN)            :: potential_parameter
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: int_eps
      INTEGER, INTENT(IN), OPTIONAL                      :: op_pos
      LOGICAL, INTENT(IN), OPTIONAL                      :: do_kpoints, desymmetrize
      INTEGER, DIMENSION(2), INTENT(IN), OPTIONAL        :: bounds_i, bounds_j, bounds_k

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_3c_integrals', &
         routineP = moduleN//':'//routineN

      INTEGER :: block_end_i, block_end_j, block_end_k, block_start_i, block_start_j, &
         block_start_k, egfi, handle, handle2, i, iatom, ibasis, ikind, ilist, imax, iset, jatom, &
         jcell, jkind, jset, katom, kcell, kkind, kset, m_max, max_ncoi, max_ncoj, max_ncok, &
         max_nset, max_nsgfi, max_nsgfj, max_nsgfk, maxli, maxlj, maxlk, natom, nbasis, ncoi, &
         ncoj, ncok, nimg, nseti, nsetj, nsetk, op_ij, op_jk, op_pos_prv, sgfi, sgfj, sgfk, unit_id
      INTEGER, DIMENSION(3)                              :: blk_size, cell_j, cell_k, &
                                                            kp_index_lbounds, kp_index_ubounds, sp
      INTEGER, DIMENSION(:), POINTER                     :: lmax_i, lmax_j, lmax_k, lmin_i, lmin_j, &
                                                            lmin_k, npgfi, npgfj, npgfk, nsgfi, &
                                                            nsgfj, nsgfk
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgf_i, first_sgf_j, first_sgf_k
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: block_not_zero, debug, desymmetrize_prv, &
                                                            do_kpoints_prv, found
      REAL(KIND=dp)                                      :: dij, dik, djk, dr_ij, dr_ik, dr_jk, &
                                                            kind_radius_i, kind_radius_j, &
                                                            kind_radius_k, prefac, sijk_ext
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: ccp_buffer, cpp_buffer, &
                                                            max_contraction_i, max_contraction_j, &
                                                            max_contraction_k
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: block_t, dummy_block_t, sijk, sijk_contr
      REAL(KIND=dp), DIMENSION(3)                        :: ri, rij, rik, rj, rjk, rk
      REAL(KIND=dp), DIMENSION(:), POINTER               :: set_radius_i, set_radius_j, set_radius_k
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: rpgf_i, rpgf_j, rpgf_k, sphi_i, sphi_j, &
                                                            sphi_k, zeti, zetj, zetk
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_2d_r_p_type), DIMENSION(:, :), POINTER     :: spi, spk, tspj
      TYPE(cp_libint_t)                                  :: lib
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_t_type)                                 :: t_3c_tmp
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gto_basis_set_type), POINTER                  :: basis_set
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(neighbor_list_3c_iterator_type)               :: nl_3c_iter
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      debug = .FALSE.

      IF (PRESENT(do_kpoints)) THEN
         do_kpoints_prv = do_kpoints
      ELSE
         do_kpoints_prv = .FALSE.
      ENDIF

      IF (PRESENT(desymmetrize)) THEN
         desymmetrize_prv = desymmetrize
      ELSE
         desymmetrize_prv = .TRUE.
      ENDIF

      op_ij = do_potential_id; op_jk = do_potential_id

      IF (PRESENT(op_pos)) THEN
         op_pos_prv = op_pos
      ELSE
         op_pos_prv = 1
      ENDIF

      SELECT CASE (op_pos_prv)
      CASE (1)
         op_ij = potential_parameter%potential_type
      CASE (2)
         op_jk = potential_parameter%potential_type
      END SELECT

      dr_ij = 0.0_dp; dr_jk = 0.0_dp; dr_ik = 0.0_dp

      IF (op_ij == do_potential_truncated .OR. op_ij == do_potential_short) THEN
         dr_ij = potential_parameter%cutoff_radius*cutoff_screen_factor
         dr_ik = potential_parameter%cutoff_radius*cutoff_screen_factor
      ELSEIF (op_ij == do_potential_coulomb) THEN
         dr_ij = 1000000.0_dp
         dr_ik = 1000000.0_dp
      ENDIF

      IF (op_jk == do_potential_truncated .OR. op_jk == do_potential_short) THEN
         dr_jk = potential_parameter%cutoff_radius*cutoff_screen_factor
         dr_ik = potential_parameter%cutoff_radius*cutoff_screen_factor
      ELSEIF (op_jk == do_potential_coulomb) THEN
         dr_jk = 1000000.0_dp
         dr_ik = 1000000.0_dp
      ENDIF

      NULLIFY (qs_kind_set, atomic_kind_set)

      CALL alloc_block_3c(t3c, nl_3c, basis_i, basis_j, basis_k, qs_env, &
                          potential_parameter, op_pos=op_pos_prv, do_kpoints=do_kpoints)

      ! get stuff
      CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, qs_kind_set=qs_kind_set, &
                      natom=natom, kpoints=kpoints, dft_control=dft_control, para_env=para_env)

      IF (do_kpoints_prv) THEN
         nimg = dft_control%nimages
         CALL get_kpoint_info(kpoints, cell_to_index=cell_to_index)
      ELSE
         nimg = 1
      END IF

      CPASSERT(ALL(SHAPE(t3c) == [nimg, nimg]))

      !Need the max l for each basis for libint and max nset, nco and nsgf for LIBXSMM contraction
      nbasis = SIZE(basis_i)
      max_nsgfi = 0
      max_ncoi = 0
      max_nset = 0
      maxli = 0
      DO ibasis = 1, nbasis
         CALL get_gto_basis_set(gto_basis_set=basis_i(ibasis)%gto_basis_set, maxl=imax, &
                                nset=iset, nsgf_set=nsgfi, npgf=npgfi)
         maxli = MAX(maxli, imax)
         max_nset = MAX(max_nset, iset)
         max_nsgfi = MAX(max_nsgfi, MAXVAL(nsgfi))
         max_ncoi = MAX(max_ncoi, MAXVAL(npgfi)*ncoset(maxli))
      END DO
      max_nsgfj = 0
      max_ncoj = 0
      maxlj = 0
      DO ibasis = 1, nbasis
         CALL get_gto_basis_set(gto_basis_set=basis_j(ibasis)%gto_basis_set, maxl=imax, &
                                nset=jset, nsgf_set=nsgfj, npgf=npgfj)
         maxlj = MAX(maxlj, imax)
         max_nset = MAX(max_nset, jset)
         max_nsgfj = MAX(max_nsgfj, MAXVAL(nsgfj))
         max_ncoj = MAX(max_ncoj, MAXVAL(npgfj)*ncoset(maxlj))
      END DO
      max_nsgfk = 0
      max_ncok = 0
      maxlk = 0
      DO ibasis = 1, nbasis
         CALL get_gto_basis_set(gto_basis_set=basis_k(ibasis)%gto_basis_set, maxl=imax, &
                                nset=kset, nsgf_set=nsgfk, npgf=npgfk)
         maxlk = MAX(maxlk, imax)
         max_nset = MAX(max_nset, kset)
         max_nsgfk = MAX(max_nsgfk, MAXVAL(nsgfk))
         max_ncok = MAX(max_ncok, MAXVAL(npgfk)*ncoset(maxlk))
      END DO
      m_max = maxli + maxlj + maxlk

      !To minimize expensive memory opsand generally optimize contraction, pre-allocate buffers and
      !contiguous sphi arrays (and transposed in the cas of sphi_i)
      ALLOCATE (cpp_buffer(max_nsgfj*max_ncok), ccp_buffer(max_nsgfj*max_nsgfk*max_ncoi))

      NULLIFY (tspj, spi, spk)
      ALLOCATE (spi(max_nset, nbasis), tspj(max_nset, nbasis), spk(max_nset, nbasis))

      DO ibasis = 1, nbasis
         DO iset = 1, max_nset
            NULLIFY (spi(iset, ibasis)%array)
            NULLIFY (tspj(iset, ibasis)%array)

            NULLIFY (spk(iset, ibasis)%array)
         END DO
      END DO

      DO ilist = 1, 3
         DO ibasis = 1, nbasis
            IF (ilist == 1) basis_set => basis_i(ibasis)%gto_basis_set
            IF (ilist == 2) basis_set => basis_j(ibasis)%gto_basis_set
            IF (ilist == 3) basis_set => basis_k(ibasis)%gto_basis_set

            DO iset = 1, basis_set%nset

               ncoi = basis_set%npgf(iset)*ncoset(basis_set%lmax(iset))
               sgfi = basis_set%first_sgf(1, iset)
               egfi = sgfi + basis_set%nsgf_set(iset) - 1

               IF (ilist == 1) THEN
                  ALLOCATE (spi(iset, ibasis)%array(ncoi, basis_set%nsgf_set(iset)))
                  spi(iset, ibasis)%array(:, :) = basis_set%sphi(1:ncoi, sgfi:egfi)

               ELSE IF (ilist == 2) THEN
                  ALLOCATE (tspj(iset, ibasis)%array(basis_set%nsgf_set(iset), ncoi))
                  tspj(iset, ibasis)%array(:, :) = TRANSPOSE(basis_set%sphi(1:ncoi, sgfi:egfi))

               ELSE
                  ALLOCATE (spk(iset, ibasis)%array(ncoi, basis_set%nsgf_set(iset)))
                  spk(iset, ibasis)%array(:, :) = basis_set%sphi(1:ncoi, sgfi:egfi)
               END IF

            END DO !iset
         END DO !ibasis
      END DO !ilist

      !Init the truncated Coulomb operator
      IF (op_ij == do_potential_truncated .OR. op_jk == do_potential_truncated) THEN

         IF (m_max > get_lmax_init()) THEN
            IF (para_env%mepos == 0) THEN
               CALL open_file(unit_number=unit_id, file_name=potential_parameter%filename)
            END IF
            CALL init(m_max, unit_id, para_env%mepos, para_env%group)
            IF (para_env%mepos == 0) THEN
               CALL close_file(unit_id)
            END IF
         END IF
      END IF

      CALL init_md_ftable(nmax=m_max)

      CALL cp_libint_init_3eri(lib, MAX(maxli, maxlj, maxlk))
      CALL cp_libint_set_contrdepth(lib, 1)

      CALL neighbor_list_3c_iterator_create(nl_3c_iter, nl_3c)
      CALL nl_3c_iter_set_bounds(nl_3c_iter, bounds_i, bounds_j, bounds_k)

      IF (do_kpoints_prv) THEN
         kp_index_lbounds = LBOUND(cell_to_index)
         kp_index_ubounds = UBOUND(cell_to_index)
      ENDIF

      DO WHILE (neighbor_list_3c_iterate(nl_3c_iter) == 0)
         CALL get_3c_iterator_info(nl_3c_iter, ikind=ikind, jkind=jkind, kkind=kkind, &
                                   iatom=iatom, jatom=jatom, katom=katom, &
                                   rij=rij, rjk=rjk, rik=rik, cell_j=cell_j, cell_k=cell_k)

         IF (do_kpoints_prv) THEN
            prefac = 0.5_dp
         ELSEIF (nl_3c%sym == symmetric_jk) THEN
            IF (jatom == katom) THEN
               ! factor 0.5 due to double-counting of diagonal blocks
               ! (we desymmetrize by adding transpose)
               prefac = 0.5_dp
            ELSE
               prefac = 1.0_dp
            ENDIF
         ELSE
            prefac = 1.0_dp
         ENDIF

         IF (do_kpoints_prv) THEN

            IF (ANY([cell_j(1), cell_j(2), cell_j(3)] < kp_index_lbounds) .OR. &
                ANY([cell_j(1), cell_j(2), cell_j(3)] > kp_index_ubounds)) CYCLE

            jcell = cell_to_index(cell_j(1), cell_j(2), cell_j(3))
            IF (jcell > nimg) CYCLE

            IF (ANY([cell_k(1), cell_k(2), cell_k(3)] < kp_index_lbounds) .OR. &
                ANY([cell_k(1), cell_k(2), cell_k(3)] > kp_index_ubounds)) CYCLE

            kcell = cell_to_index(cell_k(1), cell_k(2), cell_k(3))
            IF (kcell > nimg) CYCLE

         ELSE
            jcell = 1; kcell = 1
         END IF

         CALL get_gto_basis_set(basis_i(ikind)%gto_basis_set, first_sgf=first_sgf_i, lmax=lmax_i, lmin=lmin_i, &
                                npgf=npgfi, nset=nseti, nsgf_set=nsgfi, pgf_radius=rpgf_i, set_radius=set_radius_i, &
                                sphi=sphi_i, zet=zeti, kind_radius=kind_radius_i)

         CALL get_gto_basis_set(basis_j(jkind)%gto_basis_set, first_sgf=first_sgf_j, lmax=lmax_j, lmin=lmin_j, &
                                npgf=npgfj, nset=nsetj, nsgf_set=nsgfj, pgf_radius=rpgf_j, set_radius=set_radius_j, &
                                sphi=sphi_j, zet=zetj, kind_radius=kind_radius_j)

         CALL get_gto_basis_set(basis_k(kkind)%gto_basis_set, first_sgf=first_sgf_k, lmax=lmax_k, lmin=lmin_k, &
                                npgf=npgfk, nset=nsetk, nsgf_set=nsgfk, pgf_radius=rpgf_k, set_radius=set_radius_k, &
                                sphi=sphi_k, zet=zetk, kind_radius=kind_radius_k)

         djk = NORM2(rjk)
         dij = NORM2(rij)
         dik = NORM2(rik)

         IF (kind_radius_j + kind_radius_i + dr_ij < dij) CYCLE
         IF (kind_radius_j + kind_radius_k + dr_jk < djk) CYCLE
         IF (kind_radius_k + kind_radius_i + dr_ik < dik) CYCLE

         ALLOCATE (max_contraction_i(nseti))
         max_contraction_i = 0.0_dp
         DO iset = 1, nseti
            sgfi = first_sgf_i(1, iset)
            max_contraction_i(iset) = MAXVAL((/(SUM(ABS(sphi_i(:, i))), i=sgfi, sgfi + nsgfi(iset) - 1)/))
         ENDDO

         ALLOCATE (max_contraction_j(nsetj))
         max_contraction_j = 0.0_dp
         DO jset = 1, nsetj
            sgfj = first_sgf_j(1, jset)
            max_contraction_j(jset) = MAXVAL((/(SUM(ABS(sphi_j(:, i))), i=sgfj, sgfj + nsgfj(jset) - 1)/))
         ENDDO

         ALLOCATE (max_contraction_k(nsetk))
         max_contraction_k = 0.0_dp
         DO kset = 1, nsetk
            sgfk = first_sgf_k(1, kset)
            max_contraction_k(kset) = MAXVAL((/(SUM(ABS(sphi_k(:, i))), i=sgfk, sgfk + nsgfk(kset) - 1)/))
         ENDDO

         CALL dbcsr_t_blk_sizes(t3c(jcell, kcell), [iatom, jatom, katom], blk_size)

         ALLOCATE (block_t(blk_size(2), blk_size(3), blk_size(1)))

         block_t = 0.0_dp
         block_not_zero = .FALSE.

         DO iset = 1, nseti

            DO jset = 1, nsetj

               IF (set_radius_j(jset) + set_radius_i(iset) + dr_ij < dij) CYCLE

               DO kset = 1, nsetk

                  IF (set_radius_j(jset) + set_radius_k(kset) + dr_jk < djk) CYCLE
                  IF (set_radius_k(kset) + set_radius_i(iset) + dr_ik < dik) CYCLE

                  ncoi = npgfi(iset)*ncoset(lmax_i(iset))
                  ncoj = npgfj(jset)*ncoset(lmax_j(jset))
                  ncok = npgfk(kset)*ncoset(lmax_k(kset))

                  sgfi = first_sgf_i(1, iset)
                  sgfj = first_sgf_j(1, jset)
                  sgfk = first_sgf_k(1, kset)

                  IF (ncoj*ncok*ncoi > 0) THEN
                     ALLOCATE (sijk(ncoj, ncok, ncoi))
                     sijk(:, :, :) = 0.0_dp
                     !need positions for libint. Only relative positions are needed => set ri to 0.0
                     ri = 0.0_dp
                     rj = rij ! ri + rij
                     rk = rik ! ri + rik
                     CALL eri_3center(sijk, &
                                      lmin_j(jset), lmax_j(jset), npgfj(jset), zetj(:, jset), rpgf_j(:, jset), rj, &
                                      lmin_k(kset), lmax_k(kset), npgfk(kset), zetk(:, kset), rpgf_k(:, kset), rk, &
                                      lmin_i(iset), lmax_i(iset), npgfi(iset), zeti(:, iset), rpgf_i(:, iset), ri, &
                                      djk, dij, dik, lib, potential_parameter, int_abc_ext=sijk_ext)

                     IF (PRESENT(int_eps)) THEN
                        IF (int_eps > sijk_ext*(max_contraction_i(iset)* &
                                                max_contraction_j(jset)* &
                                                max_contraction_k(kset))) THEN
                           DEALLOCATE (sijk)
                           CYCLE
                        END IF
                     ENDIF

                     block_not_zero = .TRUE.

                     ALLOCATE (sijk_contr(nsgfj(jset), nsgfk(kset), nsgfi(iset)))
                     CALL libxsmm_abc_contract(sijk_contr, sijk, tspj(jset, jkind)%array, &
                                               spk(kset, kkind)%array, spi(iset, ikind)%array, &
                                               ncoj, ncok, ncoi, nsgfj(jset), nsgfk(kset), &
                                               nsgfi(iset), cpp_buffer, ccp_buffer)
                     DEALLOCATE (sijk)

                     block_start_j = sgfj
                     block_end_j = sgfj + nsgfj(jset) - 1
                     block_start_k = sgfk
                     block_end_k = sgfk + nsgfk(kset) - 1
                     block_start_i = sgfi
                     block_end_i = sgfi + nsgfi(iset) - 1

                     block_t(block_start_j:block_end_j, &
                             block_start_k:block_end_k, &
                             block_start_i:block_end_i) = &
                        block_t(block_start_j:block_end_j, &
                                block_start_k:block_end_k, &
                                block_start_i:block_end_i) + &
                        prefac*sijk_contr(:, :, :)
                     DEALLOCATE (sijk_contr)

                  END IF ! number of triples > 0

               END DO

            END DO

         END DO

         IF (block_not_zero) THEN
            CALL timeset(routineN//"_put_dbcsr", handle2)
            IF (debug) THEN
               CALL dbcsr_t_get_block(t3c(jcell, kcell), &
                                      [iatom, jatom, katom], dummy_block_t, found=found)
               CPASSERT(found)
            ENDIF
            sp = SHAPE(block_t)

            sp([2, 3, 1]) = sp

            CALL dbcsr_t_put_block(t3c(jcell, kcell), &
                                   [iatom, jatom, katom], sp, RESHAPE(block_t, SHAPE=sp, ORDER=[2, 3, 1]), summation=.TRUE.)

            CALL timestop(handle2)
         ENDIF

         DEALLOCATE (block_t)

         DEALLOCATE (max_contraction_i, max_contraction_j, max_contraction_k)
      END DO

      CALL cp_libint_cleanup_3eri(lib)

      CALL neighbor_list_3c_iterator_destroy(nl_3c_iter)

      IF (nl_3c%sym == symmetric_jk .OR. do_kpoints_prv) THEN
         DO jcell = 1, nimg
            DO kcell = 1, nimg
               ! need half of filter eps because afterwards we add transposed tensor
               CALL dbcsr_t_filter(t3c(jcell, kcell), filter_eps/2)
            ENDDO
         ENDDO

         IF (desymmetrize_prv) THEN
            ! add transposed of overlap integrals
            CALL dbcsr_t_create(t3c(1, 1), t_3c_tmp)
            DO jcell = 1, nimg
               DO kcell = 1, jcell
                  CALL dbcsr_t_copy(t3c(jcell, kcell), t_3c_tmp)
                  CALL dbcsr_t_copy(t_3c_tmp, t3c(kcell, jcell), order=[1, 3, 2], summation=.TRUE., move_data=.TRUE.)
                  CALL dbcsr_t_filter(t3c(kcell, jcell), filter_eps)
               ENDDO
            ENDDO
            DO jcell = 1, nimg
               DO kcell = jcell + 1, nimg
                  CALL dbcsr_t_copy(t3c(jcell, kcell), t_3c_tmp)
                  CALL dbcsr_t_copy(t_3c_tmp, t3c(kcell, jcell), order=[1, 3, 2], summation=.FALSE., move_data=.TRUE.)
                  CALL dbcsr_t_filter(t3c(kcell, jcell), filter_eps)
               ENDDO
            ENDDO
            CALL dbcsr_t_destroy(t_3c_tmp)
         ENDIF
      ELSEIF (nl_3c%sym == symmetric_none) THEN
         DO jcell = 1, nimg
            DO kcell = 1, nimg
               CALL dbcsr_t_filter(t3c(jcell, kcell), filter_eps)
            ENDDO
         ENDDO
      ELSE
         CPABORT("requested symmetric case not implemented")
      ENDIF

      DO iset = 1, max_nset
         DO ibasis = 1, nbasis
            IF (ASSOCIATED(spi(iset, ibasis)%array)) DEALLOCATE (spi(iset, ibasis)%array)
            IF (ASSOCIATED(tspj(iset, ibasis)%array)) DEALLOCATE (tspj(iset, ibasis)%array)

            IF (ASSOCIATED(spk(iset, ibasis)%array)) DEALLOCATE (spk(iset, ibasis)%array)
         END DO
      END DO

      DEALLOCATE (spi, tspj, spk)

      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param t2c empty DBCSR matrix
!> \param filter_eps Filter threshold for matrix blocks
!> \param qs_env ...
!> \param nl_2c 2-center neighborlist
!> \param basis_i ...
!> \param basis_j ...
!> \param potential_parameter ...
!> \param do_kpoints ...
!> this routine requires that libint has been static initialised somewhere else
! **************************************************************************************************
   SUBROUTINE build_2c_integrals(t2c, filter_eps, qs_env, &
                                 nl_2c, basis_i, basis_j, &
                                 potential_parameter, do_kpoints)
      TYPE(dbcsr_type), DIMENSION(:), INTENT(INOUT)      :: t2c
      REAL(KIND=dp), INTENT(IN)                          :: filter_eps
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: nl_2c
      TYPE(gto_basis_set_p_type), DIMENSION(:)           :: basis_i, basis_j
      TYPE(libint_potential_type), INTENT(IN)            :: potential_parameter
      LOGICAL, INTENT(IN), OPTIONAL                      :: do_kpoints

      CHARACTER(len=*), PARAMETER :: routineN = 'build_2c_integrals', &
         routineP = moduleN//':'//routineN

      INTEGER :: handle, iatom, ibasis, icol, ikind, imax, img, irow, iset, jatom, jkind, jset, &
         m_max, maxli, maxlj, n1, n2, natom, ncoi, ncoj, nimg, nseti, nsetj, offi, offj, op_prv, &
         sgfi, sgfj, unit_id
      INTEGER, DIMENSION(3)                              :: cell
      INTEGER, DIMENSION(:), POINTER                     :: lmax_i, lmax_j, lmin_i, lmin_j, npgfi, &
                                                            npgfj, nsgfi, nsgfj
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgf_i, first_sgf_j
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: do_kpoints_prv, do_symmetric, found, &
                                                            trans
      REAL(KIND=dp)                                      :: dab
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: sij, sij_contr, sij_rs
      REAL(KIND=dp), DIMENSION(3)                        :: ri, rij, rj
      REAL(KIND=dp), DIMENSION(:), POINTER               :: set_radius_i, set_radius_j
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: rpgf_i, rpgf_j, scon_i, scon_j, sphi_i, &
                                                            sphi_j, zeti, zetj
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(block_p_type)                                 :: block_t
      TYPE(cp_libint_t)                                  :: lib
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      IF (PRESENT(do_kpoints)) THEN
         do_kpoints_prv = do_kpoints
      ELSE
         do_kpoints_prv = .FALSE.
      ENDIF

      op_prv = potential_parameter%potential_type

      NULLIFY (qs_kind_set, atomic_kind_set, block_t%block, cell_to_index)

      ! get stuff
      CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, qs_kind_set=qs_kind_set, &
                      natom=natom, kpoints=kpoints, dft_control=dft_control, para_env=para_env)

      IF (do_kpoints_prv) THEN
         nimg = dft_control%nimages
         CALL get_kpoint_info(kpoints, cell_to_index=cell_to_index)
      ELSE
         nimg = 1
      END IF

      CPASSERT(ALL(SHAPE(t2c) == [nimg]))

      ! check for symmetry
      CPASSERT(SIZE(nl_2c) > 0)
      CALL get_neighbor_list_set_p(neighbor_list_sets=nl_2c, symmetric=do_symmetric)

      IF (do_symmetric) THEN
         DO img = 1, nimg
            CPASSERT(dbcsr_has_symmetry(t2c(img)))
         ENDDO
      ELSE
         DO img = 1, nimg
            CPASSERT(.NOT. dbcsr_has_symmetry(t2c(img)))
         ENDDO
      ENDIF

      DO img = 1, nimg
         CALL cp_dbcsr_alloc_block_from_nbl(t2c(img), nl_2c)
      ENDDO

      maxli = 0
      DO ibasis = 1, SIZE(basis_i)
         CALL get_gto_basis_set(gto_basis_set=basis_i(ibasis)%gto_basis_set, maxl=imax)
         maxli = MAX(maxli, imax)
      END DO
      maxlj = 0
      DO ibasis = 1, SIZE(basis_j)
         CALL get_gto_basis_set(gto_basis_set=basis_j(ibasis)%gto_basis_set, maxl=imax)
         maxlj = MAX(maxlj, imax)
      END DO

      m_max = maxli + maxlj

      !Init the truncated Coulomb operator
      IF (op_prv == do_potential_truncated) THEN

         IF (m_max > get_lmax_init()) THEN
            IF (para_env%mepos == 0) THEN
               CALL open_file(unit_number=unit_id, file_name=potential_parameter%filename)
            END IF
            CALL init(m_max, unit_id, para_env%mepos, para_env%group)
            IF (para_env%mepos == 0) THEN
               CALL close_file(unit_id)
            END IF
         END IF
      END IF

      CALL init_md_ftable(nmax=m_max)

      IF (op_prv /= do_potential_id) THEN
         CALL cp_libint_init_2eri(lib, MAX(maxli, maxlj))
         CALL cp_libint_set_contrdepth(lib, 1)
      ENDIF

      CALL neighbor_list_iterator_create(nl_iterator, nl_2c)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)

         CALL get_iterator_info(nl_iterator, ikind=ikind, jkind=jkind, &
                                iatom=iatom, jatom=jatom, r=rij, cell=cell)
         IF (do_kpoints_prv) THEN
            img = cell_to_index(cell(1), cell(2), cell(3))
            IF (img > nimg) CYCLE
         ELSE
            img = 1
         END IF

         CALL get_gto_basis_set(basis_i(ikind)%gto_basis_set, first_sgf=first_sgf_i, lmax=lmax_i, lmin=lmin_i, &
                                npgf=npgfi, nset=nseti, nsgf_set=nsgfi, pgf_radius=rpgf_i, set_radius=set_radius_i, &
                                sphi=sphi_i, zet=zeti, scon=scon_i)

         CALL get_gto_basis_set(basis_j(jkind)%gto_basis_set, first_sgf=first_sgf_j, lmax=lmax_j, lmin=lmin_j, &
                                npgf=npgfj, nset=nsetj, nsgf_set=nsgfj, pgf_radius=rpgf_j, set_radius=set_radius_j, &
                                sphi=sphi_j, zet=zetj, scon=scon_j)

         IF (do_symmetric) THEN
            IF (iatom <= jatom) THEN
               irow = iatom
               icol = jatom
            ELSE
               irow = jatom
               icol = iatom
            END IF
         ELSE
            irow = iatom
            icol = jatom
         END IF

         dab = NORM2(rij)

         CALL dbcsr_get_block_p(matrix=t2c(img), &
                                row=irow, col=icol, BLOCK=block_t%block, found=found)
         CPASSERT(found)
         trans = do_symmetric .AND. (iatom > jatom)

         DO iset = 1, nseti

            ncoi = npgfi(iset)*ncoset(lmax_i(iset))
            n1 = npgfi(iset)*(ncoset(lmax_i(iset)) - ncoset(lmin_i(iset) - 1))
            sgfi = first_sgf_i(1, iset)
            offi = ncoset(lmin_i(iset) - 1) + 1

            DO jset = 1, nsetj

               ncoj = npgfj(jset)*ncoset(lmax_j(jset))
               n2 = npgfj(jset)*(ncoset(lmax_j(jset)) - ncoset(lmin_j(jset) - 1))
               sgfj = first_sgf_j(1, jset)
               offj = ncoset(lmin_j(jset) - 1) + 1

               IF (ncoi*ncoj > 0) THEN
                  ALLOCATE (sij_contr(nsgfi(iset), nsgfj(jset)))
                  sij_contr(:, :) = 0.0_dp

                  IF (op_prv == do_potential_id) THEN
                     ALLOCATE (sij(n1, n2))
                     sij(:, :) = 0.0_dp

                     CALL overlap_ab(lmax_i(iset), lmin_i(iset), npgfi(iset), rpgf_i(:, iset), zeti(:, iset), &
                                     lmax_j(jset), lmin_j(jset), npgfj(jset), rpgf_j(:, jset), zetj(:, jset), &
                                     rij, sab=sij(:, :))

                     CALL ab_contract(sij_contr, sij, &
                                      scon_i(:, sgfi:), scon_j(:, sgfj:), &
                                      n1, n2, nsgfi(iset), nsgfj(jset))

                  ELSE
                     ALLOCATE (sij(ncoi, ncoj))
                     sij(:, :) = 0.0_dp

                     ri = 0.0_dp
                     rj = rij

                     CALL eri_2center(sij, lmin_i(iset), lmax_i(iset), npgfi(iset), zeti(:, iset), &
                                      rpgf_i(:, iset), ri, lmin_j(jset), lmax_j(jset), npgfj(jset), zetj(:, jset), &
                                      rpgf_j(:, jset), rj, dab, lib, potential_parameter)

                     CALL ab_contract(sij_contr, sij, &
                                      sphi_i(:, sgfi:), sphi_j(:, sgfj:), &
                                      ncoi, ncoj, nsgfi(iset), nsgfj(jset))
                  ENDIF

                  DEALLOCATE (sij)
                  IF (trans) THEN
                     ALLOCATE (sij_rs(nsgfj(jset), nsgfi(iset)))
                     sij_rs(:, :) = TRANSPOSE(sij_contr)
                  ELSE
                     ALLOCATE (sij_rs(nsgfi(iset), nsgfj(jset)))
                     sij_rs(:, :) = sij_contr
                  ENDIF

                  DEALLOCATE (sij_contr)

                  CALL block_add("IN", sij_rs, &
                                 nsgfi(iset), nsgfj(jset), block_t%block, &
                                 sgfi, sgfj, trans=trans)
                  DEALLOCATE (sij_rs)
               ENDIF
            END DO
         END DO
      ENDDO

      IF (op_prv /= do_potential_id) THEN
         CALL cp_libint_cleanup_2eri(lib)
      ENDIF

      CALL neighbor_list_iterator_release(nl_iterator)
      DO img = 1, nimg
         CALL dbcsr_finalize(t2c(img))
         CALL dbcsr_filter(t2c(img), filter_eps)
      ENDDO

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief Change process grid of tensor, a load balanced default distribution adapted to the new grid
!>        is chosen automatically.
!> \param t_3c ...
!> \param pgrid ...
!> \param nodata ...
!> \param starts_array_mc_block_1 ...
!> \param ends_array_mc_block_1 ...
!> \param starts_array_mc_block_2 ...
!> \param ends_array_mc_block_2 ...
!> \param starts_array_mc_block_3 ...
!> \param ends_array_mc_block_3 ...
!> \param blk_indices ...
!> \param compressed ...
!> \param eps ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE tensor_change_pgrid(t_3c, pgrid, nodata, &
                                  starts_array_mc_block_1, ends_array_mc_block_1, &
                                  starts_array_mc_block_2, ends_array_mc_block_2, &
                                  starts_array_mc_block_3, ends_array_mc_block_3, &
                                  blk_indices, compressed, eps, &
                                  unit_nr)
      TYPE(dbcsr_t_type), INTENT(INOUT)                  :: t_3c
      TYPE(dbcsr_t_pgrid_type), INTENT(IN)               :: pgrid
      LOGICAL, INTENT(IN), OPTIONAL                      :: nodata
      INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: starts_array_mc_block_1, &
         ends_array_mc_block_1, starts_array_mc_block_2, ends_array_mc_block_2, &
         starts_array_mc_block_3, ends_array_mc_block_3
      INTEGER, ALLOCATABLE, DIMENSION(:, :), &
         INTENT(INOUT), OPTIONAL                         :: blk_indices
      TYPE(hfx_compression_type), INTENT(INOUT), &
         OPTIONAL                                        :: compressed
      REAL(dp), INTENT(IN), OPTIONAL                     :: eps
      INTEGER, INTENT(IN), OPTIONAL                      :: unit_nr

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tensor_change_pgrid', &
         routineP = moduleN//':'//routineN

      CHARACTER(default_string_length)                   :: name
      INTEGER                                            :: handle, imem, memcut_1, memcut_2, &
                                                            memcut_3, size_cut
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: bs1, bs2, bs3, dist1, dist2, dist3
      INTEGER, DIMENSION(3)                              :: pcoord, pcoord_ref, pdims, pdims_ref, &
                                                            tdims
      TYPE(dbcsr_t_type)                                 :: t_tmp
      TYPE(dbcsr_t_distribution_type)                    :: dist
      LOGICAL                                            :: mem_aware_1, mem_aware_2, mem_aware_3
      INTEGER, DIMENSION(dbcsr_t_ndims_matrix_row(t_3c)) :: map1
      INTEGER, &
         DIMENSION(dbcsr_t_ndims_matrix_column(t_3c))    :: map2

      CALL dbcsr_t_mp_environ_pgrid(pgrid, pdims, pcoord)
      CALL dbcsr_t_mp_environ_pgrid(t_3c%pgrid, pdims_ref, pcoord_ref)

      IF (ALL(pdims == pdims_ref)) THEN
         IF (ALLOCATED(pgrid%tas_split_info) .AND. ALLOCATED(t_3c%pgrid%tas_split_info)) THEN
            IF (pgrid%tas_split_info%ngroup == t_3c%pgrid%tas_split_info%ngroup) RETURN
         ENDIF
      ENDIF

      CALL timeset(routineN, handle)

      IF (PRESENT(compressed)) THEN
         CPASSERT(PRESENT(blk_indices))
         CPASSERT(PRESENT(eps))
         IF (PRESENT(nodata)) THEN
            CPASSERT(.NOT. nodata)
         ENDIF
      ENDIF

      mem_aware_1 = PRESENT(starts_array_mc_block_1) .AND. PRESENT(ends_array_mc_block_1)
      mem_aware_2 = PRESENT(starts_array_mc_block_2) .AND. PRESENT(ends_array_mc_block_2)
      mem_aware_3 = PRESENT(starts_array_mc_block_3) .AND. PRESENT(ends_array_mc_block_3)

      IF (mem_aware_1) memcut_1 = SIZE(starts_array_mc_block_1)
      IF (mem_aware_2) memcut_2 = SIZE(starts_array_mc_block_2)
      IF (mem_aware_3) memcut_3 = SIZE(starts_array_mc_block_3)

      ALLOCATE (bs1(dbcsr_t_nblks_total(t_3c, 1)))
      ALLOCATE (bs2(dbcsr_t_nblks_total(t_3c, 2)))
      ALLOCATE (bs3(dbcsr_t_nblks_total(t_3c, 3)))
      CALL dbcsr_t_get_info(t_3c, nblks_total=tdims, blk_size_1=bs1, blk_size_2=bs2, blk_size_3=bs3, &
                            name=name)

      ALLOCATE (dist1(tdims(1)))
      ALLOCATE (dist2(tdims(2)))
      ALLOCATE (dist3(tdims(3)))

      IF (mem_aware_1) THEN
         DO imem = 1, memcut_1
            size_cut = ends_array_mc_block_1(imem) - starts_array_mc_block_1(imem) + 1
            CALL dbcsr_t_default_distvec( &
               size_cut, pdims(1), &
               bs1(starts_array_mc_block_1(imem):ends_array_mc_block_1(imem)), &
               dist1(starts_array_mc_block_1(imem):ends_array_mc_block_1(imem)))
         ENDDO
      ELSE
         CALL dbcsr_t_default_distvec(tdims(1), pdims(1), bs1, dist1)
      ENDIF

      IF (mem_aware_2) THEN
         DO imem = 1, memcut_2
            size_cut = ends_array_mc_block_2(imem) - starts_array_mc_block_2(imem) + 1
            CALL dbcsr_t_default_distvec( &
               size_cut, pdims(2), &
               bs2(starts_array_mc_block_2(imem):ends_array_mc_block_2(imem)), &
               dist2(starts_array_mc_block_2(imem):ends_array_mc_block_2(imem)))
         ENDDO
      ELSE
         CALL dbcsr_t_default_distvec(tdims(2), pdims(2), bs2, dist2)
      ENDIF

      IF (mem_aware_3) THEN
         DO imem = 1, memcut_3
            size_cut = ends_array_mc_block_3(imem) - starts_array_mc_block_3(imem) + 1
            CALL dbcsr_t_default_distvec( &
               size_cut, pdims(3), &
               bs3(starts_array_mc_block_3(imem):ends_array_mc_block_3(imem)), &
               dist3(starts_array_mc_block_3(imem):ends_array_mc_block_3(imem)))
         ENDDO
      ELSE
         CALL dbcsr_t_default_distvec(tdims(3), pdims(3), bs3, dist3)
      ENDIF

      CALL dbcsr_t_get_mapping_info(t_3c%nd_index_blk, map1_2d=map1, map2_2d=map2)
      CALL dbcsr_t_distribution_new(dist, pgrid, dist1, dist2, dist3)
      CALL dbcsr_t_create(t_tmp, name, dist, map1, map2, dbcsr_type_real_8, bs1, bs2, bs3)
      CALL dbcsr_t_distribution_destroy(dist)

      IF (PRESENT(compressed)) THEN
         CALL decompress_tensor(t_3c, blk_indices, compressed, eps)
      ENDIF

      IF (PRESENT(nodata)) THEN
         IF (.NOT. nodata) CALL dbcsr_t_copy(t_3c, t_tmp, move_data=.TRUE.)
      ELSE
         CALL dbcsr_t_copy(t_3c, t_tmp, move_data=.TRUE.)
      ENDIF

      CALL dbcsr_t_destroy(t_3c)
      t_3c = t_tmp

      IF (PRESENT(compressed)) THEN
         IF (ALLOCATED(blk_indices)) DEALLOCATE (blk_indices)
         ALLOCATE (blk_indices(dbcsr_t_get_num_blocks(t_3c), 3))
         CALL dbcsr_t_reserved_block_indices(t_3c, blk_indices)
         CALL compress_tensor(t_3c, compressed, eps)
      ENDIF

      IF (PRESENT(unit_nr)) THEN
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, "(T2,A,1X,A)") "OPTIMIZED PGRID INFO FOR", TRIM(t_3c%name)
            WRITE (unit_nr, "(T4,A,1X,3I6)") "process grid dimensions:", pdims
            CALL dbcsr_t_write_split_info(pgrid, unit_nr)
         ENDIF
      ENDIF

      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param tensor tensor with data. Data is cleared after compression.
!> \param compressed compressed tensor data
!> \param eps all entries < eps are discarded
! **************************************************************************************************
   SUBROUTINE compress_tensor(tensor, compressed, eps)
      TYPE(dbcsr_t_type), INTENT(INOUT)                  :: tensor
      TYPE(hfx_compression_type), INTENT(INOUT)          :: compressed
      REAL(dp), INTENT(IN)                               :: eps

      INTEGER                                            :: blk, buffer_left, buffer_size, &
                                                            buffer_start, i, nbits, nints, unused
      INTEGER(int_8)                                     :: estimate_to_store_int
      INTEGER, DIMENSION(3)                              :: ind
      LOGICAL                                            :: found
      REAL(dp)                                           :: spherical_estimate
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :), TARGET  :: blk_data
      REAL(dp), DIMENSION(:), POINTER                    :: blk_data_1d
      TYPE(dbcsr_t_iterator_type)                        :: iter
      TYPE(hfx_cache_type), DIMENSION(:), POINTER        :: integral_caches
      TYPE(hfx_cache_type), POINTER                      :: maxval_cache
      TYPE(hfx_container_type), DIMENSION(:), POINTER    :: integral_containers
      TYPE(hfx_container_type), POINTER                  :: maxval_container

      CALL dealloc_containers(compressed, unused)
      CALL alloc_containers(compressed, 1)

      maxval_container => compressed%maxval_container(1)
      integral_containers => compressed%integral_containers(:, 1)

      CALL hfx_init_container(maxval_container, unused, .FALSE.)
      DO i = 1, 64
         CALL hfx_init_container(integral_containers(i), unused, .FALSE.)
      END DO

      maxval_cache => compressed%maxval_cache(1)
      integral_caches => compressed%integral_caches(:, 1)

      CALL dbcsr_t_iterator_start(iter, tensor)
      DO WHILE (dbcsr_t_iterator_blocks_left(iter))
         CALL dbcsr_t_iterator_next_block(iter, ind, blk)
         CALL dbcsr_t_get_block(tensor, ind, blk_data, found)
         CPASSERT(found)
         nints = SIZE(blk_data)
         blk_data_1d(1:nints) => blk_data
         spherical_estimate = MAXVAL(ABS(blk_data_1d))
         IF (spherical_estimate == 0.0_dp) spherical_estimate = TINY(spherical_estimate)
         estimate_to_store_int = EXPONENT(spherical_estimate)
         estimate_to_store_int = MAX(estimate_to_store_int, -15_int_8)

         CALL hfx_add_single_cache_element(estimate_to_store_int, 6, &
                                           maxval_cache, maxval_container, unused, &
                                           .FALSE.)

         spherical_estimate = SET_EXPONENT(1.0_dp, estimate_to_store_int + 1)

         nbits = EXPONENT(ANINT(spherical_estimate/eps)) + 1
         IF (nbits > 64) THEN
            CALL cp_abort(__LOCATION__, &
                          "Overflow during tensor compression. Please use a larger EPS_FILTER or EPS_STORAGE_SCALING")
         ENDIF

         buffer_left = nints
         buffer_start = 1

         DO WHILE (buffer_left > 0)
            buffer_size = MIN(buffer_left, cache_size)
            CALL hfx_add_mult_cache_elements(blk_data_1d(buffer_start:), &
                                             buffer_size, nbits, &
                                             integral_caches(nbits), &
                                             integral_containers(nbits), &
                                             eps, 1.0_dp, &
                                             unused, &
                                             .FALSE.)
            buffer_left = buffer_left - buffer_size
            buffer_start = buffer_start + buffer_size
         ENDDO

         NULLIFY (blk_data_1d); DEALLOCATE (blk_data)
      ENDDO
      CALL dbcsr_t_iterator_stop(iter)

      CALL dbcsr_t_clear(tensor)

      CALL hfx_flush_last_cache(6, maxval_cache, maxval_container, unused, &
                                .FALSE.)
      DO i = 1, 64
         CALL hfx_flush_last_cache(i, integral_caches(i), integral_containers(i), &
                                   unused, .FALSE.)
      END DO

      CALL hfx_reset_cache_and_container(maxval_cache, maxval_container, unused, .FALSE.)
      DO i = 1, 64
         CALL hfx_reset_cache_and_container(integral_caches(i), integral_containers(i), &
                                            unused, .FALSE.)
      END DO

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param tensor empty tensor which is filled by decompressed data
!> \param blk_indices indices of blocks to be reserved
!> \param compressed compressed data
!> \param eps all entries < eps are discarded
! **************************************************************************************************
   SUBROUTINE decompress_tensor(tensor, blk_indices, compressed, eps)

      TYPE(dbcsr_t_type), INTENT(INOUT)                  :: tensor
      INTEGER, DIMENSION(:, :)                           :: blk_indices
      TYPE(hfx_compression_type), INTENT(INOUT)          :: compressed
      REAL(dp), INTENT(IN)                               :: eps

      INTEGER                                            :: blk, buffer_left, buffer_size, &
                                                            buffer_start, i, nbits, nints, unused
      INTEGER(int_8)                                     :: estimate_to_store_int
      INTEGER, DIMENSION(3)                              :: blk_size, ind
      REAL(dp)                                           :: spherical_estimate
      REAL(dp), ALLOCATABLE, DIMENSION(:), TARGET        :: blk_data
      REAL(dp), DIMENSION(:, :, :), POINTER              :: blk_data_3d
      TYPE(dbcsr_t_iterator_type)                        :: iter
      TYPE(hfx_cache_type), DIMENSION(:), POINTER        :: integral_caches
      TYPE(hfx_cache_type), POINTER                      :: maxval_cache
      TYPE(hfx_container_type), DIMENSION(:), POINTER    :: integral_containers
      TYPE(hfx_container_type), POINTER                  :: maxval_container

      maxval_cache => compressed%maxval_cache(1)
      maxval_container => compressed%maxval_container(1)
      integral_caches => compressed%integral_caches(:, 1)
      integral_containers => compressed%integral_containers(:, 1)

      unused = 0

      CALL hfx_decompress_first_cache(6, maxval_cache, maxval_container, unused, .FALSE.)

      DO i = 1, 64
         CALL hfx_decompress_first_cache(i, integral_caches(i), integral_containers(i), &
                                         unused, .FALSE.)
      ENDDO

      CALL dbcsr_t_reserve_blocks(tensor, blk_indices)
      CALL dbcsr_t_iterator_start(iter, tensor)
      DO WHILE (dbcsr_t_iterator_blocks_left(iter))
         CALL dbcsr_t_iterator_next_block(iter, ind, blk, blk_size=blk_size)
         nints = PRODUCT(blk_size)
         CALL hfx_get_single_cache_element( &
            estimate_to_store_int, 6, &
            maxval_cache, maxval_container, unused, &
            .FALSE.)

         spherical_estimate = SET_EXPONENT(1.0_dp, estimate_to_store_int + 1)

         nbits = EXPONENT(ANINT(spherical_estimate/eps)) + 1

         buffer_left = nints
         buffer_start = 1

         ALLOCATE (blk_data(nints))
         DO WHILE (buffer_left > 0)
            buffer_size = MIN(buffer_left, cache_size)
            CALL hfx_get_mult_cache_elements(blk_data(buffer_start), &
                                             buffer_size, nbits, &
                                             integral_caches(nbits), &
                                             integral_containers(nbits), &
                                             eps, 1.0_dp, &
                                             unused, &
                                             .FALSE.)
            buffer_left = buffer_left - buffer_size
            buffer_start = buffer_start + buffer_size
         ENDDO

         blk_data_3d(1:blk_size(1), 1:blk_size(2), 1:blk_size(3)) => blk_data
         CALL dbcsr_t_put_block(tensor, ind, blk_size, blk_data_3d)
         NULLIFY (blk_data_3d); DEALLOCATE (blk_data)
      ENDDO

      CALL dbcsr_t_iterator_stop(iter)

      CALL hfx_reset_cache_and_container(maxval_cache, maxval_container, unused, .FALSE.)
      DO i = 1, 64
         CALL hfx_reset_cache_and_container(integral_caches(i), integral_containers(i), &
                                            unused, .FALSE.)
      END DO
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param tensor ...
!> \param nze ...
!> \param occ ...
! **************************************************************************************************
   SUBROUTINE get_tensor_occupancy(tensor, nze, occ)
      TYPE(dbcsr_t_type), INTENT(IN)                     :: tensor
      INTEGER(int_8), INTENT(OUT)                        :: nze
      REAL(dp), INTENT(OUT)                              :: occ

      INTEGER, DIMENSION(dbcsr_t_ndims(tensor))          :: dims

      nze = dbcsr_t_get_nze_total(tensor)
      CALL dbcsr_t_get_info(tensor, nfull_total=dims)
      occ = REAL(nze, dp)/PRODUCT(REAL(dims, dp))

   END SUBROUTINE

END MODULE
