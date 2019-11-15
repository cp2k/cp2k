!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Utility methods to build 3-center integral tensors of various types.
! **************************************************************************************************
MODULE qs_3c_tensors
   USE ai_contraction_sphi,             ONLY: abc_contract
   USE ai_overlap3,                     ONLY: overlap3
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind_set
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_p_type
   USE cell_types,                      ONLY: cell_type
   USE cp_array_utils,                  ONLY: cp_1d_i_p_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                              cp_blacs_env_release,&
                                              cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_para_env,                     ONLY: cp_para_env_create,&
                                              cp_para_env_release
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_tensor_api,                ONLY: &
        dbcsr_t_blk_sizes, dbcsr_t_copy, dbcsr_t_create, dbcsr_t_destroy, dbcsr_t_filter, &
        dbcsr_t_get_block, dbcsr_t_get_stored_coordinates, dbcsr_t_put_block, &
        dbcsr_t_reserve_blocks, dbcsr_t_type
   USE distribution_1d_types,           ONLY: distribution_1d_type
   USE distribution_2d_types,           ONLY: distribution_2d_create_prv => distribution_2d_create,&
                                              distribution_2d_release,&
                                              distribution_2d_type
   USE input_constants,                 ONLY: do_potential_coulomb,&
                                              do_potential_id,&
                                              do_potential_short,&
                                              do_potential_truncated
   USE input_section_types,             ONLY: section_vals_val_get
   USE kinds,                           ONLY: dp
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_type
   USE message_passing,                 ONLY: mp_cart_sub,&
                                              mp_comm_free,&
                                              mp_environ
   USE molecule_types,                  ONLY: molecule_type
   USE orbital_pointers,                ONLY: ncoset
   USE particle_types,                  ONLY: particle_type
   USE qs_3c_tensors_types,             ONLY: distribution_3d_type,&
                                              neighbor_list_3c_iterator_type,&
                                              neighbor_list_3c_type,&
                                              symmetric_ij,&
                                              symmetric_ijk,&
                                              symmetric_jk,&
                                              symmetric_none,&
                                              symmetrik_ik
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_neighbor_list_types,          ONLY: deallocate_neighbor_list_set,&
                                              get_iterator_info,&
                                              neighbor_list_iterate,&
                                              neighbor_list_iterator_create,&
                                              neighbor_list_iterator_release,&
                                              neighbor_list_set_p_type,&
                                              nl_sub_iterate
   USE qs_neighbor_lists,               ONLY: atom2d_build,&
                                              atom2d_cleanup,&
                                              build_neighbor_lists,&
                                              local_atoms_type,&
                                              pair_radius_setup
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_3c_tensors'

   PUBLIC :: distribution_3d_create, distribution_3d_destroy, build_3c_neighbor_lists, &
             neighbor_list_3c_destroy, neighbor_list_3c_iterate, neighbor_list_3c_iterator_create, &
             neighbor_list_3c_iterator_destroy, get_3c_iterator_info, build_3c_integrals, &
             contiguous_tensor_dist, cyclic_tensor_dist

   TYPE one_dim_int_array
      INTEGER, DIMENSION(:), ALLOCATABLE    :: array
   END TYPE

CONTAINS
! **************************************************************************************************
!> \brief Create a 2d distribution. This mainly wraps distribution_2d_create
!>        for consistency with distribution_3d_create.
!> \param dist_2d 2d distribution object
!> \param dist1 distribution vector along 1st process grid dimension
!> \param dist2 distribution vector along 2nd process grid dimension
!> \param nkind ...
!> \param particle_set ...
!> \param mp_comm_2d MPI communicator with a 3d cartesian topology
! **************************************************************************************************
   SUBROUTINE distribution_2d_create(dist_2d, dist1, dist2, nkind, particle_set, mp_comm_2d)
      TYPE(distribution_2d_type), POINTER                :: dist_2d
      INTEGER, DIMENSION(:), INTENT(IN)                  :: dist1, dist2
      INTEGER, INTENT(IN)                                :: nkind
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      INTEGER, INTENT(IN)                                :: mp_comm_2d

      INTEGER                                            :: iatom, ikind, n, natom, nproc
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: nparticle_local_col, nparticle_local_row
      INTEGER, DIMENSION(2)                              :: mp_coor, mp_dims
      INTEGER, DIMENSION(:, :), POINTER                  :: dist1_prv, dist2_prv
      TYPE(cp_1d_i_p_type), DIMENSION(:), POINTER        :: local_particle_col, local_particle_row
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_para_env_type), POINTER                    :: para_env

      NULLIFY (blacs_env, local_particle_col, local_particle_row, para_env)

      CALL mp_environ(nproc, mp_dims, mp_coor, mp_comm_2d)
      CALL cp_para_env_create(para_env, mp_comm_2d, owns_group=.FALSE.)
      CALL cp_blacs_env_create(blacs_env, para_env, &
                               grid_2d=mp_dims)

      CPASSERT(blacs_env%mepos(1) == mp_coor(1))
      CPASSERT(blacs_env%mepos(2) == mp_coor(2))

      CALL cp_para_env_release(para_env)

      natom = SIZE(particle_set)
      ALLOCATE (dist1_prv(natom, 2), dist2_prv(natom, 2))
      dist1_prv(:, 1) = dist1
      dist2_prv(:, 1) = dist2

      ALLOCATE (local_particle_col(nkind), local_particle_row(nkind))
      ALLOCATE (nparticle_local_row(nkind), nparticle_local_col(nkind))
      nparticle_local_row = 0; nparticle_local_col = 0

      DO iatom = 1, natom
         ikind = particle_set(iatom)%atomic_kind%kind_number

         IF (dist1_prv(iatom, 1) == mp_coor(1)) nparticle_local_row(ikind) = nparticle_local_row(ikind) + 1
         IF (dist2_prv(iatom, 1) == mp_coor(2)) nparticle_local_col(ikind) = nparticle_local_col(ikind) + 1
      END DO

      DO ikind = 1, nkind
         n = nparticle_local_row(ikind)
         ALLOCATE (local_particle_row(ikind)%array(n))

         n = nparticle_local_col(ikind)
         ALLOCATE (local_particle_col(ikind)%array(n))
      END DO

      nparticle_local_row = 0; nparticle_local_col = 0
      DO iatom = 1, natom
         ikind = particle_set(iatom)%atomic_kind%kind_number

         IF (dist1_prv(iatom, 1) == mp_coor(1)) THEN
            nparticle_local_row(ikind) = nparticle_local_row(ikind) + 1
            local_particle_row(ikind)%array(nparticle_local_row(ikind)) = iatom
         END IF
         IF (dist2_prv(iatom, 1) == mp_coor(2)) THEN
            nparticle_local_col(ikind) = nparticle_local_col(ikind) + 1
            local_particle_col(ikind)%array(nparticle_local_col(ikind)) = iatom
         END IF
      END DO

      CALL distribution_2d_create_prv(dist_2d, row_distribution_ptr=dist1_prv, &
                                      col_distribution_ptr=dist2_prv, local_rows_ptr=local_particle_row, &
                                      local_cols_ptr=local_particle_col, blacs_env=blacs_env)

      CALL cp_blacs_env_release(blacs_env)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Create a 3d distribution
!> \param dist_3d 3d distribution object
!> \param dist1 distribution vector along 1st process grid dimension
!> \param dist2 distribution vector along 2nd process grid dimension
!> \param dist3 distribution vector along 3rd process grid dimension
!> \param nkind ...
!> \param particle_set ...
!> \param mp_comm_3d MPI communicator with a 3d cartesian topology
!> \param own_comm Whether mp_comm_3d should be owned by dist_3d (default false)
! **************************************************************************************************
   SUBROUTINE distribution_3d_create(dist_3d, dist1, dist2, dist3, nkind, particle_set, mp_comm_3d, own_comm)
      TYPE(distribution_3d_type)                         :: dist_3d
      INTEGER, DIMENSION(:), INTENT(IN)                  :: dist1, dist2, dist3
      INTEGER, INTENT(IN)                                :: nkind
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      INTEGER, INTENT(IN)                                :: mp_comm_3d
      LOGICAL, INTENT(IN), OPTIONAL                      :: own_comm

      CHARACTER(len=*), PARAMETER :: routineN = 'distribution_3d_create', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: comm_2d_1, comm_2d_2, handle, nproc_1, &
                                                            nproc_2
      INTEGER, DIMENSION(2)                              :: mp_coor_1, mp_coor_2, mp_dims_1, &
                                                            mp_dims_2

      CALL timeset(routineN, handle)

      IF (PRESENT(own_comm)) THEN
         IF (own_comm) dist_3d%comm_3d = mp_comm_3d
         dist_3d%owns_comm = own_comm
      ELSE
         dist_3d%owns_comm = .FALSE.
      ENDIF

      CALL mp_cart_sub(mp_comm_3d, [.TRUE., .TRUE., .FALSE.], comm_2d_1)
      CALL mp_cart_sub(mp_comm_3d, [.FALSE., .TRUE., .TRUE.], comm_2d_2)

      CALL mp_environ(nproc_1, mp_dims_1, mp_coor_1, comm_2d_1)
      CALL mp_environ(nproc_2, mp_dims_2, mp_coor_2, comm_2d_2)

      CPASSERT(mp_coor_1(2) == mp_coor_2(1))

      CALL distribution_2d_create(dist_3d%dist_2d_1, dist1, dist2, nkind, particle_set, comm_2d_1)
      CALL distribution_2d_create(dist_3d%dist_2d_2, dist2, dist3, nkind, particle_set, comm_2d_2)

      dist_3d%comm_2d_1 = comm_2d_1
      dist_3d%comm_2d_2 = comm_2d_2

      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Destroy a 3d distribution
!> \param dist ...
! **************************************************************************************************
   SUBROUTINE distribution_3d_destroy(dist)
      TYPE(distribution_3d_type)                         :: dist

      CHARACTER(len=*), PARAMETER :: routineN = 'distribution_3d_destroy', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)
      CALL distribution_2d_release(dist%dist_2d_1)
      CALL distribution_2d_release(dist%dist_2d_2)
      CALL mp_comm_free(dist%comm_2d_1)
      CALL mp_comm_free(dist%comm_2d_2)
      IF (dist%owns_comm) CALL mp_comm_free(dist%comm_3d)

      NULLIFY (dist%dist_2d_1, dist%dist_2d_2)

      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Build 2-center neighborlists adapted to different operators
!>        This mainly wraps build_neighbor_lists for consistency with build_3c_neighbor_lists
!> \param ij_list 2c neighbor list for atom pairs i, j
!> \param basis_i basis object for atoms i
!> \param basis_j basis object for atoms j
!> \param name name of 2c neighbor list
!> \param qs_env ...
!> \param sym_ij Symmetry in i, j (default .TRUE.)
!> \param molecular ...
!> \param op_ij operator
!> \param x_range in case some truncated/screened operator is used, gives its range
!> \param dist_2d optionally a custom 2d distribution
! **************************************************************************************************
   SUBROUTINE build_2c_neighbor_lists(ij_list, basis_i, basis_j, name, qs_env, &
                                      sym_ij, molecular, op_ij, x_range, dist_2d)
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: ij_list
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_i, basis_j
      CHARACTER(LEN=*), INTENT(IN)                       :: name
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN), OPTIONAL                      :: sym_ij, molecular
      INTEGER, INTENT(IN), OPTIONAL                      :: op_ij
      REAL(dp), INTENT(IN), OPTIONAL                     :: x_range
      TYPE(distribution_2d_type), OPTIONAL, POINTER      :: dist_2d

      CHARACTER(len=*), PARAMETER :: routineN = 'build_2c_neighbor_lists', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: ikind, nkind, op_ij_prv
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

      IF (PRESENT(op_ij)) THEN
         op_ij_prv = op_ij
      ELSE
         op_ij_prv = do_potential_id
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
      IF (op_ij_prv == do_potential_id) THEN

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

      ELSE IF (op_ij_prv == do_potential_coulomb) THEN

         !Coulomb operator, virtually infinite range => set j_radius to arbitrarily large number
         DO ikind = 1, nkind
            IF (ASSOCIATED(basis_j(ikind)%gto_basis_set)) THEN
               j_present(ikind) = .TRUE.
               j_radius(ikind) = 1000000.0_dp
            END IF
            IF (ASSOCIATED(basis_i(ikind)%gto_basis_set)) i_present(ikind) = .TRUE.
         END DO !ikind

      ELSE IF (op_ij_prv == do_potential_truncated .OR. op_ij_prv == do_potential_short) THEN

         !Truncated coulomb/short range: set j_radius to x_range + the kind_radii
         DO ikind = 1, nkind
            IF (ASSOCIATED(basis_i(ikind)%gto_basis_set)) THEN
               i_present(ikind) = .TRUE.
               CALL get_gto_basis_set(basis_i(ikind)%gto_basis_set, kind_radius=i_radius(ikind))
            END IF
            IF (ASSOCIATED(basis_j(ikind)%gto_basis_set)) THEN
               j_present(ikind) = .TRUE.
               CALL get_gto_basis_set(basis_j(ikind)%gto_basis_set, kind_radius=j_radius(ikind))
               j_radius(ikind) = j_radius(ikind) + x_range
            END IF
         END DO

      ELSE
         CPABORT("Operator not known")
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
!> \param name name of 3c neighbor list
!> \param qs_env ...
!> \param sym_ij Symmetry in i, j (default .FALSE.)
!> \param sym_jk Symmetry in j, k (default .FALSE.)
!> \param sym_ik Symmetry in i, k (default .FALSE.)
!> \param molecular ??? not tested
!> \param op_ij operator between i and j
!> \param op_jk operator between j and k
!> \param x_range in case some truncated/screened operator is used, gives its range
!> \param own_dist ...
! **************************************************************************************************
   SUBROUTINE build_3c_neighbor_lists(ijk_list, basis_i, basis_j, basis_k, &
                                      dist_3d, name, qs_env, sym_ij, sym_jk, sym_ik, molecular, op_ij, op_jk, &
                                      x_range, own_dist)
      TYPE(neighbor_list_3c_type), INTENT(OUT)           :: ijk_list
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_i, basis_j, basis_k
      TYPE(distribution_3d_type), INTENT(IN)             :: dist_3d
      CHARACTER(LEN=*), INTENT(IN)                       :: name
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN), OPTIONAL                      :: sym_ij, sym_jk, sym_ik, molecular
      INTEGER, INTENT(IN), OPTIONAL                      :: op_ij, op_jk
      REAL(dp), INTENT(IN), OPTIONAL                     :: x_range
      LOGICAL, INTENT(IN), OPTIONAL                      :: own_dist

      CHARACTER(len=*), PARAMETER :: routineN = 'build_3c_neighbor_lists', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, sym_level

      CALL timeset(routineN, handle)

      CALL build_2c_neighbor_lists(ijk_list%ij_list, basis_i, basis_j, TRIM(name)//"_sub_1", &
                                   qs_env, sym_ij=.FALSE., molecular=molecular, op_ij=op_ij, x_range=x_range, &
                                   dist_2d=dist_3d%dist_2d_1)

      CALL build_2c_neighbor_lists(ijk_list%jk_list, basis_j, basis_k, TRIM(name)//"_sub_2", &
                                   qs_env, sym_ij=.FALSE., molecular=molecular, op_ij=op_jk, x_range=x_range, &
                                   dist_2d=dist_3d%dist_2d_2)

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

      CHARACTER(len=*), PARAMETER :: routineN = 'neighbor_list_3c_destroy', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i

      CALL timeset(routineN, handle)
      DO i = 1, SIZE(ijk_list%ij_list)
         CALL deallocate_neighbor_list_set(ijk_list%ij_list(i)%neighbor_list_set)
      END DO
      DO i = 1, SIZE(ijk_list%jk_list)
         CALL deallocate_neighbor_list_set(ijk_list%jk_list(i)%neighbor_list_set)
      END DO
      DEALLOCATE (ijk_list%ij_list, ijk_list%jk_list)

      IF (ijk_list%owns_dist) THEN
         CALL distribution_3d_destroy(ijk_list%dist_3d)
      ENDIF

      CALL timestop(handle)
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

      CALL timestop(handle)
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

      INTEGER                                            :: handle, iatom, iter_level, jatom, &
                                                            jatom_1, jatom_2, katom
      LOGICAL                                            :: skip_this

      CALL timeset(routineN, handle)

      iter_level = iterator%iter_level

      IF (iter_level == 0) THEN
         iter_stat = neighbor_list_iterate(iterator%iter_ij)

         IF (iter_stat /= 0) THEN
            CALL timestop(handle)
            RETURN
         ENDIF
      ENDIF
      iter_stat = nl_sub_iterate(iterator%iter_jk, iterator%iter_ij)
      IF (iter_stat /= 0) THEN
         iterator%iter_level = 0
         iter_stat = neighbor_list_3c_iterate(iterator)
         CALL timestop(handle)
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

      skip_this = .TRUE.
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
         CALL timestop(handle)
         RETURN
      ENDIF

      CALL timestop(handle)
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

      INTEGER                                            :: handle
      INTEGER, DIMENSION(2)                              :: atoms_1, atoms_2, kinds_1, kinds_2
      INTEGER, DIMENSION(3)                              :: cell_1, cell_2
      REAL(KIND=dp), DIMENSION(3)                        :: r_1, r_2

      CALL timeset(routineN, handle)

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

      CALL timestop(handle)
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
!> \param do_kpoints ...
! **************************************************************************************************
   SUBROUTINE alloc_block_3c(t3c, nl_3c, basis_i, basis_j, basis_k, qs_env, do_kpoints)
      TYPE(dbcsr_t_type), DIMENSION(:, :), INTENT(INOUT) :: t3c
      TYPE(neighbor_list_3c_type), INTENT(INOUT)         :: nl_3c
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_i, basis_j, basis_k
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN), OPTIONAL                      :: do_kpoints

      CHARACTER(LEN=*), PARAMETER :: routineN = 'alloc_block_3c', routineP = moduleN//':'//routineN

      INTEGER                                            :: blk_cnt, handle, i, i_img, iatom, iblk, &
                                                            ikind, iproc, j_img, jatom, jcell, &
                                                            jkind, katom, kcell, kkind, natom, nimg
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: kind_of, tmp
      INTEGER, DIMENSION(3)                              :: cell_j, cell_k
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: do_kpoints_prv, new_block
      REAL(KIND=dp)                                      :: dij, dik, djk, kind_radius_i, &
                                                            kind_radius_j, kind_radius_k
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

      ALLOCATE (kind_of(natom))
      CALL get_atomic_kind_set(atomic_kind_set, kind_of=kind_of)

      CALL neighbor_list_3c_iterator_create(nl_3c_iter, nl_3c)
      DO WHILE (neighbor_list_3c_iterate(nl_3c_iter) == 0)
         CALL get_3c_iterator_info(nl_3c_iter, iatom=iatom, jatom=jatom, katom=katom, &
                                   rij=rij, rjk=rjk, rik=rik, cell_j=cell_j, cell_k=cell_k)

         IF (do_kpoints_prv) THEN

            jcell = cell_to_index(cell_j(1), cell_j(2), cell_j(3))
            IF (jcell > nimg) CYCLE

            kcell = cell_to_index(cell_k(1), cell_k(2), cell_k(3))
            IF (kcell > nimg) CYCLE
         ELSE
            jcell = 1; kcell = 1
         END IF

         djk = SQRT(SUM(rjk**2))
         dij = SQRT(SUM(rij**2))
         dik = SQRT(SUM(rik**2))

         ikind = kind_of(iatom)
         jkind = kind_of(jatom)
         kkind = kind_of(katom)
         CALL get_gto_basis_set(basis_i(ikind)%gto_basis_set, kind_radius=kind_radius_i)
         CALL get_gto_basis_set(basis_j(jkind)%gto_basis_set, kind_radius=kind_radius_j)
         CALL get_gto_basis_set(basis_k(kkind)%gto_basis_set, kind_radius=kind_radius_k)

         IF (kind_radius_j + kind_radius_i < dij) CYCLE
         IF (kind_radius_j + kind_radius_k < djk) CYCLE
         IF (kind_radius_k + kind_radius_i < dik) CYCLE

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
!> \param op operator, overlap integrals by default, other operators NYI
!> \param do_kpoints ...
! **************************************************************************************************
   SUBROUTINE build_3c_integrals(t3c, filter_eps, qs_env, &
                                 nl_3c, basis_i, basis_j, basis_k, op, do_kpoints)

      TYPE(dbcsr_t_type), DIMENSION(:, :), INTENT(INOUT) :: t3c
      REAL(KIND=dp), INTENT(IN)                          :: filter_eps
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(neighbor_list_3c_type), INTENT(INOUT)         :: nl_3c
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_i, basis_j, basis_k
      INTEGER, INTENT(IN), OPTIONAL                      :: op
      LOGICAL, INTENT(IN), OPTIONAL                      :: do_kpoints

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_3c_integrals', &
         routineP = moduleN//':'//routineN

      INTEGER :: block_end_i, block_end_j, block_end_k, block_start_i, block_start_j, &
         block_start_k, handle, handle2, iatom, ikind, iset, jatom, jcell, jkind, jset, katom, &
         kcell, kkind, kset, natom, ncoi, ncoj, ncok, nimg, nseti, nsetj, nsetk, op_prv, sgfi, &
         sgfj, sgfk
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: kind_of
      INTEGER, DIMENSION(3)                              :: blk_size, cell_j, cell_k
      INTEGER, DIMENSION(:), POINTER                     :: lmax_i, lmax_j, lmax_k, lmin_i, lmin_j, &
                                                            lmin_k, npgfi, npgfj, npgfk, nsgfi, &
                                                            nsgfj, nsgfk
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgf_i, first_sgf_j, first_sgf_k
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: do_kpoints_prv, found
      REAL(KIND=dp)                                      :: dij, dik, djk, kind_radius_i, &
                                                            kind_radius_j, kind_radius_k, prefac
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: block_t, dummy_block_t, sijk, sijk_contr
      REAL(KIND=dp), DIMENSION(3)                        :: rij, rik, rjk
      REAL(KIND=dp), DIMENSION(:), POINTER               :: set_radius_i, set_radius_j, set_radius_k
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: rpgf_i, rpgf_j, rpgf_k, sphi_i, sphi_j, &
                                                            sphi_k, zeti, zetj, zetk
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(dbcsr_t_type)                                 :: t_3c_tmp
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(neighbor_list_3c_iterator_type)               :: nl_3c_iter
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      IF (PRESENT(do_kpoints)) THEN
         do_kpoints_prv = do_kpoints
      ELSE
         do_kpoints_prv = .FALSE.
      ENDIF

      IF (PRESENT(op)) THEN
         op_prv = op
      ELSE
         op_prv = do_potential_id
      ENDIF

      NULLIFY (qs_kind_set, atomic_kind_set)

      CALL alloc_block_3c(t3c, nl_3c, basis_i, basis_j, basis_k, qs_env, do_kpoints=do_kpoints)

      ! get stuff
      CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, qs_kind_set=qs_kind_set, &
                      natom=natom, kpoints=kpoints, dft_control=dft_control)

      IF (do_kpoints_prv) THEN
         nimg = dft_control%nimages
         CALL get_kpoint_info(kpoints, cell_to_index=cell_to_index)
      ELSE
         nimg = 1
      END IF

      ALLOCATE (kind_of(natom))

      CALL get_atomic_kind_set(atomic_kind_set, kind_of=kind_of)

      CALL neighbor_list_3c_iterator_create(nl_3c_iter, nl_3c)
      DO WHILE (neighbor_list_3c_iterate(nl_3c_iter) == 0)
         CALL get_3c_iterator_info(nl_3c_iter, iatom=iatom, jatom=jatom, katom=katom, &
                                   rij=rij, rjk=rjk, rik=rik, cell_j=cell_j, cell_k=cell_k)

         IF (do_kpoints_prv) THEN

            jcell = cell_to_index(cell_j(1), cell_j(2), cell_j(3))
            IF (jcell > nimg) CYCLE

            kcell = cell_to_index(cell_k(1), cell_k(2), cell_k(3))
            IF (kcell > nimg) CYCLE

         ELSE
            jcell = 1; kcell = 1
         END IF

         ikind = kind_of(iatom)
         jkind = kind_of(jatom)
         kkind = kind_of(katom)
         CALL get_gto_basis_set(basis_i(ikind)%gto_basis_set, first_sgf=first_sgf_i, lmax=lmax_i, lmin=lmin_i, &
                                npgf=npgfi, nset=nseti, nsgf_set=nsgfi, pgf_radius=rpgf_i, set_radius=set_radius_i, &
                                sphi=sphi_i, zet=zeti, kind_radius=kind_radius_i)

         CALL get_gto_basis_set(basis_j(jkind)%gto_basis_set, first_sgf=first_sgf_j, lmax=lmax_j, lmin=lmin_j, &
                                npgf=npgfj, nset=nsetj, nsgf_set=nsgfj, pgf_radius=rpgf_j, set_radius=set_radius_j, &
                                sphi=sphi_j, zet=zetj, kind_radius=kind_radius_j)

         CALL get_gto_basis_set(basis_k(kkind)%gto_basis_set, first_sgf=first_sgf_k, lmax=lmax_k, lmin=lmin_k, &
                                npgf=npgfk, nset=nsetk, nsgf_set=nsgfk, pgf_radius=rpgf_k, set_radius=set_radius_k, &
                                sphi=sphi_k, zet=zetk, kind_radius=kind_radius_k)

         djk = SQRT(SUM(rjk**2))
         dij = SQRT(SUM(rij**2))
         dik = SQRT(SUM(rik**2))

         IF (kind_radius_j + kind_radius_i < dij) CYCLE
         IF (kind_radius_j + kind_radius_k < djk) CYCLE
         IF (kind_radius_k + kind_radius_i < dik) CYCLE

         CALL dbcsr_t_blk_sizes(t3c(jcell, kcell), [iatom, jatom, katom], blk_size)
         ALLOCATE (block_t(blk_size(1), blk_size(2), blk_size(3)))
         block_t = 0.0_dp

         DO iset = 1, nseti

            DO jset = 1, nsetj

               IF (set_radius_j(jset) + set_radius_i(iset) < dij) CYCLE

               DO kset = 1, nsetk

                  IF (set_radius_j(jset) + set_radius_k(kset) < djk) CYCLE
                  IF (set_radius_k(kset) + set_radius_i(iset) < dik) CYCLE

                  ncoi = npgfi(iset)*ncoset(lmax_i(iset))
                  ncoj = npgfj(jset)*ncoset(lmax_j(jset))
                  ncok = npgfk(kset)*ncoset(lmax_k(kset))

                  sgfi = first_sgf_i(1, iset)
                  sgfj = first_sgf_j(1, jset)
                  sgfk = first_sgf_k(1, kset)

                  IF (ncoj*ncok*ncoi > 0) THEN

                     ALLOCATE (sijk(ncoi, ncoj, ncok))
                     sijk(:, :, :) = 0.0_dp

                     SELECT CASE (op_prv)
                     CASE (do_potential_id)
                        CALL overlap3(lmax_i(iset), npgfi(iset), zeti(:, iset), rpgf_i(:, iset), &
                                      lmin_i(iset), &
                                      lmax_j(jset), npgfj(jset), zetj(:, jset), rpgf_j(:, jset), lmin_j(jset), &
                                      lmax_k(kset), npgfk(kset), zetk(:, kset), rpgf_k(:, kset), lmin_k(kset), &
                                      rij, dij, rik, dik, rjk, djk, sijk)
                     CASE DEFAULT
                        CPABORT("Operators other than identity not implemented")
                     END SELECT

                     ALLOCATE (sijk_contr(nsgfi(iset), nsgfj(jset), nsgfk(kset)))

                     CALL abc_contract(sijk_contr, sijk, &
                                       sphi_i(:, sgfi:), sphi_j(:, sgfj:), sphi_k(:, sgfk:), &
                                       ncoi, ncoj, ncok, nsgfi(iset), nsgfj(jset), nsgfk(kset))

                     block_start_j = sgfj
                     block_end_j = sgfj + nsgfj(jset) - 1
                     block_start_k = sgfk
                     block_end_k = sgfk + nsgfk(kset) - 1
                     block_start_i = sgfi
                     block_end_i = sgfi + nsgfi(iset) - 1

                     IF (do_kpoints_prv) THEN
                        prefac = 0.5_dp
                     ELSE
                        IF (jatom == katom) THEN
                           ! factor 0.5 due to double-counting of diagonal blocks
                           ! (we desymmetrize by adding transpose)
                           prefac = 0.5_dp
                        ELSE
                           prefac = 1.0_dp
                        ENDIF

                     ENDIF

                     block_t(block_start_i:block_end_i, &
                             block_start_j:block_end_j, &
                             block_start_k:block_end_k) = &
                        block_t(block_start_i:block_end_i, &
                                block_start_j:block_end_j, &
                                block_start_k:block_end_k) + &
                        prefac*sijk_contr(:, :, :)
                     DEALLOCATE (sijk, sijk_contr)

                  END IF ! number of triples > 0

               END DO

            END DO

         END DO
         CALL timeset(routineN//"_put_dbcsr", handle2)

         CALL dbcsr_t_get_block(t3c(jcell, kcell), &
                                [iatom, jatom, katom], dummy_block_t, found=found)
         CPASSERT(found)

         CALL dbcsr_t_put_block(t3c(jcell, kcell), &
                                [iatom, jatom, katom], SHAPE(block_t), block_t, summation=.TRUE.)

         CALL timestop(handle2)

         DEALLOCATE (block_t)

      END DO

      CALL neighbor_list_3c_iterator_destroy(nl_3c_iter)

      DO jcell = 1, nimg
         DO kcell = 1, nimg
            ! need half of filter eps because afterwards we add transposed tensor
            CALL dbcsr_t_filter(t3c(jcell, kcell), filter_eps/2)
         ENDDO
      ENDDO

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

      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief contiguous distribution of weighted elements
!> \param nel ...
!> \param nbin ...
!> \param weights ...
!> \param limits_start ...
!> \param limits_end ...
!> \param dist ...
! **************************************************************************************************
   SUBROUTINE contiguous_tensor_dist(nel, nbin, weights, limits_start, limits_end, dist)
      INTEGER, INTENT(IN)                                :: nel, nbin
      INTEGER, DIMENSION(nel), INTENT(IN)                :: weights
      INTEGER, DIMENSION(nbin), INTENT(OUT), OPTIONAL    :: limits_start, limits_end
      INTEGER, DIMENSION(nel), INTENT(OUT), OPTIONAL     :: dist

      INTEGER                                            :: el_end, el_start, end_weight, ibin, &
                                                            nel_div, nel_rem, nel_split, nel_w, &
                                                            w_partialsum

      nel_w = SUM(weights)
      nel_div = nel_w/nbin
      nel_rem = MOD(nel_w, nbin)

      w_partialsum = 0
      el_end = 0
      end_weight = 0
      DO ibin = 1, nbin
         nel_split = nel_div
         IF (ibin <= nel_rem) THEN
            nel_split = nel_split + 1
         ENDIF
         el_start = el_end + 1
         el_end = el_start
         w_partialsum = w_partialsum + weights(el_end)
         end_weight = end_weight + nel_split
         DO WHILE (w_partialsum < end_weight)
            !IF (ABS(w_partialsum + weights(el_end) - end_weight) > ABS(w_partialsum - end_weight)) EXIT
            el_end = el_end + 1
            w_partialsum = w_partialsum + weights(el_end)
         ENDDO
         IF (PRESENT(dist)) dist(el_start:el_end) = ibin - 1
         IF (PRESENT(limits_start)) limits_start(ibin) = el_start
         IF (PRESENT(limits_end)) limits_end(ibin) = el_end
      ENDDO

   END SUBROUTINE contiguous_tensor_dist

! **************************************************************************************************
!> \brief cyclic distribution of weighted elements
!> \param nel ...
!> \param nbin ...
!> \param weights ...
!> \param dist ...
! **************************************************************************************************
   SUBROUTINE cyclic_tensor_dist(nel, nbin, weights, dist)
      INTEGER, INTENT(IN)                                :: nel, nbin
      INTEGER, DIMENSION(nel), INTENT(IN)                :: weights
      INTEGER, DIMENSION(nel), INTENT(OUT)               :: dist

      INTEGER                                            :: ibin, iel, niter
      INTEGER, DIMENSION(nbin)                           :: occup

      occup(:) = 0
      ibin = 0
      DO iel = 1, nel
         niter = 0
         ibin = MOD(ibin + 1, nbin)
         DO WHILE (occup(ibin + 1) + weights(iel) .GE. MAXVAL(occup))
            IF (MINLOC(occup, DIM=1) == ibin + 1) EXIT
            ibin = MOD(ibin + 1, nbin)
            niter = niter + 1
         ENDDO
         dist(iel) = ibin
         occup(ibin + 1) = occup(ibin + 1) + weights(iel)
      ENDDO

   END SUBROUTINE cyclic_tensor_dist

END MODULE
