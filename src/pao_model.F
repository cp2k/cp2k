!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Module for equivariant PAO-ML based on PyTorch.
!> \author Ole Schuett
! **************************************************************************************************
MODULE pao_model
   USE OMP_LIB,                         ONLY: omp_init_lock,&
                                              omp_set_lock,&
                                              omp_unset_lock
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE cell_types,                      ONLY: cell_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_get_info,&
                                              dbcsr_iterator_blocks_left,&
                                              dbcsr_iterator_next_block,&
                                              dbcsr_iterator_start,&
                                              dbcsr_iterator_stop,&
                                              dbcsr_iterator_type,&
                                              dbcsr_type
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp,&
                                              int_8,&
                                              sp
   USE message_passing,                 ONLY: mp_para_env_type
   USE pao_types,                       ONLY: pao_env_type,&
                                              pao_model_type
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: angstrom
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE torch_api,                       ONLY: &
        torch_dict_create, torch_dict_get, torch_dict_insert, torch_dict_release, torch_dict_type, &
        torch_model_forward, torch_model_get_attr, torch_model_load, torch_tensor_backward, &
        torch_tensor_data_ptr, torch_tensor_from_array, torch_tensor_grad, torch_tensor_release, &
        torch_tensor_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: pao_model_load, pao_model_predict, pao_model_forces, pao_model_type

CONTAINS

! **************************************************************************************************
!> \brief Loads a PAO-ML model.
!> \param pao ...
!> \param qs_env ...
!> \param ikind ...
!> \param pao_model_file ...
!> \param model ...
! **************************************************************************************************
   SUBROUTINE pao_model_load(pao, qs_env, ikind, pao_model_file, model)
      TYPE(pao_env_type), INTENT(IN)                     :: pao
      TYPE(qs_environment_type), INTENT(IN)              :: qs_env
      INTEGER, INTENT(IN)                                :: ikind
      CHARACTER(LEN=default_path_length), INTENT(IN)     :: pao_model_file
      TYPE(pao_model_type), INTENT(OUT)                  :: model

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pao_model_load'

      CHARACTER(LEN=default_string_length)               :: kind_name
      CHARACTER(LEN=default_string_length), &
         ALLOCATABLE, DIMENSION(:)                       :: model_kind_names
      INTEGER                                            :: handle, jkind, kkind, pao_basis_size, z
      REAL(dp)                                           :: cutoff_angstrom
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(gto_basis_set_type), POINTER                  :: basis_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)
      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set, atomic_kind_set=atomic_kind_set)

      IF (pao%iw > 0) WRITE (pao%iw, '(A)') " PAO| Loading PyTorch model from: "//TRIM(pao_model_file)
      CALL torch_model_load(model%torch_model, pao_model_file)

      ! Read model attributes.
      CALL torch_model_get_attr(model%torch_model, "pao_model_version", model%version)
      CALL torch_model_get_attr(model%torch_model, "kind_name", model%kind_name)
      CALL torch_model_get_attr(model%torch_model, "atomic_number", model%atomic_number)
      CALL torch_model_get_attr(model%torch_model, "prim_basis_name", model%prim_basis_name)
      CALL torch_model_get_attr(model%torch_model, "prim_basis_size", model%prim_basis_size)
      CALL torch_model_get_attr(model%torch_model, "pao_basis_size", model%pao_basis_size)
      CALL torch_model_get_attr(model%torch_model, "num_layers", model%num_layers)
      CALL torch_model_get_attr(model%torch_model, "cutoff", cutoff_angstrom)
      CALL torch_model_get_attr(model%torch_model, "all_kind_names", model_kind_names)
      model%cutoff = cutoff_angstrom/angstrom

      ! Freeze model after all attributes have been read.
      ! TODO Re-enable once the memory leaks of torch::jit::freeze() are fixed.
      ! https://github.com/pytorch/pytorch/issues/96726
      ! CALL torch_model_freeze(model%torch_model)

      ! For each of the model's kind names lookup the corresponding atomic kind index.
      ALLOCATE (model%kinds_mapping(SIZE(atomic_kind_set)))
      model%kinds_mapping(:) = -1
      DO jkind = 1, SIZE(atomic_kind_set)
         DO kkind = 1, SIZE(model_kind_names)
            IF (TRIM(atomic_kind_set(jkind)%name) == TRIM(model_kind_names(kkind))) THEN
               model%kinds_mapping(jkind) = kkind - 1
               EXIT
            END IF
         END DO
         IF (model%kinds_mapping(jkind) < 0) THEN
            CALL cp_abort(__LOCATION__, "PAO-ML model lacks kind '"//TRIM(atomic_kind_set(jkind)%name)//"' .")
         END IF
      END DO

      ! Check compatibility
      CALL get_qs_kind(qs_kind_set(ikind), basis_set=basis_set, pao_basis_size=pao_basis_size)
      CALL get_atomic_kind(atomic_kind_set(ikind), name=kind_name, z=z)
      IF (model%version /= 2) &
         CPABORT("Model version not supported.")
      IF (TRIM(model%kind_name) .NE. TRIM(kind_name)) &
         CPABORT("Kind name does not match.")
      IF (model%atomic_number /= z) &
         CPABORT("Atomic number does not match.")
      IF (TRIM(model%prim_basis_name) .NE. TRIM(basis_set%name)) &
         CPABORT("Primary basis set name does not match.")
      IF (model%prim_basis_size /= basis_set%nsgf) &
         CPABORT("Primary basis set size does not match.")
      IF (model%pao_basis_size /= pao_basis_size) &
         CPABORT("PAO basis size does not match.")

      CALL omp_init_lock(model%lock)
      CALL timestop(handle)

   END SUBROUTINE pao_model_load

! **************************************************************************************************
!> \brief Fills pao%matrix_X based on machine learning predictions
!> \param pao ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE pao_model_predict(pao, qs_env)
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pao_model_predict'

      INTEGER                                            :: acol, arow, handle, iatom
      REAL(dp), DIMENSION(:, :), POINTER                 :: block_X
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, handle)

!$OMP PARALLEL DEFAULT(NONE) SHARED(pao,qs_env) PRIVATE(iter,arow,acol,iatom,block_X)
      CALL dbcsr_iterator_start(iter, pao%matrix_X)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, arow, acol, block_X)
         IF (SIZE(block_X) == 0) CYCLE ! pao disabled for iatom
         iatom = arow; CPASSERT(arow == acol)
         CALL predict_single_atom(pao, qs_env, iatom, block_X=block_X)
      END DO
      CALL dbcsr_iterator_stop(iter)
!$OMP END PARALLEL

      CALL timestop(handle)

   END SUBROUTINE pao_model_predict

! **************************************************************************************************
!> \brief Calculate forces contributed by machine learning
!> \param pao ...
!> \param qs_env ...
!> \param matrix_G ...
!> \param forces ...
! **************************************************************************************************
   SUBROUTINE pao_model_forces(pao, qs_env, matrix_G, forces)
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_type)                                   :: matrix_G
      REAL(dp), DIMENSION(:, :), INTENT(INOUT)           :: forces

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pao_model_forces'

      INTEGER                                            :: acol, arow, handle, iatom
      REAL(dp), DIMENSION(:, :), POINTER                 :: block_G
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, handle)

!$OMP PARALLEL DEFAULT(NONE) SHARED(pao,qs_env,matrix_G,forces) PRIVATE(iter,arow,acol,iatom,block_G)
      CALL dbcsr_iterator_start(iter, matrix_G)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, arow, acol, block_G)
         iatom = arow; CPASSERT(arow == acol)
         IF (SIZE(block_G) == 0) CYCLE ! pao disabled for iatom
         CALL predict_single_atom(pao, qs_env, iatom, block_G=block_G, forces=forces)
      END DO
      CALL dbcsr_iterator_stop(iter)
!$OMP END PARALLEL

      CALL timestop(handle)

   END SUBROUTINE pao_model_forces

! **************************************************************************************************
!> \brief Predicts a single block_X.
!> \param pao ...
!> \param qs_env ...
!> \param iatom ...
!> \param block_X ...
!> \param block_G ...
!> \param forces ...
! **************************************************************************************************
   SUBROUTINE predict_single_atom(pao, qs_env, iatom, block_X, block_G, forces)
      TYPE(pao_env_type), INTENT(IN), POINTER            :: pao
      TYPE(qs_environment_type), INTENT(IN), POINTER     :: qs_env
      INTEGER, INTENT(IN)                                :: iatom
      REAL(dp), DIMENSION(:, :), OPTIONAL                :: block_X, block_G, forces

      INTEGER                                            :: i, iedge, ikind, j, jatom, jcell, jkind, &
                                                            jneighbor, k, katom, kneighbor, m, n, &
                                                            natoms, num_edges, num_neighbors
      INTEGER(kind=int_8), ALLOCATABLE, DIMENSION(:)     :: neighbor_atom_types
      INTEGER(kind=int_8), ALLOCATABLE, DIMENSION(:, :)  :: central_edge_index, edge_index
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: neighbor_atom_index
      INTEGER, DIMENSION(:), POINTER                     :: blk_sizes_pao, blk_sizes_pri
      REAL(dp), DIMENSION(3)                             :: Ri, Rj, Rjk
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: cell_shifts, neighbor_pos
      REAL(sp), ALLOCATABLE, DIMENSION(:, :)             :: edge_vectors
      REAL(sp), ALLOCATABLE, DIMENSION(:, :, :)          :: outer_grad
      REAL(sp), DIMENSION(:, :), POINTER                 :: edge_vectors_grad
      REAL(sp), DIMENSION(:, :, :), POINTER              :: predicted_xblock
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(pao_model_type), POINTER                      :: model
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(torch_dict_type)                              :: model_inputs, model_outputs
      TYPE(torch_tensor_type) :: atom_types_tensor, central_edge_index_tensor, edge_index_tensor, &
         edge_vectors_grad_tensor, edge_vectors_tensor, outer_grad_tensor, predicted_xblock_tensor

      CALL dbcsr_get_info(pao%matrix_Y, row_blk_size=blk_sizes_pri, col_blk_size=blk_sizes_pao)
      n = blk_sizes_pri(iatom) ! size of primary basis
      m = blk_sizes_pao(iatom) ! size of pao basis

      CALL get_qs_env(qs_env, &
                      para_env=para_env, &
                      cell=cell, &
                      particle_set=particle_set, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      natom=natoms)

      CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
      Ri = particle_set(iatom)%r
      model => pao%models(ikind)
      CPASSERT(model%version > 0)
      CALL omp_set_lock(model%lock) ! TODO: might not be needed for inference.

      ! TODO: this is a quadratic algorithm, use a neighbor-list instead.

      ! Enumerate all neighboring images. TODO: should be all images within num_layers*cutoff.
      ALLOCATE (cell_shifts(27, 3))
      jcell = 0
      DO i = -1, +1
      DO j = -1, +1
      DO k = -1, +1
         jcell = jcell + 1
         cell_shifts(jcell, :) = i*cell%hmat(1, :) + j*cell%hmat(2, :) + k*cell%hmat(3, :)
      END DO
      END DO
      END DO

      ! Find neighbors, ie. atoms that are reachable within num_layers*cutoff.
      ! 1st pass to count neighbors.
      num_neighbors = 1 ! first neighbor is always the central atom
      DO jatom = 1, natoms
      DO jcell = 1, 27
         Rj = particle_set(jatom)%r + cell_shifts(jcell, :)
         IF (NORM2(Rj - Ri) < model%num_layers*model%cutoff .AND. ANY(Rj /= Ri)) THEN
            num_neighbors = num_neighbors + 1
         END IF
      END DO
      END DO

      ! 2nd pass to collect neighbors.
      ALLOCATE (neighbor_pos(num_neighbors, 3), neighbor_atom_types(num_neighbors), neighbor_atom_index(num_neighbors))
      num_neighbors = 1 ! first neighbor is always the central atom
      neighbor_pos(1, :) = Ri
      neighbor_atom_types(1) = model%kinds_mapping(ikind)
      neighbor_atom_index(1) = iatom
      DO jatom = 1, natoms
      DO jcell = 1, 27
         Rj = particle_set(jatom)%r + cell_shifts(jcell, :)
         jkind = particle_set(jatom)%atomic_kind%kind_number
         IF (NORM2(Rj - Ri) < model%num_layers*model%cutoff .AND. ANY(Rj /= Ri)) THEN
            num_neighbors = num_neighbors + 1
            neighbor_pos(num_neighbors, :) = Rj
            neighbor_atom_types(num_neighbors) = model%kinds_mapping(jkind)
            neighbor_atom_index(num_neighbors) = jatom
         END IF
      END DO
      END DO

      ! Build connectivity graph of neighbors.
      ! 1st pass to count edges.
      num_edges = 0
      DO jneighbor = 1, num_neighbors
      DO kneighbor = 1, num_neighbors
         Rjk = neighbor_pos(kneighbor, :) - neighbor_pos(jneighbor, :)
         IF (NORM2(Rjk) < model%cutoff .AND. jneighbor /= kneighbor) THEN
            num_edges = num_edges + 1
         END IF
      END DO
      END DO

      ! 2nd pass to collect edges.
      ALLOCATE (edge_index(num_edges, 2), edge_vectors(3, num_edges)) ! edge_index is transposed
      num_edges = 0
      DO jneighbor = 1, num_neighbors
      DO kneighbor = 1, num_neighbors
         Rjk = neighbor_pos(kneighbor, :) - neighbor_pos(jneighbor, :)
         IF (NORM2(Rjk) < model%cutoff .AND. jneighbor /= kneighbor) THEN
            num_edges = num_edges + 1
            edge_index(num_edges, :) = [jneighbor - 1, kneighbor - 1]
            edge_vectors(:, num_edges) = REAL(Rjk*angstrom, kind=sp)
         END IF
      END DO
      END DO

      ALLOCATE (central_edge_index(1, 2))
      central_edge_index(:, :) = 0

      ! Inference.
      CALL torch_dict_create(model_inputs)

      CALL torch_tensor_from_array(atom_types_tensor, neighbor_atom_types)
      CALL torch_dict_insert(model_inputs, "atom_types", atom_types_tensor)

      CALL torch_tensor_from_array(edge_index_tensor, edge_index)
      CALL torch_dict_insert(model_inputs, "edge_index", edge_index_tensor)

      CALL torch_tensor_from_array(edge_vectors_tensor, edge_vectors, requires_grad=PRESENT(block_G))
      CALL torch_dict_insert(model_inputs, "edge_vectors", edge_vectors_tensor)

      CALL torch_tensor_from_array(central_edge_index_tensor, central_edge_index)
      CALL torch_dict_insert(model_inputs, "central_edge_index", central_edge_index_tensor)

      CALL torch_dict_create(model_outputs)
      CALL torch_model_forward(model%torch_model, model_inputs, model_outputs)

      ! Copy predicted XBlock.
      NULLIFY (predicted_xblock)
      CALL torch_dict_get(model_outputs, "xblock", predicted_xblock_tensor)
      CALL torch_tensor_data_ptr(predicted_xblock_tensor, predicted_xblock)
      CPASSERT(SIZE(predicted_xblock, 1) == n .AND. SIZE(predicted_xblock, 2) == m .AND. SIZE(predicted_xblock, 3) == 1)
      IF (PRESENT(block_X)) THEN
         block_X = RESHAPE(predicted_xblock, [n*m, 1])
      END IF

      ! TURNING POINT (if calc forces) ------------------------------------------
      IF (PRESENT(block_G)) THEN
         ALLOCATE (outer_grad(n, m, 1))
         outer_grad(:, :, :) = REAL(RESHAPE(block_G, [n, m, 1]), kind=sp)
         CALL torch_tensor_from_array(outer_grad_tensor, outer_grad)
         CALL torch_tensor_backward(predicted_xblock_tensor, outer_grad_tensor)
         CALL torch_tensor_grad(edge_vectors_tensor, edge_vectors_grad_tensor)
         NULLIFY (edge_vectors_grad)
         CALL torch_tensor_data_ptr(edge_vectors_grad_tensor, edge_vectors_grad)
         CPASSERT(SIZE(edge_vectors_grad, 1) == 3 .AND. SIZE(edge_vectors_grad, 2) == num_edges)
         DO iedge = 1, num_edges
            jneighbor = INT(edge_index(iedge, 1) + 1)
            kneighbor = INT(edge_index(iedge, 2) + 1)
            jatom = neighbor_atom_index(jneighbor)
            katom = neighbor_atom_index(kneighbor)
            forces(jatom, :) = forces(jatom, :) + edge_vectors_grad(:, iedge)*angstrom
            forces(katom, :) = forces(katom, :) - edge_vectors_grad(:, iedge)*angstrom
         END DO
         CALL torch_tensor_release(outer_grad_tensor)
         CALL torch_tensor_release(edge_vectors_grad_tensor)
      END IF

      ! Clean up.
      CALL torch_tensor_release(atom_types_tensor)
      CALL torch_tensor_release(edge_index_tensor)
      CALL torch_tensor_release(edge_vectors_tensor)
      CALL torch_tensor_release(central_edge_index_tensor)
      CALL torch_tensor_release(predicted_xblock_tensor)
      CALL torch_dict_release(model_inputs)
      CALL torch_dict_release(model_outputs)
      CALL omp_unset_lock(model%lock)

   END SUBROUTINE predict_single_atom

END MODULE pao_model
