!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
!> \brief Native Gaussian Bloch kinetic integrals and their stationary geometry response.
MODULE qs_upf_inverse_gaussian
   USE cell_types,                      ONLY: cell_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_copy,&
                                              dbcsr_get_block_p,&
                                              dbcsr_get_info,&
                                              dbcsr_get_matrix_type,&
                                              dbcsr_p_type,&
                                              dbcsr_set,&
                                              dbcsr_type_symmetric
   USE cp_dbcsr_operations,             ONLY: dbcsr_deallocate_matrix_set
   USE kinds,                           ONLY: dp
   USE kpoint_transform,                ONLY: kp_transform_plan_create
   USE kpoint_types,                    ONLY: kpoint_type
   USE mathconstants,                   ONLY: twopi
   USE mathlib,                         ONLY: invmat
   USE particle_types,                  ONLY: particle_type
   USE qs_kinetic,                      ONLY: build_kinetic_matrix
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_env_type
   USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                              neighbor_list_iterate,&
                                              neighbor_list_iterator_create,&
                                              neighbor_list_iterator_p_type,&
                                              neighbor_list_iterator_release,&
                                              neighbor_list_set_p_type
   USE qs_upf_inverse_types,            ONLY: qs_upf_inverse_gaussian_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: qs_upf_inverse_gaussian_prepare, qs_upf_inverse_gaussian_kinetic
   PUBLIC :: qs_upf_inverse_gaussian_operator, gaussian_images
   PUBLIC :: qs_upf_inverse_gaussian_response
CONTAINS

! **************************************************************************************************
!> \brief Own native image sparsity and phases in the original-coordinate Gaussian AO basis.
!> \param ks_env Native overlap images, bases and neighbor lists at the current geometry.
!> \param map Owned image occupations and phase map, rebuilt after geometry changes.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_gaussian_prepare(ks_env, map)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      TYPE(qs_upf_inverse_gaussian_type), INTENT(OUT)    :: map

      CHARACTER(LEN=1)                                   :: symmetry
      INTEGER                                            :: a, b, ENTRY, image, info, native_cell(3)
      INTEGER, POINTER                                   :: sizes(:)
      REAL(dp)                                           :: image_error, inverse_lattice(3, 3), &
                                                            rab(3), tolerance
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_p_type), POINTER                        :: overlap(:, :)
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(neighbor_list_iterator_p_type), POINTER       :: iterator(:)
      TYPE(neighbor_list_set_p_type), POINTER            :: sab(:)
      TYPE(particle_type), POINTER                       :: atoms(:)

      CALL get_ks_env(ks_env, matrix_s_kp=overlap, kpoints=kpoints, sab_orb=sab)
      CPASSERT(ASSOCIATED(overlap) .AND. ASSOCIATED(sab))
      CALL dbcsr_get_info(overlap(1, 1)%matrix, row_blk_size=sizes, group=map%comm)
      map%nao = SUM(sizes)
      CPASSERT(map%nao > 0)
      ALLOCATE (map%density(SIZE(overlap, 2)))
      DO image = 1, SIZE(overlap, 2)
         symmetry = dbcsr_get_matrix_type(overlap(1, image)%matrix)
         CPASSERT(symmetry == dbcsr_type_symmetric)
         ALLOCATE (map%density(image)%matrix)
         CALL dbcsr_copy(map%density(image)%matrix, overlap(1, image)%matrix)
         CALL dbcsr_set(map%density(image)%matrix, 0.0_dp)
      END DO
      CPASSERT(ASSOCIATED(kpoints))
      IF (kpoints%symmetry .OR. .NOT. kpoints%full_grid) THEN
         CPABORT("USPP Gaussian integrals require an unreduced mesh")
      END IF
      CALL kp_transform_plan_create(map%plan, sab, kpoints%cell_to_index, SIZE(overlap, 2), &
                                    block_template=overlap(1, 1)%matrix)
      CPASSERT(map%plan%symmetric)
      CALL get_ks_env(ks_env, cell=cell, particle_set=atoms)
      inverse_lattice = cell%hmat
      CALL invmat(inverse_lattice, info)
      CPASSERT(info == 0)
      ! Native image matrices use folded atoms. Convert only their AO phase cells
      ! to the original coordinates used by the inverse fields; image ownership stays native.
      ENTRY = 0
      CALL neighbor_list_iterator_create(iterator, sab)
      DO WHILE (neighbor_list_iterate(iterator) == 0)
         CALL get_iterator_info(iterator, iatom=a, jatom=b, cell=native_cell, r=rab)
         image = kpoints%cell_to_index(native_cell(1), native_cell(2), native_cell(3))
         IF (image < 1 .OR. image > SIZE(overlap, 2)) CYCLE
         ENTRY = ENTRY + 1
         map%plan%cell(:, ENTRY) = NINT(MATMUL(inverse_lattice, atoms(a)%r - atoms(b)%r + rab))
         image_error = NORM2(atoms(a)%r - atoms(b)%r + rab - &
                             MATMUL(cell%hmat, REAL(map%plan%cell(:, ENTRY), dp)))
         CPASSERT(image_error < 1.e-9_dp)
      END DO
      CALL neighbor_list_iterator_release(iterator)
      CPASSERT(ENTRY == map%plan%nentry)
      ALLOCATE (map%reduced_k(3, kpoints%nkp), map%weights(kpoints%nkp))
      map%reduced_k(:, :) = kpoints%xkp
      map%weights(:) = kpoints%wkp
      CPASSERT(ALL(map%weights >= 0.0_dp))
      tolerance = 128.0_dp*EPSILON(1.0_dp)*SIZE(map%weights)
      CPASSERT(ABS(SUM(map%weights) - 1.0_dp) <= tolerance)
   END SUBROUTINE qs_upf_inverse_gaussian_prepare

! **************************************************************************************************
!> \brief Assemble native Gaussian kinetic integrals without k-point weights.
!> \param ks_env Native untransformed Gaussian image kinetic matrices.
!> \param map Original-coordinate phases and native image ownership.
!> \param matrix Replicated Hermitian Gaussian kinetic matrices indexed by global k.
!> \param k_response Optional reduced-vector derivative, yielding dT instead of T.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_gaussian_kinetic(ks_env, map, matrix, k_response)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      TYPE(qs_upf_inverse_gaussian_type), INTENT(IN)     :: map
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: matrix(:, :, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: k_response(:, :)

      TYPE(dbcsr_p_type), POINTER                        :: source(:, :)

      CALL get_ks_env(ks_env, kinetic_kp=source)
      CPASSERT(ASSOCIATED(source))
      CALL qs_upf_inverse_gaussian_operator(source(1, :), map, matrix, k_response)
   END SUBROUTINE qs_upf_inverse_gaussian_kinetic

! **************************************************************************************************
!> \brief Transform native image operators into the Gaussian Bloch basis.
!> \param source ...
!> \param map ...
!> \param matrix ...
!> \param k_response ...
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_gaussian_operator(source, map, matrix, k_response)
      TYPE(dbcsr_p_type), INTENT(IN)                     :: source(:)
      TYPE(qs_upf_inverse_gaussian_type), INTENT(IN)     :: map
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: matrix(:, :, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: k_response(:, :)

      COMPLEX(dp)                                        :: phase
      INTEGER                                            :: c, ENTRY, image, k, r
      LOGICAL                                            :: found
      REAL(dp)                                           :: angle_response
      REAL(dp), POINTER                                  :: block(:, :)

      CPASSERT(ASSOCIATED(map%density))
      CPASSERT(SIZE(source) == SIZE(map%density))
      IF (PRESENT(k_response)) THEN
         CPASSERT(ALL(SHAPE(k_response) == SHAPE(map%reduced_k)))
      END IF
      ALLOCATE (matrix(map%nao, map%nao, SIZE(map%weights)))
      matrix = CMPLX(0.0_dp, 0.0_dp, dp)
      DO ENTRY = 1, map%plan%nentry
         image = map%plan%image(ENTRY)
         CALL dbcsr_get_block_p(source(image)%matrix, map%plan%row(ENTRY), map%plan%col(ENTRY), block, found=found)
         IF (.NOT. found) CYCLE
         r = map%plan%row_offset(ENTRY)
         c = map%plan%col_offset(ENTRY)
         DO k = 1, SIZE(map%weights)
            phase = gaussian_phase(map, ENTRY, k)
            IF (PRESENT(k_response)) THEN
               angle_response = twopi*map%plan%symmetry_sign(ENTRY)*SUM(k_response(:, k)*map%plan%cell(:, ENTRY))
               phase = CMPLX(0.0_dp, angle_response, dp)*phase
            END IF
            CALL add_operator(block, phase, r, c, matrix(:, :, k))
         END DO
      END DO
      DO k = 1, SIZE(map%weights)
         CALL map%comm%sum(matrix(:, :, k))
      END DO
   END SUBROUTINE qs_upf_inverse_gaussian_operator

! **************************************************************************************************
!> \brief Contract native Gaussian kinetic derivatives at accepted Bloch AO occupations.
!>        Native density, forces, virial and operators are unchanged. Mesh weights occur once.
!> \param ks_env Native bases and neighbor lists at the current geometry.
!> \param map Owned native image occupations and original-coordinate phases.
!> \param density Accepted Hermitian occupations, without mesh weights.
!> \param atom_response Cartesian displacement of each Gaussian center.
!> \param cell_response Lattice-vector displacement.
!> \param k_response Reduced Bloch-vector displacement, zero for fixed-k forces and stress.
!> \param response Replicated Gaussian kinetic energy response.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_gaussian_response(ks_env, map, density, atom_response, cell_response, k_response, response)
      TYPE(qs_ks_env_type), INTENT(IN), TARGET           :: ks_env
      TYPE(qs_upf_inverse_gaussian_type), INTENT(INOUT)  :: map
      COMPLEX(dp), INTENT(IN)                            :: density(:, :, :)
      REAL(dp), INTENT(IN)                               :: atom_response(:, :), &
                                                            cell_response(3, 3), k_response(:, :)
      REAL(dp), INTENT(OUT)                              :: response

      COMPLEX(dp), ALLOCATABLE                           :: derivative(:, :, :)
      INTEGER                                            :: atom, image, info, k
      REAL(dp)                                           :: cell_inverse(3, 3), strain(3, 3), &
                                                            virial(3, 3)
      REAL(dp), ALLOCATABLE                              :: gradients(:, :), relative_motion(:, :)
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_p_type), POINTER                        :: occupations(:, :), operators(:, :)
      TYPE(neighbor_list_set_p_type), POINTER            :: sab(:)
      TYPE(particle_type), POINTER                       :: atoms(:)
      TYPE(qs_ks_env_type), POINTER                      :: native_ks

      native_ks => ks_env
      CALL gaussian_images(map, density)
      CALL get_ks_env(ks_env, cell=cell, particle_set=atoms, sab_orb=sab)
      CPASSERT(ASSOCIATED(sab))
      ALLOCATE (gradients(3, SIZE(atoms)), relative_motion(3, SIZE(atoms)))
      CPASSERT(ALL(SHAPE(atom_response) == SHAPE(gradients)))
      CPASSERT(ALL(SHAPE(k_response) == SHAPE(map%reduced_k)))
      cell_inverse = cell%hmat
      CALL invmat(cell_inverse, info)
      CPASSERT(info == 0)
      strain = MATMUL(cell_response, cell_inverse)
      DO atom = 1, SIZE(atoms)
         relative_motion(:, atom) = atom_response(:, atom) - MATMUL(strain, atoms(atom)%r)
      END DO
      NULLIFY (operators)
      ALLOCATE (occupations(1, SIZE(map%density)))
      DO image = 1, SIZE(map%density)
         occupations(1, image)%matrix => map%density(image)%matrix
      END DO
      CALL build_kinetic_matrix(native_ks, matrixkp_t=operators, basis_type="ORB", sab_nl=sab, &
                                calculate_forces=.TRUE., matrixkp_p=occupations, &
                                force_response=gradients, virial_response=virial)
      CALL map%comm%sum(gradients)
      CALL map%comm%sum(virial)
      response = SUM(gradients*relative_motion) - SUM(virial*strain)
      IF (ANY(k_response /= 0.0_dp)) THEN
         CALL qs_upf_inverse_gaussian_kinetic(ks_env, map, derivative, k_response)
         DO k = 1, SIZE(map%weights)
            response = response + map%weights(k)*REAL(SUM(density(:, :, k)*TRANSPOSE(derivative(:, :, k))), dp)
         END DO
      END IF
      CALL dbcsr_deallocate_matrix_set(operators)
      DEALLOCATE (occupations)
   END SUBROUTINE qs_upf_inverse_gaussian_response

! **************************************************************************************************
!> \brief Populate owned real image occupations from the same accepted Hermitian Bloch density.
!> \param map ...
!> \param matrix ...
! **************************************************************************************************
   SUBROUTINE gaussian_images(map, matrix)
      TYPE(qs_upf_inverse_gaussian_type), INTENT(INOUT)  :: map
      COMPLEX(dp), INTENT(IN)                            :: matrix(:, :, :)

      COMPLEX(dp)                                        :: phase
      INTEGER                                            :: c, ENTRY, image, k, nc, nr, r
      LOGICAL                                            :: found
      REAL(dp)                                           :: defect, tolerance
      REAL(dp), POINTER                                  :: block(:, :)

      CPASSERT(ASSOCIATED(map%density))
      CPASSERT(ALL(SHAPE(matrix) == [map%nao, map%nao, SIZE(map%weights)]))
      DO k = 1, SIZE(map%weights)
         tolerance = 128.0_dp*EPSILON(1.0_dp)*MAX(1.0_dp, MAXVAL(ABS(matrix(:, :, k))))
         defect = MAXVAL(ABS(matrix(:, :, k) - CONJG(TRANSPOSE(matrix(:, :, k)))))
         CPASSERT(defect <= tolerance)
      END DO
      DO image = 1, SIZE(map%density)
         CALL dbcsr_set(map%density(image)%matrix, 0.0_dp)
      END DO
      DO ENTRY = 1, map%plan%nentry
         image = map%plan%image(ENTRY)
         CALL dbcsr_get_block_p(map%density(image)%matrix, map%plan%row(ENTRY), map%plan%col(ENTRY), &
                                block, found=found)
         IF (.NOT. found) CYCLE
         r = map%plan%row_offset(ENTRY)
         c = map%plan%col_offset(ENTRY)
         nr = SIZE(block, 1)
         nc = SIZE(block, 2)
         DO k = 1, SIZE(map%weights)
            phase = gaussian_phase(map, ENTRY, k)
            block = block + map%weights(k)*REAL(CONJG(phase)*matrix(r:r + nr - 1, c:c + nc - 1, k), dp)
         END DO
      END DO
   END SUBROUTINE gaussian_images

! **************************************************************************************************
!> \brief Native image phase including symmetric-block orientation.
!> \param map ...
!> \param ENTRY ...
!> \param k ...
!> \return ...
! **************************************************************************************************
   FUNCTION gaussian_phase(map, ENTRY, k) RESULT(phase)
      TYPE(qs_upf_inverse_gaussian_type), INTENT(IN)     :: map
      INTEGER, INTENT(IN)                                :: entry, k
      COMPLEX(dp)                                        :: phase

      REAL(dp)                                           :: angle

      angle = twopi*SUM(map%reduced_k(:, k)*map%plan%cell(:, ENTRY))
      phase = CMPLX(COS(angle), map%plan%symmetry_sign(ENTRY)*SIN(angle), dp)
   END FUNCTION gaussian_phase

! **************************************************************************************************
!> \brief Assemble a Hermitian operator from a stored symmetric native image block.
!> \param block ...
!> \param phase ...
!> \param r ...
!> \param c ...
!> \param ao_operator ...
! **************************************************************************************************
   SUBROUTINE add_operator(block, phase, r, c, ao_operator)
      REAL(dp), INTENT(IN)                               :: block(:, :)
      COMPLEX(dp), INTENT(IN)                            :: phase
      INTEGER, INTENT(IN)                                :: r, c
      COMPLEX(dp), INTENT(INOUT)                         :: ao_operator(:, :)

      INTEGER                                            :: i, j
      REAL(dp)                                           :: factor

      factor = 1.0_dp
      IF (r == c) factor = 0.5_dp
      DO j = 1, SIZE(block, 2)
         DO i = 1, SIZE(block, 1)
            ao_operator(r + i - 1, c + j - 1) = ao_operator(r + i - 1, c + j - 1) + factor*phase*block(i, j)
            ao_operator(c + j - 1, r + i - 1) = ao_operator(c + j - 1, r + i - 1) + factor*CONJG(phase)*block(i, j)
         END DO
      END DO
   END SUBROUTINE add_operator

END MODULE qs_upf_inverse_gaussian
