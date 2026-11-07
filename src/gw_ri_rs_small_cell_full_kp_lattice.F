!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief RI lattice matrices, real-space cells and data distribution shared by χ, W, Σ and the
!>        Gamma correction of periodic RI-RS GW.
!> \author Jan Wilhelm (with help of ChatGPT, GPT-6 Sol and GPT-6 Astra)
! **************************************************************************************************
MODULE gw_ri_rs_small_cell_full_kp_lattice
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_cfm_types,                    ONLY: cp_cfm_get_info,&
                                              cp_cfm_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_create,&
                                              dbcsr_distribution_release,&
                                              dbcsr_distribution_type,&
                                              dbcsr_release,&
                                              dbcsr_type,&
                                              dbcsr_type_no_symmetry
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              cp_dbcsr_dist2d_to_dist
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_type
   USE distribution_2d_types,           ONLY: distribution_2d_type
   USE input_constants,                 ONLY: do_potential_short,&
                                              do_potential_truncated
   USE integral_library_types,          ONLY: coulomb_operator_type
   USE kinds,                           ONLY: dp,&
                                              dp_size,&
                                              int_8
   USE kpoint_types,                    ONLY: kpoint_type
   USE mathconstants,                   ONLY: pi,&
                                              twopi,&
                                              z_zero
   USE message_passing,                 ONLY: mp_para_env_type
   USE orbital_pointers,                ONLY: indco,&
                                              ncoset
   USE particle_methods,                ONLY: get_particle_set
   USE particle_types,                  ONLY: particle_type
   USE post_scf_bandstructure_types,    ONLY: post_scf_bandstructure_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                              neighbor_list_iterate,&
                                              neighbor_list_iterator_create,&
                                              neighbor_list_iterator_p_type,&
                                              neighbor_list_iterator_release,&
                                              neighbor_list_set_p_type,&
                                              release_neighbor_list_sets
   USE qs_tensors,                      ONLY: build_2c_integrals,&
                                              build_2c_neighbor_lists
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gw_ri_rs_small_cell_full_kp_lattice'

   PUBLIC :: assemble_AA_contraction, &
             block_range, &
             build_full_cell_grid, &
             build_ri_2c_lattice, &
             complex_bytes, &
             gamma_auxiliary_charges, &
             owned_mo_coefficients, &
             real_bytes, &
             ri_2c_lattice_to_kp, &
             ri_metric_head_zero

   ! Bytes per stored real and complex number.
   REAL(KIND=dp), PARAMETER :: real_bytes = REAL(dp_size, dp), complex_bytes = 2.0_dp*real_bytes

CONTAINS

! **************************************************************************************************
!> \brief Assembles the q-independent, atom-diagonal AA contraction from stored coefficients.
!> \param bs_env ...
!> \param U_Pp Atom-diagonal reference-to-optimized AA contraction.
! **************************************************************************************************
   SUBROUTINE assemble_AA_contraction(bs_env, U_Pp)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      COMPLEX(KIND=dp), ALLOCATABLE, INTENT(OUT)         :: U_Pp(:, :)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'assemble_AA_contraction'

      INTEGER                                            :: first, first_opt, first_ref, handle, &
                                                            iatom, last, nopt, nref

      CALL timeset(routineN, handle)

      IF (.NOT. bs_env%auto_ri%ready .OR. &
          .NOT. ALLOCATED(bs_env%auto_ri%sizes_ref_RI)) THEN
         CPABORT('Periodic AUTO_RI contraction is not initialized')
      END IF
      ALLOCATE (U_Pp(SUM(bs_env%auto_ri%sizes_ref_RI), &
                     SUM(bs_env%auto_ri%sizes_opt_RI)), source=z_zero)
      first_ref = 1
      first_opt = 1
      DO iatom = 1, bs_env%n_atom
         nref = bs_env%auto_ri%sizes_ref_RI(iatom)
         nopt = bs_env%auto_ri%sizes_opt_RI(iatom)
         first = bs_env%auto_ri%U_Pp_atom_offset(iatom)
         last = first + nref*nopt - 1
         U_Pp(first_ref:first_ref + nref - 1, first_opt:first_opt + nopt - 1) = &
            RESHAPE(bs_env%auto_ri%U_Pp_atom(first:last), [nref, nopt])
         first_ref = first_ref + nref
         first_opt = first_opt + nopt
      END DO
      IF (ANY(SUM(ABS(U_Pp)**2, DIM=1) <= TINY(1.0_dp))) THEN
         CPABORT('Periodic AUTO_RI contraction has an empty optimized direction')
      END IF

      CALL timestop(handle)

   END SUBROUTINE assemble_AA_contraction

! **************************************************************************************************
!> \brief Dense MO coefficients of the cyclically owned SCF k points (owner of k: MODULO(k-1, ranks)).
!> The locally stored parts of the distributed matrices are summed to the owner in batches of
!> k points: few large reductions instead of one global gather per k, and no rank stores the
!> coefficients of k points it does not own.
!> \param bs_env ...
!> \param mo_coeff distributed MO coefficients per k point and spin
!> \param coeff coeff(:, :, k_local, spin) for the owned k points
! **************************************************************************************************
   SUBROUTINE owned_mo_coefficients(bs_env, mo_coeff, coeff)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(cp_cfm_type), DIMENSION(:, :), INTENT(IN)     :: mo_coeff
      COMPLEX(KIND=dp), ALLOCATABLE, INTENT(OUT)         :: coeff(:, :, :, :)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'owned_mo_coefficients'
      REAL(KIND=dp), PARAMETER                           :: batch_bytes = 128.0E6_dp

      COMPLEX(KIND=dp), ALLOCATABLE                      :: batch(:, :, :)
      INTEGER                                            :: batch_size, first_owned, handle, &
                                                            i_batch, ikp, ispin, n_ao, n_batch, &
                                                            n_k, n_mo, n_owned, num_pe, root_rank

      CALL timeset(routineN, handle)

      n_k = SIZE(mo_coeff, 1)
      num_pe = bs_env%para_env%num_pe
      CALL cp_cfm_get_info(mo_coeff(1, 1), nrow_global=n_ao, ncol_global=n_mo)
      ALLOCATE (coeff(n_ao, n_mo, (n_k + num_pe - 1 - bs_env%para_env%mepos)/num_pe, SIZE(mo_coeff, 2)))
      batch_size = MAX(1, MIN((n_k + num_pe - 1)/num_pe, INT(batch_bytes/(complex_bytes*REAL(n_ao*n_mo, dp)))))
      ALLOCATE (batch(n_ao, n_mo, batch_size))
      DO ispin = 1, SIZE(mo_coeff, 2)
         DO root_rank = 0, MIN(num_pe, n_k) - 1
            n_owned = (n_k - root_rank - 1)/num_pe + 1
            DO first_owned = 1, n_owned, batch_size
               n_batch = MIN(batch_size, n_owned - first_owned + 1)
               DO i_batch = 1, n_batch
                  ikp = root_rank + 1 + (first_owned + i_batch - 2)*num_pe
                  CALL copy_local_cfm(mo_coeff(ikp, ispin), batch(:, :, i_batch))
               END DO
               CALL bs_env%para_env%sum(batch(:, :, 1:n_batch), root_rank)
               IF (root_rank == bs_env%para_env%mepos) THEN
                  coeff(:, :, first_owned:first_owned + n_batch - 1, ispin) = batch(:, :, 1:n_batch)
               END IF
            END DO
         END DO
      END DO
      DEALLOCATE (batch)

      CALL timestop(handle)

   END SUBROUTINE owned_mo_coefficients

! **************************************************************************************************
!> \brief Copy the locally owned elements of a distributed complex full matrix.
!> \param cfm ...
!> \param array ...
! **************************************************************************************************
   SUBROUTINE copy_local_cfm(cfm, array)
      TYPE(cp_cfm_type)                                  :: cfm
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(OUT)     :: array

      INTEGER                                            :: col, i_col, i_row, ncol_local, &
                                                            nrow_local, row
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices

      CALL cp_cfm_get_info(cfm, nrow_local=nrow_local, ncol_local=ncol_local, &
                           row_indices=row_indices, col_indices=col_indices)
      array = z_zero
      DO i_col = 1, ncol_local
         col = col_indices(i_col)
         DO i_row = 1, nrow_local
            row = row_indices(i_row)
            array(row, col) = cfm%local_data(i_row, i_col)
         END DO
      END DO
   END SUBROUTINE copy_local_cfm

! **************************************************************************************************
!> \brief Construct a complete cell grid dual to the uniform SCF k-point mesh.
!> \param bs_env ...
!> \param index_to_cell ...
! **************************************************************************************************
   SUBROUTINE build_full_cell_grid(bs_env, index_to_cell)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: index_to_cell(:, :)

      INTEGER                                            :: cell(3), i1, i2, i3, i_cell, &
                                                            i_direction, lower(3), upper(3)

      DO i_direction = 1, 3
         lower(i_direction) = -bs_env%kpoints_scf_desymm%nkp_grid(i_direction)/2
         upper(i_direction) = lower(i_direction) + &
                              bs_env%kpoints_scf_desymm%nkp_grid(i_direction) - 1
      END DO
      ALLOCATE (index_to_cell(3, PRODUCT(bs_env%kpoints_scf_desymm%nkp_grid)))
      i_cell = 0
      DO i1 = lower(1), upper(1)
         DO i2 = lower(2), upper(2)
            DO i3 = lower(3), upper(3)
               i_cell = i_cell + 1
               cell = [i1, i2, i3]
               index_to_cell(:, i_cell) = cell
            END DO
         END DO
      END DO

   END SUBROUTINE build_full_cell_grid

! **************************************************************************************************
!> \brief Zero-wave-vector Fourier value of the selected RI metric kernel.
!> \param bs_env ...
!> \return ...
! **************************************************************************************************
   FUNCTION ri_metric_head_zero(bs_env) RESULT(head)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp)                                      :: head

      IF (bs_env%ri_metric%potential_type == do_potential_short) THEN
         IF (bs_env%ri_metric%omega <= 0.0_dp) THEN
            CPABORT('RI attenuation must be positive')
         END IF
         head = pi/bs_env%ri_metric%omega**2
      ELSE
         IF (bs_env%ri_metric%potential_type /= do_potential_truncated) THEN
            CPABORT('Unsupported RI metric')
         END IF
         head = twopi*bs_env%ri_metric%cutoff_radius**2
      END IF
   END FUNCTION ri_metric_head_zero

! **************************************************************************************************
!> \brief Form distributed I_PQ(R) = (φ_P,0 | v | φ_Q,R) over every screened neighbor cell.
!>
!> The finite-range potential v and the Gaussian screening radii determine the cell set.
!> No cell is dropped because it lies outside the SCF Born-von Karman cell. Each real matrix
!> remains distributed over the existing RI matrix communicator; callers release fm_R.
!> The reciprocal interaction is I_PQ(q) = sum_R exp(i q.R) I_PQ(R).
!> \param fm_R ...
!> \param cells ...
!> \param potential ...
!> \param regularization ...
!> \param bs_env ...
!> \param qs_env ...
!> \param reference_basis ...
! **************************************************************************************************
   SUBROUTINE build_ri_2c_lattice(fm_R, cells, potential, regularization, bs_env, qs_env, reference_basis)
      TYPE(cp_fm_type), ALLOCATABLE, INTENT(OUT)         :: fm_R(:)
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: cells(:, :)
      TYPE(coulomb_operator_type), INTENT(IN)            :: potential
      REAL(KIND=dp), INTENT(IN)                          :: regularization
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN), OPTIONAL                      :: reference_basis

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_ri_2c_lattice'

      INTEGER                                            :: bounds(3), cell_index(3), handle, icell, &
                                                            ix, iy, iz, ncells
      INTEGER, ALLOCATABLE                               :: occupied(:, :, :), occupied_flat(:), &
                                                            sizes(:)
      INTEGER, POINTER                                   :: col_sizes(:), row_sizes(:)
      LOGICAL                                            :: use_reference_basis
      TYPE(cp_blacs_env_type), POINTER                   :: blacs
      TYPE(cp_fm_struct_type), POINTER                   :: reference_matrix_struct
      TYPE(dbcsr_distribution_type)                      :: distribution
      TYPE(dbcsr_type), ALLOCATABLE                      :: matrices(:)
      TYPE(distribution_2d_type), POINTER                :: dist_2d
      TYPE(kpoint_type), POINTER                         :: kpoints_ptr
      TYPE(kpoint_type), TARGET                          :: kpoints
      TYPE(mp_para_env_type), POINTER                    :: para
      TYPE(neighbor_list_iterator_p_type), POINTER       :: iterator(:)
      TYPE(neighbor_list_set_p_type), POINTER            :: neighbors(:)
      TYPE(particle_type), POINTER                       :: particles(:)
      TYPE(qs_kind_type), POINTER                        :: qs_kinds(:)

      CALL timeset(routineN, handle)

      use_reference_basis = .FALSE.
      IF (PRESENT(reference_basis)) use_reference_basis = reference_basis
      NULLIFY (blacs, para, reference_matrix_struct)
      CALL get_qs_env(qs_env, distribution_2d=dist_2d, particle_set=particles, qs_kind_set=qs_kinds, &
                      blacs_env=blacs, para_env=para)
      ALLOCATE (sizes(bs_env%n_atom))
      CALL get_particle_set(particles, qs_kinds, nsgf=sizes, basis=bs_env%basis_set_RI)
      NULLIFY (neighbors, iterator)
      CALL build_2c_neighbor_lists(neighbors, bs_env%basis_set_RI, bs_env%basis_set_RI, &
                                   potential, "full_2c_RI", qs_env, sym_ij=.FALSE., dist_2d=dist_2d)
      bounds = 0
      CALL neighbor_list_iterator_create(iterator, neighbors)
      DO WHILE (neighbor_list_iterate(iterator) == 0)
         CALL get_iterator_info(iterator, cell=cell_index)
         bounds = MAX(bounds, ABS(cell_index))
      END DO
      CALL neighbor_list_iterator_release(iterator)
      CALL bs_env%para_env%max(bounds)
      ALLOCATE (occupied(-bounds(1):bounds(1), -bounds(2):bounds(2), -bounds(3):bounds(3)))
      occupied = 0
      CALL neighbor_list_iterator_create(iterator, neighbors)
      DO WHILE (neighbor_list_iterate(iterator) == 0)
         CALL get_iterator_info(iterator, cell=cell_index)
         occupied(cell_index(1), cell_index(2), cell_index(3)) = 1
      END DO
      CALL neighbor_list_iterator_release(iterator)
      ALLOCATE (occupied_flat(SIZE(occupied)))
      occupied_flat(:) = RESHAPE(occupied, [SIZE(occupied)])
      CALL bs_env%para_env%max(occupied_flat)
      occupied(:, :, :) = RESHAPE(occupied_flat, SHAPE(occupied))
      DEALLOCATE (occupied_flat)
      ncells = SUM(occupied)
      CPASSERT(ncells > 0)
      ALLOCATE (cells(3, ncells), fm_R(ncells))
      ALLOCATE (kpoints%cell_to_index(-bounds(1):bounds(1), -bounds(2):bounds(2), -bounds(3):bounds(3)))
      kpoints%cell_to_index = 0
      icell = 0
      DO ix = -bounds(1), bounds(1)
         DO iy = -bounds(2), bounds(2)
            DO iz = -bounds(3), bounds(3)
               IF (occupied(ix, iy, iz) == 0) CYCLE
               icell = icell + 1
               cells(:, icell) = [ix, iy, iz]
               kpoints%cell_to_index(ix, iy, iz) = icell
            END DO
         END DO
      END DO
      DEALLOCATE (occupied)
      kpoints_ptr => kpoints
      ALLOCATE (matrices(ncells))
      CALL cp_dbcsr_dist2d_to_dist(dist_2d, distribution)
      ALLOCATE (row_sizes(SIZE(sizes)), col_sizes(SIZE(sizes)))
      row_sizes = sizes
      col_sizes = sizes
      CALL dbcsr_create(matrices(1), "full RI two-center", distribution, dbcsr_type_no_symmetry, row_sizes, col_sizes)
      DEALLOCATE (row_sizes, col_sizes)
      DO icell = 2, ncells
         CALL dbcsr_create(matrices(icell), template=matrices(1))
      END DO
      CALL build_2c_integrals(matrices, 0.0_dp, qs_env, neighbors, bs_env%basis_set_RI, &
                              bs_env%basis_set_RI, potential, do_kpoints=.TRUE., &
                              ext_kpoints=kpoints_ptr, regularization_RI=regularization)
      DO icell = 1, ncells
         IF (use_reference_basis) THEN
            IF (.NOT. bs_env%auto_ri%enabled) THEN
               CPABORT('Reference RI lattice requires AUTO_RI')
            END IF
            IF (.NOT. ASSOCIATED(reference_matrix_struct)) THEN
               CALL cp_fm_struct_create(reference_matrix_struct, context=blacs, nrow_global=SUM(sizes), &
                                        ncol_global=SUM(sizes), para_env=para)
            END IF
            CALL cp_fm_create(fm_R(icell), reference_matrix_struct)
         ELSE
            CALL cp_fm_create(fm_R(icell), bs_env%fm_RI_RI%matrix_struct)
         END IF
         CALL copy_dbcsr_to_fm(matrices(icell), fm_R(icell))
         CALL dbcsr_release(matrices(icell))
      END DO
      CALL dbcsr_distribution_release(distribution)
      CALL release_neighbor_list_sets(neighbors)
      IF (ASSOCIATED(reference_matrix_struct)) CALL cp_fm_struct_release(reference_matrix_struct)
      DEALLOCATE (matrices, sizes, kpoints%cell_to_index)

      CALL timestop(handle)

   END SUBROUTINE build_ri_2c_lattice

! **************************************************************************************************
!> \brief Evaluate I_PQ(q) = sum_R exp(i q.R) I_PQ(R) collectively into an existing dense matrix.
!>
!> q is in reciprocal fractional coordinates, so q.R includes 2π. The real-space matrices
!> stay distributed; no global RI-by-RI-by-cell array is replicated on the ranks.
!> \param fm_R ...
!> \param cells ...
!> \param qpoint ...
!> \param matrix_q ...
! **************************************************************************************************
   SUBROUTINE ri_2c_lattice_to_kp(fm_R, cells, qpoint, matrix_q)
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_R(:)
      INTEGER, INTENT(IN)                                :: cells(:, :)
      REAL(KIND=dp), INTENT(IN)                          :: qpoint(3)
      COMPLEX(KIND=dp), INTENT(OUT)                      :: matrix_q(:, :)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'ri_2c_lattice_to_kp'

      COMPLEX(KIND=dp)                                   :: phase
      INTEGER                                            :: handle, icell, icol, irow, ncols, nrows
      INTEGER, POINTER                                   :: cols(:), rows(:)

      CALL timeset(routineN, handle)

      CPASSERT(SIZE(fm_R) > 0)
      CPASSERT(SIZE(cells, 1) == 3)
      CPASSERT(SIZE(cells, 2) == SIZE(fm_R))
      CALL cp_fm_get_info(fm_R(1), nrow_local=nrows, ncol_local=ncols, row_indices=rows, col_indices=cols)
      matrix_q = z_zero
      DO icell = 1, SIZE(fm_R)
         phase = EXP(CMPLX(0.0_dp, twopi*DOT_PRODUCT(qpoint, REAL(cells(:, icell), dp)), dp))
         DO icol = 1, ncols
            DO irow = 1, nrows
               matrix_q(rows(irow), cols(icol)) = matrix_q(rows(irow), cols(icol)) + &
                                                  phase*fm_R(icell)%local_data(irow, icol)
            END DO
         END DO
      END DO
      CALL fm_R(1)%matrix_struct%para_env%sum(matrix_q)

      CALL timestop(handle)

   END SUBROUTINE ri_2c_lattice_to_kp

! **************************************************************************************************
!> \brief Integrate each auxiliary function for the Gamma charge-mode projection.
!> For a primitive x^(2m) exp(-α x²), I_0=sqrt(π/α), I_2m=(2m-1)I_(2m-2)/(2α).
!> Odd Cartesian moments vanish; the existing spherical contraction supplies normalization.
!> \param qs_env Atomic auxiliary bases.
!> \param charge Integrals of the reference auxiliary functions.
! **************************************************************************************************
   SUBROUTINE gamma_auxiliary_charges(qs_env, charge)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      COMPLEX(KIND=dp), INTENT(OUT)                      :: charge(:)

      INTEGER                                            :: first, iatom, ico, il, ipgf, iset, isgf, &
                                                            lmax, ncart, offset
      REAL(KIND=dp)                                      :: alpha
      REAL(KIND=dp), ALLOCATABLE                         :: moments(:), primitive(:)
      TYPE(gto_basis_set_type), POINTER                  :: basis
      TYPE(particle_type), POINTER                       :: particles(:)
      TYPE(qs_kind_type), POINTER                        :: qs_kinds(:)

      CALL get_qs_env(qs_env, particle_set=particles, qs_kind_set=qs_kinds)
      charge = z_zero
      first = 0
      DO iatom = 1, SIZE(particles)
         CALL get_qs_kind(qs_kinds(particles(iatom)%atomic_kind%kind_number), &
                          basis_set=basis, basis_type='RI_AUX')
         DO iset = 1, basis%nset
            lmax = basis%lmax(iset)
            ncart = ncoset(lmax)
            ALLOCATE (moments(0:lmax), primitive(ncart*basis%npgf(iset)))
            primitive = 0.0_dp
            DO ipgf = 1, basis%npgf(iset)
               alpha = basis%zet(ipgf, iset)
               moments = 0.0_dp
               moments(0) = SQRT(pi/alpha)
               DO il = 2, lmax, 2
                  moments(il) = REAL(il - 1, dp)*moments(il - 2)/(2.0_dp*alpha)
               END DO
               DO ico = ncoset(basis%lmin(iset) - 1) + 1, ncart
                  primitive((ipgf - 1)*ncart + ico) = PRODUCT(moments(indco(:, ico)))
               END DO
            END DO
            offset = basis%first_sgf(1, iset)
            DO isgf = offset, offset + basis%nsgf_set(iset) - 1
               charge(first + isgf) = CMPLX(SUM(basis%sphi(1:SIZE(primitive), isgf)*primitive), 0.0_dp, dp)
            END DO
            DEALLOCATE (moments, primitive)
         END DO
         first = first + basis%nsgf
      END DO
      CPASSERT(first == SIZE(charge))

   END SUBROUTINE gamma_auxiliary_charges

! **************************************************************************************************
!> \brief Return the contiguous part of 1..n owned by one rank of a block distribution.
!> \param n Number of distributed items (grid rows or ordered grid pairs).
!> \param rank ...
!> \param n_rank ...
!> \param first First owned item.
!> \param n_local Number of owned items.
! **************************************************************************************************
   PURE SUBROUTINE block_range(n, rank, n_rank, first, n_local)
      INTEGER, INTENT(IN)                                :: n, rank, n_rank
      INTEGER, INTENT(OUT)                               :: first, n_local

      first = INT(INT(n, int_8)*rank/n_rank) + 1
      n_local = INT(INT(n, int_8)*(rank + 1)/n_rank) - first + 1

   END SUBROUTINE block_range

END MODULE gw_ri_rs_small_cell_full_kp_lattice
