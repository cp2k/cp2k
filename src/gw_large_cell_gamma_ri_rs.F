!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief GW using RI-RS Approximation
!> \par History
!>      04.2026 created [Ritaj Tyagi]
! **************************************************************************************************
MODULE gw_large_cell_Gamma_ri_rs
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_add, dbcsr_binary_read, dbcsr_binary_write, dbcsr_copy, dbcsr_create, &
        dbcsr_distribution_get, dbcsr_distribution_new, dbcsr_distribution_release, &
        dbcsr_distribution_type, dbcsr_finalize, dbcsr_get_block_p, dbcsr_get_info, &
        dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, dbcsr_iterator_start, &
        dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_multiply, dbcsr_p_type, dbcsr_put_block, &
        dbcsr_release, dbcsr_scale, dbcsr_set, dbcsr_type, dbcsr_type_no_symmetry
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
                                              dbcsr_deallocate_matrix_set
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_release,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_should_output
   USE gw_integrals,                    ONLY: build_3c_integral_block
   USE gw_large_cell_gamma,             ONLY: G_occ_vir,&
                                              compute_QP_energies,&
                                              delete_unnecessary_files,&
                                              fill_fm_Sigma_c_Gamma_time,&
                                              get_W_MIC,&
                                              multiply_fm_W_MIC_time_with_Minv_Gamma
   USE gw_utils,                        ONLY: de_init_bs_env
   USE input_section_types,             ONLY: section_vals_type
   USE kinds,                           ONLY: dp
   USE machine,                         ONLY: m_walltime
   USE message_passing,                 ONLY: mp_para_env_type
   USE mp2_ri_2c,                       ONLY: RI_2c_integral_mat
   USE orbital_pointers,                ONLY: indco,&
                                              ncoset
   USE particle_types,                  ONLY: particle_type
   USE post_scf_bandstructure_types,    ONLY: post_scf_bandstructure_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gw_large_cell_Gamma_ri_rs'

   PUBLIC :: gw_calc_large_cell_Gamma_ri_rs

CONTAINS

! **************************************************************************************************
!> \brief GW calculation using RI-RS formalism
!> \param qs_env ...
!> \param bs_env ...
! **************************************************************************************************

   SUBROUTINE gw_calc_large_cell_Gamma_ri_rs(qs_env, bs_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gw_calc_large_cell_Gamma_ri_rs'

      INTEGER                                            :: handle
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_Sigma_x_Gamma, fm_W_MIC_time
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :, :)  :: fm_Sigma_c_Gamma_time

      CALL timeset(routineN, handle)

       !!========================================================================
       !! 1. Grid Generation for RI-RS
       !!    (Modified Lebedev grids from Ivan Duchemin and Xavier Blase)
       !!    Generate flattened 1D array of grid points for RI-RS.
       !!    Equation: r_g(k) = R_A + r_g(A)
       !!========================================================================
      CALL ri_rs_grid_assembler(qs_env, bs_env, bs_env%ri_rs%grid_points)

       !!========================================================================
       !! 2. Atomic Basis Evaluation
       !!    Compute values of spherical atomic basis functions at grid points.
       !!    Expression: Φ_μl = Φ_μ(r_l) (mat_phi_mu_l)
       !!========================================================================
      CALL atomic_basis_at_grid_point(qs_env, bs_env, bs_env%ri_rs%grid_points, &
                                      bs_env%ri_rs%mat_phi_mu_l)

       !!========================================================================
       !! 3. Compute RI-RS Coefficients (Z_lp)
       !!    Solve the regularized system for each atom P, where the grid domain
       !!    is restricted to r_l within a cutoff distance of atom P:
       !!    a. D_ll' = [ Σ_μ Φ_μ(r_l) Φ_μ(r_l') ]^2 (Equation 13)
       !!    b. D_lP  = Σ_{μν} Φ_μ(r_l) Φ_ν(r_l) (μν|P) (Equation 15)
       !!    c. Conditioning:
       !!       Dvec_l   = 1 / sqrt(D_ll)  (Diagonal scaling vector)
       !!       D'_ll' = Dvec_l * D_ll' * Dvec_l' + λδ_ll'
       !!       D'_lP  = Dvec_l * D_lP
       !!    d. Solve: Σ_l' D'_ll' * Z'_l'P = D'_lP (Equation 14)
       !!    e. Rescale: Z_lP = Z'_lP * Dvec_l   (Z_lP stored in mat_Z_lP)
       !!========================================================================
      CALL compute_coeff_Z_lP(qs_env, bs_env, bs_env%ri_rs%grid_points, &
                              bs_env%ri_rs%mat_phi_mu_l, bs_env%ri_rs%mat_Z_lP)

       !!========================================================================
       !! 4. Compute Independent-Particle Polarizability (χ)
       !!    G^occ_µλ(i|τ|,k=0) = sum_n^occ C_µn(k=0) e^(-|(ϵ_nk=0-ϵ_F)τ|) C_λn(k=0)
       !!    G^vir_µλ(i|τ|,k=0) = sum_n^vir C_µn(k=0) e^(-|(ϵ_nk=0-ϵ_F)τ|) C_λn(k=0)
       !!    G^occ_ll'(i|τ|,k=0) = sum_µν Φ_µ(r_l) G^occ_µν Φ_ν(r_l')
       !!    G^vir_ll'(i|τ|,k=0) = sum_µν Φ_µ(r_l) G^vir_µν Φ_ν(r_l')
       !!    χ_ll'(iτ,k=0)       = G^occ_ll'(i|τ|,k=0) * G^vir_ll'(i|τ|,k=0)
       !!    χ_PQ(iτ,k=0)        = sum_ll' Z_lP χ_ll'(iτ,k=0) Z_l'Q
       !!========================================================================
      CALL get_mat_chi_Gamma_tau(bs_env, bs_env%mat_chi_Gamma_tau, &
                                 bs_env%ri_rs%mat_phi_mu_l, bs_env%ri_rs%mat_Z_lP)

       !!========================================================================
       !! 5. Compute Screened Interaction (W^MIC)
       !!    χ_PQ(iτ,k=0) -> χ_PQ(iω,k) -> ε_PQ(iω,k) -> W_PQ(iω,k) -> W^MIC_PQ(iτ)
       !!========================================================================
      CALL get_W_MIC(bs_env, qs_env, bs_env%mat_chi_Gamma_tau, fm_W_MIC_time)

       !!========================================================================
       !! 6. Compute Exact Exchange Self-Energy (Σ^x)
       !!    D_µν          = sum_n^occ C_µn(k=0) C_νn(k=0)
       !!    D_ll'         = sum_µν Φ_µ(r_l) D_µν Φ_ν(r_l')
       !!    V^trunc_ll'   = sum_PQ Z_lP V^trunc_PQ Z_l'Q
       !!    Σ^x_ll'       = D_ll' * V^trunc_ll'
       !!    Σ^x_λσ(k=0)   = -sum_ll' Φ_λ(r_l) Σ^x_ll' Φ_σ(r_l')
       !!========================================================================
      CALL compute_Sigma_x(bs_env, qs_env, bs_env%ri_rs%mat_phi_mu_l, &
                           bs_env%ri_rs%mat_Z_lP, fm_Sigma_x_Gamma)

       !!========================================================================
       !! 7. Compute Correlation Self-Energy (Σ^c)
       !!    W^MIC_ll'(iτ,k=0) = sum_PQ Z_lP W^MIC_PQ(iτ) Z_l'Q
       !!    Σ^c_ll'(iτ,k=0)   = -G^occ_ll'(i|τ|,k=0) * W^MIC_ll'(iτ,k=0), for τ < 0
       !!    Σ^c_ll'(iτ,k=0)   =  G^vir_ll'(i|τ|,k=0) * W^MIC_ll'(iτ,k=0), for τ > 0
       !!    Σ^c_λσ(iτ,k=0)    = sum_ll' Φ_λ(r_l) Σ^c_ll'(iτ,k=0) Φ_σ(r_l')
       !!========================================================================
      CALL compute_Sigma_c(bs_env, fm_W_MIC_time, bs_env%ri_rs%mat_phi_mu_l, &
                           bs_env%ri_rs%mat_Z_lP, fm_Sigma_c_Gamma_time)

       !!========================================================================
       !! 8. Compute Quasiparticle Energies
       !!    Σ^c_λσ(iτ,k=0) -> Σ^c_nn(ϵ,k)
       !!    ϵ_nk^GW = ϵ_nk^DFT + Σ^c_nn(ϵ,k) + Σ^x_nn(k) - v^xc_nn(k)
       !!========================================================================
      CALL compute_QP_energies(bs_env, qs_env, fm_Sigma_x_Gamma, fm_Sigma_c_Gamma_time)

      CALL de_init_bs_env(bs_env)

      CALL timestop(handle)

   END SUBROUTINE gw_calc_large_cell_Gamma_ri_rs

! **************************************************************************************************
!> \brief Compute grid points for RI-RS
!>        Right now based on Ivan and Xavier implementation
!>        JCP 150, 174120 (2019), JCTC 17, 2383 (2021)
!> \param qs_env ...
!> \param bs_env ...
!> \param ri_rs_grid_points ...
! **************************************************************************************************

   SUBROUTINE ri_rs_grid_assembler(qs_env, bs_env, ri_rs_grid_points)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), ALLOCATABLE, INTENT(OUT)            :: ri_rs_grid_points(:, :)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'ri_rs_grid_assembler'

      INTEGER                                            :: atom_idx, end_idx, handle, ikind, j, &
                                                            natom, nkind, start_idx, &
                                                            total_grid_npts
      INTEGER, ALLOCATABLE                               :: atom_to_kind(:), ri_rs_grid_offsets(:)
      REAL(KIND=dp)                                      :: atomic_center(3)
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      CALL timeset(routineN, handle)

       !! Get the information about the atoms in the system
      CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, particle_set=particle_set)

      nkind = SIZE(atomic_kind_set)
      natom = SIZE(particle_set)

       !! 1. Build the grid cache
      CALL build_grid_cache(bs_env, atomic_kind_set)

       !! 2. Calculate grid counts and offsets
      ALLOCATE (ri_rs_grid_offsets(natom + 1))
      ALLOCATE (atom_to_kind(natom))

      total_grid_npts = 0
      DO ikind = 1, nkind
         DO j = 1, SIZE(atomic_kind_set(ikind)%atom_list)
            atom_idx = atomic_kind_set(ikind)%atom_list(j)
            atom_to_kind(atom_idx) = ikind
            ri_rs_grid_offsets(atom_idx) = total_grid_npts + 1
            total_grid_npts = total_grid_npts + bs_env%ri_rs%grid_cache(ikind)%npts
         END DO
      END DO

      ri_rs_grid_offsets(natom + 1) = total_grid_npts + 1

      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, FMT="(T2,A,T69,I12)") &
            'Total grid points used for RI-RS:', total_grid_npts
         WRITE (bs_env%unit_nr, "(A)") ' '
      END IF

       !! 3. Allocate the global ri_rs_grid arrays
      ALLOCATE (ri_rs_grid_points(3, total_grid_npts))

       !! 4. Parallelize the grid generation loop
      !$OMP PARALLEL DO DEFAULT(NONE) &
      !$OMP SHARED(ri_rs_grid_points, ri_rs_grid_offsets, atom_to_kind, &
      !$OMP        particle_set, bs_env, natom) &
      !$OMP PRIVATE(atom_idx, ikind, atomic_center, start_idx, end_idx) &
      !$OMP SCHEDULE(DYNAMIC, 1)
      DO atom_idx = 1, natom
         ikind = atom_to_kind(atom_idx)
         atomic_center(:) = particle_set(atom_idx)%r(:)

         start_idx = ri_rs_grid_offsets(atom_idx)
         end_idx = ri_rs_grid_offsets(atom_idx + 1) - 1

          !! Shift the cached origin grid by the atom's center
         ri_rs_grid_points(1, start_idx:end_idx) = bs_env%ri_rs%grid_cache(ikind)%raw_points(1, :) + atomic_center(1)
         ri_rs_grid_points(2, start_idx:end_idx) = bs_env%ri_rs%grid_cache(ikind)%raw_points(2, :) + atomic_center(2)
         ri_rs_grid_points(3, start_idx:end_idx) = bs_env%ri_rs%grid_cache(ikind)%raw_points(3, :) + atomic_center(3)

      END DO
      !$OMP END PARALLEL DO

       !! 5. Cleanup memory
      IF (ALLOCATED(bs_env%ri_rs%grid_cache)) THEN
         DO ikind = 1, nkind
            IF (ALLOCATED(bs_env%ri_rs%grid_cache(ikind)%raw_points)) DEALLOCATE (bs_env%ri_rs%grid_cache(ikind)%raw_points)
         END DO
         DEALLOCATE (bs_env%ri_rs%grid_cache)
      END IF

      DEALLOCATE (atom_to_kind, ri_rs_grid_offsets)
      CALL timestop(handle)

   END SUBROUTINE ri_rs_grid_assembler

! **************************************************************************************************
!> \brief Reads grids from .ion files and stores them in memory based on grid_select
!> \param bs_env ...
!> \param atomic_kind_set ...
! **************************************************************************************************

   SUBROUTINE build_grid_cache(bs_env, atomic_kind_set)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set

      CHARACTER(LEN=256)                                 :: atom_sym, filename, full_path, line, &
                                                            suffix
      INTEGER                                            :: colon_idx, i, ierr, ikind, iunit, nkind
      REAL(KIND=dp)                                      :: pt(3)

      ! Determine the file suffix based on the user's choice
      IF (bs_env%ri_rs%grid_select == 1) THEN
         suffix = "_def2-tzvp-rs.ion"
      ELSE IF (bs_env%ri_rs%grid_select == 2) THEN
         suffix = "_cc-pvtz-rs.ion"
      ELSE
         CPABORT("Unknown grid_select value. Valid options are 1 (def2-TZVPP) or 2 (cc-pVTZ).")
      END IF

      nkind = SIZE(atomic_kind_set)
      IF (.NOT. ALLOCATED(bs_env%ri_rs%grid_cache)) ALLOCATE (bs_env%ri_rs%grid_cache(nkind))

      DO ikind = 1, nkind
         atom_sym = TRIM(atomic_kind_set(ikind)%element_symbol)
         filename = TRIM(atom_sym)//TRIM(suffix)

         full_path = "ri_rs_grid/"//TRIM(filename)

         CALL open_file(file_name=TRIM(full_path), unit_number=iunit, &
                        file_action="READ", file_status="OLD")

         ! Parse the preamble for 'n points'
         bs_env%ri_rs%grid_cache(ikind)%npts = 0
         DO
            READ (iunit, '(A)', IOSTAT=ierr) line
            IF (ierr /= 0) EXIT
            IF (INDEX(line, 'n points') > 0) THEN
               colon_idx = INDEX(line, ':')
               READ (line(colon_idx + 1:), *) bs_env%ri_rs%grid_cache(ikind)%npts
               EXIT
            END IF
         END DO

         ! Allocate the cache array for this specific element
         ALLOCATE (bs_env%ri_rs%grid_cache(ikind)%raw_points(3, bs_env%ri_rs%grid_cache(ikind)%npts))

         ! Fast-forward to the <grid_points> tag
         REWIND (iunit)
         DO
            READ (iunit, '(A)', IOSTAT=ierr) line
            IF (ierr /= 0) EXIT
            IF (INDEX(line, '<grid_points>') > 0) EXIT
         END DO

         ! Read the raw grid coordinates
         DO i = 1, bs_env%ri_rs%grid_cache(ikind)%npts
            READ (iunit, *, IOSTAT=ierr) pt(1), pt(2), pt(3)
            IF (ierr /= 0) THEN
               CPABORT("Unexpected EOF in grid file ")
            END IF
            bs_env%ri_rs%grid_cache(ikind)%raw_points(:, i) = pt(:)
         END DO

         CALL close_file(unit_number=iunit)
      END DO

   END SUBROUTINE build_grid_cache

! **************************************************************************************************
!> \brief Evaluates atomic basis functions on a real-space grid and builds a sparse DBCSR matrix.
!> \param qs_env ...
!> \param bs_env ...
!> \param ri_rs_grid_points ...
!> \param mat_phi_mu_l ...
! **************************************************************************************************

   SUBROUTINE atomic_basis_at_grid_point(qs_env, bs_env, ri_rs_grid_points, mat_phi_mu_l)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), ALLOCATABLE, INTENT(INOUT)          :: ri_rs_grid_points(:, :)
      TYPE(dbcsr_type), INTENT(OUT)                      :: mat_phi_mu_l

      CHARACTER(LEN=*), PARAMETER :: routineN = 'atomic_basis_at_grid_point'

      INTEGER :: c_size, chunk_size, dimen_ORB, handle, i, i_blk, iatom, natom, npcol, nprow, &
         num_grid_chunks, r_end, r_start, total_grid_npts
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: first_sgf
      INTEGER, DIMENSION(:), POINTER                     :: c_blk_sizes, col_dist, col_dist_ks, &
                                                            r_blk_sizes, row_dist, row_dist_ks
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: atom_col_buffer
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_distribution_type)                      :: dist
      TYPE(dbcsr_distribution_type), POINTER             :: dbcsr_dist_ks
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      ! Setup Grid Blocking
      ! Right now using 256
      chunk_size = bs_env%ri_rs%chunk_size_dbcsr

      ! Extract environment variables
      CALL get_qs_env(qs_env, cell=cell, atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, particle_set=particle_set, &
                      para_env=para_env)

      natom = SIZE(particle_set)
      total_grid_npts = SIZE(ri_rs_grid_points, 2)

      ! Map the starting indices of spherical gaussian functions (SGF) for each atom
      ALLOCATE (first_sgf(natom + 1))
      CALL get_basis_offsets(particle_set, qs_kind_set, first_sgf, dimen_ORB)

      ! =========================================================================
      ! 1. SETUP DBCSR MATRIX TOPOLOGY
      ! =========================================================================

      ! A. Define Column Block Sizes (1 Block = 1 Atom's full basis set)
      ALLOCATE (c_blk_sizes(natom))
      DO iatom = 1, natom
         c_blk_sizes(iatom) = first_sgf(iatom + 1) - first_sgf(iatom)
      END DO

      ! B. Define Row Block Sizes (Grid chunks of max size 256)
      num_grid_chunks = CEILING(REAL(total_grid_npts, KIND=dp)/REAL(chunk_size, KIND=dp))
      ALLOCATE (r_blk_sizes(num_grid_chunks))
      r_blk_sizes = chunk_size
      IF (MOD(total_grid_npts, chunk_size) /= 0) THEN
         r_blk_sizes(num_grid_chunks) = MOD(total_grid_npts, chunk_size)
      END IF

      ! C. Fetch CP2K's Default Process Grid Configuration
      CALL get_qs_env(qs_env, dbcsr_dist=dbcsr_dist_ks)
      CALL dbcsr_distribution_get(dbcsr_dist_ks, row_dist=row_dist_ks, col_dist=col_dist_ks)

      ! D. Build Custom Mappings using Round-Robin across the 2D process grid
      !    (MAXVAL + 1 accounts for 0-based process indexing in DBCSR)
      nprow = MAXVAL(row_dist_ks) + 1
      npcol = MAXVAL(col_dist_ks) + 1

      ALLOCATE (row_dist(num_grid_chunks))
      DO i = 1, num_grid_chunks
         row_dist(i) = MOD(i - 1, nprow)
      END DO

      ALLOCATE (col_dist(natom))
      DO i = 1, natom
         col_dist(i) = MOD(i - 1, npcol)
      END DO

      ! E. Create the DBCSR Distribution and Initialize the Matrix
      CALL dbcsr_distribution_new(dist, template=dbcsr_dist_ks, &
                                  row_dist=row_dist, col_dist=col_dist)

      CALL dbcsr_create(mat_phi_mu_l, name="phi_val_sparse", dist=dist, &
                        matrix_type=dbcsr_type_no_symmetry, &
                        row_blk_size=r_blk_sizes, col_blk_size=c_blk_sizes)

      ! =========================================================================
      ! 2. STREAM DATA DIRECTLY INTO SPARSE MATRIX
      ! =========================================================================
      ! Iterate over the atoms assigned to this specific MPI rank
      DO iatom = para_env%mepos + 1, natom, para_env%num_pe

         c_size = c_blk_sizes(iatom)

         ! Allocate a temporary dense buffer just for this specific atom
         ALLOCATE (atom_col_buffer(total_grid_npts, c_size))
         atom_col_buffer = 0.0_dp

         ! Evaluate the basis functions on the grid
         CALL fill_phi_for_atom(atom_col_buffer, ri_rs_grid_points, total_grid_npts, &
                                iatom, particle_set, qs_kind_set, cell)

         ! Slice the dense column into chunks and insert into DBCSR
         DO i_blk = 1, num_grid_chunks
            r_start = (i_blk - 1)*chunk_size + 1
            r_end = MIN(i_blk*chunk_size, total_grid_npts)

            ! Apply dynamic sparsity filtering: Only store blocks with physical significance
            IF (MAXVAL(ABS(atom_col_buffer(r_start:r_end, 1:c_size))) > bs_env%eps_filter) THEN
               CALL dbcsr_put_block(mat_phi_mu_l, row=i_blk, col=iatom, &
                                    block=atom_col_buffer(r_start:r_end, 1:c_size))
            END IF
         END DO

         DEALLOCATE (atom_col_buffer)

      END DO

      ! Finalize triggers internal MPI communication to route blocks to their correct 2D process owners
      CALL dbcsr_finalize(mat_phi_mu_l)

      ! -------------------------------------------------------------------------
      ! CLEANUP
      ! -------------------------------------------------------------------------
      DEALLOCATE (first_sgf, r_blk_sizes, c_blk_sizes, row_dist, col_dist)
      CALL dbcsr_distribution_release(dist)

      CALL timestop(handle)

   END SUBROUTINE atomic_basis_at_grid_point

! **************************************************************************************************
!> \brief Compute value of all basis function for single atom across all grid points
!> \param phi_val ...
!> \param ri_rs_grid ...
!> \param npts ...
!> \param iatom ...
!> \param particle_set ...
!> \param qs_kind_set ...
!> \param cell ...
! **************************************************************************************************

   SUBROUTINE fill_phi_for_atom(phi_val, ri_rs_grid, npts, iatom, &
                                particle_set, qs_kind_set, cell)

      REAL(KIND=dp), INTENT(INOUT)                       :: phi_val(:, :)
      INTEGER, INTENT(IN)                                :: npts
      REAL(KIND=dp), INTENT(IN)                          :: ri_rs_grid(3, npts)
      INTEGER, INTENT(IN)                                :: iatom
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(cell_type), POINTER                           :: cell

      INTEGER                                            :: first_sgf, i_pt, ico, iend_co, ikind, &
                                                            ipgf, iset, isgf, ishell, istart_co, &
                                                            l, last_sgf, lx, ly, lz, n_cart_total, &
                                                            row_idx
      REAL(KIND=dp)                                      :: alpha, dist_vec(3), exp_val, poly, r2, &
                                                            r_atom(3), weight
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_set

      ! Get Atom Info
      ikind = particle_set(iatom)%atomic_kind%kind_number
      CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis_set, basis_type="ORB")
      IF (.NOT. ASSOCIATED(orb_basis_set)) RETURN

      r_atom = particle_set(iatom)%r

      !$OMP PARALLEL DO DEFAULT(NONE) &
      !$OMP SHARED(phi_val, ri_rs_grid, npts, orb_basis_set, r_atom, cell, ncoset, indco) &
      !$OMP PRIVATE(i_pt, dist_vec, r2, iset, n_cart_total, ishell, l, istart_co, iend_co, &
      !$OMP         first_sgf, last_sgf, ipgf, alpha, exp_val, isgf, ico, row_idx, weight, &
      !$OMP         lx, ly, lz, poly) &
      !$OMP SCHEDULE(STATIC)
      DO i_pt = 1, npts
         dist_vec = pbc(ri_rs_grid(:, i_pt) - r_atom, cell)
         r2 = DOT_PRODUCT(dist_vec, dist_vec)

         DO iset = 1, orb_basis_set%nset
            n_cart_total = ncoset(orb_basis_set%lmax(iset))

            DO ishell = 1, orb_basis_set%nshell(iset)
               l = orb_basis_set%l(ishell, iset)
               istart_co = ncoset(l - 1) + 1
               iend_co = ncoset(l)

               first_sgf = orb_basis_set%first_sgf(ishell, iset)
               last_sgf = orb_basis_set%last_sgf(ishell, iset)

               DO ipgf = 1, orb_basis_set%npgf(iset)
                  alpha = orb_basis_set%zet(ipgf, iset)
                  exp_val = EXP(-alpha*r2)

                  DO isgf = first_sgf, last_sgf
                     DO ico = istart_co, iend_co
                        row_idx = (ipgf - 1)*n_cart_total + ico
                        weight = orb_basis_set%sphi(row_idx, isgf)
                        lx = indco(1, ico)
                        ly = indco(2, ico)
                        lz = indco(3, ico)
                        poly = (dist_vec(1)**lx)*(dist_vec(2)**ly)*(dist_vec(3)**lz)

                        phi_val(i_pt, isgf) = phi_val(i_pt, isgf) + (weight*poly*exp_val)

                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
      !$OMP END PARALLEL DO

   END SUBROUTINE fill_phi_for_atom

! **************************************************************************************************
!> \brief Helper for OMP threads to fill phi_val column values
!> \param particle_set ...
!> \param qs_kind_set ...
!> \param first_sgf ...
!> \param total_sgf ...
! **************************************************************************************************

   SUBROUTINE get_basis_offsets(particle_set, qs_kind_set, first_sgf, total_sgf)

      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      INTEGER, INTENT(OUT)                               :: first_sgf(:), total_sgf

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'get_basis_offsets'

      INTEGER                                            :: handle, iatom, ikind, nsgf

      CALL timeset(routineN, handle)

      total_sgf = 0
      DO iatom = 1, SIZE(particle_set)
         first_sgf(iatom) = total_sgf + 1
         ikind = particle_set(iatom)%atomic_kind%kind_number
         CALL get_qs_kind(qs_kind_set(ikind), nsgf=nsgf, basis_type="ORB")
         total_sgf = total_sgf + nsgf
      END DO
      first_sgf(SIZE(particle_set) + 1) = total_sgf + 1

      CALL timestop(handle)

   END SUBROUTINE get_basis_offsets

! **************************************************************************************************
!> \brief Compute RI-RS Coefficients (Z_lP)
!> \param qs_env ...
!> \param bs_env ...
!> \param ri_rs_grid_points ...
!> \param mat_phi_mu_l ...
!> \param mat_Z_lP ...
! **************************************************************************************************

   SUBROUTINE compute_coeff_Z_lP(qs_env, bs_env, ri_rs_grid_points, mat_phi_mu_l, mat_Z_lP)

      ! Arguments
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), ALLOCATABLE, INTENT(INOUT)          :: ri_rs_grid_points(:, :)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: mat_phi_mu_l
      TYPE(dbcsr_type), INTENT(OUT)                      :: mat_Z_lP

      CHARACTER(LEN=*), PARAMETER :: key = 'PROPERTIES%BANDSTRUCTURE%GW%PRINT%RESTART', &
         routineN = 'compute_coeff_Z_lP'

      INTEGER :: atom_P, chunk_size, col, current_chunk_size, global_i, global_j, group_handle, &
         handle, i, i_blk, info, istart_ao, j, j_ri, l, loc_idx, max_ao_size, max_loc_ri, &
         n_ao_total, n_grid_total, n_loc_ri, n_local_grid, natom, num_grid_chunks, r_end, r_start, &
         row
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: local_grid_idx, map_global_to_local, &
                                                            row_offset
      INTEGER, DIMENSION(:), POINTER                     :: col_dist_phi, col_dist_ri, r_blk_sizes, &
                                                            ri_blk_sizes, row_dist_grid
      REAL(KIND=dp)                                      :: cutoff_ri, dist, t1
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: d_vec_local
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: D_local, d_lp_local, int_2d_local, &
                                                            phi_global, phi_local, rho_pair_local, &
                                                            Z_blk, Z_global_temp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: int_3c_local
      REAL(KIND=dp), DIMENSION(3)                        :: pos_P
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: block_ptr
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_distribution_type)                      :: dist_phi, dist_Z
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(section_vals_type), POINTER                   :: input

      CALL timeset(routineN, handle)

      ! Restrict integration domain to a localized sphere around the atom
      cutoff_ri = bs_env%ri_rs%cutoff_radius_ri_rs

      t1 = m_walltime()

      CALL get_qs_env(qs_env, para_env=para_env, particle_set=particle_set, input=input)

      natom = SIZE(bs_env%i_RI_start_from_atom)
      n_ao_total = bs_env%i_ao_end_from_atom(natom)
      n_grid_total = SIZE(ri_rs_grid_points, 2)

      ! =========================================================================
      ! 1. SETUP DBCSR TOPOLOGY & EXACT OFFSETS
      ! =========================================================================
      CALL dbcsr_get_info(mat_phi_mu_l, row_blk_size=r_blk_sizes, distribution=dist_phi)
      CALL dbcsr_distribution_get(dist_phi, row_dist=row_dist_grid, col_dist=col_dist_phi, group=group_handle)

      chunk_size = r_blk_sizes(1)
      num_grid_chunks = SIZE(r_blk_sizes)

      ALLOCATE (row_offset(num_grid_chunks))
      row_offset(1) = 0
      DO i_blk = 2, num_grid_chunks
         row_offset(i_blk) = row_offset(i_blk - 1) + r_blk_sizes(i_blk - 1)
      END DO

      ALLOCATE (ri_blk_sizes(natom), col_dist_ri(natom))
      DO atom_P = 1, natom
         ri_blk_sizes(atom_P) = bs_env%i_RI_end_from_atom(atom_P) - bs_env%i_RI_start_from_atom(atom_P) + 1
         col_dist_ri(atom_P) = MOD(atom_P - 1, MAXVAL(col_dist_phi) + 1)
      END DO

      CALL dbcsr_distribution_new(dist_Z, template=dist_phi, row_dist=row_dist_grid, col_dist=col_dist_ri)

      IF (bs_env%ri_rs%Z_lP_exists) THEN
         CALL dbcsr_binary_read(filepath=TRIM(bs_env%prefix)//"Z_lP.matrix", &
                                distribution=dist_Z, &
                                matrix_new=mat_Z_lP)
         IF (bs_env%unit_nr > 0) THEN
            WRITE (bs_env%unit_nr, '(T2,A,T57,A,F7.1,A)') &
               'Read Z_lP from file ', ' Execution time', m_walltime() - t1, ' s'
            WRITE (bs_env%unit_nr, '(A)') ' '
         END IF
      ELSE

         CALL dbcsr_create(mat_Z_lP, name="mat_Z_lP", dist=dist_Z, &
                           matrix_type=dbcsr_type_no_symmetry, &
                           row_blk_size=r_blk_sizes, col_blk_size=ri_blk_sizes)

         max_ao_size = 0
         DO j = 1, SIZE(bs_env%i_ao_start_from_atom)
            max_ao_size = MAX(max_ao_size, bs_env%i_ao_end_from_atom(j) - bs_env%i_ao_start_from_atom(j) + 1)
         END DO
         max_loc_ri = MAXVAL(ri_blk_sizes)

         ! Global mapping array: 4 bytes * N_grid_total (~4MB per million grid points)
         ALLOCATE (map_global_to_local(n_grid_total))

         ! =========================================================================
         ! 2. PRE-COMPUTE: EXTRACT FULL DBCSR TO DENSE GLOBAL ARRAY
         ! =========================================================================
         ALLOCATE (phi_global(n_grid_total, n_ao_total))
         phi_global = 0.0_dp

         CALL dbcsr_iterator_start(iter, mat_phi_mu_l)
         DO WHILE (dbcsr_iterator_blocks_left(iter))
            CALL dbcsr_iterator_next_block(iter, row, col, block_ptr)

            IF (col < 1 .OR. col > natom) CYCLE
            istart_ao = bs_env%i_ao_start_from_atom(col)

            DO j = 1, SIZE(block_ptr, 2)
               global_j = istart_ao + j - 1
               IF (global_j < 1 .OR. global_j > n_ao_total) CYCLE

               DO i = 1, SIZE(block_ptr, 1)
                  global_i = row_offset(row) + i
                  IF (global_i < 1 .OR. global_i > n_grid_total) CYCLE

                  ! Insert directly into the global dense array
                  phi_global(global_i, global_j) = block_ptr(i, j)
               END DO
            END DO
         END DO
         CALL dbcsr_iterator_stop(iter)

         CALL para_env%sum(phi_global)

         ! =========================================================================
         ! 3. MPI LOOP OVER ATOMS (Fully independent, no MPI barriers inside)
         ! This distributes the N_atoms among the available MPI processors.
         ! =========================================================================
         DO atom_P = para_env%mepos + 1, natom, para_env%num_pe

            n_loc_ri = ri_blk_sizes(atom_P)
            pos_P(:) = particle_set(atom_P)%r(:)

            ! ---------------------------------------------------------------------
            ! A. Determine Local Grid Domain based on cutoff_ri
            ! ---------------------------------------------------------------------
            n_local_grid = 0
            DO l = 1, n_grid_total
               dist = SQRT(SUM((ri_rs_grid_points(1:3, l) - pos_P(1:3))**2))
               IF (dist <= cutoff_ri) n_local_grid = n_local_grid + 1
            END DO

            ALLOCATE (local_grid_idx(n_local_grid))
            map_global_to_local = 0

            n_local_grid = 0
            DO l = 1, n_grid_total
               dist = SQRT(SUM((ri_rs_grid_points(1:3, l) - pos_P(1:3))**2))
               IF (dist <= cutoff_ri) THEN
                  n_local_grid = n_local_grid + 1
                  local_grid_idx(n_local_grid) = l
                  map_global_to_local(l) = n_local_grid
               END IF
            END DO

            ! ---------------------------------------------------------------------
            ! B. Smart Fetch of Local Phi (MPI-Safe & Bounds-Checked)
            ! ---------------------------------------------------------------------
            ALLOCATE (phi_local(n_local_grid, n_ao_total))

            !$OMP PARALLEL DO DEFAULT(NONE) &
            !$OMP SHARED(n_ao_total, n_local_grid, local_grid_idx, phi_local, phi_global) &
            !$OMP PRIVATE(global_j, loc_idx, global_i) &
            !$OMP SCHEDULE(STATIC)
            DO global_j = 1, n_ao_total
               DO loc_idx = 1, n_local_grid
                  global_i = local_grid_idx(loc_idx)
                  phi_local(loc_idx, global_j) = phi_global(global_i, global_j)
               END DO
            END DO
            !$OMP END PARALLEL DO

            ! ---------------------------------------------------------------------
            ! C. Build and Balance Local LHS Matrix (D_local)
            ! ---------------------------------------------------------------------
            ALLOCATE (D_local(n_local_grid, n_local_grid), d_vec_local(n_local_grid))
            D_local = 0.0_dp

            CALL dsyrk("L", "N", n_local_grid, n_ao_total, 1.0_dp, phi_local, n_local_grid, 0.0_dp, D_local, n_local_grid)

            !$OMP PARALLEL DO DEFAULT(NONE) &
            !$OMP SHARED(n_local_grid, D_local, d_vec_local, bs_env) &
            !$OMP PRIVATE(i) &
            !$OMP SCHEDULE(STATIC)
            DO i = 1, n_local_grid
               D_local(i, i) = D_local(i, i)**2
               d_vec_local(i) = 1.0_dp/SQRT(MAX(D_local(i, i), 1.0E-16_dp))
               D_local(i, i) = (D_local(i, i)*d_vec_local(i)**2) + bs_env%ri_rs%tikhonov
            END DO
            !$OMP END PARALLEL DO

            !$OMP PARALLEL DO DEFAULT(NONE) &
            !$OMP SHARED(n_local_grid, D_local, d_vec_local) &
            !$OMP PRIVATE(j, i) &
            !$OMP SCHEDULE(DYNAMIC)
            DO j = 1, n_local_grid
               DO i = j + 1, n_local_grid
                  D_local(i, j) = D_local(i, j)**2
                  D_local(i, j) = D_local(i, j)*d_vec_local(i)*d_vec_local(j)
                  D_local(j, i) = D_local(i, j)
               END DO
            END DO
            !$OMP END PARALLEL DO

            ! ---------------------------------------------------------------------
            ! D. Build Local RHS Matrix (d_lp_local)
            ! ---------------------------------------------------------------------
            ALLOCATE (d_lp_local(n_local_grid, n_loc_ri))
            ALLOCATE (rho_pair_local(n_local_grid, max_ao_size*max_ao_size))
            ALLOCATE (int_3c_local(max_ao_size, max_ao_size, n_loc_ri))
            ALLOCATE (int_2d_local(max_ao_size*max_ao_size, n_loc_ri))

            d_lp_local = 0.0_dp

            CALL compute_d_lp(qs_env, bs_env, phi_local, d_lp_local, n_local_grid, &
                              n_loc_ri, atom_P, rho_pair_local, int_3c_local, int_2d_local, max_ao_size)

            !$OMP PARALLEL DO DEFAULT(NONE) &
            !$OMP SHARED(n_loc_ri, n_local_grid, d_lp_local, d_vec_local) &
            !$OMP PRIVATE(j_ri, i) &
            !$OMP SCHEDULE(STATIC)
            DO j_ri = 1, n_loc_ri
               DO i = 1, n_local_grid
                  d_lp_local(i, j_ri) = d_lp_local(i, j_ri)*d_vec_local(i)
               END DO
            END DO
            !$OMP END PARALLEL DO

            ! ---------------------------------------------------------------------
            ! E. Fast Local LAPACK Solve
            ! ---------------------------------------------------------------------

            CALL dpotrf('L', n_local_grid, D_local, n_local_grid, info)

            CALL dpotrs('L', n_local_grid, n_loc_ri, D_local, n_local_grid, d_lp_local, n_local_grid, info)

            !$OMP PARALLEL DO DEFAULT(NONE) &
            !$OMP SHARED(n_loc_ri, n_local_grid, d_lp_local, d_vec_local) &
            !$OMP PRIVATE(j_ri, i) &
            !$OMP SCHEDULE(STATIC)
            DO j_ri = 1, n_loc_ri
               DO i = 1, n_local_grid
                  d_lp_local(i, j_ri) = d_lp_local(i, j_ri)*d_vec_local(i)
               END DO
            END DO
            !$OMP END PARALLEL DO

            ! ---------------------------------------------------------------------
            ! F. Scatter Local Solution Back to Global DBCSR Matrix
            ! ---------------------------------------------------------------------
            ALLOCATE (Z_global_temp(n_grid_total, n_loc_ri))
            Z_global_temp = 0.0_dp

            !$OMP PARALLEL DO DEFAULT(NONE) &
            !$OMP SHARED(n_local_grid, local_grid_idx, Z_global_temp, n_loc_ri, d_lp_local) &
            !$OMP PRIVATE(loc_idx, global_i) &
            !$OMP SCHEDULE(STATIC)
            DO loc_idx = 1, n_local_grid
               global_i = local_grid_idx(loc_idx)
               Z_global_temp(global_i, 1:n_loc_ri) = d_lp_local(loc_idx, 1:n_loc_ri)
            END DO
            !$OMP END PARALLEL DO

            ALLOCATE (Z_blk(MAXVAL(r_blk_sizes), n_loc_ri))

            DO i_blk = 1, num_grid_chunks
               r_start = row_offset(i_blk) + 1
               r_end = row_offset(i_blk) + r_blk_sizes(i_blk)
               current_chunk_size = r_blk_sizes(i_blk)

               Z_blk = 0.0_dp
               Z_blk(1:current_chunk_size, 1:n_loc_ri) = Z_global_temp(r_start:r_end, 1:n_loc_ri)

               IF (MAXVAL(ABS(Z_blk(1:current_chunk_size, 1:n_loc_ri))) > bs_env%eps_filter) THEN
                  CALL dbcsr_put_block(mat_Z_lP, row=i_blk, col=atom_P, &
                                       block=Z_blk(1:current_chunk_size, 1:n_loc_ri))
               END IF
            END DO

            DEALLOCATE (Z_global_temp, Z_blk)
            DEALLOCATE (D_local, d_vec_local, d_lp_local)
            DEALLOCATE (rho_pair_local, int_3c_local, int_2d_local)

            DEALLOCATE (local_grid_idx, phi_local)

         END DO

         ! Clean up the massive global array once the loop finishes
         DEALLOCATE (phi_global)

         CALL dbcsr_finalize(mat_Z_lP)

         IF (bs_env%unit_nr > 0) THEN
            WRITE (bs_env%unit_nr, '(T2,A,T57,A,F7.1,A)') &
               'Computed Z_lP ', ' Execution time', m_walltime() - t1, ' s'
            WRITE (bs_env%unit_nr, '(A)') ' '
         END IF

         logger => cp_get_default_logger()

         IF (BTEST(cp_print_key_should_output(logger%iter_info, input, key), cp_p_file)) THEN
            CALL dbcsr_binary_write(matrix=mat_Z_lP, filepath=TRIM(bs_env%prefix)//"Z_lP.matrix")
         END IF

         DEALLOCATE (map_global_to_local)

      END IF

      DEALLOCATE (row_offset, ri_blk_sizes, col_dist_ri)
      CALL dbcsr_distribution_release(dist_Z)

      DEALLOCATE (ri_rs_grid_points)

      CALL timestop(handle)

   END SUBROUTINE compute_coeff_Z_lP

! **************************************************************************************************
!> \brief Computes the dense localized Right-Hand Side (RHS) matrix d_lp
!> \param qs_env ...
!> \param bs_env ...
!> \param phi_val ...
!> \param d_lp ...
!> \param n_grid_total ...
!> \param n_loc_ri ...
!> \param atom_P ...
!> \param rho_pair ...
!> \param int_3c ...
!> \param int_2d ...
!> \param max_ao_size ...
! **************************************************************************************************

   SUBROUTINE compute_d_lp(qs_env, bs_env, phi_val, d_lp, n_grid_total, n_loc_ri, atom_P, &
                           rho_pair, int_3c, int_2d, max_ao_size)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: phi_val
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: d_lp
      INTEGER, INTENT(IN)                                :: n_grid_total, n_loc_ri, atom_P
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: rho_pair
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: int_3c
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: int_2d
      INTEGER, INTENT(IN)                                :: max_ao_size

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_d_lp'

      INTEGER                                            :: atom_j, atom_k, handle, j, jk_idx, &
                                                            jsize, jstart, k, ksize, kstart, l, ri

      CALL timeset(routineN, handle)

      DO atom_j = 1, SIZE(bs_env%i_ao_start_from_atom)
         jstart = bs_env%i_ao_start_from_atom(atom_j)
         jsize = bs_env%i_ao_end_from_atom(atom_j) - jstart + 1

         DO atom_k = 1, SIZE(bs_env%i_ao_start_from_atom)
            kstart = bs_env%i_ao_start_from_atom(atom_k)
            ksize = bs_env%i_ao_end_from_atom(atom_k) - kstart + 1

            int_3c(1:jsize, 1:ksize, 1:n_loc_ri) = 0.0_dp

            ! Compute B_{μν,P} = (μν|P)
            CALL build_3c_integral_block(int_3c(1:jsize, 1:ksize, 1:n_loc_ri), &
                                         qs_env, bs_env%ri_metric, &
                                         basis_j=bs_env%basis_set_AO, &
                                         basis_k=bs_env%basis_set_AO, &
                                         basis_i=bs_env%basis_set_RI, &
                                         atom_k=atom_k, atom_j=atom_j, atom_i=atom_P)

            ! Build Pair Density: ρ(l, μν) = Φ_μ(r_l) \times Φ_ν(r_l)
            !$OMP PARALLEL DO DEFAULT(NONE) &
            !$OMP SHARED(ksize, jsize, n_grid_total, rho_pair, phi_val, jstart, kstart) &
            !$OMP PRIVATE(k, j, l, jk_idx) &
            !$OMP COLLAPSE(2)
            DO k = 1, ksize
               DO j = 1, jsize
                  jk_idx = (k - 1)*jsize + j
                  DO l = 1, n_grid_total
                     rho_pair(l, jk_idx) = phi_val(l, jstart + j - 1)*phi_val(l, kstart + k - 1)
                  END DO
               END DO
            END DO
            !$OMP END PARALLEL DO

            ! Flatten 3D B_{μν, P} tensor to 2D B_{(μν), P} matrix for BLAS
            !$OMP PARALLEL DO DEFAULT(NONE) &
            !$OMP SHARED(n_loc_ri, ksize, jsize, int_2d, int_3c) &
            !$OMP PRIVATE(ri, k, j, jk_idx) &
            !$OMP COLLAPSE(2)
            DO ri = 1, n_loc_ri
               DO k = 1, ksize
                  DO j = 1, jsize
                     jk_idx = (k - 1)*jsize + j
                     int_2d(jk_idx, ri) = int_3c(j, k, ri)
                  END DO
               END DO
            END DO
            !$OMP END PARALLEL DO

            ! Perform Contraction: d_{l, P} += ρ(l, μν) \times B_{(μν), P}
            CALL dgemm("N", "N", n_grid_total, n_loc_ri, jsize*ksize, &
                       1.0_dp, rho_pair, n_grid_total, &
                       int_2d, max_ao_size*max_ao_size, &
                       1.0_dp, d_lp(1, 1), n_grid_total)
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE compute_d_lp

! **************************************************************************************************
!> \brief Computes the χ(iτ, k=0) matrix
!> \param bs_env ...
!> \param mat_chi_Gamma_tau ...
!> \param mat_phi_mu_l ...
!> \param mat_Z_lP ...
! **************************************************************************************************

   SUBROUTINE get_mat_chi_Gamma_tau(bs_env, mat_chi_Gamma_tau, mat_phi_mu_l, mat_Z_lP)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_chi_Gamma_tau
      TYPE(dbcsr_type), INTENT(INOUT)                    :: mat_phi_mu_l, mat_Z_lP

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_mat_chi_Gamma_tau'

      INTEGER                                            :: handle, i, i_t, ispin, npcol
      INTEGER, DIMENSION(:), POINTER                     :: blk_ao, blk_grid, dist_col_ao, &
                                                            dist_col_grid, dist_row_grid
      REAL(KIND=dp)                                      :: t1, tau
      TYPE(dbcsr_distribution_type)                      :: dist_grid_grid, dist_phi
      TYPE(dbcsr_type)                                   :: matrix_chi_grid, matrix_chi_grid_spin, &
                                                            matrix_G_occ_grid, matrix_G_vir_grid

      CALL timeset(routineN, handle)

      ! =========================================================================
      ! 1. SETUP CORE TOPOLOGIES
      ! =========================================================================
      CALL dbcsr_get_info(mat_phi_mu_l, distribution=dist_phi, row_blk_size=blk_grid, col_blk_size=blk_ao)
      CALL dbcsr_distribution_get(dist_phi, row_dist=dist_row_grid, col_dist=dist_col_ao)

      ! Determine number of MPI process columns
      npcol = MAXVAL(dist_col_ao) + 1

      ! Build a perfectly safe column distribution for the Grid dimension
      ALLOCATE (dist_col_grid(SIZE(blk_grid)))
      DO i = 1, SIZE(blk_grid)
         dist_col_grid(i) = MOD(i - 1, npcol)
      END DO

      CALL dbcsr_distribution_new(dist_grid_grid, template=dist_phi, &
                                  row_dist=dist_row_grid, col_dist=dist_col_grid)

      CALL dbcsr_create(matrix_G_occ_grid, "G_occ_grid", dist_grid_grid, dbcsr_type_no_symmetry, blk_grid, blk_grid)
      CALL dbcsr_create(matrix_G_vir_grid, "G_vir_grid", dist_grid_grid, dbcsr_type_no_symmetry, blk_grid, blk_grid)
      CALL dbcsr_create(matrix_chi_grid, "chi_grid", dist_grid_grid, dbcsr_type_no_symmetry, blk_grid, blk_grid)
      CALL dbcsr_create(matrix_chi_grid_spin, "chi_grid_spin", dist_grid_grid, dbcsr_type_no_symmetry, blk_grid, blk_grid)

      ! =========================================================================
      ! 2. MAIN IMAGINARY TIME LOOP
      ! =========================================================================
      DO i_t = 1, bs_env%num_time_freq_points
         t1 = m_walltime()

         tau = bs_env%imag_time_points(i_t)
         CALL dbcsr_set(matrix_chi_grid, 0.0_dp)

         ! ----------------------------------------------------------------------
         ! A. SPIN LOOP (Allocations safely encapsulated in wrappers)
         ! ----------------------------------------------------------------------
         DO ispin = 1, bs_env%n_spin

            ! G^occ_µλ(i|τ|,k=0) = sum_n^occ C_µn(k=0) e^(-|(ϵ_nk=0-ϵ_F)τ|) C_λn(k=0)
            ! G^occ_ll'(i|τ|,k=0) = sum_µν Φ_µ(r_l) G^occ_µν Φ_ν(r_l')
            CALL build_G_grid(bs_env, tau, ispin, .TRUE., .FALSE., mat_phi_mu_l, &
                              matrix_G_occ_grid, bs_env%eps_filter)

            ! G^vir_µλ(i|τ|,k=0) = sum_n^vir C_µn(k=0) e^(-|(ϵ_nk=0-ϵ_F)τ|) C_λn(k=0)
            ! G^vir_ll'(i|τ|,k=0) = sum_µν Φ_µ(r_l) G^vir_µν Φ_ν(r_l')
            CALL build_G_grid(bs_env, tau, ispin, .FALSE., .TRUE., mat_phi_mu_l, &
                              matrix_G_vir_grid, bs_env%eps_filter)

            ! -------------------------------------------------------------------
            ! B. ELEMENT-WISE HADAMARD PRODUCT
            ! -------------------------------------------------------------------
            ! χ_ll'(iτ,k=0) = G^occ_ll'(i|τ|,k=0) * G^vir_ll'(i|τ|,k=0)
            CALL hadamard_product(matrix_G_occ_grid, matrix_G_vir_grid, matrix_chi_grid_spin, bs_env%spin_degeneracy)

            ! Accumulate spin contributions
            CALL dbcsr_add(matrix_chi_grid, matrix_chi_grid_spin, 1.0_dp, 1.0_dp)

         END DO ! ispin

         ! ----------------------------------------------------------------------
         ! C. TRANSFORM TO AUXILIARY BASIS & EXPORT DIRECTLY
         ! χ_aux        = Z^T * χ_grid * Z
         ! χ_PQ(iτ,k=0) = sum_ll' Z_lP χ_ll'(iτ,k=0) Z_l'Q
         ! Result is dumped directly into the final array mat_chi_Gamma_tau!
         ! ----------------------------------------------------------------------
         CALL contract_A_B_A("T", "N", mat_Z_lP, matrix_chi_grid, &
                             mat_chi_Gamma_tau(i_t)%matrix, bs_env%eps_filter)

         IF (bs_env%unit_nr > 0) THEN
            WRITE (bs_env%unit_nr, '(T2,A,I13,A,I3,A,F7.1,A)') &
               'Computed χ(iτ,k=0) for time point', i_t, ' /', bs_env%num_time_freq_points, &
               ',    Execution time', m_walltime() - t1, ' s'
         END IF

      END DO ! i_t

      ! =========================================================================
      ! 3. FINAL CLEANUP
      ! =========================================================================
      CALL dbcsr_release(matrix_G_occ_grid)
      CALL dbcsr_release(matrix_G_vir_grid)
      CALL dbcsr_release(matrix_chi_grid)
      CALL dbcsr_release(matrix_chi_grid_spin)
      CALL dbcsr_distribution_release(dist_grid_grid)
      DEALLOCATE (dist_col_grid)

      IF (bs_env%unit_nr > 0) WRITE (bs_env%unit_nr, '(A)') ' '

      CALL timestop(handle)

   END SUBROUTINE get_mat_chi_Gamma_tau

! **************************************************************************************************
!> \brief Computes Green's Function in grid basis
!> \param bs_env ...
!> \param tau ...
!> \param ispin ...
!> \param occ ...
!> \param vir ...
!> \param mat_phi_mu_l ...
!> \param matrix_G_grid ...
!> \param eps_filter ...
! **************************************************************************************************

   SUBROUTINE build_G_grid(bs_env, tau, ispin, occ, vir, mat_phi_mu_l, matrix_G_grid, eps_filter)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), INTENT(IN)                          :: tau
      INTEGER, INTENT(IN)                                :: ispin
      LOGICAL, INTENT(IN)                                :: occ, vir
      TYPE(dbcsr_type), INTENT(INOUT)                    :: mat_phi_mu_l, matrix_G_grid
      REAL(KIND=dp), INTENT(IN)                          :: eps_filter

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'build_G_grid'

      INTEGER                                            :: handle
      INTEGER, DIMENSION(:), POINTER                     :: blk_ao, dist_row_ao
      TYPE(cp_fm_type), POINTER                          :: fm_G
      TYPE(dbcsr_distribution_type)                      :: dist_ao_ao
      TYPE(dbcsr_type)                                   :: matrix_G_ao

      CALL timeset(routineN, handle)

      ! 1. Select the correct FM matrix based on occ/vir flags
      IF (occ) THEN
         fm_G => bs_env%fm_Gocc
      ELSE
         fm_G => bs_env%fm_Gvir
      END IF

      ! 2. Compute Dense FM Green's Function
      ! G^occ/vir_µλ(i|τ|,k=0) = sum_G^occ/vir_µλn^occ/vir C_µn(k=0) e^(-|(ϵ_nk=0-ϵ_F)τ|) C_λn(k=0)
      CALL G_occ_vir(bs_env, tau, fm_G, ispin, occ=occ, vir=vir)

      ! 3. Setup AO DBCSR Topology and Create Matrix dynamically
      CALL setup_square_topology(mat_phi_mu_l, 'COL', dist_ao_ao, blk_ao, dist_row_ao)

      CALL dbcsr_create(matrix_G_ao, name="G_ao", dist=dist_ao_ao, &
                        matrix_type=dbcsr_type_no_symmetry, &
                        row_blk_size=blk_ao, col_blk_size=blk_ao)

      ! 4. Convert FM to Sparse DBCSR
      CALL copy_fm_to_dbcsr(fm_G, matrix_G_ao, keep_sparsity=.FALSE.)

      ! 5. Transform to Grid Basis: G_grid = phi * G_ao * phi^T
      ! G^occ/vir_ll'(i|τ|,k=0) = sum_µν Φ_µ(r_l) G^occ/vir_µν Φ_ν(r_l')
      CALL contract_A_B_A("N", "T", mat_phi_mu_l, matrix_G_ao, matrix_G_grid, eps_filter)

      ! 6. Release AO matrix and topology
      CALL release_dbcsr_topology_and_matrices(dist=dist_ao_ao, mapped_dist=dist_row_ao, m1=matrix_G_ao)

      CALL timestop(handle)

   END SUBROUTINE build_G_grid

! **************************************************************************************************
!> \brief Generalized routine to compute OUT = A * B * A^T  OR  OUT = A^T * B * A using DBCSR
!> \param transA_left ...
!> \param transA_right ...
!> \param matrix_A ...
!> \param matrix_B ...
!> \param matrix_out ...
!> \param eps_filter ...
! **************************************************************************************************

   SUBROUTINE contract_A_B_A(transA_left, transA_right, matrix_A, matrix_B, matrix_out, eps_filter)

      CHARACTER(LEN=1), INTENT(IN)                       :: transA_left, transA_right
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_A, matrix_B, matrix_out
      REAL(KIND=dp), INTENT(IN)                          :: eps_filter

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'contract_A_B_A'

      INTEGER                                            :: handle
      TYPE(dbcsr_type)                                   :: matrix_tmp

      CALL timeset(routineN, handle)

      CALL dbcsr_create(matrix_tmp, template=matrix_A)

      IF (transA_left == "N" .AND. transA_right == "T") THEN
         ! Path 1: Out = A * B * A^T
         CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_A, matrix_B, &
                             0.0_dp, matrix_tmp, filter_eps=eps_filter)
         CALL dbcsr_multiply("N", "T", 1.0_dp, matrix_tmp, matrix_A, &
                             0.0_dp, matrix_out, filter_eps=eps_filter)

      ELSE IF (transA_left == "T" .AND. transA_right == "N") THEN
         ! Path 2: Out = A^T * B * A
         CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_B, matrix_A, &
                             0.0_dp, matrix_tmp, filter_eps=eps_filter)
         CALL dbcsr_multiply("T", "N", 1.0_dp, matrix_A, matrix_tmp, &
                             0.0_dp, matrix_out, filter_eps=eps_filter)
      ELSE
         CPABORT("Unsupported transposition pair in contract_A_B_A")
      END IF

      CALL dbcsr_release(matrix_tmp)

      CALL timestop(handle)

   END SUBROUTINE contract_A_B_A

! **************************************************************************************************
!> \brief Computes C = A ◦ B (Element-wise Hadamard product) for sparse DBCSR matrices.
!> \param matrix_A ...
!> \param matrix_B ...
!> \param matrix_C ...
!> \param fac (Scaling factor applied to the product)
! **************************************************************************************************

   SUBROUTINE hadamard_product(matrix_A, matrix_B, matrix_C, fac)

      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_A, matrix_B, matrix_C
      REAL(KIND=dp), INTENT(IN)                          :: fac

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'hadamard_product'

      INTEGER                                            :: col, handle, row
      LOGICAL                                            :: found
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: blk_B, blk_C
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, handle)

      CALL dbcsr_copy(matrix_C, matrix_A)

      CALL dbcsr_iterator_start(iter, matrix_C)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, row, col, blk_C)

         CALL dbcsr_get_block_p(matrix_B, row, col, blk_B, found)

         IF (found) THEN
            blk_C(:, :) = fac*blk_C(:, :)*blk_B(:, :)
         ELSE
            ! If B is sparse here, the product is zero
            blk_C(:, :) = 0.0_dp
         END IF
      END DO
      CALL dbcsr_iterator_stop(iter)

      CALL timestop(handle)

   END SUBROUTINE hadamard_product

! **************************************************************************************************
!> \brief Computes the exact exchange part of the GW self-energy
!> \param bs_env ...
!> \param qs_env ...
!> \param mat_phi_mu_l ...
!> \param mat_Z_lP ...
!> \param fm_Sigma_x_Gamma ...
! **************************************************************************************************

   SUBROUTINE compute_Sigma_x(bs_env, qs_env, mat_phi_mu_l, mat_Z_lP, fm_Sigma_x_Gamma)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_type), INTENT(INOUT)                    :: mat_phi_mu_l, mat_Z_lP
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_Sigma_x_Gamma

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_Sigma_x'

      INTEGER                                            :: handle, ispin
      INTEGER, DIMENSION(:), POINTER                     :: blk_aux, blk_grid, dist_col_grid, &
                                                            dist_row_aux
      REAL(KIND=dp)                                      :: t1
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: fm_Vtr_Gamma
      TYPE(dbcsr_distribution_type)                      :: dist_aux_aux, dist_grid_grid
      TYPE(dbcsr_type)                                   :: mat_Sigma_x_Gamma, matrix_D_grid, &
                                                            matrix_Sigma_x_grid, matrix_V_aux, &
                                                            matrix_V_grid

      CALL timeset(routineN, handle)

      t1 = m_walltime()

      ALLOCATE (fm_Sigma_x_Gamma(bs_env%n_spin))
      DO ispin = 1, bs_env%n_spin
         CALL cp_fm_create(fm_Sigma_x_Gamma(ispin), bs_env%fm_s_Gamma%matrix_struct)
      END DO

      CALL dbcsr_create(mat_Sigma_x_Gamma, template=bs_env%mat_ao_ao%matrix)

      ! =========================================================================
      ! 1. SETUP CORE TOPOLOGIES
      ! =========================================================================
      CALL setup_square_topology(mat_phi_mu_l, 'ROW', dist_grid_grid, blk_grid, dist_col_grid)
      CALL setup_square_topology(mat_Z_lP, 'COL', dist_aux_aux, blk_aux, dist_row_aux)

      ! =========================================================================
      ! 2. COMPUTE V^tr_ll'
      ! =========================================================================
      CALL RI_2c_integral_mat(qs_env, fm_Vtr_Gamma, bs_env%fm_RI_RI, bs_env%n_RI, &
                              bs_env%trunc_coulomb, do_kpoints=.FALSE.)

      ! M^-1 * V^tr * M^-1 directly modifies fm_Vtr_Gamma(:, 1)
      CALL multiply_fm_W_MIC_time_with_Minv_Gamma(bs_env, qs_env, fm_Vtr_Gamma(:, 1))

      CALL dbcsr_create(matrix_V_aux, "V_aux", dist_aux_aux, dbcsr_type_no_symmetry, blk_aux, blk_aux)
      CALL dbcsr_create(matrix_V_grid, "V_grid", dist_grid_grid, dbcsr_type_no_symmetry, blk_grid, blk_grid)

      CALL copy_fm_to_dbcsr(fm_Vtr_Gamma(1, 1), matrix_V_aux, keep_sparsity=.FALSE.)

      ! V^tr_ll' = sum_PQ Z_lP V^trunc_PQ Z_l'Q
      CALL contract_A_B_A("N", "T", mat_Z_lP, matrix_V_aux, matrix_V_grid, bs_env%eps_filter)
      CALL dbcsr_release(matrix_V_aux)

      ! =========================================================================
      ! 3. SPIN LOOP FOR EXACT EXCHANGE
      ! =========================================================================
      DO ispin = 1, bs_env%n_spin

         ! Density matrix on grid is essentially G_occ at tau = 0.0
         ! D_µν  = sum_n^occ C_µn(k=0) C_νn(k=0)
         ! D_ll' = sum_µν Φ_µ(r_l) D_µν Φ_ν(r_l')
         CALL dbcsr_create(matrix_D_grid, "D_grid", dist_grid_grid, dbcsr_type_no_symmetry, blk_grid, blk_grid)
         CALL build_G_grid(bs_env, 0.0_dp, ispin, .TRUE., .FALSE., mat_phi_mu_l, matrix_D_grid, bs_env%eps_filter)

         ! Element-wise Hadamard product: Σ^x_grid = D_grid ◦ V_grid
         ! Σ^x_ll' = D_ll' * V^tr_ll'
         CALL dbcsr_create(matrix_Sigma_x_grid, template=matrix_V_grid)
         CALL hadamard_product(matrix_D_grid, matrix_V_grid, matrix_Sigma_x_grid, 1.0_dp)

         CALL dbcsr_release(matrix_D_grid)

         ! Transform back to AO basis: Σ^x_ao = -1.0 * phi^T * Σ^x_grid * phi
         ! Σ^x_λσ(k=0)   = -sum_ll' Φ_λ(r_l) Σ^x_ll' Φ_σ(r_l')
         CALL contract_A_B_A("T", "N", mat_phi_mu_l, matrix_Sigma_x_grid, mat_Sigma_x_Gamma, bs_env%eps_filter)
         CALL dbcsr_scale(mat_Sigma_x_Gamma, -1.0_dp)

         CALL dbcsr_release(matrix_Sigma_x_grid)

         ! Data I/O and Export to CP2K Full Matrices
         CALL copy_dbcsr_to_fm(mat_Sigma_x_Gamma, fm_Sigma_x_Gamma(ispin))

      END DO ! ispin

      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(T2,A,T58,A,F7.1,A)') &
            'Computed Σ^x(k=0),', ' Execution time', m_walltime() - t1, ' s'
         WRITE (bs_env%unit_nr, '(A)') ' '
      END IF

      ! =========================================================================
      ! 4. CLEANUP
      ! =========================================================================
      CALL release_dbcsr_topology_and_matrices(dist=dist_grid_grid, mapped_dist=dist_col_grid, &
                                               m1=mat_Sigma_x_Gamma, m2=matrix_V_grid)
      CALL release_dbcsr_topology_and_matrices(dist=dist_aux_aux, mapped_dist=dist_row_aux)

      CALL cp_fm_release(fm_Vtr_Gamma)

      CALL timestop(handle)

   END SUBROUTINE compute_Sigma_x

! **************************************************************************************************
!> \brief Computes the correlation part of the GW self-energy
!> \param bs_env ...
!> \param fm_W_MIC_time ...
!> \param mat_phi_mu_l ...
!> \param mat_Z_lP ...
!> \param fm_Sigma_c_Gamma_time ...
! **************************************************************************************************

   SUBROUTINE compute_Sigma_c(bs_env, fm_W_MIC_time, mat_phi_mu_l, mat_Z_lP, fm_Sigma_c_Gamma_time)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_W_MIC_time
      TYPE(dbcsr_type), INTENT(INOUT)                    :: mat_phi_mu_l, mat_Z_lP
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :, :)  :: fm_Sigma_c_Gamma_time

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_Sigma_c'

      INTEGER                                            :: handle, i_t, ispin
      INTEGER, DIMENSION(:), POINTER                     :: blk_aux, blk_grid, dist_col_grid, &
                                                            dist_row_aux
      REAL(KIND=dp)                                      :: t1, tau
      TYPE(dbcsr_distribution_type)                      :: dist_aux_aux, dist_grid_grid
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: mat_Sigma_neg_tau, mat_Sigma_pos_tau
      TYPE(dbcsr_type) :: matrix_G_occ_grid, matrix_G_vir_grid, matrix_Sigma_neg_grid, &
         matrix_Sigma_pos_grid, matrix_W_aux, matrix_W_grid

      CALL timeset(routineN, handle)

      ! =========================================================================
      ! 1. SETUP CORE TOPOLOGIES AND PRE-ALLOCATE OUTPUT ARRAYS
      ! =========================================================================
      CALL setup_square_topology(mat_phi_mu_l, 'ROW', dist_grid_grid, blk_grid, dist_col_grid)
      CALL setup_square_topology(mat_Z_lP, 'COL', dist_aux_aux, blk_aux, dist_row_aux)

      ! Pre-allocate local DBCSR matrices to act as targets for final output
      NULLIFY (mat_Sigma_neg_tau, mat_Sigma_pos_tau)
      ALLOCATE (mat_Sigma_neg_tau(bs_env%num_time_freq_points, bs_env%n_spin))
      ALLOCATE (mat_Sigma_pos_tau(bs_env%num_time_freq_points, bs_env%n_spin))

      DO i_t = 1, bs_env%num_time_freq_points
         DO ispin = 1, bs_env%n_spin
            ALLOCATE (mat_Sigma_neg_tau(i_t, ispin)%matrix)
            ALLOCATE (mat_Sigma_pos_tau(i_t, ispin)%matrix)
            CALL dbcsr_create(mat_Sigma_neg_tau(i_t, ispin)%matrix, template=bs_env%mat_ao_ao%matrix)
            CALL dbcsr_create(mat_Sigma_pos_tau(i_t, ispin)%matrix, template=bs_env%mat_ao_ao%matrix)
         END DO
      END DO

      ! =========================================================================
      ! 2. MAIN IMAGINARY TIME LOOP
      ! =========================================================================
      DO i_t = 1, bs_env%num_time_freq_points
         tau = bs_env%imag_time_points(i_t)

         ! -------------------------------------------------------------------
         ! Compute W_grid = Z * W_aux * Z^T
         ! -------------------------------------------------------------------
         CALL dbcsr_create(matrix_W_aux, "W_aux", dist_aux_aux, dbcsr_type_no_symmetry, blk_aux, blk_aux)
         CALL dbcsr_create(matrix_W_grid, "W_grid", dist_grid_grid, dbcsr_type_no_symmetry, blk_grid, blk_grid)

         CALL copy_fm_to_dbcsr(fm_W_MIC_time(i_t), matrix_W_aux, keep_sparsity=.FALSE.)

         ! W^MIC_ll'(iτ,k=0) = sum_PQ Z_lP W^MIC_PQ(iτ) Z_l'Q
         CALL contract_A_B_A("N", "T", mat_Z_lP, matrix_W_aux, matrix_W_grid, bs_env%eps_filter)

         CALL dbcsr_release(matrix_W_aux) ! Clean up aux basis immediately

         DO ispin = 1, bs_env%n_spin
            t1 = m_walltime()

            ! -------------------------------------------------------------------
            ! A. Transform Green's Functions to the Grid
            ! -------------------------------------------------------------------
            CALL dbcsr_create(matrix_G_occ_grid, "G_occ_grid", dist_grid_grid, dbcsr_type_no_symmetry, blk_grid, blk_grid)
            CALL dbcsr_create(matrix_G_vir_grid, "G_vir_grid", dist_grid_grid, dbcsr_type_no_symmetry, blk_grid, blk_grid)

            ! G^occ_µλ(i|τ|,k=0) = sum_G^occ_µλn^occ C_µn(k=0) e^(-|(ϵ_nk=0-ϵ_F)τ|) C_λn(k=0)
            ! G^occ_ll'(i|τ|,k=0) = sum_µν Φ_µ(r_l) G^occ_µν Φ_ν(r_l')
            CALL build_G_grid(bs_env, tau, ispin, .TRUE., .FALSE., mat_phi_mu_l, matrix_G_occ_grid, bs_env%eps_filter)

            ! G^vir_µλ(i|τ|,k=0) = sum_n^vir C_µn(k=0) e^(-|(ϵ_nk=0-ϵ_F)τ|) C_λn(k=0)
            ! G^vir_ll'(i|τ|,k=0) = sum_µν Φ_µ(r_l) G^vir_µν Φ_ν(r_l')
            CALL build_G_grid(bs_env, tau, ispin, .FALSE., .TRUE., mat_phi_mu_l, matrix_G_vir_grid, bs_env%eps_filter)

            ! -------------------------------------------------------------------
            ! B. Element-wise Hadamard Products for Sigma_c on Grid
            ! Σ_neg_grid = G_occ_grid ◦ W_grid
            ! Σ_pos_grid = G_vir_grid ◦ W_grid
            ! -------------------------------------------------------------------
            CALL dbcsr_create(matrix_Sigma_neg_grid, template=matrix_W_grid)
            CALL dbcsr_create(matrix_Sigma_pos_grid, template=matrix_W_grid)

            ! Σ^c_ll'(iτ,k=0)   = -G^occ_ll'(i|τ|,k=0) * W^MIC_ll'(iτ,k=0), for τ < 0
            CALL hadamard_product(matrix_G_occ_grid, matrix_W_grid, matrix_Sigma_neg_grid, 1.0_dp)

            ! Σ^c_ll'(iτ,k=0)   =  G^vir_ll'(i|τ|,k=0) * W^MIC_ll'(iτ,k=0), for τ > 0
            CALL hadamard_product(matrix_G_vir_grid, matrix_W_grid, matrix_Sigma_pos_grid, 1.0_dp)

            ! Instantly purge massive G_grid arrays to save memory
            CALL dbcsr_release(matrix_G_occ_grid)
            CALL dbcsr_release(matrix_G_vir_grid)

            ! -------------------------------------------------------------------
            ! C. Transform Sigma back to AO Basis
            ! Σ_AO = phi^T * Σ_grid * phi
            ! -------------------------------------------------------------------

            ! Σ^c_λσ(iτ,k=0)    = sum_ll' Φ_λ(r_l) Σ^c_ll'(iτ,k=0) Φ_σ(r_l'), for τ < 0
            CALL contract_A_B_A("T", "N", mat_phi_mu_l, matrix_Sigma_neg_grid, &
                                mat_Sigma_neg_tau(i_t, ispin)%matrix, bs_env%eps_filter)
            CALL dbcsr_scale(mat_Sigma_neg_tau(i_t, ispin)%matrix, -1.0_dp)

            ! Σ^c_λσ(iτ,k=0)    = sum_ll' Φ_λ(r_l) Σ^c_ll'(iτ,k=0) Φ_σ(r_l'), for τ > 0
            CALL contract_A_B_A("T", "N", mat_phi_mu_l, matrix_Sigma_pos_grid, &
                                mat_Sigma_pos_tau(i_t, ispin)%matrix, bs_env%eps_filter)

            ! Purge Grid Sigma arrays
            CALL dbcsr_release(matrix_Sigma_neg_grid)
            CALL dbcsr_release(matrix_Sigma_pos_grid)

            IF (bs_env%unit_nr > 0) THEN
               WRITE (bs_env%unit_nr, '(T2,A,I10,A,I3,A,F7.1,A)') &
                  'Computed Σ^c(iτ,k=0) for time point ', i_t, ' /', bs_env%num_time_freq_points, &
                  ',    Execution time', m_walltime() - t1, ' s'
            END IF

         END DO ! ispin

         ! Release the W_grid for this time point
         CALL dbcsr_release(matrix_W_grid)

      END DO ! i_t

      IF (bs_env%unit_nr > 0) WRITE (bs_env%unit_nr, '(A)') ' '

      ! -------------------------------------------------------------------------
      ! 3. FINALIZE AND CLEANUP
      ! -------------------------------------------------------------------------
      CALL fill_fm_Sigma_c_Gamma_time(fm_Sigma_c_Gamma_time, bs_env, &
                                      mat_Sigma_pos_tau, mat_Sigma_neg_tau)

      CALL cp_fm_release(fm_W_MIC_time)

      CALL dbcsr_deallocate_matrix_set(mat_Sigma_neg_tau)
      CALL dbcsr_deallocate_matrix_set(mat_Sigma_pos_tau)

      CALL release_dbcsr_topology_and_matrices(dist=dist_grid_grid, mapped_dist=dist_col_grid)
      CALL release_dbcsr_topology_and_matrices(dist=dist_aux_aux, mapped_dist=dist_row_aux)

      CALL delete_unnecessary_files(bs_env)
      CALL timestop(handle)

   END SUBROUTINE compute_Sigma_c

! **************************************************************************************************
!> \brief DBCSR Topology Generation
!> \param matrix_template ...
!> \param dim_type ...
!> \param square_dist ...
!> \param blk_sizes ...
!> \param mapped_dist ...
! **************************************************************************************************

   SUBROUTINE setup_square_topology(matrix_template, dim_type, square_dist, blk_sizes, mapped_dist)

      TYPE(dbcsr_type), INTENT(IN)                       :: matrix_template
      CHARACTER(LEN=*), INTENT(IN)                       :: dim_type
      TYPE(dbcsr_distribution_type), INTENT(OUT)         :: square_dist
      INTEGER, DIMENSION(:), INTENT(OUT), POINTER        :: blk_sizes, mapped_dist

      INTEGER                                            :: i, np
      INTEGER, DIMENSION(:), POINTER                     :: col_blk, col_dist, row_blk, row_dist
      TYPE(dbcsr_distribution_type)                      :: dist_template

      CALL dbcsr_get_info(matrix_template, distribution=dist_template, &
                          row_blk_size=row_blk, col_blk_size=col_blk)
      CALL dbcsr_distribution_get(dist_template, row_dist=row_dist, col_dist=col_dist)

      IF (TRIM(dim_type) == 'ROW') THEN
         ! Creates ROW x ROW (e.g., Grid x Grid from mat_phi_mu_l)
         blk_sizes => row_blk
         np = MAXVAL(col_dist) + 1  ! npcol
         ALLOCATE (mapped_dist(SIZE(blk_sizes)))
         DO i = 1, SIZE(blk_sizes)
            mapped_dist(i) = MOD(i - 1, np)
         END DO
         CALL dbcsr_distribution_new(square_dist, template=dist_template, &
                                     row_dist=row_dist, col_dist=mapped_dist)

      ELSE IF (TRIM(dim_type) == 'COL') THEN
         ! Creates COL x COL (e.g., Aux x Aux from mat_Z_lP)
         blk_sizes => col_blk
         np = MAXVAL(row_dist) + 1  ! nprow
         ALLOCATE (mapped_dist(SIZE(blk_sizes)))
         DO i = 1, SIZE(blk_sizes)
            mapped_dist(i) = MOD(i - 1, np)
         END DO
         CALL dbcsr_distribution_new(square_dist, template=dist_template, &
                                     row_dist=mapped_dist, col_dist=col_dist)
      END IF

   END SUBROUTINE setup_square_topology

! **************************************************************************************************
!> \brief DBCSR matrices deallocation
!> \param dist ...
!> \param mapped_dist    ...
!> \param m1 ...
!> \param m2 ...
!> \param m3 ...
!> \param m4 ...
! **************************************************************************************************

   SUBROUTINE release_dbcsr_topology_and_matrices(dist, mapped_dist, m1, m2, m3, m4)

      TYPE(dbcsr_distribution_type), INTENT(INOUT), &
         OPTIONAL                                        :: dist
      INTEGER, DIMENSION(:), INTENT(INOUT), OPTIONAL, &
         POINTER                                         :: mapped_dist
      TYPE(dbcsr_type), INTENT(INOUT), OPTIONAL          :: m1, m2, m3, m4

      IF (PRESENT(dist)) CALL dbcsr_distribution_release(dist)
      IF (PRESENT(mapped_dist)) THEN
         IF (ASSOCIATED(mapped_dist)) THEN
            DEALLOCATE (mapped_dist)
            NULLIFY (mapped_dist)
         END IF
      END IF
      IF (PRESENT(m1)) CALL dbcsr_release(m1)
      IF (PRESENT(m2)) CALL dbcsr_release(m2)
      IF (PRESENT(m3)) CALL dbcsr_release(m3)
      IF (PRESENT(m4)) CALL dbcsr_release(m4)

   END SUBROUTINE release_dbcsr_topology_and_matrices

END MODULE gw_large_cell_Gamma_ri_rs
