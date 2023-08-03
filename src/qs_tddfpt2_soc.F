!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
MODULE qs_tddfpt2_soc
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_cfm_diag,                     ONLY: cp_cfm_heevd
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_get_info,&
                                              cp_cfm_get_submatrix,&
                                              cp_cfm_release,&
                                              cp_cfm_type,&
                                              cp_fm_to_cfm
   USE cp_control_types,                ONLY: dft_control_type,&
                                              qs_control_type,&
                                              tddfpt2_control_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
                                              cp_dbcsr_sm_fm_multiply
   USE cp_fm_basic_linalg,              ONLY: cp_fm_scale,&
                                              cp_fm_transpose,&
                                              cp_fm_upper_to_full
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: &
        cp_fm_create, cp_fm_get_diag, cp_fm_get_info, cp_fm_get_submatrix, cp_fm_release, &
        cp_fm_set_all, cp_fm_set_element, cp_fm_to_fm_submat, cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_unit_nr
   USE dbcsr_api,                       ONLY: &
        dbcsr_copy, dbcsr_create, dbcsr_distribution_get, dbcsr_distribution_new, &
        dbcsr_distribution_release, dbcsr_distribution_type, dbcsr_get_info, dbcsr_multiply, &
        dbcsr_p_type, dbcsr_print, dbcsr_release, dbcsr_reserve_all_blocks, dbcsr_scale, &
        dbcsr_type, dbcsr_type_no_symmetry
   USE input_constants,                 ONLY: tddfpt_dipole_length
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE lebedev,                         ONLY: deallocate_lebedev_grids,&
                                              get_number_of_lebedev_grid,&
                                              init_lebedev_grids,&
                                              lebedev_grid
   USE mathlib,                         ONLY: get_diag
   USE memory_utilities,                ONLY: reallocate
   USE message_passing,                 ONLY: mp_para_env_type
   USE orbital_pointers,                ONLY: indso,&
                                              nsoset
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: evolt,&
                                              wavenumbers
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_grid_atom,                    ONLY: allocate_grid_atom,&
                                              create_grid_atom,&
                                              grid_atom_type
   USE qs_harmonics_atom,               ONLY: allocate_harmonics_atom,&
                                              create_harmonics_atom,&
                                              get_maxl_CG,&
                                              harmonics_atom_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              qs_kind_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type
   USE qs_tddfpt2_soc_types,            ONLY: soc_atom_create,&
                                              soc_atom_env_type,&
                                              soc_atom_release,&
                                              soc_env_create,&
                                              soc_env_release,&
                                              soc_env_type
   USE qs_tddfpt2_soc_utils,            ONLY: atomic_pseudopotential,&
                                              soc_contract_evect,&
                                              soc_dipole_operator
   USE qs_tddfpt2_types,                ONLY: tddfpt_ground_state_mos
   USE spherical_harmonics,             ONLY: clebsch_gordon,&
                                              clebsch_gordon_deallocate,&
                                              clebsch_gordon_init
   USE xas_tdp_atom,                    ONLY: compute_sphi_so,&
                                              integrate_soc_atoms,&
                                              truncate_radial_grid
   USE xas_tdp_utils,                   ONLY: rcs_amew_soc_elements
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_tddfpt2_soc'

   LOGICAL, PARAMETER, PRIVATE          :: debug_this_module = .FALSE.

   INTEGER, PARAMETER, PRIVATE          :: nderivs = 3
   INTEGER, PARAMETER, PRIVATE          :: maxspins = 2

   !A helper type for SOC
   TYPE dbcsr_soc_package_type
      TYPE(dbcsr_type), POINTER     :: dbcsr_sg => Null()
      TYPE(dbcsr_type), POINTER     :: dbcsr_tp => Null()
      TYPE(dbcsr_type), POINTER     :: dbcsr_sc => Null()
      TYPE(dbcsr_type), POINTER     :: dbcsr_sf => Null()
      TYPE(dbcsr_type), POINTER     :: dbcsr_prod => Null()
      TYPE(dbcsr_type), POINTER     :: dbcsr_ovlp => Null()
      TYPE(dbcsr_type), POINTER     :: dbcsr_tmp => Null()
      TYPE(dbcsr_type), POINTER     :: dbcsr_work => Null()
   END TYPE dbcsr_soc_package_type

   PUBLIC :: tddfpt_soc

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief Perform TDDFPT-SOC calculation.
!> \param qs_env  Quickstep environment
!> \param evals_a eigenvalues for the singlet states
!> \param evals_b eigenvalues for the triplet states
!> \param evects_a eigenvectors for the singlet states
!> \param evects_b eigenvectors for the triplet states
!> \param gs_mos ground state orbitlas from the TDDFPT calculation
!> \par History
!>    * 02.2023 created [Jan-Robert Vogt]
!> \note Based on tddfpt2_methods and xas_tdp_utils.
!> \note Only rcs for now, but written with the addition of os in mind
! **************************************************************************************************

   SUBROUTINE tddfpt_soc(qs_env, evals_a, evals_b, evects_a, evects_b, gs_mos)

      TYPE(qs_environment_type), INTENT(IN), POINTER     :: qs_env
      REAL(kind=dp), DIMENSION(:), INTENT(IN)            :: evals_a, evals_b
      TYPE(cp_fm_type), DIMENSION(:, :), INTENT(IN)      :: evects_a, evects_b
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         POINTER                                         :: gs_mos

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'tddfpt_soc'

      INTEGER                                            :: handle, istate, log_unit
      LOGICAL                                            :: do_os
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dft_control_type), POINTER                    :: dft_control

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      IF (logger%para_env%is_source()) THEN
         log_unit = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         log_unit = -1
      END IF

      IF (log_unit > 0) THEN
         WRITE (log_unit, "(1X,A)") "", &
            "-------------------------------------------------------------------------------", &
            "-                         START SOC CALCULATIONS                              -", &
            "-------------------------------------------------------------------------------"
      END IF

      NULLIFY (dft_control)
      CALL get_qs_env(qs_env, dft_control=dft_control)
      do_os = dft_control%uks .OR. dft_control%roks
      IF (do_os) CPABORT("SOC only implemented for closed shell.")

      IF (log_unit > 0) THEN
         WRITE (log_unit, '(A)') "Starting from TDDFPT Excited States:"
         WRITE (log_unit, '(A)') "      STATE          SINGLET/eV         TRIPLET/eV"
         DO istate = 1, SIZE(evals_a)
            WRITE (log_unit, '(6X,I3,11X,F10.5,6X,F10.5)') istate, evals_a(istate)*evolt, evals_b(istate)*evolt
         END DO
      END IF

      IF (log_unit > 0) WRITE (log_unit, '(A)') "Starting restricted closed shell:"
      CALL tddfpt_soc_rcs(qs_env, evals_a, evals_b, evects_a, evects_b, log_unit, gs_mos)

      IF (log_unit > 0) THEN
         WRITE (log_unit, '(A,/,A)') "SOC Calculation terminated", &
            "Returning to TDDFPT for Force calculation and deallocations"
      END IF

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief Will perform the soc-calculation for restricted-closed-shell systems
!> \param qs_env Quickstep Enviroment
!> \param evals_sing eigenvalues singlet states
!> \param evals_trip eigenvalues triplet states
!> \param evects_sing eigenvector singlet states
!> \param evects_trip eigenvectors triplet states
!> \param log_unit default log_unit for convinients
!> \param gs_mos ground state MOs from TDDFPT
!> \par History
!>      * created 02.2023 [Jan-Robert Vogt]
!> \note Mostly copied and modified from xas_tdp_utils include_rcs_soc
! **************************************************************************************************
   SUBROUTINE tddfpt_soc_rcs(qs_env, evals_sing, evals_trip, evects_sing, evects_trip, log_unit, gs_mos)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(kind=dp), DIMENSION(:), INTENT(in), TARGET    :: evals_sing, evals_trip
      TYPE(cp_fm_type), DIMENSION(:, :), INTENT(IN)      :: evects_sing, evects_trip
      INTEGER                                            :: log_unit
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         POINTER                                         :: gs_mos

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'tddfpt_soc_rcs'

      CHARACTER(len=1)                                   :: mult
      INTEGER                                            :: dipole_form, group, handle, iex, ii, &
                                                            isg, istate, itp, jj, nactive, nao, &
                                                            nex, npcols, nprows, nsg, ntot, ntp
      INTEGER, DIMENSION(:), POINTER                     :: col_blk_size, col_dist, row_blk_size, &
                                                            row_dist, row_dist_new
      INTEGER, DIMENSION(:, :), POINTER                  :: pgrid
      LOGICAL                                            :: print_ev, print_some, print_splitting, &
                                                            print_wn
      REAL(dp)                                           :: eps_filter, soc_gst, sqrt2, tmp, &
                                                            unit_scale
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: diag, tmp_evals
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: img_tmp, real_tmp
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :)        :: gstp_block, mo_soc_x, mo_soc_y, mo_soc_z
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_cfm_type)                                  :: evecs_cfm, hami_cfm
      TYPE(cp_fm_struct_type), POINTER                   :: full_struct, gstp_struct, prod_struct, &
                                                            vec_struct, work_struct
      TYPE(cp_fm_type)                                   :: gstp_fm, img_fm, prod_fm, real_fm, &
                                                            tmp_fm, vec_soc_x, vec_soc_y, &
                                                            vec_soc_z, work_fm
      TYPE(cp_fm_type), POINTER                          :: gs_coeffs, tp_coeffs
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_distribution_type), POINTER             :: coeffs_dist, dbcsr_dist, prod_dist
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dbcsr_soc_package_type)                       :: dbcsr_soc_package
      TYPE(dbcsr_type), POINTER                          :: dbcsr_dummy, dbcsr_ovlp, dbcsr_prod, &
                                                            dbcsr_sg, dbcsr_tmp, dbcsr_tp, &
                                                            dbcsr_work, orb_soc_x, orb_soc_y, &
                                                            orb_soc_z
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(section_vals_type), POINTER                   :: soc_print_section, soc_section, &
                                                            tddfpt_print_section
      TYPE(soc_atom_env_type), POINTER                   :: soc_atom_env
      TYPE(soc_env_type), TARGET                         :: soc_env

!! OUTPUT structures

      CALL timeset(routineN, handle)

      NULLIFY (logger, pgrid, row_dist, row_dist_new)
      NULLIFY (col_dist, col_blk_size, row_blk_size, matrix_s)
      NULLIFY (gstp_struct, tddfpt_print_section, soc_print_section, soc_section)

      logger => cp_get_default_logger()
      tddfpt_print_section => section_vals_get_subs_vals(qs_env%input, "PROPERTIES%TDDFPT%PRINT")
      soc_print_section => section_vals_get_subs_vals(qs_env%input, "PROPERTIES%TDDFPT%PRINT%SOC_PRINT")
      soc_section => section_vals_get_subs_vals(qs_env%input, "PROPERTIES%TDDFPT%SOC")

      CALL section_vals_val_get(soc_section, "EPS_FILTER", r_val=eps_filter)

      nsg = SIZE(evals_sing)   ! Number of excited singlet states
      ntp = SIZE(evals_trip)   ! Number of excited triplet states
      nex = nsg                ! Number of excited states of each multiplicity
      ntot = 1 + nsg + 3*ntp ! Number of (GS + S + T^-1 + T^0 + T^1)

      ! Initzialize Working environment
      CALL inititialize_soc(qs_env, soc_atom_env, soc_env, &
                            evects_sing, evects_trip, dipole_form)

      NULLIFY (mos)
      CALL get_qs_env(qs_env, mos=mos, para_env=para_env, blacs_env=blacs_env)
      CALL get_mo_set(mos(1), nao=nao, homo=nactive)

      ! this will create the H^SOC in an atomic basis
      CALL integrate_soc_atoms(soc_env%orb_soc, qs_env=qs_env, soc_atom_env=soc_atom_env)
      CALL soc_atom_release(soc_atom_env)

      !! Point at H^SOC and MOs for better readablity
      NULLIFY (tp_coeffs, orb_soc_x, orb_soc_y, orb_soc_z)
      tp_coeffs => soc_env%b_coeff
      soc_env%evals_a => evals_sing
      soc_env%evals_b => evals_trip
      orb_soc_x => soc_env%orb_soc(1)%matrix
      orb_soc_y => soc_env%orb_soc(2)%matrix
      orb_soc_z => soc_env%orb_soc(3)%matrix

      !! Create a matrix-structure, which links all states in this calculation
      NULLIFY (full_struct)
      CALL cp_fm_struct_create(full_struct, context=blacs_env, para_env=para_env, nrow_global=ntot, ncol_global=ntot)
      CALL cp_fm_create(real_fm, full_struct)
      CALL cp_fm_create(img_fm, full_struct)
      CALL cp_fm_set_all(real_fm, 0.0_dp)
      CALL cp_fm_set_all(img_fm, 0.0_dp)

      !  Put the excitation energies on the diagonal of the real matrix
      DO isg = 1, nsg
         CALL cp_fm_set_element(real_fm, 1 + isg, 1 + isg, evals_sing(isg))
      END DO
      DO itp = 1, ntp
         ! first T^-1, then T^0, then T^+1
         CALL cp_fm_set_element(real_fm, 1 + itp + nsg, 1 + itp + nsg, evals_trip(itp))
         CALL cp_fm_set_element(real_fm, 1 + itp + ntp + nsg, 1 + itp + ntp + nsg, evals_trip(itp))
         CALL cp_fm_set_element(real_fm, 1 + itp + 2*ntp + nsg, 1 + itp + 2*ntp + nsg, evals_trip(itp))
      END DO

      !!Create the dbcsr structures for this calculations
      CALL get_qs_env(qs_env, dbcsr_dist=dbcsr_dist)
      CALL dbcsr_distribution_get(dbcsr_dist, group=group, row_dist=row_dist, pgrid=pgrid, &
                                  npcols=npcols, nprows=nprows)
      ALLOCATE (col_dist(nex), row_dist_new(nex))                   ! Split for each excitation
      DO iex = 1, nex
         col_dist(iex) = MODULO(npcols - iex, npcols)
         row_dist_new(iex) = MODULO(nprows - iex, nprows)
      END DO
      ALLOCATE (coeffs_dist, prod_dist)
      CALL dbcsr_distribution_new(coeffs_dist, group=group, pgrid=pgrid, row_dist=row_dist, &
                                  col_dist=col_dist)
      CALL dbcsr_distribution_new(prod_dist, group=group, pgrid=pgrid, row_dist=row_dist_new, &
                                  col_dist=col_dist)

      !! Create the matrices
      ALLOCATE (col_blk_size(nex))
      col_blk_size = nactive

      CALL get_qs_env(qs_env, matrix_s=matrix_s)
      CALL dbcsr_get_info(matrix_s(1)%matrix, row_blk_size=row_blk_size)

      !! The Eigenvectors for Sg und Tp will be dived into their diffrent components again
      ALLOCATE (dbcsr_sg, dbcsr_tp, dbcsr_work, dbcsr_ovlp, dbcsr_tmp, dbcsr_prod)
      CALL dbcsr_create(matrix=dbcsr_sg, name="SINGLETS", matrix_type=dbcsr_type_no_symmetry, &
                        dist=coeffs_dist, row_blk_size=row_blk_size, col_blk_size=col_blk_size)
      CALL dbcsr_create(matrix=dbcsr_tp, name="TRIPLETS", matrix_type=dbcsr_type_no_symmetry, &
                        dist=coeffs_dist, row_blk_size=row_blk_size, col_blk_size=col_blk_size)
      CALL dbcsr_create(matrix=dbcsr_work, name="WORK", matrix_type=dbcsr_type_no_symmetry, &
                        dist=coeffs_dist, row_blk_size=row_blk_size, col_blk_size=col_blk_size)
      CALL dbcsr_create(matrix=dbcsr_prod, name="PROD", matrix_type=dbcsr_type_no_symmetry, &
                        dist=prod_dist, row_blk_size=col_blk_size, col_blk_size=col_blk_size)
      CALL dbcsr_create(matrix=dbcsr_ovlp, name="OVLP", matrix_type=dbcsr_type_no_symmetry, &
                        dist=prod_dist, row_blk_size=col_blk_size, col_blk_size=col_blk_size)

      col_blk_size = 1
      CALL dbcsr_create(matrix=dbcsr_tmp, name="TMP", matrix_type=dbcsr_type_no_symmetry, &
                        dist=prod_dist, row_blk_size=col_blk_size, col_blk_size=col_blk_size)
      CALL dbcsr_reserve_all_blocks(dbcsr_tmp)

      IF (debug_this_module) THEN
         ALLOCATE (dbcsr_dummy)
         CALL dbcsr_create(matrix=dbcsr_dummy, name="DUMMY", matrix_type=dbcsr_type_no_symmetry, &
                           dist=prod_dist, row_blk_size=col_blk_size, col_blk_size=col_blk_size)
         CALL dbcsr_reserve_all_blocks(dbcsr_dummy)
      END IF

      !! This work dbcsr matrix will be packed together for easy transfer to other subroutines
      dbcsr_soc_package%dbcsr_sg => dbcsr_sg
      dbcsr_soc_package%dbcsr_tp => dbcsr_tp
      dbcsr_soc_package%dbcsr_work => dbcsr_work
      dbcsr_soc_package%dbcsr_ovlp => dbcsr_ovlp
      dbcsr_soc_package%dbcsr_prod => dbcsr_prod
      dbcsr_soc_package%dbcsr_tmp => dbcsr_tmp

      !Filling the coeffs matrices by copying from the stored fms
      CALL copy_fm_to_dbcsr(soc_env%a_coeff, dbcsr_sg)
      CALL copy_fm_to_dbcsr(soc_env%b_coeff, dbcsr_tp)

      !Create the work and helper fms
      !CALL get_mo_set(mos(1),mo_coeff=gs_coeffs)
      gs_coeffs => gs_mos(1)%mos_occ
      CALL cp_fm_get_info(gs_coeffs, matrix_struct=vec_struct)
      CALL cp_fm_create(vec_soc_x, vec_struct)
      CALL cp_fm_create(vec_soc_y, vec_struct)
      CALL cp_fm_create(vec_soc_z, vec_struct)

      CALL cp_fm_struct_create(prod_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=nactive, ncol_global=nactive)
      CALL cp_fm_create(prod_fm, prod_struct)

      CALL cp_fm_struct_create(work_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=nex, ncol_global=nex)
      CALL cp_fm_create(work_fm, work_struct)
      CALL cp_fm_create(tmp_fm, work_struct)

      IF (log_unit > 0 .AND. debug_this_module) THEN
         WRITE (log_unit, '(A)') "Starting Precomputations"
      END IF

      !! Begin with the precomputation
      !! Prefactor due to rcs.
      !! The excitation is presented as a linear combination of
      !! alpha and beta, where dalpha=-dbeta for triplet exc.
      !! and dalpha=dbeta for the singlet case
      sqrt2 = SQRT(2.0_dp)

      !! Precompute the <phi_i^0|H^SOC|phi_j^0> matrix elements
      ALLOCATE (diag(nactive))
      ALLOCATE (mo_soc_x(nactive, nactive), mo_soc_y(nactive, nactive), mo_soc_z(nactive, nactive))

      !! This will be the GS|H|GS contribution needed for all couplings
      CALL cp_dbcsr_sm_fm_multiply(orb_soc_x, gs_coeffs, vec_soc_x, ncol=nactive)
      CALL parallel_gemm('T', 'N', nactive, &
                         nactive, nao, 1.0_dp, &
                         gs_coeffs, vec_soc_x, &
                         0.0_dp, prod_fm)
      CALL cp_fm_get_submatrix(prod_fm, mo_soc_x)

      CALL cp_dbcsr_sm_fm_multiply(orb_soc_y, gs_coeffs, vec_soc_y, ncol=nactive)
      CALL parallel_gemm('T', 'N', nactive, &
                         nactive, nao, 1.0_dp, &
                         gs_coeffs, vec_soc_y, &
                         0.0_dp, prod_fm)
      CALL cp_fm_get_submatrix(prod_fm, mo_soc_y)

      CALL cp_dbcsr_sm_fm_multiply(orb_soc_z, gs_coeffs, vec_soc_z, ncol=nactive)
      CALL parallel_gemm('T', 'N', nactive, &
                         nactive, nao, 1.0_dp, &
                         gs_coeffs, vec_soc_z, &
                         0.0_dp, prod_fm)
      CALL cp_fm_get_submatrix(prod_fm, mo_soc_z)

      !  Only have SOC between singlet-triplet triplet-triplet and ground_state-triplet, the resulting
      !  matrix is Hermitian i.e. the real part is symmetric and the imaginary part is anti-symmetric.
      !  Can only fill upper half

      IF (log_unit > 0 .AND. debug_this_module) THEN
         WRITE (log_unit, '(A)') "Starting Ground-State contributions..."
      END IF
      !Start with the ground state/triplet SOC, SOC*gs_coeffs already computed above
      !note: we are computing <0|H|T>, but have SOC*gs_coeffs instead of gs_coeffs*SOC in store.
      !Since the SOC Hamiltonian is anti-symmetric, a - signs pops up in the gemms below
      CALL cp_fm_struct_create(gstp_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=ntp*nactive, ncol_global=nactive)
      CALL cp_fm_create(gstp_fm, gstp_struct)
      ALLOCATE (gstp_block(nactive, nactive))

      !gs-triplet with Ms=+-1, imaginary part
      ! <T+-1|H_x|GS>
      ! -1 to change to <GS|H|T+-1>
      CALL parallel_gemm('T', 'N', nactive*ntp, &
                         nactive, nao, -1.0_dp, &
                         tp_coeffs, vec_soc_x, &
                         0.0_dp, gstp_fm)

      !! Seperate them into the different states again (nactive x nactive)
      DO itp = 1, ntp
         CALL cp_fm_get_submatrix(fm=gstp_fm, target_m=gstp_block, &
                                  start_row=(itp - 1)*nactive + 1, start_col=1, &
                                  n_rows=nactive, n_cols=nactive)
         diag(:) = get_diag(gstp_block)
         soc_gst = SUM(diag)
         !!  <0|H_x|T^-1>
         CALL cp_fm_set_element(img_fm, 1, 1 + nsg + itp, -1_dp*soc_gst)
         !! <0|H_x|T^+1>
         CALL cp_fm_set_element(img_fm, 1, 1 + nsg + 2*ntp + itp, soc_gst)
      END DO

      IF (debug_this_module .AND. (log_unit > 0)) THEN
         WRITE (log_unit, "(A,2F18.6)") "<0|H_x|T^-1> in Hartree / cm^-1", soc_gst, soc_gst*wavenumbers
         WRITE (log_unit, "(A,2F18.6)") "<0|H_x|T^+1> in Hartree / cm^-1", soc_gst, soc_gst*wavenumbers
      END IF

      !gs-triplet with Ms=+-1, real part
      ! <T+-1|H_y|GS>
      CALL parallel_gemm('T', 'N', nactive*ntp, &
                         nactive, nao, -1.0_dp, &
                         tp_coeffs, vec_soc_y, &
                         0.0_dp, gstp_fm)

      DO itp = 1, ntp
         CALL cp_fm_get_submatrix(fm=gstp_fm, target_m=gstp_block, start_row=(itp - 1)*nactive + 1, &
                                  start_col=1, n_rows=nactive, n_cols=nactive)
         diag(:) = get_diag(gstp_block)
         soc_gst = SUM(diag)
         ! <0|H_y|T^-1>
         CALL cp_fm_set_element(real_fm, 1, 1 + nsg + itp, -1.0_dp*soc_gst)
         ! <0|H_y|T^+1>
         CALL cp_fm_set_element(real_fm, 1, 1 + nsg + 2*ntp + itp, -1.0_dp*soc_gst)
      END DO

      IF (debug_this_module .AND. (log_unit > 0)) THEN
         WRITE (log_unit, "(A,2F18.6)") "<0|H_y|T^-1> in Hartree / cm^-1", soc_gst, soc_gst*wavenumbers
         WRITE (log_unit, "(A,2F18.6)") "<0|H_y|T^+1> in Hartree / cm^-1", soc_gst, soc_gst*wavenumbers
      END IF

      !gs-triplet with Ms=0, purely imaginary
      !< T0|H_z|GS>
      CALL parallel_gemm('T', 'N', nactive*ntp, &
                         nactive, nao, -1.0_dp, &
                         tp_coeffs, vec_soc_z, &
                         0.0_dp, gstp_fm)

      DO itp = 1, ntp
         CALL cp_fm_get_submatrix(fm=gstp_fm, &
                                  target_m=gstp_block, &
                                  start_row=(itp - 1)*nactive + 1, &
                                  start_col=1, &
                                  n_rows=nactive, &
                                  n_cols=nactive)
         diag(:) = get_diag(gstp_block)
         soc_gst = sqrt2*SUM(diag)
         CALL cp_fm_set_element(img_fm, 1, 1 + nsg + ntp + itp, soc_gst)
      END DO

      IF (debug_this_module .AND. (log_unit > 0)) THEN
         WRITE (log_unit, "(A,2F18.6)") "<0|H_z|T> in Hartree / cm^-1", soc_gst, soc_gst*wavenumbers
      END IF

      !! After all::
      !! T-1 :: -<0|H_x|T^-1> - <0|H_y|T^-1> at 1+nsg+itp
      !! T+1 :: <0|H_x|T^-1> - <0|H_y|T^-1> at 1+nsp+2tnp+itp
      !! T0  :: < T0|H_z|GS>                at 1+nsg+ntp+itp

      !gs clean-up
      CALL cp_fm_release(prod_fm)
      CALL cp_fm_release(vec_soc_x)
      CALL cp_fm_release(vec_soc_y)
      CALL cp_fm_release(vec_soc_z)
      CALL cp_fm_release(gstp_fm)
      CALL cp_fm_struct_release(gstp_struct)
      CALL cp_fm_struct_release(prod_struct)
      DEALLOCATE (gstp_block)

      IF (log_unit > 0 .AND. debug_this_module) THEN
         WRITE (log_unit, '(A)') "Starting Singlet-Triplet contributions..."
      END IF

      !Now do the singlet-triplet SOC
      !start by computing the singlet-triplet overlap
      ! <S|T>
      CALL dbcsr_multiply('N', 'N', 1.0_dp, &
                          matrix_s(1)%matrix, &
                          dbcsr_tp, 0.0_dp, &
                          dbcsr_work, &
                          filter_eps=eps_filter)
      CALL dbcsr_multiply('T', 'N', 1.0_dp, &
                          dbcsr_sg, &
                          dbcsr_work, &
                          0.0_dp, dbcsr_ovlp, &
                          filter_eps=eps_filter)

      !singlet-triplet with Ms=+-1, imaginary part
      ! First precalculate <S|H_x|T+-1>
      CALL dbcsr_multiply('N', 'N', 1.0_dp, &
                          orb_soc_x, &
                          dbcsr_tp, 0.0_dp, &
                          dbcsr_work, &
                          filter_eps=eps_filter)
      CALL dbcsr_multiply('T', 'N', 1.0_dp, &
                          dbcsr_sg, &
                          dbcsr_work, 0.0_dp, &
                          dbcsr_prod, &
                          filter_eps=eps_filter)

      !! This will lead to:
      !! -1/sqrt(2)(<S|H_x|T> - <S|T> <GS|H_x|GS>)
      CALL rcs_amew_soc_elements(dbcsr_tmp, &
                                 dbcsr_prod, &
                                 dbcsr_ovlp, &
                                 mo_soc_x, &
                                 pref_trace=-1.0_dp, &
                                 pref_overall=-0.5_dp*sqrt2)

      !<S|H_x|T^-1>
      !! Convert to fm for transfer to img_fm
      CALL copy_dbcsr_to_fm(dbcsr_tmp, tmp_fm)
      CALL cp_fm_to_fm_submat(msource=tmp_fm, &
                              mtarget=img_fm, &
                              nrow=nex, &
                              ncol=nex, &
                              s_firstrow=1, &
                              s_firstcol=1, &
                              t_firstrow=2, &
                              t_firstcol=1 + nsg + 1)

      IF (debug_this_module) THEN
         WRITE (log_unit, "(/,A)") "<S|H_x|T^-1> in Hartree / cm^-1"
         CALL dbcsr_copy(dbcsr_dummy, dbcsr_tmp)
         CALL dbcsr_scale(dbcsr_dummy, wavenumbers)
         CALL dbcsr_print(dbcsr_dummy, unit_nr=log_unit)
      END IF

      !<S|H_x|T^+1> takes a minus sign
      CALL cp_fm_scale(-1.0_dp, tmp_fm)
      CALL cp_fm_to_fm_submat(msource=tmp_fm, &
                              mtarget=img_fm, &
                              nrow=nex, &
                              ncol=nex, &
                              s_firstrow=1, &
                              s_firstcol=1, &
                              t_firstrow=2, &
                              t_firstcol=1 + nsg + 2*ntp + 1)

      IF (debug_this_module) THEN
         WRITE (log_unit, "(/,A)") "<S|H_x|T^+1> in Hartree / cm^-1"
         CALL dbcsr_copy(dbcsr_dummy, dbcsr_tmp)
         CALL dbcsr_scale(dbcsr_dummy, wavenumbers)
         CALL dbcsr_print(dbcsr_dummy, unit_nr=log_unit)
      END IF

       !!singlet-triplet with Ms=+-1, real part
       !! Precompute <S|H_y|T>
      CALL dbcsr_multiply('N', 'N', 1.0_dp, &
                          orb_soc_y, dbcsr_tp, &
                          0.0_dp, dbcsr_work, &
                          filter_eps=eps_filter)
      CALL dbcsr_multiply('T', 'N', 1.0_dp, &
                          dbcsr_sg, dbcsr_work, &
                          0.0_dp, dbcsr_prod, &
                          filter_eps=eps_filter)

      !! This will lead to -1/sqrt(2)(<S|H_y|T> - <S|T> <GS|H_y|GS>)
      CALL rcs_amew_soc_elements(dbcsr_tmp, &
                                 dbcsr_prod, &
                                 dbcsr_ovlp, &
                                 mo_soc_y, &
                                 pref_trace=-1.0_dp, &
                                 pref_overall=-0.5_dp*sqrt2)

      !<S|H_y|T^-1>
      CALL copy_dbcsr_to_fm(dbcsr_tmp, tmp_fm)
      CALL cp_fm_to_fm_submat(msource=tmp_fm, &
                              mtarget=real_fm, &
                              nrow=nex, &
                              ncol=nex, &
                              s_firstrow=1, &
                              s_firstcol=1, &
                              t_firstrow=2, &
                              t_firstcol=1 + nsg + 1)

      IF (debug_this_module) THEN
         WRITE (log_unit, "(/,A)") "<S|H_y|T^-1> in Hartree / cm^-1"
         CALL dbcsr_copy(dbcsr_dummy, dbcsr_tmp)
         CALL dbcsr_scale(dbcsr_dummy, wavenumbers)
         CALL dbcsr_print(dbcsr_tmp, unit_nr=log_unit)
      END IF

      !<S|H_y|T^+1>
      CALL cp_fm_to_fm_submat(msource=tmp_fm, &
                              mtarget=real_fm, &
                              nrow=nex, &
                              ncol=nex, &
                              s_firstrow=1, &
                              s_firstcol=1, &
                              t_firstrow=2, &
                              t_firstcol=1 + nsg + 2*ntp + 1)

      IF (debug_this_module) THEN
         WRITE (log_unit, "(/,A)") "<S|H_y|T^+1> in Hartree / cm^-1"
         CALL dbcsr_copy(dbcsr_dummy, dbcsr_tmp)
         CALL dbcsr_scale(dbcsr_dummy, wavenumbers)
         CALL dbcsr_print(dbcsr_dummy, unit_nr=log_unit)
      END IF

      !singlet-triplet with Ms=0, purely imaginary
      CALL dbcsr_multiply('N', 'N', 1.0_dp, &
                          orb_soc_z, dbcsr_tp, &
                          0.0_dp, dbcsr_work, &
                          filter_eps=eps_filter)
      CALL dbcsr_multiply('T', 'N', 1.0_dp, &
                          dbcsr_sg, dbcsr_work, &
                          0.0_dp, dbcsr_prod, &
                          filter_eps=eps_filter)

      CALL rcs_amew_soc_elements(dbcsr_tmp, &
                                 dbcsr_prod, &
                                 dbcsr_ovlp, &
                                 mo_soc_z, &
                                 pref_trace=-1.0_dp, &
                                 pref_overall=1.0_dp)

      !<S|H_z|T^0>
      CALL copy_dbcsr_to_fm(dbcsr_tmp, tmp_fm)
      CALL cp_fm_to_fm_submat(msource=tmp_fm, &
                              mtarget=img_fm, &
                              nrow=nex, &
                              ncol=nex, &
                              s_firstrow=1, &
                              s_firstcol=1, &
                              t_firstrow=2, &
                              t_firstcol=1 + nsg + ntp + 1)

      IF (debug_this_module) THEN
         WRITE (log_unit, "(/,A)") "<S|H_z|T^0> in Hartree / cm^-1"
         CALL dbcsr_copy(dbcsr_dummy, dbcsr_tmp)
         CALL dbcsr_scale(dbcsr_dummy, wavenumbers)
         CALL dbcsr_print(dbcsr_dummy, unit_nr=log_unit)
      END IF

      IF (log_unit > 0 .AND. debug_this_module) THEN
         WRITE (log_unit, '(A)') "Starting Triplet-Triplet contributions..."
      END IF

      !Now the triplet-triplet SOC
      !start by computing the overlap
      CALL dbcsr_multiply('N', 'N', 1.0_dp, &
                          matrix_s(1)%matrix, &
                          dbcsr_tp, 0.0_dp, &
                          dbcsr_work, &
                          filter_eps=eps_filter)
      CALL dbcsr_multiply('T', 'N', 1.0_dp, &
                          dbcsr_tp, dbcsr_work, &
                          0.0_dp, dbcsr_ovlp, &
                          filter_eps=eps_filter)

      !Ms=0 to Ms=+-1 SOC, imaginary part
      CALL dbcsr_multiply('N', 'N', 1.0_dp, &
                          orb_soc_x, dbcsr_tp, &
                          0.0_dp, dbcsr_work, &
                          filter_eps=eps_filter)
      CALL dbcsr_multiply('T', 'N', 1.0_dp, &
                          dbcsr_tp, dbcsr_work, &
                          0.0_dp, dbcsr_prod, &
                          filter_eps=eps_filter)

      CALL rcs_amew_soc_elements(dbcsr_tmp, &
                                 dbcsr_prod, &
                                 dbcsr_ovlp, &
                                 mo_soc_x, &
                                 pref_trace=1.0_dp, &
                                 pref_overall=-0.5_dp*sqrt2)

      !<T^0|H_x|T^+1>
      CALL copy_dbcsr_to_fm(dbcsr_tmp, tmp_fm)
      CALL cp_fm_to_fm_submat(msource=tmp_fm, &
                              mtarget=img_fm, &
                              nrow=nex, &
                              ncol=nex, &
                              s_firstrow=1, &
                              s_firstcol=1, &
                              t_firstrow=1 + nsg + ntp + 1, &
                              t_firstcol=1 + nsg + 2*ntp + 1)

      IF (debug_this_module) THEN
         WRITE (log_unit, "(/,A)") "<T^0|H_x|T^+1> in Hartree / cm^-1"
         CALL dbcsr_copy(dbcsr_dummy, dbcsr_tmp)
         CALL dbcsr_scale(dbcsr_dummy, wavenumbers)
         CALL dbcsr_print(dbcsr_dummy, unit_nr=log_unit)
      END IF

      !<T^-1|H_x|T^0>, takes a minus sign and a transpose (because computed <T^0|H_x|T^-1>)
      CALL cp_fm_transpose(tmp_fm, work_fm)
      CALL cp_fm_scale(-1.0_dp, work_fm)
      CALL cp_fm_to_fm_submat(msource=work_fm, &
                              mtarget=img_fm, &
                              nrow=nex, &
                              ncol=nex, &
                              s_firstrow=1, &
                              s_firstcol=1, &
                              t_firstrow=1 + nsg + 1, &
                              t_firstcol=1 + nsg + ntp + 1)

      IF (debug_this_module) THEN
         WRITE (log_unit, "(/,A)") "<T^-1|H_x|T^0> in Hartree / cm^-1"
         CALL dbcsr_copy(dbcsr_dummy, dbcsr_tmp)
         CALL dbcsr_scale(dbcsr_dummy, wavenumbers)
         CALL dbcsr_print(dbcsr_dummy, unit_nr=log_unit)
      END IF

      !Ms=0 to Ms=+-1 SOC, real part
      CALL dbcsr_multiply('N', 'N', 1.0_dp, &
                          orb_soc_y, dbcsr_tp, &
                          0.0_dp, dbcsr_work, &
                          filter_eps=eps_filter)
      CALL dbcsr_multiply('T', 'N', 1.0_dp, &
                          dbcsr_tp, dbcsr_work, &
                          0.0_dp, dbcsr_prod, &
                          filter_eps=eps_filter)

      CALL rcs_amew_soc_elements(dbcsr_tmp, &
                                 dbcsr_prod, &
                                 dbcsr_ovlp, &
                                 mo_soc_y, &
                                 pref_trace=1.0_dp, &
                                 pref_overall=0.5_dp*sqrt2)

      !<T^0|H_y|T^+1>
      CALL copy_dbcsr_to_fm(dbcsr_tmp, tmp_fm)
      CALL cp_fm_to_fm_submat(msource=tmp_fm, &
                              mtarget=real_fm, &
                              nrow=nex, &
                              ncol=nex, &
                              s_firstrow=1, &
                              s_firstcol=1, t_firstrow=1 + nsg + ntp + 1, &
                              t_firstcol=1 + nsg + 2*ntp + 1)

      IF (debug_this_module) THEN
         WRITE (log_unit, "(/,A)") "<T^0|H_y|T^+1> in Hartree / cm^-1"
         CALL dbcsr_copy(dbcsr_dummy, dbcsr_tmp)
         CALL dbcsr_scale(dbcsr_dummy, wavenumbers)
         CALL dbcsr_print(dbcsr_dummy, unit_nr=log_unit)
      END IF

      !<T^-1|H_y|T^0>, takes a minus sign and a transpose
      CALL cp_fm_transpose(tmp_fm, work_fm)
      CALL cp_fm_scale(-1.0_dp, work_fm)
      CALL cp_fm_to_fm_submat(msource=work_fm, &
                              mtarget=real_fm, &
                              nrow=nex, &
                              ncol=nex, &
                              s_firstrow=1, &
                              s_firstcol=1, t_firstrow=1 + nsg + 1, &
                              t_firstcol=1 + nsg + ntp + 1)

      IF (debug_this_module) THEN
         WRITE (log_unit, "(/,A)") "<T^0|H_y|T^+1> in Hartree / cm^-1"
         CALL dbcsr_copy(dbcsr_dummy, dbcsr_tmp)
         CALL dbcsr_scale(dbcsr_dummy, wavenumbers)
         CALL dbcsr_print(dbcsr_dummy, unit_nr=log_unit)
      END IF

      !Ms=1 to Ms=1 and Ms=-1 to Ms=-1 SOC, purely imaginary
      CALL dbcsr_multiply('N', 'N', 1.0_dp, &
                          orb_soc_z, dbcsr_tp, &
                          0.0_dp, dbcsr_work, &
                          filter_eps=eps_filter)
      CALL dbcsr_multiply('T', 'N', 1.0_dp, &
                          dbcsr_tp, dbcsr_work, &
                          0.0_dp, dbcsr_prod, &
                          filter_eps=eps_filter)

      CALL rcs_amew_soc_elements(dbcsr_tmp, &
                                 dbcsr_prod, &
                                 dbcsr_ovlp, &
                                 mo_soc_z, &
                                 pref_trace=1.0_dp, &
                                 pref_overall=1.0_dp)

      !<T^+1|H_z|T^+1>
      CALL copy_dbcsr_to_fm(dbcsr_tmp, tmp_fm)
      CALL cp_fm_to_fm_submat(msource=tmp_fm, &
                              mtarget=img_fm, &
                              nrow=nex, &
                              ncol=nex, &
                              s_firstrow=1, &
                              s_firstcol=1, &
                              t_firstrow=1 + nsg + 2*ntp + 1, &
                              t_firstcol=1 + nsg + 2*ntp + 1)

      IF (debug_this_module) THEN
         WRITE (log_unit, "(/,A)") "<T^+1|H_z|T^+1> in Hartree / cm^-1"
         CALL dbcsr_copy(dbcsr_dummy, dbcsr_tmp)
         CALL dbcsr_scale(dbcsr_dummy, wavenumbers)
         CALL dbcsr_print(dbcsr_dummy, unit_nr=log_unit)
      END IF

      !<T^-1|H_z|T^-1>, takes a minus sign
      CALL cp_fm_scale(-1.0_dp, tmp_fm)
      CALL cp_fm_to_fm_submat(msource=tmp_fm, &
                              mtarget=img_fm, &
                              nrow=nex, &
                              ncol=nex, &
                              s_firstrow=1, &
                              s_firstcol=1, &
                              t_firstrow=1 + nsg + 1, &
                              t_firstcol=1 + nsg + 1)

      IF (debug_this_module) THEN
         WRITE (log_unit, "(/,A)") "<T^-1|H_z|T^-1> in Hartree / cm^-1"
         CALL dbcsr_copy(dbcsr_dummy, dbcsr_tmp)
         CALL dbcsr_scale(dbcsr_dummy, wavenumbers)
         CALL dbcsr_print(dbcsr_dummy, unit_nr=log_unit)
      END IF

      !Intermediate clean-up
      CALL cp_fm_struct_release(work_struct)
      CALL cp_fm_release(work_fm)
      CALL cp_fm_release(tmp_fm)
      DEALLOCATE (diag, mo_soc_x, mo_soc_y, mo_soc_z)

      IF (log_unit > 0) THEN
         WRITE (log_unit, '(A)') "Diagonlzing SOC-Matrix"
      END IF

      !Set-up the complex hermitian perturbation matrix
      CALL cp_cfm_create(hami_cfm, full_struct)
      CALL cp_fm_to_cfm(real_fm, img_fm, hami_cfm)

      !!Optinal Output: SOC-Matrix
      IF (logger%para_env%is_source()) THEN
         log_unit = cp_print_key_unit_nr(logger, &
                                         tddfpt_print_section, &
                                         "SOC_PRINT", &
                                         extension=".socme", &
                                         file_form="FORMATTED", &
                                         file_action="WRITE", &
                                         file_status="UNKNOWN")
      ELSE
         log_unit = -1
      END IF

      !! Get the requested energy unit and optional outputs
      CALL section_vals_val_get(soc_print_section, "UNIT_eV", l_val=print_ev)
      CALL section_vals_val_get(soc_print_section, "UNIT_wn", l_val=print_wn)
      CALL section_vals_val_get(soc_print_section, "SPLITTING", l_val=print_splitting)
      CALL section_vals_val_get(soc_print_section, "SOME", l_val=print_some)

      !! save convertion unit into different varible for cleaner code
      IF (print_wn) THEN
         unit_scale = wavenumbers
      ELSE IF (print_ev) THEN
         unit_scale = evolt
      ELSE
         CPWARN("No output unit was selected, printout will be in Hartree!")
         unit_scale = 1
      END IF

      !! This needs something to deal with a parallel run
      !! This output will be needed for NEWTONX for ISC!!
      IF (print_some) THEN
         IF (log_unit > 0) THEN
            WRITE (log_unit, '(A)') "____________________________________________________"
            WRITE (log_unit, '(A)') "                   FULL SOC-Matrix                  "
            IF (print_ev) THEN
               WRITE (log_unit, '(A)') "STATE   STATE   REAL-PART[eV]       IMG-PART[eV]    "
            ELSE IF (print_wn) THEN
               WRITE (log_unit, '(A)') "STATE   STATE   REAL-PART[cm-1]       IMG-PART[cm-1]"
            ELSE
               WRITE (log_unit, '(A)') "STATE   STATE   REAL-PART[Hartree] IMG-PART[Hartree]"
            END IF
            WRITE (log_unit, '(A)') "____________________________________________________"
         END IF
         ALLOCATE (real_tmp(ntot, ntot), img_tmp(ntot, ntot))
         real_tmp = 0_dp
         img_tmp = 0_dp

         CALL cp_fm_get_submatrix(real_fm, real_tmp, 1, 1, ntot, ntot)
         CALL cp_fm_get_submatrix(img_fm, img_tmp, 1, 1, ntot, ntot)

         IF (log_unit > 0) THEN
            DO jj = 1, ntot
               DO ii = 1, ntot
                  WRITE (unit=log_unit, fmt="(I3,5X,I3,5X,f16.8,5X,f16.8)") ii, jj, &
                     real_tmp(ii, jj)*unit_scale, &
                     img_tmp(ii, jj)*unit_scale
               END DO !! jj
            END DO        !! ii
            WRITE (log_unit, '(A)') "                   END SOC-MATRIX                   "
            WRITE (log_unit, '(A)') "____________________________________________________"
            WRITE (log_unit, '(A)') "____________________________________________________"
         END IF !(log_unit)
         DEALLOCATE (real_tmp, img_tmp)
      END IF !print_some

      CALL cp_fm_release(real_fm)
      CALL cp_fm_release(img_fm)

      ! Diagonalize the Hamiltonian
      ALLOCATE (tmp_evals(ntot))
      CALL cp_cfm_create(evecs_cfm, full_struct)
      CALL cp_cfm_heevd(hami_cfm, evecs_cfm, tmp_evals)

      !  Adjust the energies so the GS has zero, and store in the donor_state (without the GS)
      ALLOCATE (soc_env%soc_evals(ntot - 1))
      soc_env%soc_evals(:) = tmp_evals(2:ntot) - tmp_evals(1)

      !! We may be interested in the ground state stabilisation energy
      IF (logger%para_env%is_source()) THEN
         log_unit = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         log_unit = -1
      END IF

      IF (log_unit > 0) THEN
         WRITE (log_unit, '(A27,6X,F18.10)') "Ground state stabilisation:", (soc_env%soc_evals(1)*unit_scale)
         IF (debug_this_module) THEN
            WRITE (log_unit, '(A)') "Calculation Dipole Moments"
         END IF
      END IF

      !Compute the dipole oscillator strengths
      soc_env%evals_a => evals_sing
      soc_env%evals_b => evals_trip
      CALL soc_dipol(qs_env, dbcsr_soc_package, soc_env, &
                     evecs_cfm, eps_filter, dipole_form, gs_mos)

      !Create final output
      !! Output unit is choosen by the keyword "UNIT_eV" and "UNIT_wn" in "SOC_PRINT" section
      IF (logger%para_env%is_source()) THEN
         log_unit = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         log_unit = -1
      END IF

      IF (log_unit > 0) THEN
         WRITE (log_unit, '(A)') "------------------------------------------------------------------------------"
         WRITE (log_unit, '(A)') "-                         SOC CORRECTED SPECTRUM                             -"
         WRITE (log_unit, '(A)') "------------------------------------------------------------------------------"
         IF (print_ev) THEN
            WRITE (log_unit, '(A)') "-     STATE   SOC-corrected exc. energies [eV]  Oscillator strengths [a.u.]  -"
         ELSE IF (print_wn) THEN
            WRITE (log_unit, '(A)') "-     STATE   SOC-corrected exc. energies [cm-1] Oscillator strengths [a.u.]-"
         ELSE
            WRITE (log_unit, '(A)') "-     STATE SOC-corrected exc.energies [Hartree] Oscillator strengths [a.u.]-"
         END IF
         WRITE (log_unit, '(A)') "------------------------------------------------------------------------------"
         DO istate = 1, SIZE(soc_env%soc_evals)
            WRITE (log_unit, '(6X,I3,11X,2F16.5)') istate, soc_env%soc_evals(istate)*unit_scale, soc_env%soc_osc(istate)
         END DO
         !! This is used for regtesting
         WRITE (log_unit, '(A,F10.6)') "TDDFT+SOC : CheckSum (eV) ", SUM(soc_env%soc_evals)*evolt
      END IF

      !! We may be interested in the differenz between the SOC and
      !! TDDFPT Eigenvalues
      !! Activated with the "SPLITTING" keyword in "SOC_PRINT" section
      IF (log_unit > 0 .AND. (print_splitting .OR. debug_this_module)) THEN
         isg = 1
         itp = 1
         iex = 1
         tmp = 0.0
         mult = ""
         WRITE (log_unit, '(A)') "------------------------------------------------------------------------"
         WRITE (log_unit, '(A)') "-                            Splittings                                -"
         WRITE (log_unit, '(A)') "------------------------------------------------------------------------"
         IF (print_ev) THEN
            WRITE (log_unit, '(A)') "-  STATE   MULT     exc.energies[eV]    splittings[eV]                 -"
         ELSE IF (print_wn) THEN
            WRITE (log_unit, '(A)') "-  STATE   MULT     exc.energies[cm-1]  splittings[cm-1]               -"
         ELSE
            WRITE (log_unit, '(A)') "-  STATE   MULT     exc.energies[Hartree]  splittings[Hartree]         -"
         END IF
         WRITE (log_unit, '(A)') "------------------------------------------------------------------------"
         DO WHILE (iex <= SIZE(soc_env%soc_evals))
            IF (itp < SIZE(evals_trip) .AND. isg < SIZE(evals_sing)) THEN
               IF (evals_sing(isg) < evals_trip(itp)) THEN
                  tmp = soc_env%soc_evals(iex) - evals_sing(isg)
                  mult = "S"
                  isg = isg + 1
               ELSE
                  tmp = soc_env%soc_evals(iex) - evals_trip(itp)
                  mult = "T"
                  itp = itp + 1
               END IF
            ELSE IF (itp > SIZE(evals_trip)) THEN
               tmp = soc_env%soc_evals(iex) - evals_sing(isg)
               mult = "S"
               isg = isg + 1
            ELSE
               tmp = soc_env%soc_evals(iex) - evals_trip(itp)
               mult = "T"
               itp = itp + 1
            END IF
            WRITE (log_unit, '(3X,I3,7X,A,7X,F16.8,7X,F16.8)') iex, mult, &
               soc_env%soc_evals(iex)*unit_scale, tmp*unit_scale
            iex = iex + 1
            IF (mult == "T") THEN
               DO ii = 1, 2
                  tmp = soc_env%soc_evals(iex) - evals_trip(itp - 1)
                  WRITE (log_unit, '(3X,I3,7X,A,7X,F16.8,7X,F16.8)') iex, mult, &
                     soc_env%soc_evals(iex)*unit_scale, tmp*unit_scale
                  iex = iex + 1
               END DO
            END IF
         END DO
         WRITE (log_unit, '(A)') "------------------------------------------------------------------------"
      END IF

      !! clean up
      CALL soc_env_release(soc_env)
      CALL cp_fm_struct_release(full_struct)
      CALL cp_cfm_release(hami_cfm)
      CALL cp_cfm_release(evecs_cfm)
      CALL dbcsr_distribution_release(coeffs_dist)
      CALL dbcsr_distribution_release(prod_dist)
      CALL dbcsr_release(dbcsr_sg)
      CALL dbcsr_release(dbcsr_tp)
      CALL dbcsr_release(dbcsr_prod)
      CALL dbcsr_release(dbcsr_ovlp)
      CALL dbcsr_release(dbcsr_tmp)
      CALL dbcsr_release(dbcsr_work)
      IF (debug_this_module) THEN
         CALL dbcsr_release(dbcsr_dummy)
         DEALLOCATE (dbcsr_dummy)
      END IF
      DEALLOCATE (dbcsr_work, dbcsr_prod, dbcsr_ovlp, dbcsr_tmp)
      DEALLOCATE (coeffs_dist, prod_dist, col_dist, col_blk_size, row_dist_new)
      DEALLOCATE (dbcsr_sg, dbcsr_tp, tmp_evals)

      CALL timestop(handle)

   END SUBROUTINE tddfpt_soc_rcs

! **************************************************************************************************
!> \brief Some additional initalizations
!> \param qs_env Quickstep environment
!> \param soc_atom_env contains information to build the ZORA-Hamiltionian
!> \param soc_env contails details and outputs for SOC Correction
!> \param evects_a singlet Eigenvectors
!> \param evects_b triplet eigenvectors
!> \param dipole_form choosen dipole-form from TDDFPT-Section
! **************************************************************************************************
   SUBROUTINE inititialize_soc(qs_env, soc_atom_env, soc_env, &
                               evects_a, evects_b, &
                               dipole_form)
      !! Different Types
      TYPE(qs_environment_type), INTENT(IN), POINTER     :: qs_env
      TYPE(soc_atom_env_type), INTENT(OUT), POINTER      :: soc_atom_env
      TYPE(soc_env_type), INTENT(OUT), TARGET            :: soc_env
      TYPE(cp_fm_type), DIMENSION(:, :), INTENT(IN)      :: evects_a, evects_b
      INTEGER                                            :: dipole_form

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'inititialize_soc'

      CHARACTER(len=default_string_length), &
         ALLOCATABLE, DIMENSION(:, :)                    :: grid_info
      CHARACTER(len=default_string_length), &
         DIMENSION(:), POINTER                           :: k_list
      INTEGER                                            :: handle, i, ikind, irep, ispin, nactive, &
                                                            nao, natom, nkind, nrep, nspins, &
                                                            nstates
      LOGICAL                                            :: all_potential_present, do_os
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: a_coeff, b_coeff
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(section_vals_type), POINTER                   :: soc_section
      TYPE(tddfpt2_control_type), POINTER                :: tddfpt_control

      CALL timeset(routineN, handle)

      NULLIFY (qs_kind_set, particle_set, blacs_env, para_env, &
               a_coeff, b_coeff, tddfpt_control)

      !! Open_shell is not implemented yet
      do_os = .FALSE.
      !IF (SIZE(evects_a,1) > 1) do_os=.TRUE.

      !! soc_env will pass most of the larger matrices to the different modules
      CALL soc_env_create(soc_env)
      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set, &
                      mos=mos, natom=natom, &
                      particle_set=particle_set, blacs_env=blacs_env, &
                      para_env=para_env, matrix_s=matrix_s, &
                      dft_control=dft_control)

      !! The Dipole form shall be consistent with the tddfpt-module!
      tddfpt_control => dft_control%tddfpt2_control
      dipole_form = tddfpt_control%dipole_form

      NULLIFY (soc_section)
      soc_section => section_vals_get_subs_vals(qs_env%input, "PROPERTIES%TDDFPT%SOC")

      !! We may want to use a bigger grid then default
      CALL section_vals_val_get(soc_section, "GRID", n_rep_val=nrep)
      ALLOCATE (grid_info(nrep, 3))
      DO irep = 1, nrep
         CALL section_vals_val_get(soc_section, &
                                   "GRID", &
                                   i_rep_val=irep, &
                                   c_vals=k_list)
         grid_info(irep, :) = k_list
      END DO

      CALL soc_dipole_operator(soc_env%dipmat, tddfpt_control, qs_env)

      !! We need to change and create structures for the exitation vectors.
      !! All excited states shall be written in one matrix, oneafter the other
      nstates = SIZE(evects_a, 2)
      nspins = SIZE(evects_a, 1)
      DO ispin = 1, nspins
         CALL cp_fm_get_info(evects_a(ispin, 1), nrow_global=nao, ncol_global=nactive)
      END DO

      a_coeff => soc_env%a_coeff
      b_coeff => soc_env%b_coeff

      NULLIFY (fm_struct)
      CALL cp_fm_struct_create(fm_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=nao, ncol_global=nspins*nactive*nstates)
      CALL cp_fm_create(a_coeff, fm_struct)
      CALL cp_fm_create(b_coeff, fm_struct)
      CALL cp_fm_struct_release(fm_struct)

      CALL soc_contract_evect(evects_a, a_coeff)
      CALL soc_contract_evect(evects_b, b_coeff)

      ALLOCATE (soc_env%orb_soc(3))
      DO i = 1, 3
         ALLOCATE (soc_env%orb_soc(i)%matrix)
      END DO

      ! Make soc_atom_env for H^SOC calculations
      ! Compute the contraction coefficients for spherical orbitals
      CALL soc_atom_create(soc_atom_env)
      IF (do_os) THEN
         soc_atom_env%nspins = 2
      ELSE
         soc_atom_env%nspins = 1
      END IF

      nkind = SIZE(qs_kind_set)
      ALLOCATE (soc_atom_env%grid_atom_set(nkind))
      ALLOCATE (soc_atom_env%harmonics_atom_set(nkind))
      ALLOCATE (soc_atom_env%orb_sphi_so(nkind))

      CALL get_qs_kind_set(qs_kind_set, all_potential_present=all_potential_present)
      IF (.NOT. all_potential_present) THEN
         CALL atomic_pseudopotential(soc_atom_env, qs_env)
      END IF

      !  IF (all_potential_present) THEN
      !! Prepare the atomic grid we need to calculate the SOC Operator in an Atomic basis
      !! copied from init_xas_atom_grid_harmo for convenience
      CALL init_atom_grid(soc_atom_env, grid_info, qs_env)

      DO ikind = 1, nkind
         NULLIFY (soc_atom_env%orb_sphi_so(ikind)%array)
         CALL compute_sphi_so(ikind, "ORB", soc_atom_env%orb_sphi_so(ikind)%array, qs_env)
      END DO !ikind
      !  ENDIF

      DEALLOCATE (grid_info)

      CALL timestop(handle)

   END SUBROUTINE inititialize_soc

! **************************************************************************************************
!> \brief Initializes the atomic grids and harmonics for the atomic calculations
!> \param soc_atom_env ...
!> \param grid_info ...
!> \param qs_env ...
!> \note Copied and modified from init_xas_atom_grid_harmo
! **************************************************************************************************
   SUBROUTINE init_atom_grid(soc_atom_env, grid_info, qs_env)

      TYPE(soc_atom_env_type)                            :: soc_atom_env
      CHARACTER(len=default_string_length), &
         ALLOCATABLE, DIMENSION(:, :)                    :: grid_info
      TYPE(qs_environment_type)                          :: qs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'init_atom_grid'

      CHARACTER(LEN=default_string_length)               :: kind_name
      INTEGER :: handle, igrid, ikind, il, iso, iso1, iso2, l1, l1l2, l2, la, lc1, lc2, lcleb, ll, &
         llmax, lp, m1, m2, max_s_harm, max_s_set, maxl, maxlgto, maxs, mm, mp, na, nr, &
         quadrature, stat
      REAL(dp)                                           :: kind_radius
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: rga
      REAL(dp), DIMENSION(:, :, :), POINTER              :: my_CG
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(grid_atom_type), POINTER                      :: grid_atom
      TYPE(gto_basis_set_type), POINTER                  :: tmp_basis
      TYPE(harmonics_atom_type), POINTER                 :: harmonics
      TYPE(qs_control_type), POINTER                     :: qs_control
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      NULLIFY (my_CG, qs_kind_set, dft_control, qs_control)
      NULLIFY (grid_atom, harmonics, tmp_basis)

      !!  Initialization of some integer for the CG coeff generation
      CALL get_qs_env(qs_env, &
                      qs_kind_set=qs_kind_set, &
                      dft_control=dft_control)

      CALL get_qs_kind_set(qs_kind_set, maxlgto=maxlgto, basis_type="ORB")
      qs_control => dft_control%qs_control

      !!To be sure:: Is the basis grid present?
      IF (maxlgto < 0) THEN
         CPABORT("tmp_basis is not Present!")
      END IF

      !maximum expansion
      llmax = 2*maxlgto
      max_s_harm = nsoset(llmax)
      max_s_set = nsoset(maxlgto)
      lcleb = llmax

      !  Allocate and compute the CG coeffs (copied from init_rho_atom)
      CALL clebsch_gordon_init(lcleb)
      CALL reallocate(my_CG, 1, max_s_set, 1, max_s_set, 1, max_s_harm)

      ALLOCATE (rga(lcleb, 2))
      DO lc1 = 0, maxlgto
         DO iso1 = nsoset(lc1 - 1) + 1, nsoset(lc1)
            l1 = indso(1, iso1)
            m1 = indso(2, iso1)
            DO lc2 = 0, maxlgto
               DO iso2 = nsoset(lc2 - 1) + 1, nsoset(lc2)
                  l2 = indso(1, iso2)
                  m2 = indso(2, iso2)
                  CALL clebsch_gordon(l1, m1, l2, m2, rga)
                  IF (l1 + l2 > llmax) THEN
                     l1l2 = llmax
                  ELSE
                     l1l2 = l1 + l2
                  END IF
                  mp = m1 + m2
                  mm = m1 - m2
                  IF (m1*m2 < 0 .OR. (m1*m2 == 0 .AND. (m1 < 0 .OR. m2 < 0))) THEN
                     mp = -ABS(mp)
                     mm = -ABS(mm)
                  ELSE
                     mp = ABS(mp)
                     mm = ABS(mm)
                  END IF
                  DO lp = MOD(l1 + l2, 2), l1l2, 2
                     il = lp/2 + 1
                     IF (ABS(mp) <= lp) THEN
                     IF (mp >= 0) THEN
                        iso = nsoset(lp - 1) + lp + 1 + mp
                     ELSE
                        iso = nsoset(lp - 1) + lp + 1 - ABS(mp)
                     END IF
                     my_CG(iso1, iso2, iso) = rga(il, 1)
                     END IF
                     IF (mp /= mm .AND. ABS(mm) <= lp) THEN
                     IF (mm >= 0) THEN
                        iso = nsoset(lp - 1) + lp + 1 + mm
                     ELSE
                        iso = nsoset(lp - 1) + lp + 1 - ABS(mm)
                     END IF
                     my_CG(iso1, iso2, iso) = rga(il, 2)
                     END IF
                  END DO
               END DO ! iso2
            END DO ! lc2
         END DO ! iso1
      END DO ! lc1
      DEALLOCATE (rga)
      CALL clebsch_gordon_deallocate()

      !  Create the Lebedev grids and compute the spherical harmonics
      CALL init_lebedev_grids()
      quadrature = qs_control%gapw_control%quadrature

      DO ikind = 1, SIZE(soc_atom_env%grid_atom_set)

         ! Allocate the grid and the harmonics for this kind
         NULLIFY (soc_atom_env%grid_atom_set(ikind)%grid_atom)
         NULLIFY (soc_atom_env%harmonics_atom_set(ikind)%harmonics_atom)
         CALL allocate_grid_atom(soc_atom_env%grid_atom_set(ikind)%grid_atom)
         CALL allocate_harmonics_atom( &
            soc_atom_env%harmonics_atom_set(ikind)%harmonics_atom)

         NULLIFY (grid_atom, harmonics)
         grid_atom => soc_atom_env%grid_atom_set(ikind)%grid_atom
         harmonics => soc_atom_env%harmonics_atom_set(ikind)%harmonics_atom

         ! Initialize some integers
         CALL get_qs_kind(qs_kind_set(ikind), &
                          ngrid_rad=nr, &
                          ngrid_ang=na, &
                          name=kind_name)

         ! take the grid dimension given as input, if none, take the GAPW ones above
         DO igrid = 1, SIZE(grid_info, 1)
            IF (grid_info(igrid, 1) == kind_name) THEN
               READ (grid_info(igrid, 2), *, iostat=stat) na
               IF (stat .NE. 0) CPABORT("The 'na' value for the GRID keyword must be an integer")
               READ (grid_info(igrid, 3), *, iostat=stat) nr
               IF (stat .NE. 0) CPABORT("The 'nr' value for the GRID keyword must be an integer")
               EXIT
            END IF
         END DO !igrid

         ll = get_number_of_lebedev_grid(n=na)
         na = lebedev_grid(ll)%n
         la = lebedev_grid(ll)%l
         grid_atom%ng_sphere = na
         grid_atom%nr = nr

         !Create the harmonics with the ORB-Basis
         CALL get_qs_kind(qs_kind_set(ikind), &
                          basis_set=tmp_basis, &
                          basis_type="ORB")
         CALL get_gto_basis_set(gto_basis_set=tmp_basis, &
                                maxl=maxl, &
                                kind_radius=kind_radius)

         CALL create_grid_atom(grid_atom, nr, na, llmax, ll, quadrature)
         CALL truncate_radial_grid(grid_atom, kind_radius)

         maxs = nsoset(maxl)
         CALL create_harmonics_atom(harmonics, &
                                    my_CG, na, &
                                    llmax, maxs, &
                                    max_s_harm, ll, &
                                    grid_atom%wa, &
                                    grid_atom%azi, &
                                    grid_atom%pol)
         CALL get_maxl_CG(harmonics, tmp_basis, llmax, max_s_harm)
      END DO !ikind

      CALL deallocate_lebedev_grids()
      DEALLOCATE (my_CG)

      CALL timestop(handle)

   END SUBROUTINE init_atom_grid

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param dbcsr_soc_package ...
!> \param soc_env ...
!> \param soc_evecs_cfm ...
!> \param eps_filter ...
!> \param dipole_form ...
!> \param gs_mos ...
! **************************************************************************************************
   SUBROUTINE soc_dipol(qs_env, dbcsr_soc_package, soc_env, &
                        soc_evecs_cfm, eps_filter, dipole_form, gs_mos)
      TYPE(qs_environment_type)                          :: qs_env
      TYPE(dbcsr_soc_package_type)                       :: dbcsr_soc_package
      TYPE(soc_env_type), TARGET                         :: soc_env
      TYPE(cp_cfm_type), INTENT(in)                      :: soc_evecs_cfm
      REAL(dp), INTENT(IN)                               :: eps_filter
      INTEGER                                            :: dipole_form
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         POINTER                                         :: gs_mos

      CHARACTER(len=*), PARAMETER                        :: routineN = 'soc_dipol'

      COMPLEX(dp), ALLOCATABLE, DIMENSION(:, :)          :: transdip
      INTEGER                                            :: handle, i, nosc, ntot
      REAL(dp), DIMENSION(:), POINTER                    :: osc_str, soc_evals
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_cfm_type)                                  :: dip_cfm, work1_cfm, work2_cfm
      TYPE(cp_fm_struct_type), POINTER                   :: dip_struct, full_struct
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: amew_dip
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      NULLIFY (para_env, blacs_env, dip_struct, full_struct, osc_str)
      NULLIFY (soc_evals)

      CALL get_qs_env(qs_env, para_env=para_env, blacs_env=blacs_env)

      soc_evals => soc_env%soc_evals
      nosc = SIZE(soc_evals)
      ntot = nosc + 1
      ALLOCATE (soc_env%soc_osc(nosc))
      osc_str => soc_env%soc_osc
      osc_str(:) = 0.0_dp

      !get some work arrays/matrix
      CALL cp_fm_struct_create(dip_struct, &
                               context=blacs_env, para_env=para_env, &
                               nrow_global=ntot, ncol_global=1)

      CALL cp_cfm_get_info(soc_evecs_cfm, matrix_struct=full_struct)
      CALL cp_cfm_create(dip_cfm, dip_struct)
      CALL cp_cfm_create(work1_cfm, full_struct)
      CALL cp_cfm_create(work2_cfm, full_struct)

      ALLOCATE (transdip(ntot, 1))

      CALL get_rcs_amew_op(amew_dip, soc_env, dbcsr_soc_package, &
                           eps_filter, qs_env, gs_mos)

      DO i = 1, 3 !cartesian coord x, y, z

         !Convert the real dipole into the cfm format for calculations
         CALL cp_fm_to_cfm(msourcer=amew_dip(i), mtarget=work1_cfm)

         !compute amew_coeffs^dagger * amew_dip * amew_gs to get the transition moments
         CALL parallel_gemm('C', 'N', ntot, &
                            ntot, ntot, &
                            (1.0_dp, 0.0_dp), &
                            soc_evecs_cfm, &
                            work1_cfm, &
                            (0.0_dp, 0.0_dp), &
                            work2_cfm)
         CALL parallel_gemm('N', 'N', ntot, &
                            1, ntot, &
                            (1.0_dp, 0.0_dp), &
                            work2_cfm, &
                            soc_evecs_cfm, &
                            (0.0_dp, 0.0_dp), &
                            dip_cfm)

         CALL cp_cfm_get_submatrix(dip_cfm, transdip)

         !transition dipoles are real numbers
         osc_str(:) = osc_str(:) + REAL(transdip(2:ntot, 1))**2 + AIMAG(transdip(2:ntot, 1))**2

      END DO !i

      !multiply with appropriate prefac depending in the rep
      IF (dipole_form == tddfpt_dipole_length) THEN
         osc_str(:) = 2.0_dp/3.0_dp*soc_evals(:)*osc_str(:)
      ELSE
         osc_str(:) = 2.0_dp/3.0_dp/soc_evals(:)*osc_str(:)
      END IF

      !clean-up
      CALL cp_fm_struct_release(dip_struct)
      CALL cp_cfm_release(work1_cfm)
      CALL cp_cfm_release(work2_cfm)
      CALL cp_cfm_release(dip_cfm)
      DO i = 1, 3
         CALL cp_fm_release(amew_dip(i))
      END DO
      DEALLOCATE (amew_dip, transdip)

      CALL timestop(handle)

   END SUBROUTINE soc_dipol

!**********************************************************************************
!> \brief: This will create the Dipole operator within the amew basis
!> \param amew_op Output dipole operator
!> \param soc_env ...
!> \param dbcsr_soc_package Some work matrices
!> \param eps_filter ...
!> \param qs_env ...
!> \param gs_mos ...
! **************************************************************************************************
   SUBROUTINE get_rcs_amew_op(amew_op, soc_env, dbcsr_soc_package, eps_filter, qs_env, gs_mos)

      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:), &
         INTENT(OUT)                                     :: amew_op
      TYPE(soc_env_type), TARGET                         :: soc_env
      TYPE(dbcsr_soc_package_type)                       :: dbcsr_soc_package
      REAL(dp), INTENT(IN)                               :: eps_filter
      TYPE(qs_environment_type)                          :: qs_env
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         POINTER                                         :: gs_mos

      CHARACTER(len=*), PARAMETER                        :: routineN = 'get_rcs_amew_op'

      INTEGER                                            :: dim_op, handle, i, isg, nactive, nao, &
                                                            nex, nsg, ntot, ntp
      REAL(dp)                                           :: op, sqrt2
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: diag, gs_diag, gsgs_op
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: mo_op, sggs_block
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: full_struct, gsgs_struct, sggs_struct, &
                                                            std_struct, tmp_struct
      TYPE(cp_fm_type)                                   :: gs_fm, sggs_fm, tmp_fm, vec_op, work_fm
      TYPE(cp_fm_type), POINTER                          :: gs_coeffs
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: ao_op, matrix_s
      TYPE(dbcsr_type), POINTER                          :: ao_op_i, dbcsr_ovlp, dbcsr_prod, &
                                                            dbcsr_sg, dbcsr_tmp, dbcsr_tp, &
                                                            dbcsr_work
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      NULLIFY (mos, gs_coeffs)
      NULLIFY (matrix_s)
      NULLIFY (para_env, blacs_env)
      NULLIFY (full_struct, gsgs_struct, std_struct, tmp_struct, sggs_struct)
      NULLIFY (ao_op_i, ao_op)
      NULLIFY (dbcsr_tp, dbcsr_sg, dbcsr_ovlp, dbcsr_work, dbcsr_tmp, dbcsr_prod)

      ! Initialization
      ao_op => soc_env%dipmat
      nsg = SIZE(soc_env%evals_a)
      ntp = nsg; nex = nsg !all the same by construction, keep them separate for clarity
      ntot = 1 + nsg + 3*ntp

      CALL get_qs_env(qs_env, matrix_s=matrix_s, mos=mos, para_env=para_env, blacs_env=blacs_env)
      sqrt2 = SQRT(2.0_dp)
      dim_op = SIZE(ao_op)

      dbcsr_sg => dbcsr_soc_package%dbcsr_sg
      dbcsr_tp => dbcsr_soc_package%dbcsr_tp
      dbcsr_work => dbcsr_soc_package%dbcsr_work
      dbcsr_prod => dbcsr_soc_package%dbcsr_prod
      dbcsr_ovlp => dbcsr_soc_package%dbcsr_ovlp
      dbcsr_tmp => dbcsr_soc_package%dbcsr_tmp

      !  Create the amew_op matrix
      CALL cp_fm_struct_create(full_struct, &
                               context=blacs_env, para_env=para_env, &
                               nrow_global=ntot, ncol_global=ntot)
      ALLOCATE (amew_op(dim_op))
      DO i = 1, dim_op
         CALL cp_fm_create(amew_op(i), full_struct)
      END DO !i

      !  Deal with the GS-GS contribution <0|0> = 2*sum_j <phi_j|op|phi_j>
      CALL get_mo_set(mos(1), nao=nao, homo=nactive)
      gs_coeffs => gs_mos(1)%mos_occ
      CALL cp_fm_get_info(gs_coeffs, matrix_struct=std_struct)

      CALL cp_fm_struct_create(gsgs_struct, &
                               context=blacs_env, para_env=para_env, &
                               nrow_global=nactive, ncol_global=nactive)

      CALL cp_fm_create(gs_fm, gsgs_struct) !nactive nactive
      CALL cp_fm_create(work_fm, std_struct) !nao nactive

      ALLOCATE (gsgs_op(dim_op))
      ALLOCATE (gs_diag(nactive))

      DO i = 1, dim_op

         ao_op_i => ao_op(i)%matrix

         CALL cp_dbcsr_sm_fm_multiply(ao_op_i, gs_coeffs, work_fm, ncol=nactive)
         CALL parallel_gemm('T', 'N', nactive, nactive, nao, &
                            1.0_dp, gs_coeffs, work_fm, 0.0_dp, gs_fm)
         CALL cp_fm_get_diag(gs_fm, gs_diag)
         gsgs_op(i) = 2.0_dp*SUM(gs_diag)

      END DO !i
      DEALLOCATE (gs_diag)

      CALL cp_fm_release(work_fm)

      !  Create the work and helper fms
      CALL cp_fm_create(vec_op, std_struct) !nao nactive
      CALL cp_fm_struct_release(gsgs_struct)

      CALL cp_fm_struct_create(tmp_struct, &
                               context=blacs_env, para_env=para_env, &
                               nrow_global=nex, ncol_global=nex)
      CALL cp_fm_struct_create(sggs_struct, &
                               context=blacs_env, para_env=para_env, &
                               nrow_global=nactive*nsg, ncol_global=nactive)

      CALL cp_fm_create(tmp_fm, tmp_struct)
      CALL cp_fm_create(work_fm, full_struct)
      CALL cp_fm_create(sggs_fm, sggs_struct)

      ALLOCATE (diag(nactive))
      ALLOCATE (mo_op(nactive, nactive))
      ALLOCATE (sggs_block(nactive, nactive))

      ! Iterate over the dimensions of the operator
      ! Note: operator matrices are asusmed symmetric, can only do upper half
      DO i = 1, dim_op

         ao_op_i => ao_op(i)%matrix

         ! The GS-GS contribution
         CALL cp_fm_set_element(amew_op(i), 1, 1, gsgs_op(i))

         ! Compute the operator in the MO basis
         CALL cp_dbcsr_sm_fm_multiply(ao_op_i, gs_coeffs, vec_op, ncol=nactive)
         CALL cp_fm_get_submatrix(gs_fm, mo_op)  ! instead of prod_fm

         ! Compute the ground-state/singlet components. ao_op*gs_coeffs already stored in vec_op
         CALL parallel_gemm('T', 'N', nactive*nsg, nactive, nao, &
                            1.0_dp, soc_env%a_coeff, vec_op, 0.0_dp, sggs_fm)
         DO isg = 1, nsg
            CALL cp_fm_get_submatrix(fm=sggs_fm, &
                                     target_m=sggs_block, &
                                     start_row=(isg - 1)*nactive + 1, &
                                     start_col=1, &
                                     n_rows=nactive, &
                                     n_cols=nactive)
            diag(:) = get_diag(sggs_block)
            op = sqrt2*SUM(diag)
            CALL cp_fm_set_element(amew_op(i), 1, 1 + isg, op)
         END DO !isg

         ! do the singlet-singlet components
         !start with the overlap
         CALL dbcsr_multiply('N', 'N', 1.0_dp, &
                             matrix_s(1)%matrix, dbcsr_sg, &
                             0.0_dp, dbcsr_work, &
                             filter_eps=eps_filter)
         CALL dbcsr_multiply('T', 'N', 1.0_dp, &
                             dbcsr_sg, dbcsr_work, &
                             0.0_dp, dbcsr_ovlp, &
                             filter_eps=eps_filter)

         !then the operator in the LR orbital basis
         CALL dbcsr_multiply('N', 'N', 1.0_dp, &
                             ao_op_i, dbcsr_sg, &
                             0.0_dp, dbcsr_work, &
                             filter_eps=eps_filter)
         CALL dbcsr_multiply('T', 'N', 1.0_dp, &
                             dbcsr_sg, dbcsr_work, &
                             0.0_dp, dbcsr_prod, &
                             filter_eps=eps_filter)

         !use the soc routine, it is compatible
         CALL rcs_amew_soc_elements(dbcsr_tmp, &
                                    dbcsr_prod, &
                                    dbcsr_ovlp, &
                                    mo_op, &
                                    pref_trace=-1.0_dp, &
                                    pref_overall=1.0_dp, &
                                    pref_diags=gsgs_op(i), &
                                    symmetric=.TRUE.)  !!FOR DEBUG

         CALL copy_dbcsr_to_fm(dbcsr_tmp, tmp_fm)
         CALL cp_fm_to_fm_submat(msource=tmp_fm, &
                                 mtarget=amew_op(i), &
                                 nrow=nex, &
                                 ncol=nex, &
                                 s_firstrow=1, &
                                 s_firstcol=1, &
                                 t_firstrow=2, &
                                 t_firstcol=2)

         ! compute the triplet-triplet components
         !the overlap
         CALL dbcsr_multiply('N', 'N', 1.0_dp, &
                             matrix_s(1)%matrix, &
                             dbcsr_tp, 0.0_dp, &
                             dbcsr_work, &
                             filter_eps=eps_filter)
         CALL dbcsr_multiply('T', 'N', 1.0_dp, &
                             dbcsr_tp, dbcsr_work, &
                             0.0_dp, dbcsr_ovlp, &
                             filter_eps=eps_filter)

         !the operator in the LR orbital basis
         CALL dbcsr_multiply('N', 'N', 1.0_dp, &
                             ao_op_i, dbcsr_sg, &
                             0.0_dp, dbcsr_work, &
                             filter_eps=eps_filter)
         CALL dbcsr_multiply('T', 'N', 1.0_dp, &
                             dbcsr_sg, dbcsr_work, &
                             0.0_dp, dbcsr_prod, &
                             filter_eps=eps_filter)

         CALL rcs_amew_soc_elements(dbcsr_tmp, &
                                    dbcsr_prod, &
                                    dbcsr_ovlp, &
                                    mo_op, &
                                    pref_trace=-1.0_dp, &
                                    pref_overall=1.0_dp, &
                                    pref_diags=gsgs_op(i), &
                                    symmetric=.TRUE.)  !!FOR DEBUG

         CALL copy_dbcsr_to_fm(dbcsr_tmp, tmp_fm)
         !<T^-1|op|T^-1>
         CALL cp_fm_to_fm_submat(msource=tmp_fm, &
                                 mtarget=amew_op(i), &
                                 nrow=nex, &
                                 ncol=nex, &
                                 s_firstrow=1, &
                                 s_firstcol=1, &
                                 t_firstrow=1 + nsg + 1, &
                                 t_firstcol=1 + nsg + 1)

         !<T^0|op|T^0>
         CALL cp_fm_to_fm_submat(msource=tmp_fm, &
                                 mtarget=amew_op(i), &
                                 nrow=nex, &
                                 ncol=nex, &
                                 s_firstrow=1, &
                                 s_firstcol=1, &
                                 t_firstrow=1 + nsg + ntp + 1, &
                                 t_firstcol=1 + nsg + ntp + 1)
         !<T^+1|op|T^+1>
         CALL cp_fm_to_fm_submat(msource=tmp_fm, &
                                 mtarget=amew_op(i), &
                                 nrow=nex, &
                                 ncol=nex, &
                                 s_firstrow=1, &
                                 s_firstcol=1, &
                                 t_firstrow=1 + nsg + 2*ntp + 1, &
                                 t_firstcol=1 + nsg + 2*ntp + 1)

         ! Symmetrize the matrix (only upper triangle built)
         CALL cp_fm_upper_to_full(amew_op(i), work_fm)

      END DO !i

      ! Clean-up
      CALL cp_fm_release(gs_fm)
      CALL cp_fm_release(work_fm)
      CALL cp_fm_release(tmp_fm)
      CALL cp_fm_release(vec_op)
      CALL cp_fm_release(sggs_fm)
      CALL cp_fm_struct_release(full_struct)
      CALL cp_fm_struct_release(tmp_struct)
      CALL cp_fm_struct_release(sggs_struct)

      DEALLOCATE (sggs_block)
      DEALLOCATE (diag, gsgs_op, mo_op)
      CALL timestop(handle)

   END SUBROUTINE get_rcs_amew_op

END MODULE qs_tddfpt2_soc
