!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

MODULE qs_fb_env_methods

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_p_type,&
                                             gto_basis_set_type
  USE cell_types,                      ONLY: cell_type
  USE cp_blacs_env,                    ONLY: cp_blacs_env_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_allocate_matrix_set, cp_dbcsr_col_block_sizes, &
       cp_dbcsr_create, cp_dbcsr_deallocate_matrix_set, &
       cp_dbcsr_desymmetrize, cp_dbcsr_distribution, cp_dbcsr_finalize, &
       cp_dbcsr_get_info, cp_dbcsr_get_matrix_type, cp_dbcsr_init, &
       cp_dbcsr_iterator, cp_dbcsr_iterator_blocks_left, &
       cp_dbcsr_iterator_next_block, cp_dbcsr_iterator_start, &
       cp_dbcsr_iterator_stop, cp_dbcsr_multiply, cp_dbcsr_nfullcols_total, &
       cp_dbcsr_nfullrows_total, cp_dbcsr_p_type, cp_dbcsr_release, &
       cp_dbcsr_reserve_blocks, cp_dbcsr_row_block_sizes, cp_dbcsr_set, &
       cp_dbcsr_type, dbcsr_type_no_symmetry
  USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm
  USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm,&
                                             cp_fm_symm,&
                                             cp_fm_triangular_invert,&
                                             cp_fm_triangular_multiply,&
                                             cp_fm_upper_to_full
  USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose
  USE cp_fm_diag,                      ONLY: choose_eigv_solver,&
                                             cp_fm_power
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_release,&
                                             cp_fm_set_all,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_gemm_interface,               ONLY: cp_gemm
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_units,                        ONLY: cp_unit_from_cp2k
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE message_passing,                 ONLY: mp_max
  USE orbital_pointers,                ONLY: nco,&
                                             ncoset
  USE particle_types,                  ONLY: particle_type
  USE qs_diis,                         ONLY: qs_diis_b_step
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_fb_atomic_halo_types,         ONLY: &
       fb_atomic_halo_build_halo_atoms, fb_atomic_halo_cost, &
       fb_atomic_halo_create, fb_atomic_halo_list_create, &
       fb_atomic_halo_list_nullify, fb_atomic_halo_list_obj, &
       fb_atomic_halo_list_release, fb_atomic_halo_list_set, &
       fb_atomic_halo_list_write_info, fb_atomic_halo_nelectrons_estimate_Z, &
       fb_atomic_halo_nullify, fb_atomic_halo_obj, fb_atomic_halo_set, &
       fb_atomic_halo_sort, fb_build_pair_radii
  USE qs_fb_env_types,                 ONLY: fb_env_get,&
                                             fb_env_has_data,&
                                             fb_env_obj,&
                                             fb_env_set
  USE qs_fb_filter_matrix_methods,     ONLY: fb_fltrmat_build,&
                                             fb_fltrmat_build_2
  USE qs_fb_trial_fns_types,           ONLY: fb_trial_fns_create,&
                                             fb_trial_fns_nullify,&
                                             fb_trial_fns_obj,&
                                             fb_trial_fns_release,&
                                             fb_trial_fns_set
  USE qs_integral_utils,               ONLY: basis_set_list_setup
  USE qs_kind_types,                   ONLY: get_qs_kind,&
                                             qs_kind_type
  USE qs_matrix_pools,                 ONLY: mpools_create,&
                                             mpools_rebuild_fm_pools,&
                                             mpools_release,&
                                             qs_matrix_pools_type
  USE qs_mo_methods,                   ONLY: calculate_density_matrix
  USE qs_mo_occupation,                ONLY: set_mo_occupation
  USE qs_mo_types,                     ONLY: allocate_mo_set,&
                                             deallocate_mo_set,&
                                             get_mo_set,&
                                             init_mo_set,&
                                             mo_set_p_type,&
                                             mo_set_type,&
                                             set_mo_set
  USE qs_scf_types,                    ONLY: qs_scf_env_type
  USE scf_control_types,               ONLY: scf_control_type
  USE string_utilities,                ONLY: compress,&
                                             uppercase
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_fb_env_methods'

  PUBLIC :: fb_env_do_diag,&
            fb_env_read_input,&
            fb_env_build_rcut_auto,&
            fb_env_build_atomic_halos,&
            fb_env_write_info


CONTAINS


! *****************************************************************************
!> \brief Do filtered matrix method diagonalisation
!> \param fb_env : the filter matrix environment
!> \param qs_env : quickstep environment
!> \param matrix_ks : DBCSR system (unfiltered) input KS matrix
!> \param matrix_s  : DBCSR system (unfiltered) input overlap matrix
!> \param scf_section : SCF input section
!> \param diis_step : whether we are doing a DIIS step
!> \param error : cp2k error container
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_env_do_diag(fb_env, &
                            qs_env, &
                            matrix_ks, &
                            matrix_s, &
                            scf_section, &
                            diis_step, &
                            error)
    TYPE(fb_env_obj), INTENT(INOUT)          :: fb_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_s
    TYPE(section_vals_type), POINTER         :: scf_section
    LOGICAL, INTENT(INOUT)                   :: diis_step
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'fb_env_do_diag', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=2)                         :: spin_string
    CHARACTER(len=default_string_length)     :: name
    INTEGER :: filtered_nfullrowsORcols_total, handle, homo_filtered, ispin, &
      lfomo_filtered, my_nmo, ndep, nelectron, nmo, nmo_filtered, nspin, &
      original_nfullrowsORcols_total, stat
    INTEGER, DIMENSION(:), POINTER           :: filtered_rowORcol_block_sizes,&
                                                original_rowORcol_block_sizes
    LOGICAL                                  :: collective_com, failure
    REAL(kind=dp) :: diis_error, eps_diis, eps_eigval, fermi_level, &
      filter_temp, flexible_electron_count, KTS_filtered, maxocc, mu_filtered
    REAL(KIND=dp), DIMENSION(:), POINTER     :: eigenvalues, &
                                                eigenvalues_filtered, occ, &
                                                occ_filtered
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    TYPE(cp_dbcsr_p_type)                    :: matrix_ks_desymm, &
                                                matrix_s_desymm
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_filter
    TYPE(cp_dbcsr_type)                      :: matrix_filtered_ks, &
                                                matrix_filtered_s, matrix_tmp
    TYPE(cp_dbcsr_type), POINTER             :: matrix_filtered_p
    TYPE(cp_fm_struct_type), POINTER         :: filter_fm_struct, fm_struct
    TYPE(cp_fm_type), POINTER :: fm_matrix_filter, fm_matrix_filtered_ks, &
      fm_matrix_filtered_s, fm_matrix_ortho, fm_matrix_work, mo_coeff, &
      mo_coeff_filtered
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(fb_atomic_halo_list_obj)            :: atomic_halos
    TYPE(fb_trial_fns_obj)                   :: trial_fns
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos, mos_filtered
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_matrix_pools_type), POINTER      :: my_mpools
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(scf_control_type), POINTER          :: scf_control

! TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER :: sab_orb

    CALL timeset(routineN, handle)

    failure = .FALSE.

    NULLIFY(scf_env, scf_control, para_env, blacs_env, particle_set)
    NULLIFY(eigenvalues, eigenvalues_filtered, occ, occ_filtered)
    NULLIFY(mos, mos_filtered)
    NULLIFY(my_mpools)
    NULLIFY(matrix_filter, matrix_filtered_p)
    NULLIFY(fm_struct, filter_fm_struct)
    NULLIFY(fm_matrix_filter, fm_matrix_filtered_s, &
            fm_matrix_filtered_ks, fm_matrix_work, &
            fm_matrix_ortho, mo_coeff_filtered, mo_coeff)
    ! NULLIFY(sab_orb)
    CALL fb_atomic_halo_list_nullify(atomic_halos)
    CALL fb_trial_fns_nullify(trial_fns)
    NULLIFY(original_rowORcol_block_sizes,filtered_rowORcol_block_sizes)

    ! get qs_env information
    CALL get_qs_env(qs_env=qs_env, &
                    scf_env=scf_env, &
                    scf_control=scf_control, &
                    para_env=para_env, &
                    blacs_env=blacs_env, &
                    particle_set=particle_set, &
                    mos=mos, &
                    error=error)

    nspin = SIZE(matrix_ks)

    ! ----------------------------------------------------------------------
    ! DIIS step - based on non-filtered matrices and MOs
    ! ----------------------------------------------------------------------

    DO ispin = 1, nspin
       CALL copy_dbcsr_to_fm(matrix_ks(ispin)%matrix, &
                             scf_env%scf_work1(ispin)%matrix, &
                             error=error)
    END DO

    eps_diis = scf_control%eps_diis
    eps_eigval = scf_control%eps_eigval

    IF (scf_env%iter_count > 1 .AND. .NOT. scf_env%skip_diis) THEN
       CALL qs_diis_b_step(scf_env%scf_diis_buffer, mos, scf_env%scf_work1, &
                           scf_env%scf_work2, scf_env%iter_delta, &
                           diis_error, diis_step, eps_diis, scf_control%nmixing, &
                           s_matrix=matrix_s, scf_section=scf_section, &
                           error=error)
    ELSE
       diis_step = .FALSE.
    END IF

    IF (diis_step) THEN
       scf_env%iter_param = diis_error
       scf_env%iter_method = "DIIS/Filter"
    ELSE
       IF (scf_env%mixing_method == 0) THEN
          scf_env%iter_method = "NoMix/Filter"
       ELSE IF(scf_env%mixing_method == 1) THEN
          scf_env%iter_param = scf_env%p_mix_alpha
          scf_env%iter_method = "P_Mix/Filter"
       ELSE IF(scf_env%mixing_method > 1) THEN
          scf_env%iter_param = scf_env%mixing_store%alpha
          scf_env%iter_method = TRIM(scf_env%mixing_store%iter_method)//"/Filter"
       END IF
    END IF

    ! ----------------------------------------------------------------------
    ! Construct Filter Matrix
    ! ----------------------------------------------------------------------

    CALL fb_env_get(fb_env=fb_env, &
                    filter_temperature=filter_temp, &
                    atomic_halos=atomic_halos, &
                    error=error)

    ! construct trial functions
    CALL get_mo_set(mo_set=mos(1)%mo_set, maxocc=maxocc)
    CALL fb_env_build_trial_fns_auto(fb_env, qs_env, maxocc, error)
    CALL fb_env_get(fb_env=fb_env, &
                    trial_fns=trial_fns, &
                    error=error)

    ! allocate filter matrix (matrix_filter(ispin)%matrix are
    ! nullified by cp_dbcsr_allocate_matrix_set)
    CALL cp_dbcsr_allocate_matrix_set(matrix_filter, nspin, error)
    DO ispin = 1, nspin
       ! get system-wide fermi energy and occupancy, we use this to
       ! define the filter function used for the filter matrix
       CALL get_mo_set(mo_set=mos(ispin)%mo_set, &
                       mu=fermi_level, &
                       maxocc=maxocc)

       ! ----------------------------------------------------------------------
       ! desymmetrise matrix first to avoid the complex and hard to
       ! understand DBCSR symmetric layout first.
       !
       !TODO:LT:2015/02/04: This is inefficient, need to use symmetry!
       !
       NULLIFY(matrix_ks_desymm%matrix)
       ALLOCATE(matrix_ks_desymm%matrix)
       CALL cp_dbcsr_init(matrix_ks_desymm%matrix, error)
       CALL cp_dbcsr_desymmetrize(matrix_ks(ispin)%matrix, matrix_ks_desymm%matrix, error)
       NULLIFY(matrix_s_desymm%matrix)
       ALLOCATE(matrix_s_desymm%matrix)
       CALL cp_dbcsr_init(matrix_s_desymm%matrix, error)
       CALL cp_dbcsr_desymmetrize(matrix_s(1)%matrix, matrix_s_desymm%matrix, error)
       ! ----------------------------------------------------------------------

       ! get filter matrix name
       WRITE (spin_string, FMT="(I1)") ispin
       name = TRIM("FILTER MATRIX SPIN "//spin_string)
       CALL compress(name)
       CALL uppercase(name)
       ! calculate filter matrix (matrix_s(1) is the overlap, the rest
       ! in the array are its derivatives)
!BEG:ORIG:LT:2015/01/2
       ! CALL fb_fltrmat_build(H_mat=matrix_ks(ispin)%matrix, &
       !                       S_mat=matrix_s(1)%matrix, &
       !                       atomic_halos=atomic_halos, &
       !                       trial_fns=trial_fns, &
       !                       para_env=para_env, &
       !                       particle_set=particle_set, &
       !                       fermi_level=fermi_level, &
       !                       filter_temp=filter_temp, &
       !                       name=name, &
       !                       filter_mat=matrix_filter(ispin)%matrix, &
       !                       error=error)
!END:ORIG:LT:2015/01/2
!BEG:REPL:LT:2015/01/2
       CALL fb_env_get(fb_env=fb_env, &
                       collective_com=collective_com, &
                       error=error)
       IF (collective_com) THEN
          CALL fb_fltrmat_build_2(H_mat=matrix_ks_desymm%matrix, &
                                  S_mat=matrix_s_desymm%matrix, &
                                  atomic_halos=atomic_halos, &
                                  trial_fns=trial_fns, &
                                  para_env=para_env, &
                                  particle_set=particle_set, &
                                  fermi_level=fermi_level, &
                                  filter_temp=filter_temp, &
                                  name=name, &
                                  filter_mat=matrix_filter(ispin)%matrix, &
                                  error=error)
       ELSE
          CALL fb_fltrmat_build(H_mat=matrix_ks_desymm%matrix, &
                                S_mat=matrix_s_desymm%matrix, &
                                atomic_halos=atomic_halos, &
                                trial_fns=trial_fns, &
                                para_env=para_env, &
                                particle_set=particle_set, &
                                fermi_level=fermi_level, &
                                filter_temp=filter_temp, &
                                name=name, &
                                filter_mat=matrix_filter(ispin)%matrix, &
                                error=error)
       END IF
!END:REPL:LT:2015/01/2

       ! ----------------------------------------------------------------------
       ! Deallocate the desymmetrised matrices
       CALL cp_dbcsr_release(matrix_ks_desymm%matrix, error)
       CALL cp_dbcsr_release(matrix_s_desymm%matrix, error)
       DEALLOCATE(matrix_ks_desymm%matrix)
       DEALLOCATE(matrix_s_desymm%matrix)
       ! ----------------------------------------------------------------------

    END DO ! ispin

    ! ----------------------------------------------------------------------
    ! Do Filtered Diagonalisation
    ! ----------------------------------------------------------------------

    ! Obtain matrix dimensions. KS and S matrices are symmetric, so
    ! row_block_sizes and col_block_sizes should be identical. The
    ! same applies to the filtered block sizes. Note that filter
    ! matrix will have row_block_sizes equal to that of the original,
    ! and col_block_sizes equal to that of the filtered.  We assume
    ! also that the matrix dimensions are identical for both spin
    ! channels.
    original_rowORcol_block_sizes => cp_dbcsr_row_block_sizes(matrix_ks(1)%matrix)
    filtered_rowORcol_block_sizes => cp_dbcsr_col_block_sizes(matrix_filter(1)%matrix)
    original_nfullrowsORcols_total = cp_dbcsr_nfullrows_total(matrix_ks(1)%matrix)
    filtered_nfullrowsORcols_total = cp_dbcsr_nfullcols_total(matrix_filter(1)%matrix)

    ! filter diagonalisation works on a smaller basis set, and thus
    ! requires a new mo_set (molecular orbitals | eigenvectors) and
    ! the corresponding matrix pools for the eigenvector coefficients
    ALLOCATE(mos_filtered(nspin), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    DO ispin = 1, nspin
       CALL get_mo_set(mo_set=mos(ispin)%mo_set, &
                       maxocc=maxocc, &
                       nelectron=nelectron, &
                       flexible_electron_count=flexible_electron_count)
       NULLIFY(mos_filtered(ispin)%mo_set)
       CALL allocate_mo_set(mo_set=mos_filtered(ispin)%mo_set, &
                            nao=filtered_nfullrowsORcols_total, &
                            nmo=filtered_nfullrowsORcols_total, &
                            nelectron=nelectron, &
                            n_el_f=REAL(nelectron, dp), &
                            maxocc=maxocc, &
                            flexible_electron_count=flexible_electron_count, &
                            error=error)
    END DO ! ispin
    CALL mpools_create(mpools=my_mpools, error=error)
    CALL mpools_rebuild_fm_pools(mpools=my_mpools,&
                                 mos=mos_filtered,&
                                 blacs_env=blacs_env,&
                                 para_env=para_env,&
                                 error=error)

    ! create DBCSR filtered KS matrix, this is reused for each spin
    ! channel
    CALL cp_dbcsr_init(matrix_filtered_ks, error=error)
    ! both row_blk_size and col_blk_size should be that of
    ! col_blk_size of the filter matrix
    CALL cp_dbcsr_create(matrix=matrix_filtered_ks, &
                         name=TRIM("FILTERED_KS_MATRIX"), &
                         dist=cp_dbcsr_distribution(matrix_ks(1)%matrix), &
                         matrix_type=dbcsr_type_no_symmetry, &
                         row_blk_size=filtered_rowORcol_block_sizes, &
                         col_blk_size=filtered_rowORcol_block_sizes, &
                         nblks=0, &
                         nze=0, &
                         error=error)
    CALL cp_dbcsr_finalize(matrix_filtered_ks, error=error)

    ! create DBCSR filtered S (overlap) matrix. Note that
    ! matrix_s(1)%matrix is the orginal overlap matrix---the rest in
    ! the array are derivatives, and it should not depend on
    ! spin. HOWEVER, since the filter matrix is constructed from KS
    ! matrix, and does depend on spin, the filtered S also becomes
    ! spin dependent. Nevertheless this matrix is reused for each spin
    ! channel
    CALL cp_dbcsr_init(matrix_filtered_s, error=error)
    ! both row_blk_size and col_blk_size should be that of
    ! col_blk_size of the filter matrix
    CALL cp_dbcsr_create(matrix=matrix_filtered_s, &
                         name=TRIM("FILTERED_S_MATRIX"), &
                         dist=cp_dbcsr_distribution(matrix_s(1)%matrix), &
                         matrix_type=dbcsr_type_no_symmetry, &
                         row_blk_size=filtered_rowORcol_block_sizes, &
                         col_blk_size=filtered_rowORcol_block_sizes, &
                         nblks=0, &
                         nze=0, &
                         error=error)
    CALL cp_dbcsr_finalize(matrix_filtered_s, error=error)

    ! create temporary matrix for constructing filtered KS and S
    CALL cp_dbcsr_init(matrix_tmp, error=error)
    ! the temporary matrix won't be square
    CALL cp_dbcsr_create(matrix=matrix_tmp, &
                         name=TRIM("TEMPORARY_MATRIX"), &
                         dist=cp_dbcsr_distribution(matrix_s(1)%matrix), &
                         matrix_type=dbcsr_type_no_symmetry, &
                         row_blk_size=original_rowORcol_block_sizes, &
                         col_blk_size=filtered_rowORcol_block_sizes, &
                         nblks=0, &
                         nze=0, &
                         error=error)
    CALL cp_dbcsr_finalize(matrix_tmp, error=error)

    ! create fm format matrices used for diagonalisation
    CALL cp_fm_struct_create(fmstruct=fm_struct, &
                             para_env=para_env, &
                             context=blacs_env, &
                             nrow_global=filtered_nfullrowsORcols_total, &
                             ncol_global=filtered_nfullrowsORcols_total, &
                             error=error)
    ! both fm_matrix_filtered_s and fm_matrix_filtered_ks are reused
    ! for each spin channel
    CALL cp_fm_create(fm_matrix_filtered_s, &
                      fm_struct, &
                      name="FM_MATRIX_FILTERED_S", &
                      error=error)
    CALL cp_fm_create(fm_matrix_filtered_ks, &
                      fm_struct, &
                      name="FM_MATRIX_FILTERED_KS", &
                      error=error)
    ! creaate work matrix
    CALL cp_fm_create(fm_matrix_work, fm_struct, name="FM_MATRIX_WORK", error=error)
    CALL cp_fm_create(fm_matrix_ortho, fm_struct, name="FM_MATRIX_ORTHO", error=error)
    ! all fm matrices are created, so can release fm_struct
    CALL cp_fm_struct_release(fm_struct, error=error)

    ! construct filtered KS, S matrix and diagonalise
    DO ispin = 1, nspin

       ! construct filtered KS matrix
       CALL cp_dbcsr_multiply("N", "N", 1.0_dp, &
                              matrix_ks(ispin)%matrix, matrix_filter(ispin)%matrix, &
                              0.0_dp, matrix_tmp, error=error)
       CALL cp_dbcsr_multiply("T", "N", 1.0_dp, &
                              matrix_filter(ispin)%matrix, matrix_tmp, &
                              0.0_dp, matrix_filtered_ks, error=error)
       ! construct filtered S_matrix
       CALL cp_dbcsr_multiply("N", "N", 1.0_dp, &
                              matrix_s(1)%matrix, matrix_filter(ispin)%matrix, &
                              0.0_dp, matrix_tmp, error=error)
       CALL cp_dbcsr_multiply("T", "N", 1.0_dp, &
                              matrix_filter(ispin)%matrix, matrix_tmp, &
                              0.0_dp, matrix_filtered_s, error=error)

       ! now that we have the filtered KS and S matrices for this spin
       ! channel, perform ordinary diagonalisation

       ! convert DBCSR matrices to fm format
       CALL copy_dbcsr_to_fm(matrix_filtered_s, fm_matrix_filtered_s, error=error)
       CALL copy_dbcsr_to_fm(matrix_filtered_ks, fm_matrix_filtered_ks, error=error)

       ! setup matrix pools for the molecular orbitals
       CALL init_mo_set(mos_filtered(ispin)%mo_set,&
                        my_mpools%ao_mo_fm_pools(ispin)%pool,&
                        name="FILTERED_MOS",&
                        error=error)

       ! now diagonalise
       CALL fb_env_eigensolver(fm_matrix_filtered_ks, &
                               fm_matrix_filtered_s, &
                               mos_filtered(ispin)%mo_set, &
                               fm_matrix_ortho, &
                               fm_matrix_work, &
                               eps_eigval, &
                               ndep, &
                               error)
    END DO ! ispin

    ! release temporary matrices
    CALL cp_dbcsr_release(matrix_filtered_s, error)
    CALL cp_dbcsr_release(matrix_filtered_ks, error)
    CALL cp_fm_release(fm_matrix_filtered_s, error=error)
    CALL cp_fm_release(fm_matrix_filtered_ks, error=error)
    CALL cp_fm_release(fm_matrix_work, error=error)
    CALL cp_fm_release(fm_matrix_ortho, error=error)

    ! ----------------------------------------------------------------------
    ! Construct New Density Matrix
    ! ----------------------------------------------------------------------

    ! calculate filtered molecular orbital occupation numbers and fermi
    ! level etc
    CALL set_mo_occupation(mo_array=mos_filtered, &
                           smear=scf_control%smear, &
                           error=error)

    ! get the filtered density matrix and then convert back to the
    ! full basis version in scf_env ready to be used outside this
    ! subroutine
    ALLOCATE(matrix_filtered_p, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    CALL cp_dbcsr_init(matrix_filtered_p, error=error)
    ! the filtered density matrix should have the same sparse
    ! structure as the original density matrix, we must copy the
    ! sparse structure here, since construction of the density matrix
    ! preserves its sparse form, and therefore matrix_filtered_p must
    ! have its blocks allocated here now. We assume the original
    ! density matrix scf_env%p_mix_new has the same sparse structure
    ! in both spin channels.
    CALL cp_dbcsr_create(matrix=matrix_filtered_p, &
                         name=TRIM("FILTERED_MATRIX_P"), &
                         dist=cp_dbcsr_distribution(scf_env%p_mix_new(1,1)%matrix), &
                         matrix_type=cp_dbcsr_get_matrix_type(scf_env%p_mix_new(1,1)%matrix), &
                         row_blk_size=filtered_rowORcol_block_sizes, &
                         col_blk_size=filtered_rowORcol_block_sizes, &
                         nblks=0, &
                         nze=0, &
                         error=error)
    CALL cp_dbcsr_finalize(matrix_filtered_p, error=error)
    CALL fb_dbcsr_copy_sparse_struct(matrix_filtered_p, &
                                     scf_env%p_mix_new(1,1)%matrix, &
                                     error)
    ! old implementation, using sab_orb to allocate the blocks in matrix_filtered_p
    ! CALL get_qs_env(qs_env=qs_env, sab_orb=sab_orb, error=error)
    ! CALL cp_dbcsr_alloc_block_from_nbl(matrix_filtered_p, sab_orb, error=error)
    CALL cp_dbcsr_set(matrix_filtered_p, 0.0_dp, error=error)

    DO ispin = 1, nspin
      ! calculate matrix_filtered_p
       CALL calculate_density_matrix(mos_filtered(ispin)%mo_set, &
                                     matrix_filtered_p, &
                                     error=error)
       ! convert back to full basis p
       CALL cp_dbcsr_multiply("N", "N", 1.0_dp, &
                              matrix_filter(ispin)%matrix, matrix_filtered_p, &
                              0.0_dp, matrix_tmp, error=error)
       CALL cp_dbcsr_multiply("N", "T", 1.0_dp, &
                              matrix_tmp, matrix_filter(ispin)%matrix, &
                              0.0_dp, scf_env%p_mix_new(ispin,1)%matrix, &
                              retain_sparsity=.TRUE., error=error)
       ! note that we want to retain the sparse structure of
       ! scf_env%p_mix_new
    END DO ! ispin

    ! release temporary matrices
    CALL cp_dbcsr_release(matrix_tmp, error)
    CALL cp_dbcsr_release(matrix_filtered_p, error)
    DEALLOCATE(matrix_filtered_p, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

    ! ----------------------------------------------------------------------
    ! Update MOs
    ! ----------------------------------------------------------------------

    ! we still need to convert mos_filtered back to the full basis
    ! version (mos) for this, we need to update mo_coeff (and/or
    ! mo_coeff_b --- the DBCSR version, if used) of mos

    ! note also that mo_eigenvalues cannot be fully updated, given
    ! that the eigenvalues are computed in a smaller basis, and thus
    ! do not give the full spectron. Printing of molecular states
    ! (molecular DOS) at each SCF step is therefore not recommended
    ! when using this method. The idea is that if one wants a full
    ! molecular DOS, then one should perform a full diagonalisation
    ! without the filters once the SCF has been achieved.

    ! NOTE: from reading the source code, it appears that mo_coeff_b
    ! is actually never used by default (DOUBLE CHECK?!). Even
    ! subroutine eigensolver_dbcsr updates mo_coeff, and not
    ! mo_coeff_b.

    ! create FM format filter matrix
    CALL cp_fm_struct_create(fmstruct=filter_fm_struct, &
                             para_env=para_env, &
                             context=blacs_env, &
                             nrow_global=original_nfullrowsORcols_total, &
                             ncol_global=filtered_nfullrowsORcols_total, &
                             error=error)
    CALL cp_fm_create(fm_matrix_filter, &
                      filter_fm_struct, &
                      name="FM_MATRIX_FILTER", &
                      error=error)
    CALL cp_fm_struct_release(filter_fm_struct, error=error)

    DO ispin = 1, nspin
       ! now the full basis mo_set should only contain the reduced
       ! number of eigenvectors and eigenvalues
       CALL get_mo_set(mo_set=mos_filtered(ispin)%mo_set, &
                       homo=homo_filtered, &
                       lfomo=lfomo_filtered, &
                       nmo=nmo_filtered, &
                       eigenvalues=eigenvalues_filtered, &
                       occupation_numbers=occ_filtered, &
                       mo_coeff=mo_coeff_filtered, &
                       kTS=kTS_filtered, &
                       mu=mu_filtered)
       ! first set all the relevent scalars
       CALL set_mo_set(mo_set=mos(ispin)%mo_set, &
                       homo=homo_filtered, &
                       lfomo=lfomo_filtered, &
                       kTS=kTS_filtered, &
                       mu=mu_filtered, &
                       error=error)
       ! now set the arrays and fm_matrices
       CALL get_mo_set(mo_set=mos(ispin)%mo_set, &
                       nmo=nmo, &
                       occupation_numbers=occ, &
                       eigenvalues=eigenvalues, &
                       mo_coeff=mo_coeff)
       ! number of mos in original mo_set may sometimes be less than
       ! nmo_filtered, so we must make sure we do not go out of bounds
       my_nmo = MIN(nmo, nmo_filtered)
       eigenvalues(:) = 0.0_dp
       eigenvalues(1:my_nmo) = eigenvalues_filtered(1:my_nmo)
       occ(:) = 0.0_dp
       occ(1:my_nmo) = occ_filtered(1:my_nmo)
       ! convert mo_coeff_filtered back to original basis
       CALL cp_fm_set_all(matrix=mo_coeff, alpha=0.0_dp, error=error)
       CALL copy_dbcsr_to_fm(matrix_filter(ispin)%matrix, fm_matrix_filter, error=error)
       CALL cp_fm_gemm("N", "N", &
                       original_nfullrowsORcols_total, &
                       my_nmo, &
                       filtered_nfullrowsORcols_total, &
                       1.0_dp, fm_matrix_filter, mo_coeff_filtered, &
                       0.0_dp, mo_coeff, error=error)

    END DO ! ispin

    ! release temporary matrices
    CALL cp_fm_release(fm_matrix_filter, error=error)

    ! ----------------------------------------------------------------------
    ! Final Clean Up
    ! ----------------------------------------------------------------------

    CALL mpools_release(mpools=my_mpools, error=error)
    DO ispin = 1, nspin
       CALL deallocate_mo_set(mo_set=mos_filtered(ispin)%mo_set, &
                              error=error)
    END DO
    DEALLOCATE(mos_filtered, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    CALL cp_dbcsr_deallocate_matrix_set(matrix_filter, error)

    CALL timestop(handle)

  END SUBROUTINE fb_env_do_diag


! *****************************************************************************
!> \brief The main parallel eigensolver engine for filter matrix diagonalisation
!> \param fm_KS : the BLACS distributed Kohn-Sham matrix, input only
!> \param fm_S  : the BLACS distributed overlap matrix, input only
!> \param mo_set : upon output contains the molecular orbitals (eigenvectors)
!>                 and eigenvalues
!> \param fm_ortho : one of the work matrices, on output, the BLACS distributed
!>                   matrix for orthogalising the eigen problem. E.g. if using
!>                   Cholesky inversse, then the upper triangle part contains
!>                   the inverse of Cholesky U; if not using Cholesky, then it
!>                   contains the S^-1/2.
!> \param fm_work : work matrix used by eigen solver
!> \param eps_eigval : used for quenching the small numbers when computing S^-1/2
!>                     any values less than eps_eigval is truncated to zero.
!> \param ndep : if the overlap is not positive definite, then ndep > 0,
!>               and equals to the number of linear dependent basis functions
!>               in the filtered basis set
!> \param error : cp2k error container
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_env_eigensolver(fm_KS, fm_S, mo_set, fm_ortho, &
                                fm_work, eps_eigval, ndep, error)
    TYPE(cp_fm_type), POINTER                :: fm_KS, fm_S
    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(cp_fm_type), POINTER                :: fm_ortho, fm_work
    REAL(KIND=dp), INTENT(IN)                :: eps_eigval
    INTEGER, INTENT(OUT)                     :: ndep
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_env_eigensolver', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, info, nao, nmo
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mo_eigenvalues
    TYPE(cp_fm_type), POINTER                :: mo_coeff

    CALL timeset(routineN, handle)

    CALL get_mo_set(mo_set=mo_set, &
                    nao=nao, &
                    nmo=nmo, &
                    eigenvalues=mo_eigenvalues, &
                    mo_coeff=mo_coeff)
    ndep = 0

    ! first, obtain orthogonalisation (ortho) matrix
    CALL cp_fm_to_fm(fm_S, fm_ortho, error=error)
    CALL cp_fm_cholesky_decompose(fm_ortho, info_out=info, error=error)
    IF (info == 0) THEN
       ! we are able to use cholesky inverse method
       CALL cp_fm_triangular_invert(fm_ortho, error=error)
       ! solve eigen equation using inverse of cholesky U
       CALL cp_fm_upper_to_full(fm_KS, fm_work, error=error)
       CALL cp_fm_triangular_multiply(fm_ortho, &
                                      fm_KS, &
                                      side="R", &
                                      transpose_tr=.FALSE., &
                                      invert_tr=.FALSE., &
                                      uplo_tr="U", &
                                      n_rows=nao, &
                                      n_cols=nao, &
                                      alpha=1.0_dp, &
                                      error=error)
       CALL cp_fm_triangular_multiply(fm_ortho, &
                                      fm_KS, &
                                      side="L", &
                                      transpose_tr=.TRUE., &
                                      invert_tr=.FALSE., &
                                      uplo_tr="U", &
                                      n_rows=nao, &
                                      n_cols=nao, &
                                      alpha=1.0_dp, &
                                      error=error)
       CALL choose_eigv_solver(fm_KS, fm_work, mo_eigenvalues, error=error)
       CALL cp_fm_triangular_multiply(fm_ortho, &
                                      fm_work, &
                                      side="L", &
                                      transpose_tr=.FALSE., &
                                      invert_tr=.FALSE., &
                                      uplo_tr="U", &
                                      n_rows=nao, &
                                      n_cols=nmo, &
                                      alpha=1.0_dp, &
                                      error=error)
       CALL cp_fm_to_fm(fm_work, mo_coeff, nmo, 1, 1)
    ELSE
       ! cholesky decomposition has failed, revert to calculating
       ! ortho as S^-1/2 using diagonalisation of S, and solve
       ! accordingly
       CALL cp_fm_to_fm(fm_S, fm_ortho, error=error)
       CALL cp_fm_power(fm_ortho, fm_work, -0.5_dp, &
                        eps_eigval, ndep, error=error)
       ! solve eigen equatoin using S^-1/2
       CALL cp_fm_symm("L", "U", nao, nao, 1.0_dp, fm_KS, fm_ortho, &
                       0.0_dp, fm_work, error=error)
       CALL cp_gemm("T", "N", nao, nao, nao, 1.0_dp, fm_ortho, &
                    fm_work, 0.0_dp, fm_KS, error=error)
       CALL choose_eigv_solver(fm_KS, fm_work, mo_eigenvalues, error=error)
       CALL cp_gemm("N", "N", nao, nmo, nao, 1.0_dp, fm_ortho, &
                    fm_work, 0.0_dp, mo_coeff, error=error)
    END IF

    CALL timestop(handle)

  END SUBROUTINE fb_env_eigensolver


! *****************************************************************************
!> \brief Read input sections for filter matrix method
!> \param fb_env : the filter matrix environment
!> \param scf_section : SCF input section
!> \param error : cp2k error container
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_env_read_input(fb_env, scf_section, error)

    TYPE(fb_env_obj), INTENT(INOUT)          :: fb_env
    TYPE(section_vals_type), POINTER         :: scf_section
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_env_read_input', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure, l_val
    REAL(KIND=dp)                            :: r_val
    TYPE(section_vals_type), POINTER         :: fb_section

    CALL timeset(routineN, handle)

    failure = .FALSE.

    NULLIFY(fb_section)
    fb_section => section_vals_get_subs_vals(scf_section, &
                                             "DIAGONALIZATION%FILTER_MATRIX", &
                                             error=error)
    ! filter_temperature
    CALL section_vals_val_get(fb_section, "FILTER_TEMPERATURE", &
                              r_val=r_val, &
                              error=error)
    CALL fb_env_set(fb_env=fb_env, &
                    filter_temperature=r_val, &
                    error=error)
    ! auto_cutoff_scale
    CALL section_vals_val_get(fb_section, "AUTO_CUTOFF_SCALE", &
                              r_val=r_val, &
                              error=error)
    CALL fb_env_set(fb_env=fb_env, &
                    auto_cutoff_scale=r_val, &
                    error=error)
    ! communication model
    CALL section_vals_val_get(fb_section, "COLLECTIVE_COMMUNICATION", &
                              l_val=l_val, &
                              error=error)
    CALL fb_env_set(fb_env=fb_env, &
                    collective_com=l_val, &
                    error=error)

    CALL timestop(handle)

  END SUBROUTINE fb_env_read_input


! *****************************************************************************
!> \brief Automatically generate the cutoff radii of atoms used for
!>        constructing the atomic halos, based on basis set cutoff
!>        ranges for each kind
!> \param fb_env : the filter matrix environment
!> \param qs_env : quickstep environment
!> \param error : cp2k error container
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_env_build_rcut_auto(fb_env, qs_env, error)
    TYPE(fb_env_obj), INTENT(INOUT)          :: fb_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_env_build_rcut_auto', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ikind, nkinds, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: auto_cutoff_scale, kind_radius
    REAL(KIND=dp), DIMENSION(:), POINTER     :: rcut
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(gto_basis_set_p_type), &
      DIMENSION(:), POINTER                  :: basis_set_list
    TYPE(gto_basis_set_type), POINTER        :: basis_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set

    CALL timeset(routineN, handle)

    failure = .FALSE.

    NULLIFY(rcut, qs_kind_set, dft_control)

    CALL get_qs_env(qs_env=qs_env, &
                    qs_kind_set=qs_kind_set, &
                    dft_control=dft_control, &
                    error=error)
    CALL fb_env_get(fb_env=fb_env, &
                    auto_cutoff_scale=auto_cutoff_scale, &
                    error=error)

    nkinds = SIZE(qs_kind_set)
    ALLOCATE(rcut(nkinds), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

    ! reading from the other parts of the code, it seemed that
    ! aux_fit_basis_set is only used when do_admm is TRUE. This can be
    ! seen from the calls to generate_qs_task_list subroutine in
    ! qs_create_task_list, found in qs_environment_methods.F:
    ! basis_type is only set as input parameter for do_admm
    ! calculations, and if not set, the task list is generated using
    ! the default basis_set="ORB".
    ALLOCATE(basis_set_list(nkinds), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (dft_control%do_admm) THEN
       CALL basis_set_list_setup(basis_set_list, "AUX_FIT", qs_kind_set, error=error)
    ELSE
       CALL basis_set_list_setup(basis_set_list, "ORB", qs_kind_set, error=error)
    END IF

    DO ikind = 1, nkinds
       basis_set => basis_set_list(ikind)%gto_basis_set
       CALL get_gto_basis_set(gto_basis_set=basis_set,kind_radius=kind_radius)
       rcut(ikind) = kind_radius * auto_cutoff_scale
    END DO

    CALL fb_env_set(fb_env=fb_env, &
                    rcut=rcut, &
                    error=error)

    ! cleanup
    DEALLOCATE(basis_set_list, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

    CALL timestop(handle)

  END SUBROUTINE fb_env_build_rcut_auto


! *****************************************************************************
!> \brief Builds an fb_atomic_halo_list object using information
!>        from fb_env
!> \param fb_env the fb_env object
!> \param qs_env : quickstep environment (need this to access particle)
!>                 positions and their kinds as well as which particles
!>                 are local to this process
!> \param scf_section : SCF input section, for printing output
!> \param error : cp2k error container
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_env_build_atomic_halos(fb_env, qs_env, scf_section, error)
    TYPE(fb_env_obj), INTENT(INOUT)          :: fb_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(section_vals_type), POINTER         :: scf_section
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_env_build_atomic_halos', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, iatom, ihalo, max_natoms_local, natoms_global, &
      natoms_local, nelectrons, nhalo_atoms, nkinds_global, owner_id_in_halo, &
      stat
    INTEGER, DIMENSION(:), POINTER           :: halo_atoms, local_atoms
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: cost
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: pair_radii
    REAL(KIND=dp), DIMENSION(:), POINTER     :: rcut
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(fb_atomic_halo_list_obj)            :: atomic_halos
    TYPE(fb_atomic_halo_obj), DIMENSION(:), &
      POINTER                                :: halos
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set

    CALL timeset(routineN, handle)

    failure = .FALSE.

    CPPrecondition(fb_env_has_data(fb_env), cp_failure_level, routineP, error, failure)

    NULLIFY(cell, halos, halo_atoms, rcut, particle_set, para_env, &
            qs_kind_set, local_atoms)
    CALL fb_atomic_halo_list_nullify(atomic_halos)

    ! get relevant data from fb_env
    CALL fb_env_get(fb_env=fb_env, &
                    rcut=rcut, &
                    local_atoms=local_atoms, &
                    nlocal_atoms=natoms_local, &
                    error=error)

    ! create atomic_halos
    CALL fb_atomic_halo_list_create(atomic_halos, error)

    ! get the number of atoms and kinds:
    CALL get_qs_env(qs_env=qs_env, &
                    natom=natoms_global, &
                    particle_set=particle_set, &
                    qs_kind_set=qs_kind_set, &
                    nkind=nkinds_global, &
                    para_env=para_env, &
                    cell=cell, &
                    error=error)

    ! get the maximum number of local atoms across the procs.
    max_natoms_local = natoms_local
    CALL mp_max(max_natoms_local, para_env%group)

    ! create the halos, one for each local atom
    ALLOCATE(halos(natoms_local), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    DO ihalo = 1, natoms_local
       CALL fb_atomic_halo_nullify(halos(ihalo))
       CALL fb_atomic_halo_create(halos(ihalo), error)
    END DO
    CALL fb_atomic_halo_list_set(atomic_halos=atomic_halos, &
                                 nhalos=natoms_local, &
                                 max_nhalos=max_natoms_local, &
                                 error=error)
    ! build halos
    ALLOCATE(pair_radii(nkinds_global,nkinds_global), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    CALL fb_build_pair_radii(rcut, nkinds_global, pair_radii)
    ihalo = 0
    DO iatom = 1, natoms_local
       ihalo = ihalo + 1
       CALL fb_atomic_halo_build_halo_atoms(local_atoms(iatom), &
                                            particle_set, &
                                            cell, &
                                            pair_radii, &
                                            halo_atoms, &
                                            nhalo_atoms, &
                                            owner_id_in_halo, &
                                            error)
       CALL fb_atomic_halo_set(atomic_halo=halos(ihalo), &
                               owner_atom=local_atoms(iatom), &
                               owner_id_in_halo=owner_id_in_halo, &
                               natoms=nhalo_atoms, &
                               halo_atoms=halo_atoms, &
                               error=error)
       ! prepare halo_atoms for another halo, do not deallocate, as
       ! original data is being pointed at by the atomic halo data
       ! structure
       NULLIFY(halo_atoms)
       ! calculate the number of electrons in each halo
       nelectrons = fb_atomic_halo_nelectrons_estimate_Z(halos(ihalo), &
                                                         particle_set)
       ! calculate cost
       cost = fb_atomic_halo_cost(halos(ihalo), particle_set, qs_kind_set, error=error)
       CALL fb_atomic_halo_set(atomic_halo=halos(ihalo), &
                               nelectrons=nelectrons, &
                               cost=cost, &
                               error=error)
       ! sort atomic halo
       CALL fb_atomic_halo_sort(halos(ihalo), error)
    END DO ! iatom
    DEALLOCATE(pair_radii, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

    ! finalise
    CALL fb_atomic_halo_list_set(atomic_halos=atomic_halos, &
                                 halos=halos, &
                                 error=error)
    CALL fb_env_set(fb_env=fb_env, &
                    atomic_halos=atomic_halos, &
                    error=error)
    CALL fb_atomic_halo_list_release(atomic_halos, error)

    ! print info
    CALL fb_atomic_halo_list_write_info(atomic_halos, &
                                        para_env, &
                                        scf_section, &
                                        error)

    CALL timestop(handle)

  END SUBROUTINE fb_env_build_atomic_halos


! *****************************************************************************
!> \brief Automatically construct the trial functiosn used for generating
!>        the filter matrix. It tries to use the single zeta subset from
!>        the system GTO basis set as the trial functions
!> \param fb_env : the filter matrix environment
!> \param qs_env : quickstep environment
!> \param maxocc : maximum occupancy for an orbital
!> \param error : cp2k error container
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_env_build_trial_fns_auto(fb_env, qs_env, maxocc, error)

    TYPE(fb_env_obj), INTENT(INOUT)          :: fb_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    REAL(KIND=dp), INTENT(IN)                :: maxocc
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_env_build_trial_fns_auto', &
      routineP = moduleN//':'//routineN

    INTEGER :: counter, handle, icgf, ico, ikind, iset, ishell, itrial, &
      lshell, max_n_trial, nkinds, nset, old_lshell, stat
    INTEGER, DIMENSION(:), POINTER           :: lmax, nfunctions, nshell
    INTEGER, DIMENSION(:, :), POINTER        :: functions
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: zeff
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(fb_trial_fns_obj)                   :: trial_fns
    TYPE(gto_basis_set_p_type), &
      DIMENSION(:), POINTER                  :: basis_set_list
    TYPE(gto_basis_set_type), POINTER        :: basis_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set

    CALL timeset(routineN, handle)

    failure = .FALSE.
    CPPrecondition(fb_env_has_data(fb_env), cp_failure_level, routineP, error, failure)
    NULLIFY(nfunctions, functions, basis_set, basis_set_list, qs_kind_set, dft_control)
    CALL fb_trial_fns_nullify(trial_fns)

    ! create a new trial_fn object
    CALL fb_trial_fns_create(trial_fns, error)

    CALL get_qs_env(qs_env=qs_env, &
                    qs_kind_set=qs_kind_set, &
                    dft_control=dft_control, &
                    error=error)

    nkinds = SIZE(qs_kind_set)

    ! reading from the other parts of the code, it seemed that
    ! aux_fit_basis_set is only used when do_admm is TRUE. This can be
    ! seen from the calls to generate_qs_task_list subroutine in
    ! qs_create_task_list, found in qs_environment_methods.F:
    ! basis_type is only set as input parameter for do_admm
    ! calculations, and if not set, the task list is generated using
    ! the default basis_set="ORB".
    ALLOCATE(basis_set_list(nkinds), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (dft_control%do_admm) THEN
       CALL basis_set_list_setup(basis_set_list, "AUX_FIT", qs_kind_set, error=error)
    ELSE
       CALL basis_set_list_setup(basis_set_list, "ORB", qs_kind_set, error=error)
    END IF

    ALLOCATE(nfunctions(nkinds), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    nfunctions = 0

    DO ikind = 1, nkinds
       ! "gto = gaussian type orbital"
       basis_set => basis_set_list(ikind)%gto_basis_set
       CALL get_gto_basis_set(gto_basis_set=basis_set, &
                              nset=nset, &
                              lmax=lmax, &
                              nshell=nshell)
       CALL get_qs_kind(qs_kind=qs_kind_set(ikind), &
                        zeff=zeff,error=error)

       bset1: DO iset = 1, nset
!          old_lshell = lmax(iset)
          old_lshell = -1
          DO ishell = 1, nshell(iset)
             lshell = basis_set%l(ishell,iset)
             counter = 0
             ! loop over orbitals within the same l
             DO ico = ncoset(lshell-1)+1, ncoset(lshell)
                counter = counter + 1
                ! only include the first zeta orbitals
                IF ((lshell .GT. old_lshell) .AND. (counter .LE. nco(lshell))) THEN
                   nfunctions(ikind) = nfunctions(ikind) + 1
                END IF
             END DO
             ! we have got enough trial functions when we have enough
             ! basis functions to accomodate the number of electrons,
             ! AND that that we have included all the first zeta
             ! orbitals of an angular momentum quantum number l
             IF (((lshell .GT. old_lshell) .OR. (lshell .EQ. lmax(iset))) .AND. &
                 (maxocc * REAL(nfunctions(ikind),dp) .GE. zeff)) THEN
                EXIT bset1
             END IF
             old_lshell = lshell
          END DO
       END DO bset1
    END DO ! ikind

    ! now that we have the number of trial functions get the trial
    ! functions
    max_n_trial = MAXVAL(nfunctions)
    ALLOCATE(functions(max_n_trial,nkinds), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    functions(:,:) = 0
    ! redo the loops to get the trial function indices within the basis set
    DO ikind = 1, nkinds
       ! "gto = gaussian type orbital"
       basis_set => basis_set_list(ikind)%gto_basis_set
       CALL get_gto_basis_set(gto_basis_set=basis_set, &
                              nset=nset, &
                              lmax=lmax, &
                              nshell=nshell)
       CALL get_qs_kind(qs_kind=qs_kind_set(ikind), &
                        zeff=zeff,error=error)
       icgf = 0
       itrial = 0
       bset2: DO iset = 1, nset
          old_lshell = -1
          DO ishell = 1, nshell(iset)
             lshell = basis_set%l(ishell,iset)
             counter = 0
             ! loop over orbitals within the same l
             DO ico = ncoset(lshell-1)+1, ncoset(lshell)
                icgf = icgf + 1
                counter = counter + 1
                ! only include the first zeta orbitals
                IF ((lshell .GT. old_lshell) .AND. (counter .LE. nco(lshell))) THEN
                   itrial = itrial + 1
                   functions(itrial,ikind) = icgf
                END IF
             END DO
             ! we have got enough trial functions when we have more
             ! basis functions than the number of electrons (obtained
             ! from atomic z), AND that that we have included all the
             ! first zeta orbitals of an angular momentum quantum
             ! number l
             IF (((lshell .GT. old_lshell) .OR. (lshell .EQ. lmax(iset))) .AND. &
                 (maxocc * REAL(itrial,dp) .GE. zeff)) THEN
                EXIT bset2
             END IF
             old_lshell = lshell
          END DO
       END DO bset2
    END DO ! ikind

    ! set trial_functions
    CALL fb_trial_fns_set(trial_fns=trial_fns, &
                          nfunctions=nfunctions, &
                          functions=functions, &
                          error=error)
    ! set fb_env
    CALL fb_env_set(fb_env=fb_env, &
                    trial_fns=trial_fns, &
                    error=error)
    CALL fb_trial_fns_release(trial_fns, error)

    ! cleanup
    DEALLOCATE(basis_set_list, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE fb_env_build_trial_fns_auto


! *****************************************************************************
!> \brief Copy the sparse structure of a DBCSR matrix to another, this
!>        means the other matrix will have the same number of blocks
!>        and their corresponding logical locations allocated, although
!>        the blocks does not have to be the same size as the original
!> \param matrix_out : DBCSR matrix whose blocks are to be allocated
!> \param matrix_in  : DBCSR matrix with exising sparse structure that
!>                     is to be copied
!> \param error : cp2k error container
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_dbcsr_copy_sparse_struct(matrix_out, matrix_in, error)

    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_out
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix_in
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_dbcsr_copy_sparse_struct', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iatom, iblk, jatom, &
                                                nblkcols_total, &
                                                nblkrows_total, nblks, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: cols, rows
    LOGICAL                                  :: failure
    REAL(dp), DIMENSION(:, :), POINTER       :: mat_block
    TYPE(cp_dbcsr_iterator)                  :: iter

    failure = .FALSE.

    CALL cp_dbcsr_get_info(matrix=matrix_in, &
                           nblkrows_total=nblkrows_total, &
                           nblkcols_total=nblkcols_total)

    nblks = nblkrows_total * nblkcols_total
    ALLOCATE(rows(nblks), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    ALLOCATE(cols(nblks), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    rows(:) = 0
    cols(:) = 0
    iblk = 0
    nblks = 0
    CALL cp_dbcsr_iterator_start(iter, matrix_in)
    DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
       CALL cp_dbcsr_iterator_next_block(iter, iatom, jatom, mat_block, iblk)
       rows(iblk) = iatom
       cols(iblk) = jatom
       nblks = nblks + 1
    END DO
    CALL cp_dbcsr_iterator_stop(iter)
    CALL cp_dbcsr_reserve_blocks(matrix_out, rows(1:nblks), cols(1:nblks), error=error)
    CALL cp_dbcsr_finalize(matrix_out, error=error)

    ! cleanup
    DEALLOCATE(rows, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    DEALLOCATE(cols, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

  END SUBROUTINE fb_dbcsr_copy_sparse_struct


! *****************************************************************************
!> \brief Write out parameters used for the filter matrix method to
!>        output
!> \param fb_env : the filter matrix environment
!> \param qs_env : quickstep environment
!> \param scf_section : SCF input section
!> \param error : cp2k error container
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_env_write_info(fb_env, qs_env, scf_section, error)
    TYPE(fb_env_obj), INTENT(IN)             :: fb_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(section_vals_type), POINTER         :: scf_section
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_env_write_info', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=2)                         :: element_symbol
    INTEGER                                  :: handle, ikind, nkinds, unit_nr
    LOGICAL                                  :: collective_com, failure
    REAL(KIND=dp)                            :: auto_cutoff_scale, &
                                                filter_temperature
    REAL(KIND=dp), DIMENSION(:), POINTER     :: rcut
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN, handle)

    failure = .FALSE.
    NULLIFY(rcut, atomic_kind_set, logger)

    CALL get_qs_env(qs_env=qs_env, &
                    atomic_kind_set=atomic_kind_set, &
                    error=error)
    CALL fb_env_get(fb_env=fb_env, &
                    filter_temperature=filter_temperature, &
                    auto_cutoff_scale=auto_cutoff_scale, &
                    rcut=rcut, &
                    collective_com=collective_com, &
                    error=error)

    nkinds = SIZE(atomic_kind_set)

    logger => cp_error_get_logger(error=error)
    unit_nr = cp_print_key_unit_nr(logger, scf_section, &
                                   "PRINT%FILTER_MATRIX", &
                                   extension="", &
                                   error=error)
    IF (unit_nr > 0) THEN
       IF (collective_com) THEN
          WRITE (UNIT=unit_nr, FMT="(/,A,T71,A)") &
               " FILTER_MAT_DIAG| MPI communication method:", &
               "Collective"
       ELSE
          WRITE (UNIT=unit_nr, FMT="(/,A,T71,A)") &
               " FILTER_MAT_DIAG| MPI communication method:", &
               "At each step"
       END IF
       WRITE (UNIT=unit_nr, FMT="(A,T71,g10.4)") &
             " FILTER_MAT_DIAG| Filter temperature [K]:", &
             cp_unit_from_cp2k(filter_temperature, "K", error=error)
       WRITE (UNIT=unit_nr, FMT="(A,T71,f10.4)") &
             " FILTER_MAT_DIAG| Filter temperature [a.u.]:", &
             filter_temperature
       WRITE (UNIT=unit_nr, FMT="(A,T71,f10.4)") &
             " FILTER_MAT_DIAG| Auto atomic cutoff radius scale:", &
             auto_cutoff_scale
       WRITE (UNIT=unit_nr, FMT="(A)") &
             " FILTER_MAT_DIAG| atomic cutoff radii [a.u.]"
       DO ikind = 1, nkinds
          CALL get_atomic_kind(atomic_kind=atomic_kind_set(ikind), &
                               element_symbol=element_symbol)
          WRITE (UNIT=unit_nr, FMT="(A,A,T71,f10.4)") &
                " FILTER_MAT_DIAG|   ", element_symbol, rcut(ikind)
       END DO ! ikind
    END IF
    CALL cp_print_key_finished_output(unit_nr, logger, scf_section, &
                                      "PRINT%FILTER_MATRIX", error=error)

    CALL timestop(handle)

  END SUBROUTINE fb_env_write_info


END MODULE qs_fb_env_methods
