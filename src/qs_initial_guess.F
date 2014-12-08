!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Routines to somehow generate an intial guess
!> \par History
!>       2006.03 Moved here from qs_scf.F [Joost VandeVondele]
! *****************************************************************************
MODULE qs_initial_guess
  USE atom_kind_orbitals,              ONLY: calculate_atomic_orbitals
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_interface,              ONLY: &
       array_i1d_obj, array_release, cp_create_bl_distribution, &
       cp_dbcsr_checksum, cp_dbcsr_copy, cp_dbcsr_distribution, &
       cp_dbcsr_distribution_release, cp_dbcsr_filter, cp_dbcsr_get_diag, &
       cp_dbcsr_get_num_blocks, cp_dbcsr_get_occupation, cp_dbcsr_init, &
       cp_dbcsr_iterator, cp_dbcsr_iterator_blocks_left, &
       cp_dbcsr_iterator_next_block, cp_dbcsr_iterator_start, &
       cp_dbcsr_iterator_stop, cp_dbcsr_multiply, cp_dbcsr_nfullrows_total, &
       cp_dbcsr_p_type, cp_dbcsr_release, cp_dbcsr_row_block_sizes, &
       cp_dbcsr_scale, cp_dbcsr_set, cp_dbcsr_set_diag, cp_dbcsr_trace, &
       cp_dbcsr_type, cp_dbcsr_verify_matrix, dbcsr_distribution_mp, &
       dbcsr_distribution_new, dbcsr_distribution_obj, &
       dbcsr_distribution_row_dist, dbcsr_mp_npcols
  USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                             copy_fm_to_dbcsr,&
                                             cp_dbcsr_from_fm,&
                                             cp_dbcsr_sm_fm_multiply
  USE cp_fm_types,                     ONLY: &
       cp_fm_create, cp_fm_get_submatrix, cp_fm_init_random, cp_fm_p_type, &
       cp_fm_release, cp_fm_set_all, cp_fm_set_submatrix, cp_fm_to_fm, &
       cp_fm_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE external_potential_types,        ONLY: all_potential_type,&
                                             gth_potential_type
  USE input_constants,                 ONLY: &
       atomic_guess, core_guess, densities_guess, history_guess, mopac_guess, &
       no_guess, random_guess, restart_guess, sparse_guess
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type
  USE kinds,                           ONLY: default_path_length,&
                                             dp
  USE kpoint_types,                    ONLY: kpoint_type
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_sum
  USE particle_methods,                ONLY: get_particle_set
  USE particle_types,                  ONLY: particle_type
  USE qs_collocate_density,            ONLY: collocate_atomic_charge_density
  USE qs_dftb_utils,                   ONLY: get_dftb_atom_param
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_kind_types,                   ONLY: get_qs_kind,&
                                             get_qs_kind_set,&
                                             qs_kind_type
  USE qs_ks_methods,                   ONLY: qs_ks_update_qs_env
  USE qs_ks_types,                     ONLY: qs_ks_did_change
  USE qs_mo_io,                        ONLY: read_mo_set,&
                                             wfn_restart_file_name
  USE qs_mo_methods,                   ONLY: calculate_density_matrix,&
                                             make_basis_lowdin,&
                                             make_basis_simple,&
                                             make_basis_sm
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type,&
                                             mo_set_restrict,&
                                             set_mo_occupation
  USE qs_rho_methods,                  ONLY: qs_rho_update_rho
  USE qs_rho_types,                    ONLY: qs_rho_get,&
                                             qs_rho_type
  USE qs_scf_methods,                  ONLY: eigensolver,&
                                             eigensolver_simple
  USE qs_scf_types,                    ONLY: block_davidson_diag_method_nr,&
                                             block_krylov_diag_method_nr,&
                                             ot_diag_method_nr,&
                                             qs_scf_env_type
  USE qs_wf_history_methods,           ONLY: wfi_update
  USE scf_control_types,               ONLY: scf_control_type
  USE scp_environment_types,           ONLY: scp_environment_type
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: sort
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_initial_guess'

  PUBLIC ::  calculate_first_density_matrix, calculate_atomic_block_dm, calculate_mopac_dm

  TYPE atom_matrix_type
    REAL(KIND=dp), DIMENSION(:,:), POINTER   :: mat
  END TYPE atom_matrix_type

CONTAINS

! *****************************************************************************
!> \brief can use a variety of methods to come up with an initial
!>      density matrix and optionally an initial wavefunction
!> \param scf_env  SCF environment information
!> \param qs_env   QS environment
!> \param error    CP2K error handling
!> \par History
!>      03.2006 moved here from qs_scf [Joost VandeVondele]
!>      06.2007 allow to skip the initial guess [jgh]
!>      08.2014 kpoints [JGH]
!> \note
!>      badly needs to be split in subroutines each doing one of the possible
!>      schemes
! *****************************************************************************
  SUBROUTINE calculate_first_density_matrix(scf_env,qs_env,error)

    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'calculate_first_density_matrix', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: file_name, filename
    INTEGER :: atom_a, blk, density_guess, group, handle, homo, i, iatom, ic, &
      icol, id_nr, ikind, irow, iseed(4), ispin, istart_col, istart_row, &
      istat, j, last_read, n, n_cols, n_rows, nao, natom, natoms, natoms_tmp, &
      nelectron, nmo, nmo_tmp, not_read, nsgf, nspin, nvec, output_unit, &
      qs_env_id, safe_density_guess, size_atomic_kind_set, z
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: first_sgf, kind_of, last_sgf
    INTEGER, DIMENSION(:), POINTER           :: atom_list, elec_conf, &
                                                nelec_kind, sort_kind
    LOGICAL :: did_guess, do_kpoints, do_std_diag, exist, failure, &
      has_unit_metric, natom_mismatch, ofgpw
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: buff, buff2
    REAL(dp), DIMENSION(:, :), POINTER       :: DATA
    REAL(KIND=dp)                            :: checksum, eps, maxocc, occ, &
                                                scale, total_rho, trps1, zeff
    REAL(KIND=dp), DIMENSION(0:3)            :: edftb
    TYPE(array_i1d_obj)                      :: col_blk_size, col_dist
    TYPE(atom_matrix_type), DIMENSION(:), &
      POINTER                                :: pmat
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_dbcsr_iterator)                  :: iter
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: h_core_sparse, matrix_ks, &
                                                p_rmpv, s_sparse
    TYPE(cp_dbcsr_p_type), DIMENSION(:, :), &
      POINTER                                :: matrix_h_kp, matrix_ks_kp, &
                                                matrix_s_kp, rho_ao_kp
    TYPE(cp_dbcsr_type)                      :: mo_dbcsr, mo_tmp_dbcsr
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: work1
    TYPE(cp_fm_type), POINTER                :: mo_coeff, moa, mob, ortho, &
                                                sv, work2
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dbcsr_distribution_obj)             :: dist
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(kpoint_type), POINTER               :: kpoints
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mo_array
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(qs_kind_type), POINTER              :: qs_kind
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(scp_environment_type), POINTER      :: scp_env
    TYPE(section_vals_type), POINTER         :: dft_section, input, &
                                                subsys_section

    logger => cp_error_get_logger(error)
    failure = .FALSE.
    NULLIFY(atomic_kind, qs_kind, mo_coeff, sv, orb_basis_set, atomic_kind_set,&
         qs_kind_set, particle_set, ortho, work2, work1, mo_array, s_sparse, &
         scf_control, dft_control, p_rmpv, ortho, work2, work1, para_env,&
         s_sparse, scf_control, dft_control, h_core_sparse, matrix_ks, rho)
    NULLIFY(dft_section, input, subsys_section)
    NULLIFY(matrix_s_kp, matrix_h_kp, matrix_ks_kp, rho_ao_kp)
    NULLIFY(moa,mob)
    NULLIFY (atom_list, elec_conf, scp_env, kpoints)
    edftb = 0.0_dp

    CALL timeset(routineN,handle)

    CALL get_qs_env(qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    qs_kind_set=qs_kind_set,&
                    particle_set=particle_set,&
                    mos=mo_array,&
                    matrix_s_kp=matrix_s_kp,&
                    matrix_h_kp=matrix_h_kp,&
                    matrix_ks_kp=matrix_ks_kp,&
                    input=input,&
                    scp_env=scp_env,&
                    scf_control=scf_control,&
                    id_nr=qs_env_id,&
                    dft_control=dft_control,&
                    has_unit_metric=has_unit_metric,&
                    do_kpoints=do_kpoints,&
                    kpoints=kpoints,&
                    rho=rho,&
                    para_env=para_env,&
                    error=error)

    CALL qs_rho_get(rho, rho_ao_kp=rho_ao_kp, error=error)

    ! just initialize the first image, the other density are set to zero
    DO ispin=1, dft_control%nspins
       DO ic=1,SIZE(rho_ao_kp,2)
          CALL cp_dbcsr_set(rho_ao_kp(ispin,ic)%matrix, 0.0_dp,error=error)
       END DO
    END DO
    s_sparse => matrix_s_kp(:,1)
    h_core_sparse => matrix_h_kp(:,1)
    matrix_ks => matrix_ks_kp(:,1)
    p_rmpv => rho_ao_kp(:,1)

    work1 => scf_env%scf_work1
    work2 => scf_env%scf_work2
    ortho => scf_env%ortho

    dft_section =>  section_vals_get_subs_vals(input,"DFT",error=error)

    nspin = dft_control%nspins
    ofgpw = dft_control%qs_control%ofgpw
    density_guess=scf_control%density_guess
    do_std_diag = .FALSE.

    IF(do_kpoints .AND. density_guess==restart_guess) &
       STOP "calculate_first_density_matrix: restart not implemented for k-points"

    safe_density_guess = atomic_guess
    IF ( dft_control%qs_control%semi_empirical .OR. dft_control%qs_control%dftb .OR. &
         dft_control%qs_control%scptb ) THEN
       IF (density_guess == atomic_guess) density_guess = mopac_guess
       ! in case we need to bail to a safe restart type later on
       safe_density_guess = mopac_guess
    END IF

    IF (scf_control%use_ot.AND.&
        (.NOT.((density_guess == random_guess).OR.&
               (density_guess == atomic_guess).OR.&
               (density_guess == mopac_guess).OR.&
               (density_guess == sparse_guess).OR.&
               (((density_guess == restart_guess).OR.&
                (density_guess == history_guess)).AND.&
                (scf_control%level_shift == 0.0_dp))))) THEN
       CALL stop_program(routineN,moduleN,__LINE__,&
            "OT needs GUESS ATOMIC / RANDOM / SPARSE / RESTART / HISTORY RESTART: other options NYI")
    END IF

    ! if a restart was requested, check that the file exists,
    ! if not we fall back to an atomic guess. No kidding, the file name should remain
    ! in sync with read_mo_set
    id_nr=0
    IF (density_guess == restart_guess) THEN
        ! only check existence on I/O node, otherwise if file exists there but
        ! not on compute nodes, everything goes crazy even though only I/O
        ! node actually reads the file
        IF (para_env%ionode) &
           CALL wfn_restart_file_name(file_name,exist,dft_section,logger,error=error)
        CALL mp_bcast(exist, para_env%source, para_env%group)
        CALL mp_bcast(file_name, para_env%source, para_env%group)
        IF (.NOT.exist) THEN
           CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
                "User requested to restart the wavefunction from the file named: "//&
                TRIM(file_name)//". This file does not exist. Please check the existence of"//&
                " the file or change properly the value of the keyword WFN_RESTART_FILE_NAME."//&
                " Calculation continues using ATOMIC GUESS. "//&
CPSourceFileRef,&
                only_ionode=.TRUE.)
           density_guess = safe_density_guess
        END IF
    ELSE IF (density_guess == history_guess) THEN
       IF (para_env%ionode) &
          CALL wfn_restart_file_name(file_name,exist,dft_section,logger,error=error)
       CALL mp_bcast(exist, para_env%source, para_env%group)
       CALL mp_bcast(file_name, para_env%source, para_env%group)
       nvec = qs_env%wf_history%memory_depth
       not_read = nvec+1
       ! At this level we read the saved backup RESTART files..
       DO i=1,nvec
          j = i - 1
          filename = TRIM(file_name)
          IF (j/=0) filename = TRIM(file_name)//".bak-"//ADJUSTL(cp_to_string(j))
          IF (para_env%ionode) &
             INQUIRE(FILE=filename,exist=exist)
          CALL mp_bcast(exist, para_env%source, para_env%group)
          IF ((.NOT. exist) .AND. (i < not_read)) THEN
             not_read = i
          END IF
       END DO
       IF (not_read == 1) THEN
          density_guess = restart_guess
          filename = TRIM(file_name)
          IF (para_env%ionode) &
             INQUIRE(FILE=filename,exist=exist)
          CALL mp_bcast(exist, para_env%source, para_env%group)
          IF (.NOT. exist) THEN
             CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
                  "User requested to restart the wavefunction from a series of restart files named: "//&
                  TRIM(file_name)//" with extensions (.bak-n). These files do not exist."//&
                  " Even trying to switch to a plain restart wave-function failes because the"//&
                  " file named: "//TRIM(file_name)//" does not exist. Please check the existence of"//&
                  " the file or change properly the value of the keyword WFN_RESTART_FILE_NAME. "//&
                  " Calculation continues using ATOMIC GUESS. "//&
CPSourceFileRef,&
                  only_ionode=.TRUE.)
             density_guess = safe_density_guess
          END IF
       END IF
       last_read = not_read - 1
    END IF

    did_guess = .FALSE.

    IF (density_guess == restart_guess) THEN

       CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set)
       CALL read_mo_set(mo_array,atomic_kind_set,qs_kind_set,particle_set,para_env,&
            id_nr=id_nr,multiplicity=dft_control%multiplicity,dft_section=dft_section,&
            natom_mismatch=natom_mismatch, kpoints=kpoints, error=error)

       IF (natom_mismatch) THEN
         density_guess = safe_density_guess
       ELSE
          DO ispin=1,nspin
             IF (scf_control%level_shift /= 0.0_dp) THEN
                CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,&
                                mo_coeff=mo_coeff)
                CALL cp_fm_to_fm(mo_coeff,ortho,error=error)
             END IF

             ! make all nmo vectors present orthonormal
             CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,&
                  mo_coeff=mo_coeff, nmo=nmo, homo=homo)

             IF(has_unit_metric) THEN
               CALL make_basis_simple(mo_coeff,nmo,error=error)
             ELSEIF(dft_control%smear)THEN
               CALL make_basis_lowdin(vmatrix=mo_coeff,ncol=nmo,&
                    matrix_s=s_sparse(1)%matrix,error=error)
             ELSE
               ! ortho so that one can restart for different positions (basis sets?)
               CALL make_basis_sm(mo_coeff,homo,s_sparse(1)%matrix,error=error)
             ENDIF
            ! only alpha spin is kept for restricted
             IF (dft_control%restricted) EXIT
          ENDDO
          IF (dft_control%restricted) CALL mo_set_restrict(mo_array,error=error)

          CALL set_mo_occupation(mo_array,smear=qs_env%scf_control%smear,error=error)

          DO ispin=1,nspin

             IF(scf_control%use_ot .OR. scf_env%method==ot_diag_method_nr) THEN!fm->dbcsr
                CALL copy_fm_to_dbcsr(mo_array(ispin)%mo_set%mo_coeff,mo_array(ispin)%mo_set%mo_coeff_b,error=error)!fm->dbcsr
             ENDIF!fm->dbcsr

             CALL calculate_density_matrix(mo_array(ispin)%mo_set,&
                  p_rmpv(ispin)%matrix,error=error)
          ENDDO
       ENDIF ! natom_mismatch

       did_guess = .TRUE.
    END IF
    IF (density_guess == history_guess) THEN
       IF (not_read > 1) THEN
          DO i=1, last_read
             j = last_read -i
             CALL read_mo_set(mo_array,atomic_kind_set,qs_kind_set,particle_set,para_env,&
                  id_nr=j,multiplicity=dft_control%multiplicity,&
                  dft_section=dft_section, kpoints=kpoints, error=error)

             DO ispin=1,nspin
                IF (scf_control%level_shift /= 0.0_dp) THEN
                   CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,&
                                   mo_coeff=mo_coeff)
                   CALL cp_fm_to_fm(mo_coeff,ortho,error=error)
                END IF

                ! make all nmo vectors present orthonormal
                CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,mo_coeff=mo_coeff, nmo=nmo, homo=homo)

                IF(has_unit_metric) THEN
                   CALL make_basis_simple(mo_coeff,nmo,error=error)
                ELSE
                   ! ortho so that one can restart for different positions (basis sets?)
                   CALL make_basis_sm(mo_coeff,homo,s_sparse(1)%matrix,error=error)
                ENDIF
                ! only alpha spin is kept for restricted
                IF (dft_control%restricted) EXIT
             END DO
             IF (dft_control%restricted) CALL mo_set_restrict(mo_array,error=error)

             DO ispin=1,nspin
                CALL set_mo_occupation(mo_set=mo_array(ispin)%mo_set,&
                                       smear=qs_env%scf_control%smear,&
                                       error=error)
             ENDDO

             DO ispin=1,nspin
          IF(scf_control%use_ot .OR. scf_env%method==ot_diag_method_nr) THEN!fm->dbcsr
CALL copy_fm_to_dbcsr(mo_array(ispin)%mo_set%mo_coeff,mo_array(ispin)%mo_set%mo_coeff_b,error=error)!fm->dbcsr
          ENDIF!fm->dbcsr
                CALL calculate_density_matrix(mo_array(ispin)%mo_set, &
                                              p_rmpv(ispin)%matrix,error=error)
             ENDDO

             ! Write to extrapolation pipeline
             CALL wfi_update(wf_history=qs_env%wf_history, qs_env=qs_env, dt=1.0_dp, error=error)
          END DO
       END IF

       did_guess = .TRUE.
    END IF
    IF (density_guess == random_guess) THEN

       DO ispin=1,nspin
          CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,&
               mo_coeff=mo_coeff, nmo=nmo)
          CALL cp_fm_init_random(mo_coeff,nmo,error=error)
          IF(has_unit_metric) THEN
            CALL make_basis_simple(mo_coeff,nmo,error=error)
          ELSE
            CALL make_basis_sm(mo_coeff,nmo,s_sparse(1)%matrix,error=error)
          ENDIF
          ! only alpha spin is kept for restricted
          IF (dft_control%restricted) EXIT
       ENDDO
       IF (dft_control%restricted) CALL mo_set_restrict(mo_array,error=error)

       DO ispin=1,nspin
          CALL set_mo_occupation(mo_set=mo_array(ispin)%mo_set,&
                                 smear=qs_env%scf_control%smear,&
                                 error=error)
       ENDDO

       DO ispin=1,nspin

          IF(scf_control%use_ot .OR. scf_env%method==ot_diag_method_nr) THEN!fm->dbcsr
              CALL copy_fm_to_dbcsr(mo_array(ispin)%mo_set%mo_coeff,&
                                   mo_array(ispin)%mo_set%mo_coeff_b,error=error)!fm->dbcsr
          ENDIF!fm->dbcsr

          CALL calculate_density_matrix(mo_array(ispin)%mo_set,&
               p_rmpv(ispin)%matrix,error=error)
       ENDDO

       did_guess = .TRUE.
    END IF
    IF (density_guess == core_guess) THEN

       ispin = 1

       ! Load core Hamiltonian into work matrix

       CALL copy_dbcsr_to_fm(h_core_sparse(1)%matrix,work1(ispin)%matrix,error=error)

       ! Diagonalize the core Hamiltonian matrix and retrieve a first set of
       ! molecular orbitals (MOs)

       IF (has_unit_metric) THEN
          CALL eigensolver_simple(matrix_ks=work1(ispin)%matrix,&
                                  mo_set=mo_array(ispin)%mo_set,&
                                  work=work2,&
                                  do_level_shift=.FALSE.,&
                                  level_shift=0.0_dp,&
                                  use_jacobi=.FALSE.,jacobi_threshold=0._dp,&
                                  error=error)
       ELSE
          CALL eigensolver(matrix_ks_fm=work1(ispin)%matrix,&
                           mo_set=mo_array(ispin)%mo_set,&
                           ortho=ortho,&
                           work=work2,&
                           cholesky_method=scf_env%cholesky_method,&
                           use_jacobi=.FALSE.,&
                           error=error)
       END IF

       ! Open shell case: copy alpha MOs to beta MOs

       IF (nspin == 2) THEN
          CALL get_mo_set(mo_set=mo_array(1)%mo_set,mo_coeff=moa)
          CALL get_mo_set(mo_set=mo_array(2)%mo_set,mo_coeff=mob,nmo=nmo)
          CALL cp_fm_to_fm(moa,mob,nmo)
       END IF

       ! Build an initial density matrix (for each spin in the case of
       ! an open shell calculation) from the first MOs set

       DO ispin=1,nspin
         CALL set_mo_occupation(mo_set=mo_array(ispin)%mo_set,&
                                smear=scf_control%smear,&
                                error=error)
         CALL calculate_density_matrix(mo_array(ispin)%mo_set,&
                                       p_rmpv(ispin)%matrix,&
                                       error=error)
       END DO

       did_guess = .TRUE.
    END IF
    IF (density_guess == atomic_guess) THEN

       subsys_section => section_vals_get_subs_vals(input,"SUBSYS",error=error)
       output_unit = cp_print_key_unit_nr(logger,subsys_section,"PRINT%KINDS",extension=".Log",&
                                          error=error)
       IF (output_unit > 0) THEN
         WRITE (UNIT=output_unit,FMT="(/,(T2,A))")&
           "Atomic guess: The first density matrix is obtained in terms of atomic orbitals",&
           "              and electronic configurations assigned to each atomic kind"
       END IF

       DO ispin=1,nspin

          CALL get_mo_set(mo_set=mo_array(ispin)%mo_set, nelectron=nelectron)

          CALL calculate_atomic_block_dm(p_rmpv(ispin)%matrix,s_sparse(1)%matrix, &
                                         particle_set, atomic_kind_set, qs_kind_set, &
                                         ispin, nspin, nelectron, output_unit, error)

          ! The orbital transformation method (OT) requires not only an
          ! initial density matrix, but also an initial wavefunction (MO set)
          IF (ofgpw .AND. (scf_control%use_ot .OR. scf_env%method==ot_diag_method_nr)) THEN
             ! get orbitals later
          ELSE
             IF(ASSOCIATED(scf_env%krylov_space)) do_std_diag = (scf_env%krylov_space%eps_std_diag > 0.0_dp)
             IF (scf_control%use_ot .OR. scf_env%method==ot_diag_method_nr .OR. &
                  (scf_env%method==block_krylov_diag_method_nr .AND. .NOT. do_std_diag) &
                  .OR. dft_control%do_admm .OR. scf_env%method==block_davidson_diag_method_nr) THEN
                IF (dft_control%restricted.AND.(ispin == 2)) THEN
                   CALL mo_set_restrict(mo_array,error=error)
                ELSE
                   CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,&
                                   mo_coeff=mo_coeff,&
                                   nmo=nmo, nao=nao, homo=homo)

                   CALL cp_fm_set_all(mo_coeff,0.0_dp,error=error)
                   CALL cp_fm_init_random(mo_coeff,nmo,error=error)

                   CALL cp_fm_create(sv,mo_coeff%matrix_struct,"SV",error=error)
                   ! multiply times PS
                   IF (has_unit_metric) THEN
                      CALL cp_fm_to_fm(mo_coeff,sv,error=error)
                   ELSE
                      ! PS*C(:,1:nomo)+C(:,nomo+1:nmo) (nomo=NINT(nelectron/maxocc))
                      CALL cp_dbcsr_sm_fm_multiply(s_sparse(1)%matrix,mo_coeff,sv,nmo,error=error)
                   END IF
                   CALL cp_dbcsr_sm_fm_multiply(p_rmpv(ispin)%matrix,sv,mo_coeff,homo,error=error)

                   CALL cp_fm_release(sv,error=error)
                   ! and ortho the result
                   IF (has_unit_metric) THEN
                      CALL make_basis_simple(mo_coeff,nmo,error=error)
                   ELSE
                      CALL make_basis_sm(mo_coeff,nmo,s_sparse(1)%matrix,error=error)
                   END IF
                END IF

                CALL set_mo_occupation(mo_set=mo_array(ispin)%mo_set,&
                                       smear=qs_env%scf_control%smear,&
                                       error=error)

                CALL copy_fm_to_dbcsr(mo_array(ispin)%mo_set%mo_coeff,&
                                      mo_array(ispin)%mo_set%mo_coeff_b,error=error)!fm->dbcsr

                CALL calculate_density_matrix(mo_array(ispin)%mo_set,&
                                              p_rmpv(ispin)%matrix,&
                                              error=error)
             END IF
          END IF

       END DO

       IF (ofgpw .AND. (scf_control%use_ot .OR. scf_env%method==ot_diag_method_nr)) THEN
          ! We fit a function to the square root of the density
          CALL qs_rho_update_rho(rho,qs_env,error=error)
          CPPostcondition(1==0,cp_failure_level,routineP,error,failure)
!         CALL cp_fm_create(sv,mo_coeff%matrix_struct,"SV",error=error)
!         DO ispin=1,nspin
!           CALL integrate_ppl_rspace(qs%rho%rho_r(ispin),qs_env,error=error)
!           CALL cp_cfm_solve(overlap,mos,error)
!           CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,&
!                           mo_coeff=mo_coeff, nmo=nmo, nao=nao)
!           CALL cp_fm_init_random(mo_coeff,nmo,error=error)
!         END DO
!         CALL cp_fm_release(sv,error=error)
       END IF

       CALL cp_print_key_finished_output(output_unit,logger,subsys_section,&
                                         "PRINT%KINDS",error=error)


       did_guess = .TRUE.
    END IF
    IF (density_guess == sparse_guess) THEN

       IF (ofgpw) CALL stop_program(routineN,moduleN,__LINE__,&
                                    "SPARSE_GUESS not implemented for OFGPW")
       IF(.NOT.scf_control%use_ot) CALL stop_program(routineN,moduleN,__LINE__,&
                                                     "OT needed!")

       eps = 1.0E-5_dp

       output_unit= cp_logger_get_default_io_unit(logger)
       group = para_env%group
       natoms = SIZE(particle_set)
       ALLOCATE (kind_of(natoms),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE (first_sgf(natoms),last_sgf(natoms),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

       checksum = cp_dbcsr_checksum(s_sparse(1)%matrix, error=error)
       i = cp_dbcsr_get_num_blocks (s_sparse(1)%matrix); CALL mp_sum(i,group)
       IF(output_unit>0)WRITE(output_unit,*) 'S nblks',i,' checksum',checksum
       CALL cp_dbcsr_filter(s_sparse(1)%matrix, eps, error=error)
       checksum = cp_dbcsr_checksum(s_sparse(1)%matrix, error=error)
       i = cp_dbcsr_get_num_blocks (s_sparse(1)%matrix); CALL mp_sum(i,group)
       IF(output_unit>0)WRITE(output_unit,*) 'S nblks',i,' checksum',checksum

       CALL get_particle_set(particle_set,qs_kind_set,first_sgf=first_sgf,&
                             last_sgf=last_sgf,error=error)
       CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,kind_of=kind_of)

       ALLOCATE (pmat(SIZE(atomic_kind_set)),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

       DO ispin=1,nspin
          scale = 1._dp
          IF (nspin==2) scale=0.5_dp
          DO ikind=1,SIZE(atomic_kind_set)
            atomic_kind => atomic_kind_set(ikind)
            qs_kind     => qs_kind_set(ikind)
            NULLIFY(pmat(ikind)%mat)
            CALL calculate_atomic_orbitals(atomic_kind,qs_kind, pmat=pmat(ikind)%mat,ispin=ispin,error=error)
            NULLIFY(atomic_kind)
          END DO
          CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,&
                          maxocc=maxocc,&
                          nelectron=nelectron)
          !
          CALL cp_dbcsr_iterator_start(iter, p_rmpv(ispin)%matrix)
          DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
             CALL cp_dbcsr_iterator_next_block(iter, irow, icol, DATA, blk)
             ikind = kind_of(irow)
             IF(icol.EQ.irow) DATA(:,:) =  pmat(ikind)%mat(:,:)*scale
          ENDDO
          CALL cp_dbcsr_iterator_stop(iter)

          !CALL cp_dbcsr_verify_matrix(p_rmpv(ispin)%matrix, error=error)
          checksum = cp_dbcsr_checksum(p_rmpv(ispin)%matrix, error=error)
          occ = cp_dbcsr_get_occupation(p_rmpv(ispin)%matrix)
          IF(output_unit>0)WRITE(output_unit,*) 'P_init occ',occ,' checksum',checksum
          ! so far p needs to have the same sparsity as S
          !CALL cp_dbcsr_filter(p_rmpv(ispin)%matrix, eps, error=error)
          !CALL cp_dbcsr_verify_matrix(p_rmpv(ispin)%matrix, error=error)
          checksum = cp_dbcsr_checksum(p_rmpv(ispin)%matrix, error=error)
          occ = cp_dbcsr_get_occupation(p_rmpv(ispin)%matrix)
          IF(output_unit>0)WRITE(output_unit,*) 'P_init occ',occ,' checksum',checksum

          CALL cp_dbcsr_trace(p_rmpv(ispin)%matrix, s_sparse(1)%matrix, trps1, error=error)
          scale=REAL(nelectron,dp)/trps1
          CALL cp_dbcsr_scale(p_rmpv(ispin)%matrix, scale, error=error)

          !CALL cp_dbcsr_verify_matrix(p_rmpv(ispin)%matrix, error=error)
          checksum = cp_dbcsr_checksum(p_rmpv(ispin)%matrix, error=error)
          occ = cp_dbcsr_get_occupation(p_rmpv(ispin)%matrix)
          IF(output_unit>0)WRITE(output_unit,*) 'P occ',occ,' checksum',checksum
          !
          ! The orbital transformation method (OT) requires not only an
          ! initial density matrix, but also an initial wavefunction (MO set)
          IF (dft_control%restricted.AND.(ispin == 2)) THEN
             CALL mo_set_restrict(mo_array,error=error)
          ELSE
             CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,&
                             mo_coeff=mo_coeff,&
                             nmo=nmo, nao=nao, homo=homo)
             CALL cp_fm_set_all(mo_coeff,0.0_dp,error=error)

             n = MAXVAL(last_sgf-first_sgf)+1
             size_atomic_kind_set = SIZE(atomic_kind_set)

             ALLOCATE(buff(n,n),sort_kind(size_atomic_kind_set),&
                      nelec_kind(size_atomic_kind_set),STAT=istat)
             CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
             !
             ! sort kind vs nbr electron
             DO ikind = 1,size_atomic_kind_set
                atomic_kind => atomic_kind_set(ikind)
                qs_kind     => qs_kind_set(ikind)
                CALL get_atomic_kind(atomic_kind=atomic_kind,&
                                     natom=natom,&
                                     atom_list=atom_list,&
                                     z=z)
                CALL get_qs_kind(qs_kind, nsgf=nsgf,elec_conf=elec_conf,&
                                 orb_basis_set=orb_basis_set, zeff=zeff)
                nelec_kind(ikind) = SUM(elec_conf)
             ENDDO
             CALL sort(nelec_kind,size_atomic_kind_set,sort_kind)
             !
             ! a -very- naive sparse guess
             nmo_tmp = nmo
             natoms_tmp = natoms
             istart_col = 1
             iseed(1)=4;iseed(2)=3;iseed(3)=2;iseed(4)=1! set the seed for dlarnv
             DO i = 1,size_atomic_kind_set
                ikind = sort_kind(i)
                atomic_kind => atomic_kind_set(ikind)
                CALL get_atomic_kind(atomic_kind=atomic_kind,&
                                     natom=natom,atom_list=atom_list)
                DO iatom = 1,natom
                   !
                   atom_a = atom_list(iatom)
                   istart_row = first_sgf(atom_a)
                   n_rows = last_sgf(atom_a)-first_sgf(atom_a)+1
                   !
                   ! compute the "potential" nbr of states for this atom
                   n_cols = MAX(INT(REAL(nmo_tmp,dp)/REAL(natoms_tmp,dp)),1)
                   IF(n_cols.GT.n_rows) n_cols = n_rows
                   !
                   nmo_tmp = nmo_tmp - n_cols
                   natoms_tmp = natoms_tmp - 1
                   IF (nmo_tmp.LT.0.OR.natoms_tmp.LT.0) THEN
                      CALL stop_program(routineN,moduleN,__LINE__,"Wrong1!")
                   END IF
                   DO j = 1,n_cols
                      CALL dlarnv(1,iseed,n_rows,buff(1,j))
                   ENDDO
                   CALL cp_fm_set_submatrix(mo_coeff,buff,istart_row,istart_col,&
                                            n_rows,n_cols,error=error)
                   istart_col = istart_col + n_cols
                ENDDO
             ENDDO

             IF (istart_col.LE.nmo) THEN
                CALL stop_program(routineN,moduleN,__LINE__,"Wrong2!")
             END IF

             DEALLOCATE(buff,nelec_kind,sort_kind,STAT=istat)
             CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

             IF(.FALSE.) THEN
                ALLOCATE(buff(nao,1),buff2(nao,1))
                DO i = 1,nmo
                   CALL cp_fm_get_submatrix(mo_coeff,buff,1,i,nao,1,error=error)
                   IF(SUM(buff**2).LT.1E-10_dp) THEN
                      WRITE(*,*) 'wrong',i,SUM(buff**2)
                   ENDIF
                   buff=buff/SQRT(DOT_PRODUCT(buff(:,1),buff(:,1)))
                   DO j = i+1,nmo
                      CALL cp_fm_get_submatrix(mo_coeff,buff2,1,j,nao,1,error=error)
                      buff2=buff2/SQRT(DOT_PRODUCT(buff2(:,1),buff2(:,1)))
                      IF(ABS(DOT_PRODUCT(buff(:,1),buff2(:,1))-1.0_dp).LT.1E-10_dp) THEN
                         WRITE(*,*) 'wrong2',i,j,DOT_PRODUCT(buff(:,1),buff2(:,1))
                         DO ikind=1,nao
                            IF(ABS(mo_coeff%local_data(ikind,i)).gt.1e-10_dp) THEN
                               WRITE(*,*) 'c1',ikind,mo_coeff%local_data(ikind,i)
                            ENDIF
                            IF(ABS(mo_coeff%local_data(ikind,j)).gt.1e-10_dp) THEN
                               WRITE(*,*) 'c2',ikind,mo_coeff%local_data(ikind,j)
                            ENDIF
                         ENDDO
                         STOP
                      ENDIF
                   ENDDO
                ENDDO
                DEALLOCATE(buff,buff2)

             ENDIF
             !
             CALL cp_create_bl_distribution (col_dist, col_blk_size, nmo, &
                  dbcsr_mp_npcols(dbcsr_distribution_mp(cp_dbcsr_distribution( s_sparse(1)%matrix ))))
             CALL dbcsr_distribution_new (dist, dbcsr_distribution_mp (cp_dbcsr_distribution( s_sparse(1)%matrix )),&
                  dbcsr_distribution_row_dist(cp_dbcsr_distribution( s_sparse(1)%matrix )), col_dist)
             !
             CALL cp_dbcsr_init(mo_dbcsr, error=error)
             CALL cp_dbcsr_from_fm(mo_dbcsr, mo_coeff, 0.0_dp, dist,&
                                row_blk_size=cp_dbcsr_row_block_sizes (s_sparse(1)%matrix),&
                                col_blk_size=col_blk_size,&
                                error=error)
             !CALL cp_dbcsr_verify_matrix(mo_dbcsr, error=error)
             checksum = cp_dbcsr_checksum(mo_dbcsr, error=error)

             occ = cp_dbcsr_get_occupation(mo_dbcsr)
             IF(output_unit>0)WRITE(output_unit,*) 'C occ',occ,' checksum',checksum
             CALL cp_dbcsr_filter(mo_dbcsr, eps, error=error)
             !CALL cp_dbcsr_verify_matrix(mo_dbcsr, error=error)
             occ = cp_dbcsr_get_occupation(mo_dbcsr)
             checksum = cp_dbcsr_checksum(mo_dbcsr, error=error)
             IF(output_unit>0)WRITE(output_unit,*) 'C occ',occ,' checksum',checksum
             !
             ! multiply times PS
             IF (has_unit_metric) THEN
                CALL stop_program(routineN,moduleN,__LINE__,&
                                  "has_unit_metric will be removed soon")
             END IF
             !
             ! S*C
             CALL cp_dbcsr_init(mo_tmp_dbcsr, error)
             CALL cp_dbcsr_copy(mo_tmp_dbcsr, mo_dbcsr, name="mo_tmp", error=error)
             CALL cp_dbcsr_multiply("N", "N", 1.0_dp, s_sparse(1)%matrix, mo_dbcsr,&
                  0.0_dp, mo_tmp_dbcsr,&
                  retain_sparsity=.TRUE., error=error)
             !CALL cp_dbcsr_verify_matrix(mo_tmp_dbcsr, error=error)
             checksum = cp_dbcsr_checksum(mo_tmp_dbcsr, error=error)
             occ = cp_dbcsr_get_occupation(mo_tmp_dbcsr)
             IF(output_unit>0)WRITE(output_unit,*) 'S*C occ',occ,' checksum',checksum
             CALL cp_dbcsr_filter(mo_tmp_dbcsr, eps, error=error)
             !CALL cp_dbcsr_verify_matrix(mo_tmp_dbcsr, error=error)
             checksum = cp_dbcsr_checksum(mo_tmp_dbcsr, error=error)
             occ = cp_dbcsr_get_occupation(mo_tmp_dbcsr)
             IF(output_unit>0)WRITE(output_unit,*) 'S*C occ',occ,' checksum',checksum
             !
             ! P*SC
             ! the destroy is needed for the moment to avoid memory leaks !
             ! This one is not needed because _destroy takes care of zeroing.
             CALL cp_dbcsr_multiply("N", "N", 1.0_dp, p_rmpv(ispin)%matrix,&
                  mo_tmp_dbcsr, 0.0_dp, mo_dbcsr, error=error)
             IF(.FALSE.)CALL cp_dbcsr_verify_matrix(mo_dbcsr, error=error)
             checksum = cp_dbcsr_checksum(mo_dbcsr, error=error)
             occ = cp_dbcsr_get_occupation(mo_dbcsr)
             IF(output_unit>0)WRITE(output_unit,*) 'P*SC occ',occ,' checksum',checksum
             CALL cp_dbcsr_filter(mo_dbcsr, eps, error=error)
             !CALL cp_dbcsr_verify_matrix(mo_dbcsr, error=error)
             checksum = cp_dbcsr_checksum(mo_dbcsr, error=error)
             occ = cp_dbcsr_get_occupation(mo_dbcsr)
             IF(output_unit>0)WRITE(output_unit,*) 'P*SC occ',occ,' checksum',checksum
             !
             CALL copy_dbcsr_to_fm(mo_dbcsr, mo_coeff,error=error)

             CALL cp_dbcsr_release(mo_dbcsr, error=error)
             CALL cp_dbcsr_release(mo_tmp_dbcsr, error=error)
             CALL cp_dbcsr_distribution_release(dist)
             CALL array_release(col_blk_size)

             ! and ortho the result
             CALL make_basis_sm(mo_coeff,nmo,s_sparse(1)%matrix,error=error)
          END IF

          CALL set_mo_occupation(mo_set=mo_array(ispin)%mo_set,&
                                 smear=qs_env%scf_control%smear,&
                                 error=error)

          CALL copy_fm_to_dbcsr(mo_array(ispin)%mo_set%mo_coeff,&
                                mo_array(ispin)%mo_set%mo_coeff_b,error=error)!fm->dbcsr

          CALL calculate_density_matrix(mo_array(ispin)%mo_set,&
                                        p_rmpv(ispin)%matrix,&
                                        error=error)
          DO ikind=1,SIZE(atomic_kind_set)
            IF(ASSOCIATED(pmat(ikind)%mat)) THEN
              DEALLOCATE (pmat(ikind)%mat,STAT=istat)
              CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
            END IF
          END DO
       END DO

       DEALLOCATE (pmat,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

       DEALLOCATE (kind_of,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

       DEALLOCATE (first_sgf,last_sgf,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

       did_guess = .TRUE.
    END IF
    IF (density_guess == mopac_guess) THEN

       DO ispin=1,nspin

          CALL get_mo_set(mo_set=mo_array(ispin)%mo_set, nelectron=nelectron)

          CALL calculate_mopac_dm(p_rmpv(ispin)%matrix,s_sparse(1)%matrix, has_unit_metric, dft_control,&
                                  particle_set, atomic_kind_set, qs_kind_set, ispin, nspin, nelectron, para_env, error)

          ! The orbital transformation method (OT) requires not only an
          ! initial density matrix, but also an initial wavefunction (MO set)
          IF(ASSOCIATED(scf_env%krylov_space)) do_std_diag = (scf_env%krylov_space%eps_std_diag > 0.0_dp)
          IF (scf_control%use_ot .OR. scf_env%method==ot_diag_method_nr .OR. &
               (scf_env%method==block_krylov_diag_method_nr .AND. .NOT.do_std_diag)) THEN
             IF (dft_control%restricted.AND.(ispin == 2)) THEN
                CALL mo_set_restrict(mo_array,error=error)
             ELSE
                CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,&
                                mo_coeff=mo_coeff,&
                                nmo=nmo, homo=homo)
                CALL cp_fm_init_random(mo_coeff,nmo,error=error)
                CALL cp_fm_create(sv,mo_coeff%matrix_struct,"SV",error=error)
                ! multiply times PS
                IF (has_unit_metric) THEN
                   CALL cp_fm_to_fm(mo_coeff,sv,error=error)
                ELSE
                   CALL cp_dbcsr_sm_fm_multiply(s_sparse(1)%matrix,mo_coeff,sv,nmo,error=error)
                END IF
                ! here we could easily multiply with the diag that we actually have replicated already
                CALL cp_dbcsr_sm_fm_multiply(p_rmpv(ispin)%matrix,sv,mo_coeff,homo,error=error)
                CALL cp_fm_release(sv,error=error)
                ! and ortho the result
                IF (has_unit_metric) THEN
                   CALL make_basis_simple(mo_coeff,nmo,error=error)
                ELSE
                   CALL make_basis_sm(mo_coeff,nmo,s_sparse(1)%matrix,error=error)
                END IF
             END IF

             CALL set_mo_occupation(mo_set=mo_array(ispin)%mo_set,&
                                    smear=qs_env%scf_control%smear,&
                                    error=error)
             CALL copy_fm_to_dbcsr(mo_array(ispin)%mo_set%mo_coeff,&
                                   mo_array(ispin)%mo_set%mo_coeff_b,error=error)

             CALL calculate_density_matrix(mo_array(ispin)%mo_set,&
                                           p_rmpv(ispin)%matrix,&
                                           error=error)
          END IF
       END DO

       did_guess = .TRUE.
    END IF
    IF (density_guess == densities_guess) THEN

       ! Collocation of the density into the PW-grid
       CALL collocate_atomic_charge_density(total_rho=total_rho, qs_env=qs_env, error=error)

       ! do some assertions here on these matrices having the same structure,
       ! as is currently required
       ! end sparsity check

       qs_env%scf_env%iter_count = 1

       CALL qs_ks_did_change(ks_env=qs_env%ks_env, rho_changed=.TRUE., error=error)

       CALL qs_ks_update_qs_env(qs_env, calculate_forces=.TRUE.,&
                                just_energy=.FALSE., error=error)

       ! diagonalization
       DO ispin=1, dft_control%nspins
          CALL copy_dbcsr_to_fm(matrix_ks(ispin)%matrix,&
                             qs_env%scf_env%scf_work1(ispin)%matrix,&
                             error=error)
       END DO

       qs_env%scf_env%iter_method = "Mixing/Diag"
       qs_env%scf_env%iter_delta = 0.0_dp

       DO ispin=1, dft_control%nspins
         CALL eigensolver(matrix_ks_fm=qs_env%scf_env%scf_work1(ispin)%matrix, &
                          mo_set=mo_array(ispin)%mo_set,&
                          ortho=qs_env%scf_env%ortho, &
                          work=qs_env%scf_env%scf_work2,&
                          cholesky_method=scf_env%cholesky_method,&
                          use_jacobi=.FALSE.,&
                          error=error)
         CALL set_mo_occupation(mo_set=mo_array(ispin)%mo_set,&
                                smear=scf_control%smear,&
                                error=error)
         CALL calculate_density_matrix(mo_array(ispin)%mo_set,&
                                       p_rmpv(ispin)%matrix,&
                                       error=error)
       END DO

       did_guess = .TRUE.
    END IF

    IF (density_guess == no_guess) THEN
       did_guess = .TRUE.
    END IF

    IF (.NOT. did_guess) THEN
       CALL stop_program(routineN,moduleN,__LINE__,&
             "An invalid keyword for the initial density guess was specified")
    END IF

    CALL timestop(handle)

  END SUBROUTINE calculate_first_density_matrix

! *****************************************************************************
!> \brief returns a block diagonal density matrix. Blocks correspond to the atomic densities.
!> \param matrix_p ...
!> \param matrix_s ...
!> \param particle_set ...
!> \param atomic_kind_set ...
!> \param qs_kind_set ...
!> \param ispin ...
!> \param nspin ...
!> \param nelectron_spin ...
!> \param output_unit ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE calculate_atomic_block_dm(matrix_p,matrix_s, particle_set, atomic_kind_set, &
                                       qs_kind_set, ispin, nspin, nelectron_spin, output_unit, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_p, matrix_s
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    INTEGER, INTENT(IN)                      :: ispin, nspin, nelectron_spin, &
                                                output_unit
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_atomic_block_dm', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, handle, icol, ikind, &
                                                irow, istat, natom
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: kind_of
    LOGICAL                                  :: failure
    REAL(dp), DIMENSION(:, :), POINTER       :: DATA
    REAL(KIND=dp)                            :: scale, trps1
    TYPE(atom_matrix_type), ALLOCATABLE, &
      DIMENSION(:)                           :: pmat
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_dbcsr_iterator)                  :: iter
    TYPE(qs_kind_type), POINTER              :: qs_kind

    failure = .FALSE.

    CALL timeset(routineN,handle)

    IF ((output_unit > 0).AND.(nspin > 1)) THEN
      WRITE (UNIT=output_unit,FMT="(/,T2,A,I0)") "Spin ", ispin
    END IF

    CALL cp_dbcsr_set(matrix_p,0.0_dp,error=error)

    natom = SIZE(particle_set)

    ALLOCATE (kind_of(natom),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,kind_of=kind_of)

    ALLOCATE (pmat(SIZE(atomic_kind_set)),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    ! precompute the atomic blocks corresponding to spherical atoms
    DO ikind=1,SIZE(atomic_kind_set)
      atomic_kind => atomic_kind_set(ikind)
      qs_kind     => qs_kind_set(ikind)
      NULLIFY(pmat(ikind)%mat)
      IF (output_unit > 0) THEN
        WRITE (UNIT=output_unit,FMT="(/,T2,A)")&
          "Guess for atomic kind: "//TRIM(atomic_kind%name)
      END IF
      CALL calculate_atomic_orbitals(atomic_kind,qs_kind,iunit=output_unit,pmat=pmat(ikind)%mat,&
                                     ispin=ispin,error=error)
    END DO

    scale = 1.0_dp
    IF (nspin == 2) scale=0.5_dp

    CALL cp_dbcsr_iterator_start(iter, matrix_p)
    DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
       CALL cp_dbcsr_iterator_next_block(iter, irow, icol, DATA, blk)
       ikind = kind_of(irow)
       IF(icol.EQ.irow) DATA(:,:) =  pmat(ikind)%mat(:,:)*scale
    ENDDO
    CALL cp_dbcsr_iterator_stop(iter)

    CALL cp_dbcsr_trace(matrix_p, matrix_s, trps1, error=error)
    scale=REAL(nelectron_spin,dp)/trps1
    CALL cp_dbcsr_scale(matrix_p, scale, error=error)

    IF (output_unit > 0) THEN
      IF (nspin > 1) THEN
        WRITE (UNIT=output_unit,FMT="(T2,A,I1)")&
          "Re-scaling the density matrix to get the right number of electrons for spin ",ispin
      ELSE
        WRITE (UNIT=output_unit,FMT="(T2,A)")&
          "Re-scaling the density matrix to get the right number of electrons"
      END IF
      WRITE (output_unit,'(T19,A,T44,A,T67,A)') "# Electrons","Trace(P)","Scaling factor"
      WRITE (output_unit,'(T20,I10,T40,F12.3,T67,F14.3)') nelectron_spin,trps1,scale
    END IF

    DO ikind=1,SIZE(atomic_kind_set)
      IF(ASSOCIATED(pmat(ikind)%mat)) THEN
        DEALLOCATE (pmat(ikind)%mat,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF
    END DO

    DEALLOCATE (pmat,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    DEALLOCATE (kind_of,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE calculate_atomic_block_dm

! *****************************************************************************
!> \brief returns a block diagonal density matrix. Blocks correspond to the mopac initial guess.
!> \param matrix_p ...
!> \param matrix_s ...
!> \param has_unit_metric ...
!> \param dft_control ...
!> \param particle_set ...
!> \param atomic_kind_set ...
!> \param qs_kind_set ...
!> \param ispin ...
!> \param nspin ...
!> \param nelectron_spin ...
!> \param para_env ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE calculate_mopac_dm(matrix_p,matrix_s, has_unit_metric, &
                                dft_control, particle_set, atomic_kind_set, qs_kind_set, &
                                ispin, nspin, nelectron_spin, para_env, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_p, matrix_s
    LOGICAL                                  :: has_unit_metric
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    INTEGER, INTENT(IN)                      :: ispin, nspin, nelectron_spin
    TYPE(cp_para_env_type)                   :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_mopac_dm', &
      routineP = moduleN//':'//routineN

    INTEGER :: atom_a, group, handle, iatom, ikind, iset, isgf, isgfa, &
      ishell, istat, la, maxl, maxll, nao, natom, ncount, nset, nsgf, z
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: first_sgf
    INTEGER, DIMENSION(:), POINTER           :: atom_list, elec_conf, nshell
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, l, last_sgfa
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: maxocc, my_sum, nelec, paa, &
                                                scale, trps1, trps2, yy, zeff
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: econf, pdiag, sdiag
    REAL(KIND=dp), DIMENSION(0:3)            :: edftb
    TYPE(all_potential_type), POINTER        :: all_potential
    TYPE(gth_potential_type), POINTER        :: gth_potential
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set

    CALL timeset(routineN,handle)

    group = para_env%group
    natom = SIZE(particle_set)
    nao = cp_dbcsr_nfullrows_total(matrix_p)
    IF (nspin==1) THEN
       maxocc=2.0_dp
    ELSE
       maxocc=1.0_dp
    ENDIF

    ALLOCATE (first_sgf(natom),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL get_particle_set(particle_set,qs_kind_set,first_sgf=first_sgf,error=error)
    CALL get_qs_kind_set(qs_kind_set, maxlgto=maxl)

    ALLOCATE (econf(0:maxl),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE (pdiag(nao),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    pdiag(:) = 0.0_dp

    ALLOCATE (sdiag(nao),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    sdiag(:) = 0.0_dp
    IF (has_unit_metric) THEN
       sdiag(:) = 1.0_dp
    ELSE
       CALL cp_dbcsr_get_diag(matrix_s, sdiag, error=error)
       CALL mp_sum(sdiag,group)
    END IF

    ncount   = 0
    trps1    = 0.0_dp
    trps2    = 0.0_dp
    pdiag(:) = 0.0_dp

    IF (nelectron_spin /= 0) THEN
       DO ikind=1,SIZE(atomic_kind_set)

          CALL get_atomic_kind(atomic_kind_set(ikind), natom=natom, atom_list=atom_list)
          CALL get_qs_kind(qs_kind_set(ikind), orb_basis_set=orb_basis_set,&
                               all_potential=all_potential,&
                               gth_potential=gth_potential)

          IF ( dft_control%qs_control%dftb ) THEN
             CALL get_dftb_atom_param(qs_kind_set(ikind)%dftb_parameter,&
                  lmax=maxll,occupation=edftb)
             maxll = MIN(maxll,maxl)
             econf(0:maxl)=edftb(0:maxl)
          ELSEIF (ASSOCIATED(all_potential)) THEN
             CALL get_atomic_kind(atomic_kind_set(ikind), z=z)
             CALL get_qs_kind(qs_kind_set(ikind), nsgf=nsgf, elec_conf=elec_conf, zeff=zeff)
             maxll = MIN(SIZE(elec_conf) - 1,maxl)
             econf(:) = 0.0_dp
             econf(0:maxll) = 0.5_dp*maxocc*REAL(elec_conf(0:maxll),dp)
          ELSE IF (ASSOCIATED(gth_potential)) THEN
             CALL get_atomic_kind(atomic_kind_set(ikind), z=z)
             CALL get_qs_kind(qs_kind_set(ikind), nsgf=nsgf, elec_conf=elec_conf, zeff=zeff)
             maxll = MIN(SIZE(elec_conf) - 1,maxl)
             econf(:) = 0.0_dp
             econf(0:maxll) = 0.5_dp*maxocc*REAL(elec_conf(0:maxll),dp)
          ELSE
             CYCLE
          END IF

          ! MOPAC TYEP GUESS
          IF (dft_control%qs_control%dftb) THEN
             DO iatom=1,natom
                atom_a = atom_list(iatom)
                isgfa=first_sgf(atom_a)
                DO la=0,maxll
                   SELECT CASE (la)
                   CASE (0)
                      pdiag(isgfa) = econf(0)
                   CASE (1)
                      pdiag(isgfa+1) = econf(1)/3._dp
                      pdiag(isgfa+2) = econf(1)/3._dp
                      pdiag(isgfa+3) = econf(1)/3._dp
                   CASE (2)
                      pdiag(isgfa+4) = econf(2)/5._dp
                      pdiag(isgfa+5) = econf(2)/5._dp
                      pdiag(isgfa+6) = econf(2)/5._dp
                      pdiag(isgfa+7) = econf(2)/5._dp
                      pdiag(isgfa+8) = econf(2)/5._dp
                   CASE (3)
                      pdiag(isgfa+ 9) = econf(3)/7._dp
                      pdiag(isgfa+10) = econf(3)/7._dp
                      pdiag(isgfa+11) = econf(3)/7._dp
                      pdiag(isgfa+12) = econf(3)/7._dp
                      pdiag(isgfa+13) = econf(3)/7._dp
                      pdiag(isgfa+14) = econf(3)/7._dp
                      pdiag(isgfa+15) = econf(3)/7._dp
                   CASE DEFAULT
                      CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
                   END SELECT
                END DO
             END DO
          ELSEIF (dft_control%qs_control%semi_empirical) THEN
            yy = REAL(dft_control%charge,KIND=dp)/REAL(nao,KIND=dp)
            DO iatom=1,natom
              atom_a = atom_list(iatom)
              isgfa=first_sgf(atom_a)
              SELECT CASE (nsgf)
              CASE (1) ! s-basis
                pdiag(isgfa   ) = (zeff         - yy)*0.5_dp*maxocc
              CASE (4) ! sp-basis
                IF (z == 1 ) THEN
                  ! special case: hydrogen with sp basis
                  pdiag(isgfa   ) = (zeff         - yy)*0.5_dp*maxocc
                  pdiag(isgfa+ 1) = 0._dp
                  pdiag(isgfa+ 2) = 0._dp
                  pdiag(isgfa+ 3) = 0._dp
                ELSE
                  pdiag(isgfa   ) = (zeff*0.25_dp - yy)*0.5_dp*maxocc
                  pdiag(isgfa+ 1) = (zeff*0.25_dp - yy)*0.5_dp*maxocc
                  pdiag(isgfa+ 2) = (zeff*0.25_dp - yy)*0.5_dp*maxocc
                  pdiag(isgfa+ 3) = (zeff*0.25_dp - yy)*0.5_dp*maxocc
                END IF
              CASE (9) ! spd-basis
                IF (z < 21 .OR. z > 30 .AND. z < 39 .OR. z > 48 .AND. z < 57) THEN
                   !   Main Group Element:  The "d" shell is formally empty.
                   pdiag(isgfa   ) = (zeff*0.25_dp - yy)*0.5_dp*maxocc
                   pdiag(isgfa+ 1) = (zeff*0.25_dp - yy)*0.5_dp*maxocc
                   pdiag(isgfa+ 2) = (zeff*0.25_dp - yy)*0.5_dp*maxocc
                   pdiag(isgfa+ 3) = (zeff*0.25_dp - yy)*0.5_dp*maxocc
                   pdiag(isgfa+ 4) = (             - yy)*0.5_dp*maxocc
                   pdiag(isgfa+ 5) = (             - yy)*0.5_dp*maxocc
                   pdiag(isgfa+ 6) = (             - yy)*0.5_dp*maxocc
                   pdiag(isgfa+ 7) = (             - yy)*0.5_dp*maxocc
                   pdiag(isgfa+ 8) = (             - yy)*0.5_dp*maxocc
                ELSE IF (z < 99) THEN
                   my_sum = zeff - 9.0_dp*yy
                   !   First, put 2 electrons in the 's' shell
                   pdiag(isgfa   ) = (MAX(0.0_dp, MIN(my_sum, 2.0_dp)))*0.5_dp*maxocc
                   my_sum = my_sum - 2.0_dp
                   IF (my_sum > 0.0_dp) THEN
                      !   Now put as many electrons as possible into the 'd' shell
                      pdiag(isgfa+ 4) = (MAX(0.0_dp, MIN(my_sum*0.2_dp, 2.0_dp)))*0.5_dp*maxocc
                      pdiag(isgfa+ 5) = (MAX(0.0_dp, MIN(my_sum*0.2_dp, 2.0_dp)))*0.5_dp*maxocc
                      pdiag(isgfa+ 6) = (MAX(0.0_dp, MIN(my_sum*0.2_dp, 2.0_dp)))*0.5_dp*maxocc
                      pdiag(isgfa+ 7) = (MAX(0.0_dp, MIN(my_sum*0.2_dp, 2.0_dp)))*0.5_dp*maxocc
                      pdiag(isgfa+ 8) = (MAX(0.0_dp, MIN(my_sum*0.2_dp, 2.0_dp)))*0.5_dp*maxocc
                      my_sum = MAX(0.0_dp, my_sum-10.0_dp)
                      !   Put the remaining electrons in the 'p' shell
                      pdiag(isgfa+ 1) = (my_sum/3.0_dp)*0.5_dp*maxocc
                      pdiag(isgfa+ 2) = (my_sum/3.0_dp)*0.5_dp*maxocc
                      pdiag(isgfa+ 3) = (my_sum/3.0_dp)*0.5_dp*maxocc
                   END IF
                END IF
              CASE DEFAULT
                CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
              END SELECT
            END DO
          ELSE
            CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                                   nset=nset,&
                                   nshell=nshell,&
                                   l=l,&
                                   first_sgf=first_sgfa,&
                                   last_sgf=last_sgfa)

            DO iset=1,nset
               DO ishell=1,nshell(iset)
                  la = l(ishell,iset)
                  nelec = maxocc*REAL(2*la + 1,dp)
                  IF (econf(la) > 0.0_dp) THEN
                     IF (econf(la) >= nelec) THEN
                        paa = maxocc
                        econf(la) = econf(la) - nelec
                     ELSE
                        paa = maxocc*econf(la)/nelec
                        econf(la) = 0.0_dp
                        ncount = ncount + NINT(nelec/maxocc)
                     END IF
                     DO isgfa=first_sgfa(ishell,iset),last_sgfa(ishell,iset)
                        DO iatom=1,natom
                           atom_a = atom_list(iatom)
                           isgf = first_sgf(atom_a) + isgfa - 1
                           pdiag(isgf) = paa
                           IF (paa == maxocc) THEN
                              trps1 = trps1 + paa*sdiag(isgf)
                           ELSE
                              trps2 = trps2 + paa*sdiag(isgf)
                           END IF
                        END DO
                     END DO
                  END IF
               END DO ! ishell
            END DO ! iset
          END IF
       END DO ! ikind

       IF (trps2 == 0.0_dp) THEN
          DO isgf=1,nao
             IF (sdiag(isgf) > 0.0_dp) pdiag(isgf) = pdiag(isgf)/sdiag(isgf)
          END DO
       ELSE
          scale = (REAL(nelectron_spin,dp) - trps1)/trps2
          DO isgf=1,nao
             IF (pdiag(isgf) < maxocc) pdiag(isgf) = scale*pdiag(isgf)
          END DO
       END IF
    END IF

    CALL cp_dbcsr_set_diag(matrix_p, pdiag, error=error)

    DEALLOCATE (econf,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    DEALLOCATE (first_sgf,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    DEALLOCATE (pdiag,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    DEALLOCATE (sdiag,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE calculate_mopac_dm

END MODULE qs_initial_guess
