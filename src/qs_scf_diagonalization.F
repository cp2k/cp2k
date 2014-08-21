!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Different diagonalization schemes that can be used
!>        for the iterative solution of the eigenvalue problem
!> \par History
!>      started from routines previously located in the qs_scf module
!>      05.2009
! *****************************************************************************
MODULE qs_scf_diagonalization

  USE cp_array_r_utils,                ONLY: cp_1d_r_p_type
  USE cp_cfm_basic_linalg,             ONLY: cp_cfm_add_fm
  USE cp_cfm_diag,                     ONLY: cp_cfm_geeig
  USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                             cp_cfm_release,&
                                             cp_cfm_to_fm,&
                                             cp_cfm_type
  USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_allocate_matrix_set, cp_dbcsr_col_block_sizes, cp_dbcsr_copy, &
       cp_dbcsr_create, cp_dbcsr_deallocate_matrix, cp_dbcsr_desymmetrize, &
       cp_dbcsr_distribution, cp_dbcsr_init, cp_dbcsr_p_type, &
       cp_dbcsr_row_block_sizes, cp_dbcsr_set, cp_dbcsr_type, &
       dbcsr_type_no_symmetry, dbcsr_type_symmetric
  USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                             copy_fm_to_dbcsr,&
                                             cp_dbcsr_sm_fm_multiply
  USE cp_fm_basic_linalg,              ONLY: cp_fm_symm,&
                                             cp_fm_upper_to_full
  USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_reduce,&
                                             cp_fm_cholesky_restore
  USE cp_fm_diag,                      ONLY: choose_eigv_solver,&
                                             cp_fm_geeig
  USE cp_fm_pool_types,                ONLY: cp_fm_pool_p_type,&
                                             fm_pool_create_fm,&
                                             fm_pool_give_back_fm
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_add_to_element,&
                                             cp_fm_copy_general,&
                                             cp_fm_create,&
                                             cp_fm_get_info,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_gemm_interface,               ONLY: cp_gemm
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE input_constants,                 ONLY: &
       cholesky_dbcsr, cholesky_inverse, cholesky_off, cholesky_reduce, &
       cholesky_restore, core_guess, general_roks, high_spin_roks, &
       restart_guess
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: dp
  USE kpoint_methods,                  ONLY: kpoint_density_matrices,&
                                             kpoint_density_transform,&
                                             kpoint_set_mo_occupation,&
                                             rskp_transform
  USE kpoint_types,                    ONLY: get_kpoint_info,&
                                             kpoint_env_type,&
                                             kpoint_type
  USE machine,                         ONLY: m_flush,&
                                             m_walltime
  USE preconditioner,                  ONLY: prepare_preconditioner,&
                                             restart_preconditioner
  USE qs_density_mixing_types,         ONLY: direct_mixing_nr,&
                                             gspace_mixing_nr
  USE qs_diis,                         ONLY: qs_diis_b_step
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_gspace_mixing,                ONLY: gspace_mixing,&
                                             mixing_allocate,&
                                             mixing_init,&
                                             self_consistency_check
  USE qs_ks_methods,                   ONLY: qs_ks_did_change,&
                                             qs_ks_update_qs_env
  USE qs_ks_types,                     ONLY: qs_ks_env_type
  USE qs_matrix_pools,                 ONLY: mpools_get,&
                                             qs_matrix_pools_type
  USE qs_mo_methods,                   ONLY: calculate_density_matrix,&
                                             calculate_subspace_eigenvalues
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type,&
                                             set_mo_occupation
  USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
  USE qs_ot_eigensolver,               ONLY: ot_eigensolver
  USE qs_rho_atom_types,               ONLY: rho_atom_type
  USE qs_rho_methods,                  ONLY: qs_rho_update_rho
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE qs_scf_block_davidson,           ONLY: generate_extended_space,&
                                             generate_extended_space_sparse
  USE qs_scf_lanczos,                  ONLY: lanczos_refinement,&
                                             lanczos_refinement_2v
  USE qs_scf_methods,                  ONLY: combine_ks_matrices,&
                                             eigensolver,&
                                             eigensolver_dbcsr,&
                                             eigensolver_simple,&
                                             eigensolver_symm,&
                                             scf_env_density_mixing
  USE qs_scf_types,                    ONLY: qs_scf_env_type,&
                                             subspace_env_type
  USE scf_control_types,               ONLY: scf_control_type
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xas_env_types,                   ONLY: xas_environment_type
  USE xas_restart,                     ONLY: find_excited_core_orbital
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_scf_diagonalization'

  PUBLIC :: do_general_diag, do_general_diag_kp, do_roks_diag, &
            do_special_diag, do_ot_diag, do_block_davidson_diag, &
            do_block_krylov_diag, do_scf_diag_subspace, diag_subspace_allocate


CONTAINS

! *****************************************************************************
!> \brief the inner loop of scf, specific to diagonalization with S matrix
!>       basically, in goes the ks matrix out goes a new p matrix
!> \param scf_env ...
!> \param mos ...
!> \param matrix_ks ...
!> \param matrix_s ...
!> \param scf_control ...
!> \param scf_section ...
!> \param diis_step ...
!> \param xas_env ...
!> \param error ...
!> \par History
!>      03.2006 created [Joost VandeVondele]
! *****************************************************************************

  SUBROUTINE do_general_diag(scf_env,mos,matrix_ks,&
                             matrix_s,scf_control,scf_section,&
                             diis_step,xas_env,error)

    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_s
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(section_vals_type), POINTER         :: scf_section
    LOGICAL, INTENT(INOUT)                   :: diis_step
    TYPE(xas_environment_type), OPTIONAL, &
      POINTER                                :: xas_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'do_general_diag', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ispin, nspin
    LOGICAL                                  :: do_level_shift, use_jacobi
    REAL(kind=dp)                            :: diis_error, eps_diis
    TYPE(cp_dbcsr_type), POINTER             :: ortho_dbcsr
    TYPE(cp_fm_type), POINTER                :: ortho

    nspin = SIZE(matrix_ks)
    NULLIFY(ortho,ortho_dbcsr)

    DO ispin=1,nspin
      CALL copy_dbcsr_to_fm(matrix_ks(ispin)%matrix,&
                         scf_env%scf_work1(ispin)%matrix,&
                         error=error)
    END DO

    IF(PRESENT(xas_env) ) THEN
       eps_diis = xas_env%eps_diis
    ELSE
       eps_diis = scf_control%eps_diis
    END IF

    IF (scf_env%iter_count > 1 .AND. .NOT. scf_env%skip_diis) THEN
      CALL qs_diis_b_step(scf_env%scf_diis_buffer,mos,scf_env%scf_work1,&
                          scf_env%scf_work2,scf_env%iter_delta,diis_error,diis_step,&
                          eps_diis,scf_control%nmixing,&
                          s_matrix=matrix_s,&
                          scf_section=scf_section,error=error)
    ELSE
      diis_step = .FALSE.
    END IF

    do_level_shift = ((scf_control%level_shift /= 0.0_dp).AND.&
                      ((scf_control%density_guess == core_guess).OR.&
                       (scf_env%iter_count > 1)))

    IF ((scf_env%iter_count > 1).AND.&
        (scf_env%iter_delta < scf_control%diagonalization%eps_jacobi)) THEN
      use_jacobi = .TRUE.
    ELSE
      use_jacobi = .FALSE.
    END IF

    IF (diis_step) THEN
      scf_env%iter_param = diis_error
      IF (use_jacobi) THEN
        scf_env%iter_method = "DIIS/Jacobi"
      ELSE
        scf_env%iter_method = "DIIS/Diag."
      END IF
    ELSE
      IF(scf_env%mixing_method ==0) THEN
        scf_env%iter_method = "NoMix/Diag."
      ELSE IF(scf_env%mixing_method ==1) THEN
        scf_env%iter_param = scf_env%p_mix_alpha
        IF (use_jacobi) THEN
          scf_env%iter_method = "P_Mix/Jacobi"
        ELSE
          scf_env%iter_method = "P_Mix/Diag."
        END IF
      ELSEIF(scf_env%mixing_method>1) THEN
        scf_env%iter_param = scf_env%mixing_store%alpha
        IF (use_jacobi) THEN
          scf_env%iter_method = TRIM(scf_env%mixing_store%iter_method)//"/Jacobi"
        ELSE
          scf_env%iter_method = TRIM(scf_env%mixing_store%iter_method)//"/Diag."
        END IF
      END IF
    END IF

    IF(scf_env%cholesky_method==cholesky_dbcsr) THEN
       ortho_dbcsr => scf_env%ortho_dbcsr
       DO ispin=1,nspin
         CALL eigensolver_dbcsr(matrix_ks=matrix_ks(ispin)%matrix,matrix_ks_fm=scf_env%scf_work1(ispin)%matrix,&
              mo_set=mos(ispin)%mo_set,&
              ortho_dbcsr=ortho_dbcsr,&
              ksbuf1=scf_env%buf1_dbcsr,ksbuf2=scf_env%buf2_dbcsr,&
              work=scf_env%scf_work2, error=error)
       END DO

    ELSE IF (scf_env%cholesky_method>cholesky_off) THEN
      DO ispin=1,nspin
        IF(scf_env%cholesky_method==cholesky_inverse) THEN
          ortho => scf_env%ortho_m1
        ELSE
          ortho => scf_env%ortho
        END IF
        CALL eigensolver(matrix_ks_fm=scf_env%scf_work1(ispin)%matrix,&
                       mo_set=mos(ispin)%mo_set,&
                       ortho=ortho,&
                       work=scf_env%scf_work2,&
                       cholesky_method=scf_env%cholesky_method,&
                       use_jacobi=use_jacobi,&
                       error=error)
      END DO
    ELSE
      ortho => scf_env%ortho
      DO ispin=1,nspin
        CALL eigensolver_symm(matrix_ks_fm=scf_env%scf_work1(ispin)%matrix,&
                       mo_set=mos(ispin)%mo_set,&
                       ortho=ortho,&
                       work=scf_env%scf_work2,&
                       do_level_shift=do_level_shift,&
                       level_shift=scf_control%level_shift,&
                       use_jacobi=use_jacobi,&
                       jacobi_threshold=scf_control%diagonalization%jacobi_threshold,&
                       error=error)
      END DO
    END IF

    IF (PRESENT(xas_env)) THEN
      CALL find_excited_core_orbital(xas_env,mos,matrix_s,error=error)
      CALL set_mo_occupation(mo_set=mos(1)%mo_set,&
                             smear=xas_env%smear,&
                             xas_env=xas_env,&
                             error=error)
      DO ispin=2,nspin
         CALL set_mo_occupation(mo_set=mos(ispin)%mo_set,&
                                smear=scf_control%smear,&
                                error=error)
      END DO
    ELSE
      CALL set_mo_occupation(mo_array=mos,&
                             smear=scf_control%smear,&
                             error=error)
    END IF

    DO ispin=1,nspin
      CALL calculate_density_matrix(mos(ispin)%mo_set,&
                                    scf_env%p_mix_new(ispin)%matrix,&
                                    error=error)
    END DO

  END SUBROUTINE do_general_diag

! *****************************************************************************
!> \brief Kpoint diagonalization routine 
!>        Transforms matrices to kpoint, distributes kpoint groups, performs
!>        general diagonalization (no storgae of overlap decomposition), stores
!>        MOs, calculates occupation numbers, calculates density matrices 
!>        in kpoint representation, transforms density matrices to real space
!> \param scf_env      SCF environment
!> \param matrix_ks    Kohn-sham matrices (RS indices, global)
!> \param matrix_s     Overlap matrices (RS indices, global)
!> \param kpoints      Kpoint environment
!> \param scf_control  SCF control variables
!> \param scf_section  SCF section (input)
!> \param diis_step    Have we done a DIIS step
!> \param error        CP2K error handling
!> \par History
!>      08.2014 created [JGH]
! *****************************************************************************

  SUBROUTINE do_general_diag_kp(scf_env,matrix_ks,matrix_s,kpoints,&
                                scf_control,scf_section,diis_step,error)

    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:, :), &
      POINTER                                :: matrix_ks, matrix_s
    TYPE(kpoint_type), POINTER               :: kpoints
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(section_vals_type), POINTER         :: scf_section
    LOGICAL, INTENT(INOUT)                   :: diis_step
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'do_general_diag_kp', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp) :: cone = CMPLX(1.0_dp,0.0_dp,KIND=dp), &
      czero = CMPLX(0.0_dp,0.0_dp,KIND=dp), &
      ione = CMPLX(0.0_dp,1.0_dp,KIND=dp)
    INTEGER                                  :: ik, ikp, ikpgr, ispin, &
                                                kplocal, nkp, nkp_groups, &
                                                nspin, stat
    INTEGER, DIMENSION(2)                    :: kp_range
    INTEGER, DIMENSION(:, :), POINTER        :: kp_dist
    INTEGER, DIMENSION(:, :, :), POINTER     :: cell_to_index
    LOGICAL                                  :: do_level_shift, failure, &
                                                use_real_wfn
    REAL(kind=dp), DIMENSION(:), POINTER     :: eigenvalues
    REAL(kind=dp), DIMENSION(:, :), POINTER  :: xkp
    TYPE(cp_cfm_type), POINTER               :: cksmat, cmos, csmat, cwork
    TYPE(cp_dbcsr_type), POINTER             :: cmatrix, rmatrix, tmpmat
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: fmwork
    TYPE(cp_fm_pool_p_type), DIMENSION(:), &
      POINTER                                :: ao_ao_fm_pools
    TYPE(cp_fm_struct_type), POINTER         :: matrix_struct, mo_struct
    TYPE(cp_fm_type), POINTER                :: fmlocal, imos, mo_coeff, &
                                                rksmat, rmos, rsmat
    TYPE(kpoint_env_type), POINTER           :: kp
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_nl
    TYPE(qs_matrix_pools_type), POINTER      :: mpools

    failure = .FALSE.
    diis_step = .FALSE.
    nspin = SIZE(matrix_ks,1)

    do_level_shift = ((scf_control%level_shift /= 0.0_dp).AND.&
                      ((scf_control%density_guess == core_guess).OR.&
                       (scf_env%iter_count > 1)))

    IF(scf_env%mixing_method ==0) THEN
       scf_env%iter_method = "NoMix/Diag."
    ELSE IF(scf_env%mixing_method ==1) THEN
       scf_env%iter_param = scf_env%p_mix_alpha
       scf_env%iter_method = "P_Mix/Diag."
    ELSEIF(scf_env%mixing_method>1) THEN
       scf_env%iter_param = scf_env%mixing_store%alpha
       scf_env%iter_method = TRIM(scf_env%mixing_store%iter_method)//"/Diag."
    END IF

    ! We always do Cholesky
    CPAssert(scf_env%cholesky_method>cholesky_off,cp_failure_level,routineP,error,failure)

    NULLIFY(sab_nl)
    CALL get_kpoint_info(kpoints,nkp=nkp,xkp=xkp,use_real_wfn=use_real_wfn,kp_range=kp_range,&
         nkp_groups=nkp_groups,kp_dist=kp_dist,sab_nl=sab_nl,&
         cell_to_index=cell_to_index,error=error)
    CPAssert(ASSOCIATED(sab_nl),cp_failure_level,routineP,error,failure)
    kplocal = kp_range(2) - kp_range(1) + 1

    ! allocate some work matrices
    ALLOCATE(rmatrix,cmatrix,tmpmat,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL cp_dbcsr_init(rmatrix,error=error)
    CALL cp_dbcsr_init(cmatrix,error=error)
    CALL cp_dbcsr_init(tmpmat,error=error)
    CALL cp_dbcsr_create(rmatrix,template=matrix_ks(1,1)%matrix,error=error)
    CALL cp_dbcsr_create(cmatrix,template=matrix_ks(1,1)%matrix,error=error)
    CALL cp_dbcsr_create(tmpmat,template=matrix_ks(1,1)%matrix,&
         matrix_type=dbcsr_type_no_symmetry,error=error)

    ! init all blocks !?
    CALL cp_dbcsr_copy(rmatrix, matrix_ks(1,1)%matrix, error=error)
    CALL cp_dbcsr_copy(cmatrix, matrix_ks(1,1)%matrix, error=error)

    fmwork => scf_env%scf_work1

    ! fm pools to be used within a kpoint group
    CALL get_kpoint_info(kpoints,mpools=mpools,error=error)
    CALL mpools_get(mpools,ao_ao_fm_pools=ao_ao_fm_pools,error=error)

    CALL fm_pool_create_fm(ao_ao_fm_pools(1)%pool,fmlocal,error=error)
    CALL cp_fm_get_info(fmlocal,matrix_struct=matrix_struct,error=error)
    IF(use_real_wfn) THEN
       CALL cp_fm_create(rksmat,matrix_struct,error=error)
       CALL cp_fm_create(rsmat,matrix_struct,error=error)
    ELSE
       CALL cp_cfm_create(cksmat,matrix_struct,error=error)
       CALL cp_cfm_create(csmat,matrix_struct,error=error)
       CALL cp_cfm_create(cwork,matrix_struct,error=error)
       kp => kpoints%kp_env(1)%kpoint_env
       CALL get_mo_set(kp%mos(1,1)%mo_set,mo_coeff=mo_coeff)
       CALL cp_fm_get_info(mo_coeff,matrix_struct=mo_struct,error=error)
       CALL cp_cfm_create(cmos,mo_struct,error=error)
    END IF

    DO ikpgr = 1,kplocal
       DO ispin=1,nspin
          DO ikp = 1,nkp_groups
             ! number of current kpoint
             ik = kp_dist(1,ikp)+ikpgr-1
             IF(use_real_wfn) THEN
                ! FT of matrices KS and S, then transfer to FM type
                CALL rskp_transform(rmatrix=rmatrix,rsmat=matrix_ks,ispin=ispin,&
                     xkp=xkp(1:3,ik),cell_to_index=cell_to_index,sab_nl=sab_nl,error=error)
                CALL cp_dbcsr_desymmetrize(rmatrix,tmpmat,error=error)
                CALL copy_dbcsr_to_fm(tmpmat,fmwork(1)%matrix,error=error)
                CALL rskp_transform(rmatrix=rmatrix,rsmat=matrix_s,ispin=ispin,&
                     xkp=xkp(1:3,ik),cell_to_index=cell_to_index,sab_nl=sab_nl,error=error)
                CALL cp_dbcsr_desymmetrize(rmatrix,tmpmat,error=error)
                CALL copy_dbcsr_to_fm(tmpmat,fmwork(3)%matrix,error=error)
             ELSE
                ! FT of matrices KS and S, then transfer to FM type
                CALL rskp_transform(rmatrix=rmatrix,cmatrix=cmatrix,rsmat=matrix_ks,ispin=ispin,&
                     xkp=xkp(1:3,ik),cell_to_index=cell_to_index,sab_nl=sab_nl,error=error)
                CALL cp_dbcsr_desymmetrize(rmatrix,tmpmat,error=error)
                CALL copy_dbcsr_to_fm(tmpmat,fmwork(1)%matrix,error=error)
                CALL cp_dbcsr_desymmetrize(cmatrix,tmpmat,error=error)
                CALL copy_dbcsr_to_fm(tmpmat,fmwork(2)%matrix,error=error)
                CALL rskp_transform(rmatrix=rmatrix,cmatrix=cmatrix,rsmat=matrix_s,ispin=ispin,&
                     xkp=xkp(1:3,ik),cell_to_index=cell_to_index,sab_nl=sab_nl,error=error)
                CALL cp_dbcsr_desymmetrize(rmatrix,tmpmat,error=error)
                CALL copy_dbcsr_to_fm(tmpmat,fmwork(3)%matrix,error=error)
                CALL cp_dbcsr_desymmetrize(cmatrix,tmpmat,error=error)
                CALL copy_dbcsr_to_fm(tmpmat,fmwork(4)%matrix,error=error)
             END IF
             ! transfer to kpoint group
             ! redistribution of matrices, new blacs environment
             IF(nkp_groups == 1) THEN
                ! we need still to change the blacs environment
                ! fmwork -> fmlocal -> rksmat
                IF(use_real_wfn) THEN
                   CALL cp_fm_copy_general(fmwork(1)%matrix,rksmat,kpoints%para_env%group,error)
                   CALL cp_fm_copy_general(fmwork(3)%matrix,rsmat,kpoints%para_env%group,error)
                ELSE
                   ! Hermitian ??????
                   CALL cp_fm_copy_general(fmwork(1)%matrix,fmlocal,kpoints%para_env%group,error)
                   CALL cp_cfm_add_fm(czero,cksmat,cone,fmlocal,error)
                   CALL cp_fm_copy_general(fmwork(2)%matrix,fmlocal,kpoints%para_env%group,error)
                   CALL cp_cfm_add_fm(cone,cksmat,ione,fmlocal,error)
                   CALL cp_fm_copy_general(fmwork(3)%matrix,fmlocal,kpoints%para_env%group,error)
                   CALL cp_cfm_add_fm(czero,csmat,cone,fmlocal,error)
                   CALL cp_fm_copy_general(fmwork(4)%matrix,fmlocal,kpoints%para_env%group,error)
                   CALL cp_cfm_add_fm(cone,csmat,ione,fmlocal,error)
                END IF
             ELSE
                ! this requires more complex reshuffling of data
                ! CALL cp_fm_copy_general(fmwork(1)%matrix,rksmat,kpoints%para_env%group,error)
                CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
             END IF
          END DO
          ! Each kpoint group has now information on a kpoint to be diagonalized
          ! General eigensolver Hermite or Symmetric
          kp => kpoints%kp_env(ikpgr)%kpoint_env
          IF(use_real_wfn) THEN
             CALL get_mo_set(kp%mos(1,ispin)%mo_set,mo_coeff=mo_coeff,eigenvalues=eigenvalues)
             CALL cp_fm_geeig(rksmat,rsmat,mo_coeff,eigenvalues,fmlocal,error)
          ELSE
             CALL get_mo_set(kp%mos(1,ispin)%mo_set,mo_coeff=rmos,eigenvalues=eigenvalues)
             CALL get_mo_set(kp%mos(2,ispin)%mo_set,mo_coeff=imos)
             CALL cp_cfm_geeig(cksmat,csmat,cmos,eigenvalues,cwork,error)
             ! copy eigenvalues to imag set (keep them in sync)
             kp%mos(2,ispin)%mo_set%eigenvalues=eigenvalues
             ! split real and imaginary part of mos
             CALL cp_cfm_to_fm(cmos,rmos,imos,error)
          END IF
       END DO
    END DO
    ! MO occupations
    CALL kpoint_set_mo_occupation(kpoints,scf_control%smear,error)

    ! density matrices
    CALL kpoint_density_matrices(kpoints,error)
    ! density matrices in real space
    CALL kpoint_density_transform(kpoints,scf_env%p_mix_kp,fmwork,error)

    CALL cp_dbcsr_deallocate_matrix(rmatrix,error=error)
    CALL cp_dbcsr_deallocate_matrix(cmatrix,error=error)
    CALL cp_dbcsr_deallocate_matrix(tmpmat,error=error)

    IF(use_real_wfn) THEN
       CALL cp_fm_release(rksmat,error=error)
       CALL cp_fm_release(rsmat,error=error)
    ELSE
       CALL cp_cfm_release(cksmat,error=error)
       CALL cp_cfm_release(csmat,error=error)
       CALL cp_cfm_release(cwork,error=error)
       CALL cp_cfm_release(cmos,error=error)
    END IF
    CALL fm_pool_give_back_fm(ao_ao_fm_pools(1)%pool, fmlocal, error=error)

  END SUBROUTINE do_general_diag_kp

! *****************************************************************************
!> \brief inner loop within MOS subspace, to refine occupation and density,
!>        before next diagonalization of the Hamiltonian
!> \param qs_env ...
!> \param scf_env ...
!> \param subspace_env ...
!> \param mos ...
!> \param rho ...
!> \param ks_env ...
!> \param scf_section ...
!> \param scf_control ...
!> \param error ...
!> \par History
!>      09.2009 created [MI]
!> \note  it is assumed that when diagonalization is used, also some mixing procedure is active
! *****************************************************************************
  SUBROUTINE do_scf_diag_subspace(qs_env,scf_env,subspace_env,mos,rho,&
                  ks_env,scf_section,scf_control,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(subspace_env_type), POINTER         :: subspace_env
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(section_vals_type), POINTER         :: scf_section
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'do_scf_diag_subspace', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: rone = 1.0_dp, rzero = 0.0_dp

    INTEGER                                  :: handle, i, iloop, ispin, &
                                                istat, nao, nmo, nspin, &
                                                output_unit
    LOGICAL                                  :: converged, failure
    REAL(dp)                                 :: ene_diff, ene_old, &
                                                iter_delta, max_val, &
                                                sum_band, sum_val, t1, t2
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mo_eigenvalues, mo_occupations
    TYPE(cp_1d_r_p_type), ALLOCATABLE, &
      DIMENSION(:)                           :: eval_first, occ_first
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_s
    TYPE(cp_fm_type), POINTER                :: c0, chc, evec, mo_coeff, work
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(rho_atom_type), DIMENSION(:), &
      POINTER                                :: rho_atom

    failure = .FALSE.
    CALL timeset(routineN,handle)
    NULLIFY(c0, chc, energy, evec, matrix_ks, mo_coeff, mo_eigenvalues, &
         mo_occupations, work)

    logger => cp_error_get_logger(error)
    output_unit = cp_print_key_unit_nr(logger,scf_section,"PRINT%DIAG_SUB_SCF",&
            extension=".scfLog",error=error)

    !Extra loop keeping mos unchanged and refining the subspace occupation
    work  => scf_env%scf_work2

    nspin = SIZE(mos)

    ALLOCATE( eval_first(nspin), STAT=istat)
      CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE( occ_first(nspin),STAT=istat)
      CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    DO ispin=1,nspin
       CALL get_mo_set(mo_set=mos(ispin)%mo_set,&
                      nmo=nmo,&
                      eigenvalues=mo_eigenvalues,&
                      occupation_numbers=mo_occupations)
       ALLOCATE( eval_first(ispin)%array(nmo),STAT=istat)
      CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE( occ_first(ispin)%array(nmo),STAT=istat)
      CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
       eval_first(ispin)%array(1:nmo) = mo_eigenvalues(1:nmo)
       occ_first(ispin)%array(1:nmo) = mo_occupations(1:nmo)
    END DO

    DO ispin=1,nspin
       CALL cp_dbcsr_copy(subspace_env%p_matrix_store(ispin)%matrix, rho%rho_ao(ispin)%matrix,error=error )
       CALL cp_dbcsr_copy( rho%rho_ao(ispin)%matrix, scf_env%p_mix_new(ispin)%matrix,error=error )
    END DO

    subspace_env%p_matrix_mix => scf_env%p_mix_new

    NULLIFY(matrix_ks, energy, para_env, matrix_s)
    CALL get_qs_env(qs_env=qs_env,matrix_ks=matrix_ks,energy=energy, &
          matrix_s=matrix_s,&
          para_env=para_env, error=error)

    ! mixing storage allocation
    IF (subspace_env%mixing_method>=gspace_mixing_nr) THEN
         CALL mixing_allocate(qs_env,subspace_env%mixing_method,subspace_env%mixing_store,&
              nspins=nspin,error=error)
         IF(qs_env % dft_control % qs_control%gapw) THEN
            CALL get_qs_env(qs_env=qs_env,&
                 rho_atom_set=rho_atom,error=error)
            CALL mixing_init(subspace_env%mixing_method,rho,subspace_env%mixing_store,&
                           para_env,rho_atom=rho_atom,error=error)
         ELSE
            CALL mixing_init(subspace_env%mixing_method,rho,subspace_env%mixing_store,&
                           para_env,error=error)
         END IF
    END IF

    ene_old = 0.0_dp
    ene_diff = 0.0_dp
    IF (output_unit > 0) THEN
          WRITE(output_unit,"(/T19,A)") '<<<<<<<<<   SUBSPACE ROTATION    <<<<<<<<<<'
          WRITE(output_unit,"(T4,A,T13,A,T21,A,T38,A,T51,A,T65,A/,T4,A)")   &
           "In-step","Time","Convergence", "Band ene.","Total ene.","Energy diff.",REPEAT("-",74)
    END IF
   ! update of density
    CALL qs_rho_update_rho(rho, qs_env=qs_env, error=error)

    DO iloop = 1,subspace_env%max_iter
      t1 = m_walltime()
      converged=.FALSE.
      ene_old =  energy%total

      CALL qs_ks_did_change(ks_env,rho_changed=.TRUE.,error=error)
      CALL qs_ks_update_qs_env(ks_env,qs_env=qs_env,calculate_forces=.FALSE.,&
             just_energy=.FALSE.,print_active=.FALSE.,error=error)

      max_val = 0.0_dp
      sum_val = 0.0_dp
      sum_band = 0.0_dp
      DO ispin = 1,SIZE(matrix_ks)
         CALL get_mo_set(mo_set=mos(ispin)%mo_set,&
                      nao=nao,&
                      nmo=nmo,&
                      eigenvalues=mo_eigenvalues,&
                      occupation_numbers=mo_occupations,&
                      mo_coeff=mo_coeff)

        !compute C'HC
         chc => subspace_env%chc_mat(ispin)%matrix
         evec => subspace_env%c_vec(ispin)%matrix
         c0 => subspace_env%c0(ispin)%matrix
         CALL cp_fm_to_fm(mo_coeff,c0,error=error)
         CALL cp_dbcsr_sm_fm_multiply(matrix_ks(ispin)%matrix,c0,work,nmo,error=error)
         CALL cp_gemm('T','N',nmo,nmo,nao,rone,c0,work,rzero,chc,error=error)
        !diagonalize C'HC
         CALL choose_eigv_solver(chc,evec,mo_eigenvalues,error=error)

        !rotate the mos by the eigenvectors of C'HC
         CALL cp_gemm('N','N',nao,nmo,nmo,rone,c0,evec,rzero,mo_coeff,error=error)

         CALL set_mo_occupation(mo_set=mos(ispin)%mo_set,&
                           smear=scf_control%smear,&
                           error=error)

         CALL calculate_density_matrix(mos(ispin)%mo_set,&
                                  subspace_env%p_matrix_mix(ispin)%matrix,&
                                 error=error)

         DO i = 1, nmo
           sum_band = sum_band + mo_eigenvalues(i)*mo_occupations(i)
         END DO

        !check for self consistency
      END DO

      IF(subspace_env%mixing_method==direct_mixing_nr) THEN
         CALL scf_env_density_mixing(subspace_env%p_matrix_mix, &
              scf_env%mixing_store, rho%rho_ao, qs_env%para_env, iter_delta, iloop, &
              error=error)
      ELSE
          CALL self_consistency_check(rho%rho_ao,scf_env%p_delta,para_env,&
               subspace_env%p_matrix_mix,delta=iter_delta,error=error)
      END IF

      DO ispin=1,nspin
         CALL cp_dbcsr_copy( rho%rho_ao(ispin)%matrix,subspace_env%p_matrix_mix(ispin)%matrix ,error=error )
      END DO
      ! update of density
      CALL qs_rho_update_rho(rho, qs_env=qs_env, error=error)
      ! Mixing in reciprocal space
      IF(subspace_env%mixing_method>=gspace_mixing_nr) THEN
         CALL gspace_mixing(qs_env, scf_env, subspace_env%mixing_store, &
            rho, qs_env%para_env, error=error)
      END IF

      ene_diff =  energy%total-ene_old
      converged=(ABS(ene_diff)<subspace_env%eps_ene .AND. &
               iter_delta<subspace_env%eps_adapt*scf_env%iter_delta)
      t2 = m_walltime()
      IF (output_unit > 0) THEN
           WRITE(output_unit,"(T4,I5,T11,F8.3,T18,E14.4,T34,F12.5,T46,F16.8,T62,E14.4)")   &
                iloop, t2-t1, iter_delta, sum_band, energy%total, ene_diff
          CALL m_flush(output_unit)
      END IF
      IF(converged) THEN
        IF (output_unit > 0) WRITE(output_unit,"(T10,A,I6,A,/)")&
                   " Reached convergence in ",iloop," iterations "
        EXIT
      END IF

    END DO ! iloop

    NULLIFY (subspace_env%p_matrix_mix)
    DO ispin=1,nspin
        CALL cp_dbcsr_copy( scf_env%p_mix_new(ispin)%matrix,rho%rho_ao(ispin)%matrix,error=error )
        CALL cp_dbcsr_copy( rho%rho_ao(ispin)%matrix,  subspace_env%p_matrix_store(ispin)%matrix ,error=error )

        DEALLOCATE(eval_first(ispin)%array,occ_first(ispin)%array)
    END DO
    DEALLOCATE(eval_first,occ_first)

    CALL timestop(handle)

END SUBROUTINE  do_scf_diag_subspace


! *****************************************************************************
!> \brief ...
!> \param subspace_env ...
!> \param qs_env ...
!> \param mos ...
!> \param error ...
! *****************************************************************************
SUBROUTINE diag_subspace_allocate(subspace_env,qs_env,mos,error)

    TYPE(subspace_env_type), POINTER         :: subspace_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'diag_subspace_allocate', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ispin, istat, nmo, &
                                                nspin
    LOGICAL                                  :: failure
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb

    failure = .FALSE.

    CALL timeset(routineN,handle)

    NULLIFY(sab_orb, matrix_s)
    CALL get_qs_env(qs_env=qs_env, sab_orb=sab_orb, &
         matrix_s=matrix_s, &
         error=error)

    nspin = SIZE(mos)
!   *** allocate p_atrix_store ***
    IF (.NOT.ASSOCIATED(subspace_env%p_matrix_store)) THEN
       CALL cp_dbcsr_allocate_matrix_set(subspace_env%p_matrix_store,nspin,error=error)

       DO i=1,nspin
          ALLOCATE(subspace_env%p_matrix_store(i)%matrix)
          CALL cp_dbcsr_init(subspace_env%p_matrix_store(i)%matrix,error=error)
          CALL cp_dbcsr_create(matrix=subspace_env%p_matrix_store(i)%matrix, &
               name="DENSITY_STORE", &
               dist=cp_dbcsr_distribution(matrix_s(1)%matrix), matrix_type=dbcsr_type_symmetric,&
               row_blk_size=cp_dbcsr_row_block_sizes(matrix_s(1)%matrix), &
               col_blk_size=cp_dbcsr_col_block_sizes(matrix_s(1)%matrix), &
               nblks=0, nze=0, error=error)
          CALL cp_dbcsr_alloc_block_from_nbl(subspace_env%p_matrix_store(i)%matrix,&
               sab_orb, error=error)
          CALL cp_dbcsr_set(subspace_env%p_matrix_store(i)%matrix,0.0_dp,error=error)
       ENDDO

    END IF

    ALLOCATE(subspace_env%chc_mat(nspin),STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(subspace_env%c_vec(nspin),STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(subspace_env%c0(nspin),STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    DO ispin = 1,nspin
       CALL get_mo_set(mos(ispin)%mo_set,mo_coeff=mo_coeff, nmo=nmo)
       CALL cp_fm_create(subspace_env%c0(ispin)%matrix,mo_coeff%matrix_struct,error=error)
       NULLIFY(fm_struct_tmp)
       CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=nmo, ncol_global=nmo,&
                                   para_env=mo_coeff%matrix_struct%para_env, &
                                   context=mo_coeff%matrix_struct%context,error=error)
       CALL cp_fm_create(subspace_env%chc_mat(ispin)%matrix,fm_struct_tmp,"chc",error=error)
       CALL cp_fm_create(subspace_env%c_vec(ispin)%matrix,fm_struct_tmp,"vec",error=error)
       CALL cp_fm_struct_release(fm_struct_tmp,error=error)
    END DO

    CALL timestop(handle)

END SUBROUTINE diag_subspace_allocate

! *****************************************************************************
!> \brief the inner loop of scf, specific to diagonalization without S matrix
!>       basically, in goes the ks matrix out goes a new p matrix
!> \param scf_env ...
!> \param mos ...
!> \param matrix_ks ...
!> \param scf_control ...
!> \param scf_section ...
!> \param diis_step ...
!> \param error ...
!> \par History
!>      03.2006 created [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE do_special_diag(scf_env,mos,matrix_ks,scf_control,&
                                         scf_section,diis_step,error)

    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(section_vals_type), POINTER         :: scf_section
    LOGICAL, INTENT(INOUT)                   :: diis_step
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: ispin, nspin
    LOGICAL                                  :: do_level_shift, use_jacobi
    REAL(kind=dp)                            :: diis_error

    nspin = SIZE(matrix_ks)

    DO ispin=1,nspin
      CALL copy_dbcsr_to_fm(matrix_ks(ispin)%matrix,scf_env%scf_work1(ispin)%matrix,error=error)
    END DO
    IF (scf_env%iter_count > 1 .AND. .NOT. scf_env%skip_diis) THEN
      CALL qs_diis_b_step(scf_env%scf_diis_buffer,mos,scf_env%scf_work1,&
                          scf_env%scf_work2,scf_env%iter_delta,diis_error,diis_step,&
                          scf_control%eps_diis,scf_control%nmixing,&
                          scf_section=scf_section,&
                          error=error)
    ELSE
      diis_step = .FALSE.
    END IF

    IF ((scf_env%iter_count > 1).AND.(scf_env%iter_delta < scf_control%diagonalization%eps_jacobi)) THEN
      use_jacobi = .TRUE.
    ELSE
      use_jacobi = .FALSE.
    END IF

    do_level_shift = ((scf_control%level_shift /= 0.0_dp).AND.&
                      ((scf_control%density_guess == core_guess).OR.(scf_env%iter_count > 1)))
    IF (diis_step) THEN
      scf_env%iter_param = diis_error
      IF (use_jacobi) THEN
        scf_env%iter_method = "DIIS/Jacobi"
      ELSE
        scf_env%iter_method = "DIIS/Diag."
      END IF
    ELSE
      IF(scf_env%mixing_method ==1) THEN
        scf_env%iter_param = scf_env%p_mix_alpha
        IF (use_jacobi) THEN
          scf_env%iter_method = "P_Mix/Jacobi"
        ELSE
          scf_env%iter_method = "P_Mix/Diag."
        END IF
      ELSEIF(scf_env%mixing_method>1) THEN
        scf_env%iter_param = scf_env%mixing_store%alpha
        IF (use_jacobi) THEN
          scf_env%iter_method = TRIM(scf_env%mixing_store%iter_method)//"/Jacobi"
        ELSE
          scf_env%iter_method = TRIM(scf_env%mixing_store%iter_method)//"/Diag."
        END IF
      END IF
    END IF
    scf_env%iter_delta = 0.0_dp

    DO ispin=1,nspin
      CALL eigensolver_simple(matrix_ks=scf_env%scf_work1(ispin)%matrix,&
                              mo_set=mos(ispin)%mo_set,&
                              work=scf_env%scf_work2,&
                              do_level_shift=do_level_shift,&
                              level_shift=scf_control%level_shift,&
                              use_jacobi=use_jacobi,&
                              jacobi_threshold=scf_control%diagonalization%jacobi_threshold,&
                              error=error)
    END DO

    CALL set_mo_occupation(mo_array=mos,&
                           smear=scf_control%smear,&
                           error=error)

    DO ispin=1,nspin
      CALL calculate_density_matrix(mos(ispin)%mo_set,&
                                    scf_env%p_mix_new(ispin)%matrix,&
                                    error=error)
    END DO

  END SUBROUTINE do_special_diag

! *****************************************************************************
!> \brief the inner loop of scf, specific to iterative diagonalization using OT
!>        with S matrix; basically, in goes the ks matrix out goes a new p matrix
!> \param scf_env ...
!> \param mos ...
!> \param matrix_ks ...
!> \param matrix_s ...
!> \param scf_control ...
!> \param scf_section ...
!> \param diis_step ...
!> \param error ...
!> \par History
!>      10.2008 created [JGH]
! *****************************************************************************
  SUBROUTINE do_ot_diag(scf_env,mos,matrix_ks,matrix_s,&
               scf_control,scf_section,diis_step,error)

    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_s
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(section_vals_type), POINTER         :: scf_section
    LOGICAL, INTENT(INOUT)                   :: diis_step
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'do_ot_diag', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: homo, ispin, nmo, nspin
    LOGICAL                                  :: do_level_shift
    REAL(kind=dp)                            :: diis_error, eps_iter
    REAL(KIND=dp), DIMENSION(:), POINTER     :: eigenvalues
    TYPE(cp_fm_type), POINTER                :: mo_coeff

    NULLIFY (eigenvalues)

    nspin = SIZE(matrix_ks)

    DO ispin=1,nspin
      CALL copy_dbcsr_to_fm(matrix_ks(ispin)%matrix,&
                            scf_env%scf_work1(ispin)%matrix,&
                            error=error)
    END DO

    IF ((scf_env%iter_count > 1).AND.(.NOT.scf_env%skip_diis)) THEN
      CALL qs_diis_b_step(scf_env%scf_diis_buffer,mos,scf_env%scf_work1,&
                          scf_env%scf_work2,scf_env%iter_delta,diis_error,diis_step,&
                          scf_control%eps_diis,scf_control%nmixing,&
                          s_matrix=matrix_s,&
                          scf_section=scf_section,error=error)
    ELSE
      diis_step = .FALSE.
    END IF

    do_level_shift = ((scf_control%level_shift /= 0.0_dp).AND.&
                      ((scf_control%density_guess == core_guess).OR.&
                       (scf_env%iter_count > 1)))

    eps_iter = scf_control%diagonalization%eps_iter
    IF (diis_step) THEN
      scf_env%iter_param = diis_error
      scf_env%iter_method = "DIIS/OTdiag"
      DO ispin=1,nspin
         CALL copy_fm_to_dbcsr(scf_env%scf_work1(ispin)%matrix,&
                               matrix_ks(ispin)%matrix,keep_sparsity=.TRUE.,&
                               error=error)
      END DO
      eps_iter = MAX(eps_iter,scf_control%diagonalization%eps_adapt*diis_error)
    ELSE
      IF(scf_env%mixing_method ==1) THEN
        scf_env%iter_param = scf_env%p_mix_alpha
        scf_env%iter_method = "P_Mix/OTdiag."
      ELSEIF(scf_env%mixing_method>1) THEN
        scf_env%iter_param = scf_env%mixing_store%alpha
        scf_env%iter_method = TRIM(scf_env%mixing_store%iter_method)//"/OTdiag."
      END IF
    END IF

    scf_env%iter_delta = 0.0_dp

    DO ispin=1,nspin
      CALL get_mo_set(mos(ispin)%mo_set,&
                      mo_coeff=mo_coeff,&
                      eigenvalues=eigenvalues,&
                      nmo=nmo,&
                      homo=homo)
      CALL ot_eigensolver(matrix_h=matrix_ks(ispin)%matrix,&
                          matrix_s=matrix_s(1)%matrix,&
                          matrix_c_fm=mo_coeff,&
                          preconditioner=scf_env%ot_preconditioner(1)%preconditioner,&
                          eps_gradient=eps_iter,&
                          iter_max=scf_control%diagonalization%max_iter,&
                          silent=.TRUE.,&
                          ot_settings=scf_control%diagonalization%ot_settings,&
                          error=error)
      CALL calculate_subspace_eigenvalues(mo_coeff,matrix_ks(ispin)%matrix,&
                                          evals_arg=eigenvalues,&
                                          do_rotation=.TRUE.,&
                                          error=error)
      !MK WRITE(*,*) routinen//' copy_dbcsr_to_fm'
      CALL copy_fm_to_dbcsr(mos(ispin)%mo_set%mo_coeff,&
                            mos(ispin)%mo_set%mo_coeff_b,&
                            error=error)!fm->dbcsr
    END DO

    CALL set_mo_occupation(mo_array=mos,&
                           smear=scf_control%smear,&
                           error=error)

    DO ispin=1,nspin
      CALL calculate_density_matrix(mos(ispin)%mo_set,&
                                    scf_env%p_mix_new(ispin)%matrix,&
                                    error=error)
    END DO

  END SUBROUTINE do_ot_diag

! *****************************************************************************
!> \brief Solve a set restricted open Kohn-Sham (ROKS) equations based on the
!>         alpha and beta Kohn-Sham matrices from unrestricted Kohn-Sham.
!> \param scf_env ...
!> \param mos ...
!> \param matrix_ks ...
!> \param matrix_s ...
!> \param scf_control ...
!> \param scf_section ...
!> \param diis_step ...
!> \param orthogonal_basis ...
!> \param error ...
!> \par History
!>      04.2006 created [MK]
!>      Revised (01.05.06,MK)
!> \note
!>         this is only a high-spin ROKS.
! *****************************************************************************
  SUBROUTINE do_roks_diag(scf_env,mos,matrix_ks,matrix_s,&
                                      scf_control,scf_section,diis_step,&
                                      orthogonal_basis,error)

    ! Literature: - C. C. J. Roothaan, Rev. Mod. Phys. 32, 179 (1960)
    !             - M. F. Guest and V. R. Saunders, Mol. Phys. 28(3), 819 (1974)
    !             - M. Filatov and S. Shaik, Chem. Phys. Lett. 288, 689 (1998)

    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_s
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(section_vals_type), POINTER         :: scf_section
    LOGICAL, INTENT(INOUT)                   :: diis_step
    LOGICAL, INTENT(IN)                      :: orthogonal_basis
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'do_roks_diag', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, homoa, homob, imo, &
                                                nalpha, nao, nbeta, nmo
    REAL(KIND=dp)                            :: diis_error, level_shift_loc
    REAL(KIND=dp), DIMENSION(:), POINTER     :: eiga, eigb, occa, occb
    TYPE(cp_fm_type), POINTER                :: ksa, ksb, mo2ao, moa, mob, &
                                                ortho, work

! -------------------------------------------------------------------------

    CALL timeset(routineN,handle)

    IF(scf_env%cholesky_method==cholesky_inverse) THEN
      ortho => scf_env%ortho_m1
    ELSE
      ortho => scf_env%ortho
    END IF
    work  => scf_env%scf_work2

    ksa => scf_env%scf_work1(1)%matrix
    ksb => scf_env%scf_work1(2)%matrix

    CALL copy_dbcsr_to_fm(matrix_ks(1)%matrix,ksa,error=error)
    CALL copy_dbcsr_to_fm(matrix_ks(2)%matrix,ksb,error=error)

    ! Get MO information

    CALL get_mo_set(mo_set=mos(1)%mo_set,&
                    nao=nao,&
                    nmo=nmo,&
                    nelectron=nalpha,&
                    homo=homoa,&
                    eigenvalues=eiga,&
                    occupation_numbers=occa,&
                    mo_coeff=moa)

    CALL get_mo_set(mo_set=mos(2)%mo_set,&
                    nelectron=nbeta,&
                    homo=homob,&
                    eigenvalues=eigb,&
                    occupation_numbers=occb,&
                    mo_coeff=mob)

    ! Define the amount of level-shifting

    IF ((scf_control%level_shift /= 0.0_dp).AND.&
        ((scf_control%density_guess == core_guess).OR.&
         (scf_control%density_guess == restart_guess).OR.&
         (scf_env%iter_count > 1))) THEN
      level_shift_loc = scf_control%level_shift
    ELSE
      level_shift_loc = 0.0_dp
    END IF

    IF ((scf_env%iter_count > 1).OR.&
        (scf_control%density_guess == core_guess).OR.&
        (scf_control%density_guess == restart_guess)) THEN

      ! Transform the spin unrestricted alpha and beta Kohn-Sham matrices
      ! from AO basis to MO basis: K(MO) = C(T)*K(AO)*C

      CALL cp_fm_symm("L","U",nao,nao,1.0_dp,ksa,moa,0.0_dp,work,error=error)
      CALL cp_gemm("T","N",nao,nao,nao,1.0_dp,moa,work,0.0_dp,ksa,error=error)

      CALL cp_fm_symm("L","U",nao,nao,1.0_dp,ksb,moa,0.0_dp,work,error=error)
      CALL cp_gemm("T","N",nao,nao,nao,1.0_dp,moa,work,0.0_dp,ksb,error=error)

      ! Combine the spin unrestricted alpha and beta Kohn-Sham matrices
      ! in the MO basis

      IF (scf_control%roks_scheme == general_roks) THEN
        CALL combine_ks_matrices(ksa,ksb,occa,occb,scf_control%roks_f,&
                                 nalpha,nbeta,error=error)
      ELSE IF (scf_control%roks_scheme == high_spin_roks) THEN
        CALL combine_ks_matrices(ksa,ksb,occa,occb,scf_control%roks_parameter,&
                                 error=error)
      ELSE
        CALL stop_program(routineN,moduleN,__LINE__,&
                          "Unknown ROKS scheme requested")
      END IF

      ! Back-transform the restricted open Kohn-Sham matrix from MO basis
      ! to AO basis

      IF (orthogonal_basis) THEN
        ! Q = C
        mo2ao => moa
      ELSE
        ! Q = S*C
        mo2ao => mob
!MK     CALL copy_sm_to_fm(matrix_s(1)%matrix,work)
!MK     CALL cp_fm_symm("L","U",nao,nao,1.0_dp,work,moa,0.0_dp,mo2ao)
        CALL cp_dbcsr_sm_fm_multiply(matrix_s(1)%matrix,moa,mo2ao,nao,error=error)
      END IF

      ! K(AO) = Q*K(MO)*Q(T)

      CALL cp_gemm("N","T",nao,nao,nao,1.0_dp,ksa,mo2ao,0.0_dp,work,error=error)
      CALL cp_gemm("N","N",nao,nao,nao,1.0_dp,mo2ao,work,0.0_dp,ksa,error=error)

    ELSE

      ! No transformation matrix available, yet. The closed shell part,
      ! i.e. the beta Kohn-Sham matrix in AO basis, is taken.
      ! There might be better choices, anyhow.

      CALL cp_fm_to_fm(ksb,ksa,error=error)

    END IF

    ! Update DIIS buffer and possibly perform DIIS extrapolation step

    IF (scf_env%iter_count > 1) THEN
      IF (orthogonal_basis) THEN
        CALL qs_diis_b_step(diis_buffer=scf_env%scf_diis_buffer,&
                            mo_array=mos,&
                            kc=scf_env%scf_work1,&
                            sc=work,&
                            delta=scf_env%iter_delta,&
                            error_max=diis_error,&
                            diis_step=diis_step,&
                            eps_diis=scf_control%eps_diis,&
                            scf_section=scf_section,&
                            roks=.TRUE.,&
                            error=error)
      ELSE
        CALL qs_diis_b_step(diis_buffer=scf_env%scf_diis_buffer,&
                            mo_array=mos,&
                            kc=scf_env%scf_work1,&
                            sc=work,&
                            delta=scf_env%iter_delta,&
                            error_max=diis_error,&
                            diis_step=diis_step,&
                            eps_diis=scf_control%eps_diis,&
                            scf_section=scf_section,&
                            s_matrix=matrix_s,&
                            roks=.TRUE.,&
                            error=error)
      END IF
    END IF

    IF (diis_step) THEN
      scf_env%iter_param = diis_error
      scf_env%iter_method = "DIIS/Diag."
    ELSE
      IF(scf_env%mixing_method ==1) THEN
        scf_env%iter_param = scf_env%p_mix_alpha
        scf_env%iter_method = "P_Mix/Diag."
      ELSEIF(scf_env%mixing_method>1) THEN
        scf_env%iter_param = scf_env%mixing_store%alpha
        scf_env%iter_method = TRIM(scf_env%mixing_store%iter_method)//"/Diag."
      END IF
    END IF

    scf_env%iter_delta = 0.0_dp

    IF (level_shift_loc /= 0.0_dp) THEN

      ! Transform the current Kohn-Sham matrix from AO to MO basis
      ! for level-shifting using the current MO set

      CALL cp_fm_symm("L","U",nao,nao,1.0_dp,ksa,moa,0.0_dp,work,error=error)
      CALL cp_gemm("T","N",nao,nao,nao,1.0_dp,moa,work,0.0_dp,ksa,error=error)

      ! Apply level-shifting using 50:50 split of the shift (could be relaxed)

      DO imo=homob+1,homoa
        CALL cp_fm_add_to_element(ksa,imo,imo,0.5_dp*level_shift_loc,error)
      END DO
      DO imo=homoa+1,nmo
        CALL cp_fm_add_to_element(ksa,imo,imo,level_shift_loc,error)
      END DO

    ELSE IF (.NOT.orthogonal_basis) THEN

      ! Transform the current Kohn-Sham matrix to an orthogonal basis
      SELECT CASE(scf_env%cholesky_method)
      CASE(cholesky_reduce)
        CALL cp_fm_cholesky_reduce(ksa,ortho,error=error)
      CASE(cholesky_restore)
         CALL cp_fm_upper_to_full(ksa,work,error=error)
         CALL cp_fm_cholesky_restore(ksa,nao,ortho,work,&
                "SOLVE",pos="RIGHT",error=error)
         CALL cp_fm_cholesky_restore(work,nao,ortho,ksa,&
               "SOLVE",pos="LEFT",transa="T",error=error)
      CASE(cholesky_inverse)
         CALL cp_fm_upper_to_full(ksa,work,error=error)
         CALL cp_fm_cholesky_restore(ksa,nao,ortho,work,&
                "MULTIPLY",pos="RIGHT",error=error)
         CALL cp_fm_cholesky_restore(work,nao,ortho,ksa,&
               "MULTIPLY",pos="LEFT",transa="T",error=error)
      CASE(cholesky_off)
        CALL cp_fm_symm("L","U",nao,nao,1.0_dp,ksa,ortho,0.0_dp,work,error=error)
        CALL cp_gemm("N","N",nao,nao,nao,1.0_dp,ortho,work,0.0_dp,ksa,error=error)
      END SELECT

    END IF

    ! Diagonalization of the ROKS operator matrix

     CALL choose_eigv_solver(ksa,work,eiga,error=error)

    ! Back-transformation of the orthonormal eigenvectors if needed

    IF (level_shift_loc /= 0.0_dp) THEN
      ! Use old MO set for back-transformation if level-shifting was applied
      CALL cp_fm_to_fm(moa,ortho,error=error)
      CALL cp_gemm("N","N",nao,nmo,nao,1.0_dp,ortho,work,0.0_dp,moa,error=error)
    ELSE
      IF (orthogonal_basis) THEN
        CALL cp_fm_to_fm(work,moa,error=error)
      ELSE
        SELECT CASE(scf_env%cholesky_method)
        CASE(cholesky_reduce,cholesky_restore)
          CALL cp_fm_cholesky_restore(work,nmo,ortho,moa,"SOLVE",error=error)
        CASE(cholesky_inverse)
          CALL cp_fm_cholesky_restore(work,nmo,ortho,moa,"MULTIPLY",error=error)
        CASE(cholesky_off)
          CALL cp_gemm("N","N",nao,nmo,nao,1.0_dp,ortho,work,0.0_dp,moa,error=error)
        END SELECT
      END IF
    END IF

    ! Correct MO eigenvalues, if level-shifting was applied

    IF (level_shift_loc /= 0.0_dp) THEN
      DO imo=homob+1,homoa
        eiga(imo) = eiga(imo) - 0.5_dp*level_shift_loc
      END DO
      DO imo=homoa+1,nmo
        eiga(imo) = eiga(imo) - level_shift_loc
      END DO
    END IF

    ! Update also the beta MO set

    eigb(:) = eiga(:)
    CALL cp_fm_to_fm(moa,mob,error=error)

    ! Calculate the new alpha and beta density matrix

    CALL calculate_density_matrix(mos(1)%mo_set,scf_env%p_mix_new(1)%matrix,&
                                  error=error)
    CALL calculate_density_matrix(mos(2)%mo_set,scf_env%p_mix_new(2)%matrix,&
                                  error=error)

    CALL timestop(handle)

  END SUBROUTINE do_roks_diag

! *****************************************************************************
!> \brief iterative diagonalization using the block Krylov-space approach
!> \param scf_env ...
!> \param mos ...
!> \param matrix_ks ...
!> \param scf_control ...
!> \param scf_section ...
!> \param check_moconv_only ...
!> \param error ...
!> \param
!> \par History
!>      05.2009 created [MI]
! *****************************************************************************

  SUBROUTINE do_block_krylov_diag(scf_env,mos,matrix_ks,&
                     scf_control, scf_section, check_moconv_only, error)

    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(section_vals_type), POINTER         :: scf_section
    LOGICAL, INTENT(IN), OPTIONAL            :: check_moconv_only
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'do_block_krylov_diag', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: rone = 1.0_dp, rzero = 0.0_dp

    INTEGER                                  :: handle, homo, ispin, iter, &
                                                nao, nmo, output_unit
    LOGICAL                                  :: converged, failure, &
                                                my_check_moconv_only
    REAL(dp)                                 :: eps_iter, t1, t2
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mo_eigenvalues
    TYPE(cp_fm_type), POINTER                :: c0, c1, chc, evec, ks, &
                                                mo_coeff, ortho, work
    TYPE(cp_logger_type), POINTER            :: logger

    failure = .FALSE.
    logger => cp_error_get_logger(error)
    CALL timeset(routineN,handle)

    output_unit = cp_print_key_unit_nr(logger,scf_section,"PRINT%LANCZOS",&
            extension=".scfLog",error=error)


    my_check_moconv_only = .FALSE.
    IF(PRESENT(check_moconv_only)) my_check_moconv_only=check_moconv_only

    NULLIFY (mo_coeff,ortho,work, ks)
    NULLIFY (mo_eigenvalues)
    NULLIFY (c0, c1)

    IF(scf_env%cholesky_method==cholesky_inverse) THEN
      ortho => scf_env%ortho_m1
    ELSE
      ortho => scf_env%ortho
    END IF
    work  => scf_env%scf_work2

    DO ispin=1,SIZE(matrix_ks)
      CALL copy_dbcsr_to_fm(matrix_ks(ispin)%matrix,&
                         scf_env%scf_work1(ispin)%matrix,&
                         error=error)
    END DO

    IF(scf_env%mixing_method ==1) THEN
        scf_env%iter_param = scf_env%p_mix_alpha
        scf_env%iter_method = "P_Mix/Lanczos"
    ELSE
!        scf_env%iter_param = scf_env%mixing_store%alpha
        scf_env%iter_method = TRIM(scf_env%mixing_store%iter_method)//"/Lanc."
    END IF

    DO ispin=1,SIZE(matrix_ks)

      ks => scf_env%scf_work1(ispin)%matrix
      CALL cp_fm_upper_to_full(ks,work,error=error)

      CALL get_mo_set(mo_set=mos(ispin)%mo_set,&
                      nao=nao,&
                      nmo=nmo,&
                      homo=homo,&
                      eigenvalues=mo_eigenvalues,&
                      mo_coeff=mo_coeff)

      NULLIFY(c0,c1)
      c0 => scf_env%krylov_space%mo_conv(ispin)%matrix
      c1 => scf_env%krylov_space%mo_refine(ispin)%matrix
      SELECT CASE(scf_env%cholesky_method)
      CASE(cholesky_reduce)
         CALL cp_fm_cholesky_reduce(ks,ortho,error=error)
         CALL cp_fm_upper_to_full(ks,work,error=error)
         CALL cp_fm_cholesky_restore(mo_coeff,nmo,ortho,c0,"MULTIPLY",error=error)
      CASE(cholesky_restore)
         CALL cp_fm_cholesky_restore(ks,nao,ortho,work,&
                "SOLVE",pos="RIGHT",error=error)
         CALL cp_fm_cholesky_restore(work,nao,ortho,ks,&
               "SOLVE",pos="LEFT",transa="T",error=error)
         CALL cp_fm_cholesky_restore(mo_coeff,nmo,ortho,c0,"MULTIPLY",error=error)
      CASE(cholesky_inverse)
         CALL cp_fm_cholesky_restore(ks,nao,ortho,work,&
                "MULTIPLY",pos="RIGHT",error=error)
         CALL cp_fm_cholesky_restore(work,nao,ortho,ks,&
               "MULTIPLY",pos="LEFT",transa="T",error=error)
         CALL cp_fm_cholesky_restore(mo_coeff,nmo,ortho,c0,"SOLVE",error=error)
      END SELECT

      scf_env%krylov_space%nmo_nc  = nmo
      scf_env%krylov_space%nmo_conv  = 0

      t1 = m_walltime()
      IF (output_unit > 0) THEN
          WRITE(output_unit,"(/T15,A)") '<<<<<<<<<   LANCZOS REFINEMENT    <<<<<<<<<<'
          WRITE(output_unit,"(T8,A,T15,A,T23,A,T36,A,T49,A,T60,A,/,T8,A)")   &
               " Spin ", " Cycle ", &
               " conv. MOS ", " B2MAX ", " B2MIN ", " Time",  REPEAT("-",60)
      END IF
      eps_iter = MAX(scf_env%krylov_space%eps_conv, scf_env%krylov_space%eps_adapt*scf_env%iter_delta)
      iter = 0
      converged = .FALSE.
      !Check convergence of MOS
      IF(my_check_moconv_only) THEN

         CALL lanczos_refinement(scf_env%krylov_space, ks, c0, c1, mo_eigenvalues,&
              nao,  eps_iter, ispin, check_moconv_only=my_check_moconv_only, error=error)
         t2 = m_walltime()
         IF(output_unit > 0)&
               WRITE(output_unit,'(T8,I3,T16,I5,T24,I6,T33,E12.4,2x,E12.4,T60,F8.3)')  &
               ispin,  iter, scf_env%krylov_space%nmo_conv, &
               scf_env%krylov_space%max_res_norm, scf_env%krylov_space%min_res_norm, t2-t1

         CYCLE
      ELSE
         !Block Lanczos refinement
         DO iter = 1,scf_env%krylov_space%max_iter
           CALL lanczos_refinement_2v(scf_env%krylov_space, ks, c0, c1, mo_eigenvalues,&
                 nao, eps_iter, ispin,  error=error)
            t2 = m_walltime()
            IF (output_unit > 0) THEN
               WRITE(output_unit,'(T8,I3,T16,I5,T24,I6,T33,E12.4,2x,E12.4,T60,F8.3)')&
                      ispin, iter, scf_env%krylov_space%nmo_conv, &
                      scf_env%krylov_space%max_res_norm, scf_env%krylov_space%min_res_norm, t2-t1
            END IF
            t1=m_walltime()
            IF(scf_env%krylov_space%max_res_norm < eps_iter) THEN
              converged = .TRUE.
              IF (output_unit > 0) WRITE(output_unit,*)&
                   " Reached convergence in ",iter," iterations "
              EXIT
            END IF
         END DO

         IF(.NOT. converged .AND. output_unit > 0) THEN
           WRITE(output_unit,"(T4, A)") " WARNING Lanczos refinement could "//&
                "not converge all the mos:"
           WRITE(output_unit,"(T40,A,T70,I10)")    " number of not converged mos ",&
                 scf_env%krylov_space%nmo_nc
           WRITE(output_unit,"(T40,A,T70,E10.2)")  " max norm of the residual " , &
                 scf_env%krylov_space%max_res_norm

         END IF

         ! For the moment skip the re-orthogonalization
         IF(.FALSE.) THEN
           !Re-orthogonalization
           NULLIFY( chc, evec)
           chc => scf_env%krylov_space%chc_mat(ispin)%matrix
           evec => scf_env%krylov_space%c_vec(ispin)%matrix
           CALL cp_gemm('N','N',nao,nmo,nao,rone,ks,c0,rzero,work,error=error)
           CALL cp_gemm('T','N',nmo,nmo,nao,rone,c0,work,rzero,chc,error=error)
           !Diagonalize  (C^t)HC
           CALL choose_eigv_solver(chc,evec,mo_eigenvalues,error=error)
           !Rotate the C vectors
           CALL cp_gemm('N','N',nao,nmo,nmo,rone,c0,evec,rzero,c1,error=error)
           c0 => scf_env%krylov_space%mo_refine(ispin)%matrix
         END IF

         IF(scf_env%cholesky_method==cholesky_inverse) THEN
           CALL cp_fm_cholesky_restore(c0,nmo,ortho,mo_coeff,"MULTIPLY",error=error)
         ELSE
           CALL cp_fm_cholesky_restore(c0,nmo,ortho,mo_coeff,"SOLVE",error=error)
         END IF

         CALL set_mo_occupation(mo_set=mos(ispin)%mo_set,&
                                smear=scf_control%smear,&
                                error=error)

         CALL calculate_density_matrix(mos(ispin)%mo_set,&
                                       scf_env%p_mix_new(ispin)%matrix,&
                                       error=error)
      END  IF
    END DO ! ispin

    IF (output_unit > 0) THEN
          WRITE(output_unit,"(T15,A/)") '<<<<<<<<< END LANCZOS REFINEMENT  <<<<<<<<<<'
    END IF

    CALL cp_print_key_finished_output(output_unit,logger,scf_section,&
            "PRINT%LANCZOS", error=error)

    CALL timestop(handle)

  END SUBROUTINE do_block_krylov_diag

! *****************************************************************************
!> \brief iterative diagonalization using the block davidson space approach
!> \param qs_env ...
!> \param scf_env ...
!> \param mos ...
!> \param matrix_ks ...
!> \param matrix_s ...
!> \param scf_control ...
!> \param scf_section ...
!> \param check_moconv_only ...
!> \param error ...
!> \param
!> \par History
!>      05.2011 created [MI]
! *****************************************************************************

  SUBROUTINE do_block_davidson_diag (qs_env,scf_env,mos,matrix_ks,matrix_s,&
                     scf_control,scf_section, check_moconv_only, error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_s
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(section_vals_type), POINTER         :: scf_section
    LOGICAL, INTENT(IN), OPTIONAL            :: check_moconv_only
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'do_block_davidson_diag', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: rone = 1.0_dp, rzero = 0.0_dp

    INTEGER                                  :: handle, ispin, nspins, &
                                                output_unit
    LOGICAL                                  :: do_prec, failure, &
                                                my_check_moconv_only
    TYPE(cp_logger_type), POINTER            :: logger

    failure = .FALSE.
    logger => cp_error_get_logger(error)
    CALL timeset(routineN,handle)

    output_unit = cp_print_key_unit_nr(logger,scf_section,"PRINT%DAVIDSON",&
            extension=".scfLog",error=error)

    IF(output_unit > 0) &
        WRITE(output_unit ,"(/T15,A)") '<<<<<<<<<  DAVIDSON ITERATIONS   <<<<<<<<<<'

    IF(scf_env%mixing_method ==1) THEN
        scf_env%iter_param = scf_env%p_mix_alpha
        scf_env%iter_method = "P_Mix/Dav."
    ELSE
        scf_env%iter_param = scf_env%mixing_store%alpha
        scf_env%iter_method = TRIM(scf_env%mixing_store%iter_method)//"/Dav."
    END IF


    my_check_moconv_only = .FALSE.
    IF(PRESENT(check_moconv_only)) my_check_moconv_only=check_moconv_only
    do_prec=.FALSE.
    IF( scf_env%block_davidson_env(1)%prec_type/=0 .AND.&
        scf_env%iter_count>=scf_env%block_davidson_env(1)%first_prec)THEN
       do_prec=.TRUE.
    END IF

    nspins = SIZE(matrix_ks)

    IF( do_prec .AND. (scf_env%iter_count==scf_env%block_davidson_env(1)%first_prec  .OR. &
          MODULO(scf_env%iter_count,scf_env%block_davidson_env(1)%niter_new_prec)==0)) THEN
!dbg
!  write(*,*) '++++  RECOMPUTE PRECONDITIONER ++++++++'
!dbg
       CALL restart_preconditioner(qs_env,scf_env%ot_preconditioner,&
            prec_type=scf_env%block_davidson_env(1)%prec_type,nspins=nspins,error=error)
       CALL prepare_preconditioner(qs_env,mos,matrix_ks,matrix_s,scf_env%ot_preconditioner,&
            scf_env%block_davidson_env(1)%prec_type, &
            scf_env%block_davidson_env(1)%solver_type,&
            scf_env%block_davidson_env(1)%energy_gap,nspins,error=error)
    END IF

    DO ispin = 1,nspins
      IF(scf_env%block_davidson_env(ispin)%use_sparse_mos) THEN
        IF(.NOT. do_prec) THEN
          CALL generate_extended_space_sparse(scf_env%block_davidson_env(ispin),mos(ispin)%mo_set,&
             matrix_ks(ispin)%matrix,matrix_s(1)%matrix, output_unit,&
             error=error)
        ELSE
          CALL generate_extended_space_sparse(scf_env%block_davidson_env(ispin),mos(ispin)%mo_set,&
             matrix_ks(ispin)%matrix,matrix_s(1)%matrix,output_unit,&
             scf_env%ot_preconditioner(ispin)%preconditioner,error=error)
        END IF

      ELSE
        IF(.NOT. do_prec) THEN
          CALL generate_extended_space(scf_env%block_davidson_env(ispin),mos(ispin)%mo_set,&
             matrix_ks(ispin)%matrix,matrix_s(1)%matrix, output_unit,&
             error=error)
        ELSE
          CALL generate_extended_space(scf_env%block_davidson_env(ispin),mos(ispin)%mo_set,&
             matrix_ks(ispin)%matrix,matrix_s(1)%matrix,output_unit,&
             scf_env%ot_preconditioner(ispin)%preconditioner,error=error)
        END IF
      END IF
    END DO !ispin

    CALL set_mo_occupation(mo_array=mos,&
                             smear=scf_control%smear,&
                             error=error)

    DO ispin=1,nspins
      CALL calculate_density_matrix(mos(ispin)%mo_set,&
                                    scf_env%p_mix_new(ispin)%matrix,&
                                    error=error)
    END DO


    IF (output_unit > 0) THEN
          WRITE(output_unit,"(T15,A/)") '<<<<<<<<< END DAVIDSON ITERATION  <<<<<<<<<<'
    END IF

    CALL cp_print_key_finished_output(output_unit,logger,scf_section,&
            "PRINT%DAVIDSON", error=error)

    CALL timestop(handle)

  END SUBROUTINE  do_block_davidson_diag

END MODULE qs_scf_diagonalization
