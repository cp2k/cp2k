!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2001 - 2002  CP2K developers group
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****** cp2k/qs_scf [1.0] *
!!
!!   NAME
!!     qs_scf
!!
!!   FUNCTION
!!     Routines for the Quickstep SCF run.
!!
!!   AUTHOR
!!     Matthias Krack (30.04.2001)
!!
!!   MODIFICATION HISTORY
!!     - Joost VandeVondele (02.2002)
!!          added code for: incremental (pab and gvg) update
!!                           initialisation (init_cube, l_info)
!!     - Joost VandeVondele (02.2002)
!!          called the poisson code of the classical part
!!          this takes into account the spherical cutoff and allows for
!!          isolated systems
!!     - Joost VandeVondele (02.2002)
!!          added multiple grid feature
!!          changed to spherical cutoff consistently (?)
!!          therefore removed the gradient correct functionals
!!     - updated with the new QS data structures (10.04.02,MK)
!!     - copy_matrix replaced by transfer_matrix (11.04.02,MK)
!!     - nrebuild_rho and nrebuild_gvg unified (12.04.02,MK)
!!     - set_mo_occupation for smearing of the MO occupation numbers
!!       (17.04.02,MK)
!!     - MO level shifting added (22.04.02,MK)
!!     - Usage of TYPE mo_set_p_type
!!     - Joost VandeVondele (05.2002)
!!           added cholesky based diagonalisation
!!     - 05.2002 added pao method [fawzi]
!!     - parallel FFT (JGH 22.05.2002)
!!     - 06.2002 moved KS matrix construction to qs_build_KS_matrix.F [fawzi]
!!     - started to include more LSD (01.2003,Joost VandeVondele)
!!
!!   SOURCE
!******************************************************************************

MODULE qs_scf
! *****************************************************************************
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE checkpoint_handler,              ONLY: write_checkpoint_information
  USE coefficient_types,               ONLY: coeff_allocate,&
                                             coeff_copy,&
                                             coeff_sumup,&
                                             coeff_transform_space,&
                                             coeff_type,&
                                             coeff_zero
  USE cp_array_utils,                  ONLY: cp_1d_r_output
  USE cp_block_matrix,                 ONLY: cp_block_matrix_p_type,&
                                             cp_block_matrix_type,&
                                             cp_sp_multiplies_blacs
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_debug,&
                                             cp_error_dealloc_ref,&
                                             cp_error_get_logger,&
                                             cp_error_init,&
                                             cp_error_message,&
                                             cp_error_type,&
                                             cp_internal_error
  USE cp_fm_pool_types,                ONLY: cp_fm_pool_p_type,&
                                             cp_fm_pool_type,&
                                             fm_pool_get_mstruct
  USE cp_fm_struct,                    ONLY: cp_fm_struct_get,&
                                             cp_fm_struct_type
  USE cp_lbfgs_optimizer_gerd_f,       ONLY: cp_lbfgs_opt_gerd_f_type,&
                                             cp_opt_gerd_f_dealloc_ref,&
                                             cp_opt_gerd_f_init,&
                                             cp_opt_gerd_f_next
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_log,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_matrix_utils,                 ONLY: cp_sm_output
  USE cube_utils,                      ONLY: cube_info_type,&
                                             destroy_cube_info,&
                                             init_cube_info
  USE dft_types,                       ONLY: dft_control_type
  USE external_potential_types,        ONLY: all_potential_type,&
                                             get_potential,&
                                             gth_potential_type
  USE global_types,                    ONLY: DEBUG,&
                                             global_environment_type
  USE hartree,                         ONLY: calculate_hartree
  USE kinds,                           ONLY: int_size,&
                                             wp => dp
  USE l_utils,                         ONLY: destroy_l_info,&
                                             init_l_info,&
                                             l_info_type
  USE machine,                         ONLY: m_cputime,&
                                             m_flush
  USE mathconstants,                   ONLY: fourpi,&
                                             twopi
  USE message_passing,                 ONLY: mp_allgather,&
                                             mp_max,&
                                             mp_range,&
                                             mp_sum,&
                                             mp_sync
  USE pao_glob_angles_methods,         ONLY: pao_g_ang_get,&
                                             pao_g_ang_set
  USE pao_l_angles_methods,            ONLY: pao_guarantee_owned_NUi
  USE pao_obj_function,                ONLY: pao_gerd_f_dealloc_ref,&
                                             pao_gerd_f_init,&
                                             pao_gerd_functional_type
  USE pao_proj_methods,                ONLY: sm_matrix_transf_p_to_m,&
                                             sm_matrix_transf_to_f
  USE pao_qs_env_methods,              ONLY: pao_qs_env_did_change
  USE pao_types,                       ONLY: pao_env_get,&
                                             pao_env_type,&
                                             pao_glob_angles_type
  USE particle_types,                  ONLY: particle_type
  USE preconditioner,                  ONLY: destroy_preconditioner,&
                                             init_preconditioner,&
                                             make_preconditioner,&
                                             preconditioner_type
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_pool_types,                   ONLY: pw_pool_give_back_coeff,&
                                             pw_pool_init_coeff,&
                                             pw_pool_p_type,&
                                             pw_pool_type,&
                                             pw_pools_give_back_coeffs,&
                                             pw_pools_init_coeffs
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             COMPLEXDATA3D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_type
  USE qs_arpack,                       ONLY: arpack_diag_sy
  USE qs_arpack_full,                  ONLY: op_symv_type,&
                                             setup_op_symv
  USE qs_blacs,                        ONLY: blacs_add_to_element,&
                                             blacs_block_jacobi_classic,&
                                             blacs_get_element,&
                                             blacs_make_basis,&
                                             blacs_syevx,&
                                             copy_blacs_to_blacs_matrix,&
                                             copy_full_matrix_to_real_matrix,&
                                             copy_real_matrix_to_full_matrix,&
                                             cp_fm_add,&
                                             cp_fm_cholesky_decompose,&
                                             cp_fm_cholesky_reduce,&
                                             cp_fm_cholesky_restore,&
                                             cp_fm_create2,&
                                             cp_fm_gemm,&
                                             cp_fm_get_info,&
                                             cp_fm_init_random,&
                                             cp_fm_release,&
                                             cp_fm_set_all,&
                                             cp_fm_set_element,&
                                             cp_fm_symm,&
                                             cp_fm_trace,&
                                             cp_full_matrix_p_type,&
                                             cp_full_matrix_type,&
                                             cp_sm_fm_multiply,&
                                             power_blacs_matrix,&
                                             replicate_blacs_matrix,&
                                             symmetrise_blacs_matrix,&
                                             w_blacs_matrix => write_blacs_matrix
  USE qs_charges_types,                ONLY: qs_charges_type
  USE qs_collocate_density,            ONLY: calculate_rho_core,&
                                             calculate_total_rho,&
                                             calculate_wavefunction
  USE qs_core_energies,                ONLY: calculate_ecore,&
                                             calculate_ecore_overlap,&
                                             calculate_ecore_self
  USE qs_core_hamiltonian,             ONLY: build_core_hamiltonian_matrix
  USE qs_diis,                         ONLY: scf_diis
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_ks_methods,                   ONLY: qs_ks_did_change,&
                                             qs_ks_update_qs_env
  USE qs_ks_types,                     ONLY: qs_ks_env_p_type,&
                                             qs_ks_env_type,&
                                             qs_ks_release,&
                                             qs_ks_retain
  USE qs_mo_types,                     ONLY: allocate_mo_set,&
                                             calculate_density_matrix,&
                                             calculate_subspace_eigenvalues,&
                                             correct_mo_eigenvalues,&
                                             get_mo_set,&
                                             init_mo_set,&
                                             mo_set_p_type,&
                                             mo_set_type,&
                                             read_mo_set,&
                                             set_mo_occupation,&
                                             write_mo_set
  USE qs_ot,                           ONLY: qs_ot_new_preconditioner,&
                                             qs_ot_type
  USE qs_ot_eigensolver,               ONLY: ot_eigensolver
  USE qs_ot_scf,                       ONLY: ot_scf_destroy,&
                                             ot_scf_init,&
                                             ot_scf_mini,&
                                             ot_scf_read_input
  USE qs_overlap,                      ONLY: write_blacs_matrix,&
                                             write_sparse_matrix
  USE qs_parser,                       ONLY: close_file,&
                                             open_file
  USE qs_rho_types,                    ONLY: qs_rho_did_change,&
                                             qs_rho_get,&
                                             qs_rho_type
  USE realspace_grid_types,            ONLY: rs_pw_to_cube
  USE scf_control_types,               ONLY: scf_control_type
  USE simulation_cell,                 ONLY: cell_type,&
                                             get_cell
  USE sparse_matrix_types,             ONLY: add_all_real_matrix_blocks,&
                                             add_matrices,&
                                             allocate_matrix,&
                                             allocate_matrix_set,&
                                             deallocate_matrix,&
                                             first_block_node,&
                                             get_block_node,&
                                             get_matrix_info,&
                                             next_block_node,&
                                             real_block_node_type,&
                                             real_matrix_p_type,&
                                             real_matrix_type,&
                                             replicate_matrix,&
                                             replicate_matrix_structure,&
                                             set_matrix,&
                                             symmetrise_diagonal_blocks,&
                                             transfer_matrix,&
                                             write_blocks_maxabsval
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(LEN=*), PARAMETER :: module_name = "qs_scf"

! *** Global variables ***

  TYPE(cp_full_matrix_p_type), DIMENSION(:), POINTER :: scf_work1
  TYPE(cp_full_matrix_type), POINTER  :: ortho,scf_work2

  TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set
  TYPE(particle_type), DIMENSION(:), POINTER    :: particle_set

  LOGICAL, PARAMETER, PRIVATE :: debug_this_module=.FALSE.

! *** Public subroutines ***

  PUBLIC :: scf

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE scf(qs_env,globenv,error)

    TYPE(global_environment_type), INTENT(IN)  :: globenv
    TYPE(qs_environment_type), INTENT(INOUT)   :: qs_env
    TYPE(cp_error_type),INTENT(inout),OPTIONAL :: error

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "scf"
    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE "//routine_name//" (MODULE "//module_name//")"

!   *** Local variables ***

    TYPE(dft_control_type), POINTER   :: dft_control
    TYPE(mo_set_type), POINTER :: mo_set
    TYPE(real_matrix_p_type), DIMENSION(:), POINTER :: h,s
    REAL(wp) :: delta,diis_error,method_par,pmix,t1,t2,tmp
    INTEGER  :: handle,handle2,homo,igrid_level,imo,iscf,istat,&
                output_unit,n,k,nrow_block,ncol_block, type_switch, &
                fileunit
    LOGICAL  :: diis_step,do_level_shift,ionode, gth_potential_present,&
                use_jacobi,use_cholesky
    TYPE(cp_full_matrix_type), POINTER :: mo_eigenvectors

    LOGICAL :: pao, failure
    TYPE(pao_gerd_functional_type), POINTER :: pao_obj_f
    TYPE(cp_lbfgs_opt_gerd_f_type), POINTER :: pao_optimizer
    TYPE(pao_env_type),POINTER :: pao_env
    TYPE(cp_block_matrix_p_type),DIMENSION(:),POINTER :: min_density_m
    TYPE(cp_full_matrix_type), POINTER :: min_mo,full_mo
    REAL(kind=wp)  ::pao_f,pao_last_f,pao_actual_grad
    
    TYPE(real_matrix_p_type), DIMENSION(:), POINTER :: ks_rmpv, p_rmpv, kinetic, p_new
    TYPE(real_matrix_p_type), POINTER :: ks,p

    TYPE(qs_rho_type), POINTER :: rho
    TYPE(real_matrix_type), POINTER :: pfull
    TYPE(qs_energy_type), POINTER :: energy
    LOGICAL :: energy_only, rebuild, ot_initialized
    TYPE(cp_full_matrix_type), POINTER :: unoccupied_orbs
    TYPE(scf_control_type), POINTER   :: scf_control
    TYPE(cp_logger_type), POINTER :: logger
    TYPE(pw_env_type), POINTER :: pw_env
    TYPE(qs_charges_type), POINTER :: qs_charges
    TYPE(coeff_type) :: wf_r, wf_g
    TYPE(coeff_type), DIMENSION(:), POINTER :: mgrid_rspace, mgrid_gspace
    TYPE(pw_pool_type), POINTER :: auxbas_pw_pool
    TYPE(pw_pool_p_type), DIMENSION(:), POINTER :: pw_pools
    TYPE(qs_ks_env_type), POINTER :: ks_env

    CHARACTER(LEN=6)   :: extension
    CHARACTER(LEN=200) :: file_name
    INTEGER            :: file_unit,ivector,nao,nlumos
    CHARACTER(LEN=60)  :: method
    CHARACTER(LEN=80)  :: filename
    TYPE(preconditioner_type), TARGET :: ot_preconditioner
    TYPE(qs_ot_type)                  :: qs_ot_env
    TYPE(cp_full_matrix_p_type), DIMENSION(:), POINTER :: ot_orbs_to_rotate
    TYPE(cell_type), POINTER          :: cell ! is this one used here 
    INTEGER :: ispin
    TYPE(mo_set_p_type), DIMENSION(:), POINTER    :: c
    INTEGER :: nelectron

!   ---------------------------------------------------------------------------

    NULLIFY(mo_eigenvectors,pao_obj_f,pao_optimizer,pao_env,min_density_m,&
         min_mo,full_mo, ks, p, ks_rmpv, p_rmpv,qs_charges, &
         pfull, energy, unoccupied_orbs,scf_control,kinetic,h,s, ks_env)
    t1 = m_cputime()
    ot_initialized=.FALSE.
    failure=.FALSE.
    delta = 0.0_wp

    CALL get_qs_env(qs_env=qs_env,&
         atomic_kind_set=atomic_kind_set,&
         cell=cell,&
         dft_control=dft_control,&
         h=h,&
         s=s,energy=energy,&
         kinetic=kinetic, &
         particle_set=particle_set,&
         pao_env=pao_env,&
         pw_env=pw_env,&
         qs_charges=qs_charges,&
         ks_env=ks_env,&
         scf_control=scf_control)

!   *** Quick return, if no SCF iteration is requested ***

    IF (scf_control%max_scf < 1) RETURN

    CALL write_checkpoint_information("entering "//routine_name,globenv)

    CALL timeset(routine_name,"I","",handle)

    ionode = globenv%ionode
    output_unit = globenv%scr
    energy_only = .FALSE.
        
    pao=dft_control%qs_control%pao.AND.ASSOCIATED(pao_env)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             gth_potential_present=gth_potential_present,&
                             nelectron=nelectron)

    nelectron = nelectron - dft_control%charge ! total number of electrons

    CALL get_cell(cell=cell)

    IF (ionode.AND.globenv%print%scf) THEN
      WRITE (UNIT=output_unit,FMT="(/,/,T2,A)")&
        "SCF WAVEFUNCTION OPTIMIZATION"
    END IF


    CALL init_scf_run(qs_env=qs_env,&
                      globenv=globenv,&
                      error=error)

    CALL qs_ks_did_change(ks_env,qs_env,p_changed=.TRUE.,error=error)

    CALL pao_opt_create(pao_optimizer,pao_obj_f,qs_env,globenv,error)

    CALL get_qs_env(qs_env,p=p_rmpv, k=ks_rmpv,rho=rho,c=c)

    NULLIFY(p_new)
    ALLOCATE(p_new(dft_control%nspins))
    DO ispin=1,dft_control%nspins
       NULLIFY (p_new(ispin)%matrix)
       CALL replicate_matrix(p_rmpv(ispin)%matrix, p_new(ispin)%matrix,"NEW DENSITY MATRIX")
    ENDDO

    iscf = 0
    diis_step = .FALSE.
    rebuild = .TRUE.
    use_jacobi = .FALSE.

    CALL init_preconditioner(ot_preconditioner,para_env=globenv%para_env,&
         blacs_env=globenv%blacs_env)

    IF (scf_control%use_ot) THEN
       CALL ot_scf_read_input(qs_ot_env,globenv)
    ENDIF


    IF (ionode.AND.globenv%print%scf) THEN
      WRITE (UNIT=output_unit,&
             FMT="(/,T3,A,T9,A,T36,A,T49,A,T68,A,/,T3,A)")&
        "Step","Update method","Time","Convergence","Total energy",&
        REPEAT("-",77)
    END IF

!   *** SCF loop ***


    scf_loop: DO

      CALL timeset("scf_iter","I","",handle2)

      IF (ionode) CALL m_flush(output_unit)

      iscf = iscf + 1

      IF (.NOT. globenv%print%timecumul) t1 = m_cputime()


      CALL qs_ks_update_qs_env(ks_env,qs_env=qs_env,&
           rebuild=rebuild,error=error,&
           calculate_forces=.FALSE.,just_energy=energy_only)

      DO ispin=1,dft_control%nspins

        IF (globenv%print%density_matrix.AND.globenv%print%each_scf_step) THEN
            CALL write_sparse_matrix(p_rmpv(ispin)%matrix,4,6,qs_env,globenv)
        END IF

        IF (globenv%print%kohn_sham_matrix.AND.globenv%print%each_scf_step) THEN
            CALL write_sparse_matrix(ks_rmpv(ispin)%matrix,4,6,qs_env,globenv)
        END IF

      ENDDO

      rebuild = (MODULO(iscf,scf_control%nrebuild) == 0)

                  type_switch=1 ! default with diagonalisation
      IF (pao)    type_switch=0
      IF (scf_control%use_ot) type_switch=2

      ispin=1
      ks => ks_rmpv(ispin)
      p => p_rmpv(ispin)
      mo_set => c(ispin)%mo_set

      SELECT CASE (type_switch)
      CASE (0) ! pao

        IF (dft_control%nspins.ne.1) CALL stop_program(routine,"NYI")

        CALL copy_real_matrix_to_full_matrix(ks%matrix,scf_work1(1)%matrix)
        IF (iscf > 1) THEN
          CALL scf_diis(c,scf_work1,scf_work2,delta,diis_error,diis_step,&
                        qs_env,globenv)
        END IF
         IF (cp_debug .AND. debug_this_module) THEN
            logger => cp_error_get_logger(error)
            CALL cp_sm_output(logger=logger, &
                 outputName="scf_ks", &
                 fromWhere=routine, iter=iscf, matrix=ks%matrix,&
                 para_env=globenv%para_env,error=error)
            CALL cp_sm_output(logger=logger,&
                 outputName="scf_p", &
                 fromWhere=routine, iter=iscf, matrix=p%matrix,&
                 para_env=globenv%para_env,error=error)
         END IF
         CALL pao_qs_env_did_change(qs_env,global_env=globenv,h_changed=.TRUE.)
         DO ! optimize angles
            IF (.NOT.cp_opt_gerd_f_next(pao_optimizer,f=pao_f,&
                 last_f=pao_last_f,&
                 projected_gradient=pao_actual_grad, error=error)) THEN
               CALL cp_assert(.FALSE.,cp_warning_level,&
                    error_nr=cp_internal_error,fromWhere=routine,&
                    message="pao optimizer reached convergence",&
                    error=error)
               EXIT
            END IF
            IF (pao_actual_grad < pao_env%wanted_subspace_pos_precision) THEN
               EXIT
            END IF
         END DO
         CALL pao_env_get(pao_env, min_density_m=min_density_m,error=error)
         DO ispin=1,dft_control%nspins
            CPPrecondition(ispin==1,cp_failure_level,routine,error,failure)
            ! pao_env%angles_att should be spin dependent
            CALL sm_matrix_transf_to_f(source_m=min_density_m(ispin)%matrix%sm,&
                 target_m=p_new(ispin)%matrix, glob_angles=pao_env%angles_att,&
                 qs_env=qs_env, global_env=globenv,error=error)
            CALL pao_guarantee_owned_NUi(pao_env%angles_att, qs_env=qs_env,&
                 global_env=globenv,error=error)
            CALL get_mo_set(pao_obj_f%mo_set(ispin)%mo_set,eigenvectors=min_mo)
            CALL get_mo_set(mo_set,eigenvectors=full_mo)
            CALL cp_sp_multiplies_blacs(pao_env%angles_att%NUi_injection,&
                 min_mo,full_mo)
         END DO
         IF (cp_debug .AND. debug_this_module) THEN
            logger => cp_error_get_logger(error)
            CALL cp_1d_r_output(logger=logger,&
                 outputName="scf_x",fromWhere=routine,&
                 iter=iscf,array=pao_env%angles_att%angles,&
                 error=error)
            WRITE (unit=output_unit, fmt='(a)') "min_mo"
            CALL w_blacs_matrix(min_mo,output_unit,globenv,formatted=.TRUE.)
            WRITE (unit=output_unit, fmt='(a)') "full_mo"
            CALL w_blacs_matrix(full_mo,output_unit,globenv,formatted=.TRUE.)
            CALL cp_sm_output(logger=logger, &
                 outputName="scf_new_p", &
                 fromWhere=routine, iter=iscf, matrix=p_new(1)%matrix,&
                 para_env=globenv%para_env,error=error)
         END IF

         t2 = m_cputime()

         IF (diis_step) THEN
            delta=0.0_wp
            DO ispin=1,dft_control%nspins
               CALL density_mixing(p_new(ispin),p_rmpv(ispin),1.0_wp,tmp,globenv)
               delta=MAX(delta,tmp)
            END DO
            IF (ionode.AND.globenv%print%scf) THEN
               WRITE (UNIT=output_unit,&
                    FMT="(T2,I5,2X,A,T15,E10.2,T30,F8.2,T40,2F20.10)")&
                    iscf,"DIIS",diis_error,t2 - t1,delta,energy%total
            END IF
         ELSE
           delta=0.0_wp
           DO ispin=1,dft_control%nspins
               CALL density_mixing(p_new(ispin),p_rmpv(ispin),scf_control%p_mix,tmp,globenv)
               delta=MAX(delta,tmp)
           END DO
           IF (ionode.AND.globenv%print%scf) THEN
             WRITE (UNIT=output_unit,&
                    FMT="(T2,I5,2X,A,T15,F6.2,T30,F8.2,T40,2F20.10)")&
               iscf,"Mixing",scf_control%p_mix,t2 - t1,delta,energy%total
             WRITE (UNIT=output_unit, FMT="(A,E10.2)") "diis error=",&
                  diis_error 
           END IF
         END IF

      CASE(1) ! diagonalisation (default)

        DO ispin=1,dft_control%nspins 
           CALL copy_real_matrix_to_full_matrix(ks_rmpv(ispin)%matrix,scf_work1(ispin)%matrix)
        ENDDO

        IF (iscf > 1) THEN
            CALL scf_diis(c,scf_work1,scf_work2,delta,diis_error,diis_step,&
                          qs_env,globenv)
        END IF

        do_level_shift = ((scf_control%level_shift /= 0.0_wp).AND.&
                         ((scf_control%density_guess == "CORE").OR.(iscf > 1)))

        IF (diis_step) THEN
           pmix = 1.0_wp
           method_par = diis_error
           IF (use_jacobi) THEN
             method = "DIIS/Jacobi"
           ELSE
             method = "DIIS/SYEVX"
           END IF
        ELSE
           pmix = scf_control%p_mix
           method_par = pmix
           IF (use_jacobi) THEN
             method = "Mixing/Jacobi"
           ELSE
             method = "Mixing/SYEVX"
           END IF
        END IF

        IF ((iscf > 1).AND.(delta < scf_control%eps_jacobi)) THEN
          use_jacobi = .TRUE.
        ELSE
          use_cholesky = scf_control%use_cholesky
          use_jacobi = .FALSE.
        END IF

        delta = 0.0_wp

        DO ispin=1,dft_control%nspins

           CALL eigensolver(scf_work1(ispin)%matrix,c(ispin)%mo_set,ortho,scf_work2,&
                            do_level_shift,scf_control%level_shift,globenv,&
                            use_arpack=scf_control%use_arpack,&
                            use_cholesky=use_cholesky,&
                            work_syevx=scf_control%work_syevx,&
                            use_jacobi=use_jacobi,&
                            eps_jacobi=scf_control%eps_jacobi,&
                            jacobi_threshold=scf_control%jacobi_threshold)

           IF (scf_control%smear /= 0.0_wp) THEN
               CALL set_mo_occupation(c(ispin)%mo_set,scf_control%smear)
           END IF

           IF (globenv%print%each_scf_step) THEN
              CALL write_mo_set(c(ispin)%mo_set,atomic_kind_set,particle_set,4,6,globenv)
           END IF

           CALL calculate_density_matrix(c(ispin)%mo_set,p_new(ispin)%matrix,globenv)

           CALL density_mixing(p_new(ispin),p_rmpv(ispin),pmix,tmp,globenv)
           
           delta=MAX(delta,tmp)

        ENDDO

        t2 = m_cputime()

        IF (ionode.AND.globenv%print%scf) THEN
           WRITE (UNIT=output_unit,&
                  FMT="(T2,I5,2X,A,T22,E10.2,T32,F8.2,T40,2F20.10)")&
                  iscf,TRIM(method),method_par,t2 - t1,delta,energy%total
        END IF

      CASE(2) ! orbital transforms

        IF (dft_control%nspins.ne.1) CALL stop_program(routine,"NYI")

        CALL get_mo_set(mo_set=mo_set,eigenvectors=mo_eigenvectors)

        ! might need the KS matrix to init properly, that's why we're here
        IF (.NOT. ot_initialized) THEN
           ot_initialized=.TRUE.

           IF (qs_ot_env%settings%preconditioner_type .NE. "NONE") &
           CALL make_preconditioner(ot_preconditioner, &
                                    qs_ot_env%settings%preconditioner_type, &
                                    ks%matrix,s(1)%matrix,kinetic(1)%matrix, &
                                    mo_eigenvectors,0.2_wp) 

           CALL ot_scf_init(mo_eigenvectors,s(1)%matrix,qs_ot_env,globenv)

           IF (qs_ot_env%settings%preconditioner_type .NE. "NONE") &
           CALL qs_ot_new_preconditioner(qs_ot_env,ot_preconditioner) 
        ENDIF

        CALL ot_scf_mini(mo_eigenvectors,ks%matrix,s(1)%matrix,energy%total, &
                                          energy_only,delta,qs_ot_env,globenv)

        CALL calculate_density_matrix(mo_set,p%matrix,globenv)
        t2 = m_cputime()
        IF (ionode.AND.globenv%print%scf) THEN
              WRITE (UNIT=output_unit,&
                 FMT="(T2,I5,2X,A,T20,F6.2,T30,F8.2,T40,2F20.10)")&
                         iscf,qs_ot_env%OT_METHOD_FULL,qs_ot_env%ds_min, &
                                                t2 - t1,delta,energy%total
        END IF


      END SELECT

      CALL qs_ks_did_change(ks_env,qs_env,p_changed=.TRUE.,error=error)

      CALL timestop(0.0_wp,handle2)

      IF (delta < scf_control%eps_scf) THEN
        IF (ionode.AND.globenv%print%scf) THEN
          WRITE(UNIT=output_unit,FMT="(/,T3,A,/)")&
            "*** SCF run converged ***"
        END IF
        EXIT scf_loop
      ELSE IF (iscf == scf_control%max_scf) THEN
        IF (ionode.AND.globenv%print%scf) THEN
          WRITE(UNIT=output_unit,FMT="(/,T3,A,/)")&
            "*** SCF run NOT converged ***"
        END IF
        EXIT scf_loop
      END IF

    END DO scf_loop

    IF (scf_control%use_ot) CALL ot_scf_destroy(qs_ot_env,globenv)

    IF (ionode.AND.globenv%print%scf) THEN
      WRITE (UNIT=output_unit,FMT="(/,(T3,A,T40,2F20.10))")&
        "Total electronic density (r-space): ",&
        SUM(rho%tot_rho_r),&
        SUM(rho%tot_rho_r)+ REAL(nelectron,wp),&
        "Total core charge density (r-space):",&
        qs_charges%total_rho_core_rspace,&
        qs_charges%total_rho_core_rspace - REAL(nelectron+dft_control%charge,wp)
      WRITE (UNIT=output_unit,FMT="(T3,A,T60,F20.10)")&
        "Total charge density (r-space):     ",&
        SUM(rho%tot_rho_r)+&
        qs_charges%total_rho_core_rspace,&
        "Total charge density (g-space):     ",qs_charges%total_rho_gspace
      WRITE (UNIT=output_unit,FMT="(/,(T3,A,T60,F20.10))")&
        "Overlap energy of the core charge distribution:",energy%core_overlap,&
        "Self energy of the core charge distribution:   ",energy%core_self,&
        "Core Hamiltonian energy:                       ",energy%core,&
        "Hartree energy:                                ",energy%hartree,&
        "Exchange-correlation energy:                   ",energy%exc,&
        "Total energy:                                  ",energy%total
    END IF

!   *** Undo mixing of the density matrix (restore original density) ***
 
    DO ispin=1,dft_control%nspins 
     SELECT CASE (type_switch)
     CASE (0) ! pao
       IF (.NOT.diis_step) CALL transfer_matrix(p_new(ispin)%matrix,p_rmpv(ispin)%matrix)
     CASE (1) ! diagonalisation
       IF (.NOT.diis_step) CALL transfer_matrix(p_new(ispin)%matrix,p_rmpv(ispin)%matrix)
     CASE (2)
     END SELECT
    ENDDO

!   *** Save last density matrix in redundant form for the force calculation ***
    CALL calculate_density_matrix(mo_set,rho%rho_ao(1)%matrix,globenv)
!debCALL qs_rho_did_change(rho,rho_ao_valid=.TRUE.,error=error)
    rho%rho_ao_valid=.TRUE.
!deb hack for SGI

!   *** Write last wavefunction to screen ??????? ***

    do ispin=1,dft_control%nspins
       CALL write_mo_set(c(ispin)%mo_set,atomic_kind_set,particle_set,4,6,globenv)
    ENDDO

!   *** write additional information about the converged system
    IF (globenv%print%e_density_cube) THEN
      ! write the electron density
      IF (globenv%ionode) OPEN(UNIT=21,FILE="E_DENSITY.cube")
      CALL rs_pw_to_cube(rho%rho_r(1)%pw,21,globenv%ionode,"ELECTRON DENSITY")
      IF (globenv%ionode) CLOSE(21)
    ENDIF

    IF (globenv%print%v_hartree_cube) THEN
      ! write the hartree potential
      IF (globenv%ionode) OPEN(UNIT=21,FILE="V_HARTREE.cube")
      CALL rs_pw_to_cube(ks_env%v_hartree_rspace%pw,21,globenv%ionode,"HARTREE POTENTIAL")
      IF (globenv%ionode) CLOSE(21)
    ENDIF

    IF (globenv%print%homos .OR. globenv%print%lumos) THEN
       ! gets the tmp grids
       CALL get_qs_env(qs_env=qs_env,pw_env=pw_env)
       CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool,&
            pw_pools=pw_pools)
       CALL pw_pool_init_coeff(auxbas_pw_pool,wf_r,&
            use_data = REALDATA3D,&
            in_space = REALSPACE, error=error)
       CALL pw_pool_init_coeff(auxbas_pw_pool,wf_g,&
            use_data = COMPLEXDATA1D,&
            in_space = RECIPROCALSPACE, error=error)
       ALLOCATE(mgrid_rspace(SIZE(pw_pools)), mgrid_gspace(SIZE(pw_pools)),&
            stat=istat)
       CPPostcondition(istat==0,cp_failure_level,routine,error,failure)
       CALL pw_pools_init_coeffs(pw_pools,mgrid_gspace,&
            use_data = COMPLEXDATA1D,&
            in_space = RECIPROCALSPACE, error=error)
       CALL pw_pools_init_coeffs(pw_pools,mgrid_rspace,&
            use_data = REALDATA3D,&
            in_space = REALSPACE, error=error)
    END IF

    IF (globenv%print%homos) THEN
      DO ispin=1,dft_control%nspins
         ! write eigenvalues and a set of cubes
         CALL get_mo_set(mo_set=c(ispin)%mo_set,eigenvectors=mo_eigenvectors,homo=homo)
         IF (globenv%ionode) WRITE(globenv%scr,*) " "
         IF (globenv%ionode) WRITE(globenv%scr,*) " Eigenvalues of the occupied subspace spin ",ispin
         IF (globenv%ionode) WRITE(globenv%scr,*) " ---------------------------------------------"
         CALL calculate_subspace_eigenvalues(mo_eigenvectors,ks_rmpv(ispin)%matrix, &
                 para_env=globenv%para_env,ctxt=globenv%blacs_env,scr=globenv%scr, &
                 ionode=globenv%ionode)

         DO ivector=MAX(1,homo-globenv%print%nhomos+1),homo
            CALL calculate_wavefunction(mo_eigenvectors,ivector,wf_r, &
                                        wf_g,mgrid_rspace,&
                                        mgrid_gspace,&
                                        pw_env%gridlevel_info,pw_env%cube_info, &
                                        pw_env%l_info,&
                                        qs_env,globenv)

            WRITE(filename,'(a13,I5.5,a1,I1.1,a5)') "WAVEFUNCTION_",ivector,"_",ispin,".cube"
            IF (globenv%ionode) OPEN(UNIT=21,FILE=filename)
            WRITE(filename,*) "WAVEFUNCTION ",ivector," spin ",ispin," i.e. HOMO - ",ivector-homo
            CALL rs_pw_to_cube(wf_r%pw,21,globenv%ionode,filename)
            IF (globenv%ionode) CLOSE(21)
         ENDDO
      ENDDO
    ENDIF

    IF (globenv%print%lumos) THEN
      DO ispin=1,dft_control%nspins
         ! write eigenvalues and a set of cubes
         IF (globenv%ionode) WRITE(globenv%scr,*) " "
         IF (globenv%ionode) WRITE(globenv%scr,*) " Lowest Eigenvalues of the unoccupied subspace spin ",ispin
         IF (globenv%ionode) WRITE(globenv%scr,*) " -----------------------------------------------------"
         CALL get_mo_set(mo_set=c(ispin)%mo_set,eigenvectors=mo_eigenvectors,homo=homo,nao=nao)
         CALL cp_fm_get_info(mo_eigenvectors, nrow_global=n, ncol_global=k)
         nlumos=MAX(1,MIN(globenv%print%nlumos,nao-homo))
         CALL cp_fm_create2(new_matrix=unoccupied_orbs, &
              nrow_block=scf_control%nrow_block, &
              ncol_block=scf_control%ncol_block, &
              nrow_global=n, &
              ncol_global=nlumos, &
              name="unoccupied_orbs", &
              globenv=globenv)
         CALL cp_fm_init_random(unoccupied_orbs,nlumos)

         CALL ot_eigensolver(matrix_h=ks_rmpv(ispin)%matrix,matrix_s=s(1)%matrix, &
                             matrix_c=unoccupied_orbs, &
                             matrix_orthogonal_space=mo_eigenvectors, &
                             eps_gradient=scf_control%eps_lumos, &
                             preconditioner=ot_preconditioner, &
                             iter_max=scf_control%max_iter_lumos,globenv=globenv)
  
         CALL calculate_subspace_eigenvalues(unoccupied_orbs,ks_rmpv(ispin)%matrix,&
                 para_env=globenv%para_env,ctxt=globenv%blacs_env,scr=globenv%scr, &
                 ionode=globenv%ionode)

         DO ivector=1,nlumos
            CALL calculate_wavefunction(mo_vectors=unoccupied_orbs,&
                 ivector=ivector,rho=wf_r, &
                 rho_gspace=wf_g,&
                 mgrid_rspace=mgrid_rspace,&
                 mgrid_gspace=mgrid_gspace,&
                 gridlevel_info=pw_env%gridlevel_info,&
                 cube_info=pw_env%cube_info, &
                 l_info=pw_env%l_info,&
                 qs_env=qs_env,globenv=globenv)

            WRITE(filename,'(a13,I5.5,a1,I1.1,a5)') "WAVEFUNCTION_",homo+ivector,"_",ispin,".cube"
            IF (globenv%ionode) OPEN(UNIT=21,FILE=filename)
            WRITE(filename,*) "WAVEFUNCTION ",homo+ivector," spin ",ispin," i.e. LUMO + ",ivector-1
            CALL rs_pw_to_cube(wf_r%pw,21,globenv%ionode,filename)
            IF (globenv%ionode) CLOSE(21)

         ENDDO
      ENDDO
   ENDIF

   IF (globenv%print%lumos .OR. globenv%print%homos) THEN
      ! put back the grids
      CALL pw_pool_give_back_coeff(auxbas_pw_pool,wf_r,&
           error=error)
      CALL pw_pool_give_back_coeff(auxbas_pw_pool,wf_g,&
           error=error)
      CALL pw_pools_give_back_coeffs(pw_pools,mgrid_gspace,&
           error=error)
      CALL pw_pools_give_back_coeffs(pw_pools,mgrid_rspace,&
           error=error)
      DEALLOCATE(mgrid_rspace, mgrid_gspace, stat=istat)
      CPPostconditionNoFail(istat==0,cp_warning_level,routine,error)
      
   END IF
    IF (scf_control%use_ot) THEN
       CALL destroy_preconditioner(ot_preconditioner)
    ENDIF

    IF (globenv%print%density_matrix) THEN
      CALL write_sparse_matrix(p%matrix,4,6,qs_env,globenv)
    END IF

    IF (Globenv%print%density_matrix_magnitude) THEN
        IF (Globenv%num_pe > 1) THEN
           WRITE (Unit=Extension,Fmt="(I6)") Globenv%mepos
           File_name = "P_magnitude"//"."//ADJUSTL(Extension)
        ELSE
           File_name = "P_magnitude"
        END IF
        CALL Open_file(File_name=File_name,&
                       File_action="Write",&
                       File_form="Formatted",&
                       File_status="Replace",&
                       Unit_number=File_unit)
      
       NULLIFY(Pfull) 
       CALL Replicate_matrix(S(1)%matrix,Pfull,"Full Density Matrix",.FALSE.)
       CALL Add_all_real_matrix_blocks(Pfull,Globenv)
       CALL Calculate_density_matrix(Mo_set,Pfull,Globenv)
       CALL Write_blocks_maxabsval(Pfull,File_unit)
       CALL Deallocate_matrix(Pfull)
       CALL Close_file(Unit_number=File_unit)
    ENDIF

    IF (Globenv%print%kohn_sham_matrix) THEN
      CALL Write_sparse_matrix(Ks%matrix,4,6,Qs_env,Globenv)
    END IF

!   *** Release SCF work storage ***

    IF (ASSOCIATED(scf_work1)) THEN
       DO ispin=1,SIZE(scf_work1)
         CALL cp_fm_release(scf_work1(ispin)%matrix)
       ENDDO
       DEALLOCATE(scf_work1)
    ENDIF
    IF (ASSOCIATED(scf_work2)) CALL cp_fm_release(scf_work2)
    IF (ASSOCIATED(ortho)) CALL cp_fm_release(ortho)

    DO ispin=1,dft_control%nspins
      CALL deallocate_matrix(p_new(ispin)%matrix)
    ENDDO
    DEALLOCATE(p_new)

!   *** Write restart file ***

    CALL write_mo_set(c,atomic_kind_set,particle_set,globenv,&
         id_nr=qs_env%id_nr)

!    call pao_opt_dealloc(pao_optimizer,pao_obj_f,qs_env,globenv,error)

    CALL timestop(0.0_wp,handle)

    CALL write_checkpoint_information("leaving "//routine_name,globenv)

  END SUBROUTINE scf

! *****************************************************************************

  SUBROUTINE init_scf_run(qs_env,globenv,error)

!   Purpose: Initialise a SCF run.

!   History: - Creation (30.04.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(qs_environment_type), INTENT(INOUT)  :: qs_env
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error
 
!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "init_scf_run"
    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE "//routine_name//" (MODULE "//module_name//")"

!   *** Local variables ***

    TYPE(cp_full_matrix_type), POINTER :: mo_eigenvectors
    TYPE(mo_set_type), POINTER         :: mo_set

    INTEGER :: handle,homo,ikind,istat,nao,ndep,nmo,npcol,nprow,&
               output_unit, nrow_block, ncol_block,i
    LOGICAL :: ionode, pao, gth_potential_present, failure
    REAL(kind=wp):: ecore_self, ecore_overlap

    INTEGER, DIMENSION(0:globenv%num_pe-1) :: nblock_pe,nelement_pe
    INTEGER :: ispin,nelectron_spin
    TYPE(real_matrix_p_type), DIMENSION(:), POINTER :: k_rmpv, p_rmpv
    TYPE(scf_control_type), POINTER :: scf_control
    TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set
    TYPE(cp_fm_pool_p_type),DIMENSION(:), POINTER :: ao_mo_fm_pools
    TYPE(real_matrix_p_type), DIMENSION(:), POINTER :: h,s
    TYPE(cp_fm_struct_type), POINTER :: ao_mo_mstruct
    TYPE(cp_full_matrix_type), POINTER :: eigenvectors
    TYPE(dft_control_type), POINTER   :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), POINTER    :: c

!   ---------------------------------------------------------------------------

    CALL timeset("init_scf_run","I","",handle)

    NULLIFY(k_rmpv, p_rmpv, scf_control, atomic_kind_set,&
         ao_mo_fm_pools,h,s)
    failure=.FALSE.
    CALL get_qs_env(qs_env=qs_env,&
                    scf_control=scf_control,&
                    dft_control=dft_control,&
                    atomic_kind_set=atomic_kind_set,&
                    c=c,&
                    k=k_rmpv,&
                    p=p_rmpv,&
                    ao_mo_fm_pools=ao_mo_fm_pools,&
                    h=h,s=s)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             gth_potential_present=gth_potential_present)

    nrow_block = scf_control%nrow_block
    ncol_block = scf_control%ncol_block

    ionode = globenv%ionode
    output_unit = globenv%scr

!   *** Get the dimension of the full SCF matrices, ***
!   *** i.e. the total number of atomic orbitals    ***

!   *** Initializes the MOs ***
    CPPrecondition(ASSOCIATED(c),cp_failure_level,routine,error,failure)
    IF (.NOT.failure) THEN
       DO i=1,SIZE(c)
          CALL get_mo_set(c(i)%mo_set,eigenvectors=eigenvectors)
          IF (.NOT.ASSOCIATED(eigenvectors)) THEN
             CALL init_mo_set(c(i)%mo_set, ao_mo_fm_pools(i)%pool,&
                  name="qs_env"//TRIM(ADJUSTL(cp_to_string(qs_env%id_nr)))//&
                  "%mo"//TRIM(ADJUSTL(cp_to_string(i))),&
                  error=error)
          END IF
       END DO
    END IF

!   *** Get BLACS block size of the MO eigenvector matrix      ***
!   *** which has to fit to the other distributed SCF matrices ***

    mo_set => c(1)%mo_set
    CALL get_mo_set(mo_set=mo_set,&
                    eigenvectors=mo_eigenvectors,&
                    nao=nao)

    CALL cp_fm_get_info(matrix=mo_eigenvectors,&
                               nrow_block=nrow_block,&
                               ncol_block=ncol_block)

!   *** Allocate the distributed SCF matrices ***
  
    NULLIFY(scf_work1)
    ALLOCATE(scf_work1(dft_control%nspins))

    CALL cp_fm_create2(new_matrix=scf_work1(1)%matrix,&
                               nrow_global=nao,&
                               ncol_global=nao,&
                               nrow_block=nrow_block,&
                               ncol_block=ncol_block,&
                               name="SCF WORK MATRIX 1",&
                               globenv=globenv)

    IF (.NOT.scf_control%use_ot) THEN
      DO ispin=2,dft_control%nspins
       CALL cp_fm_create2(new_matrix=scf_work1(ispin)%matrix,&
                               nrow_global=nao,&
                               ncol_global=nao,&
                               nrow_block=nrow_block,&
                               ncol_block=ncol_block,&
                               name="SCF WORK MATRIX 1",&
                               globenv=globenv)
      ENDDO
      CALL cp_fm_create2(new_matrix=ortho,&
                                 nrow_global=nao,&
                                 ncol_global=nao,&
                                 nrow_block=nrow_block,&
                                 ncol_block=ncol_block,&
                                 name="ORTHOGONALIZATION MATRIX",&
                                 globenv=globenv)
      CALL cp_fm_create2(new_matrix=scf_work2,&
                                 nrow_global=nao,&
                                 ncol_global=nao,&
                                 nrow_block=nrow_block,&
                                 ncol_block=ncol_block,&
                                 name="SCF WORK MATRIX 2",&
                                 globenv=globenv)
    END IF

!   *** Allocate k and p and put them in the QS environment ***

    CALL allocate_matrix_set(k_rmpv,dft_control%nspins)
    CALL allocate_matrix_set(p_rmpv,dft_control%nspins)
    CALL set_qs_env(qs_env=qs_env,&
                    k=k_rmpv,&
                    p=p_rmpv)

    DO ispin=1,dft_control%nspins
       NULLIFY (k_rmpv(ispin)%matrix)
       CALL replicate_matrix_structure(s(1)%matrix,k_rmpv(ispin)%matrix,&
                                       "KOHN-SHAM MATRIX")
       NULLIFY (p_rmpv(ispin)%matrix)
       CALL replicate_matrix_structure(s(1)%matrix,p_rmpv(ispin)%matrix,&
                                       "DENSITY MATRIX")
    ENDDO


!   *** Redistribute the core Hamiltonian matrix in ***
!   *** order to eliminate redundant atomic blocks  ***

    CALL copy_real_matrix_to_full_matrix(h(1)%matrix,scf_work1(1)%matrix)

    IF ((globenv%num_pe > 1).AND.gth_potential_present) THEN
      CALL deallocate_matrix(h(1)%matrix)
      CALL replicate_matrix_structure(s(1)%matrix,h(1)%matrix,&
                                      "CORE HAMILTONIAN MATRIX")
      CALL copy_full_matrix_to_real_matrix(scf_work1(1)%matrix,h(1)%matrix)
      CALL set_qs_env(qs_env,h=h)
    END IF

    IF (scf_control%use_ot) THEN 
        CALL cp_fm_release(scf_work1(1)%matrix)
        DEALLOCATE(scf_work1)
        NULLIFY(scf_work1)
    ENDIF

!   *** done redistribute ***

    CALL calculate_ecore_self(qs_env)
    CALL calculate_ecore_overlap(qs_env,globenv,.FALSE.)

    ecore_overlap = qs_env%energy%core_overlap
    ecore_self = qs_env%energy%core_self

    IF (ionode.AND.globenv%print%scf_energies) THEN
      WRITE (UNIT=output_unit,FMT="(/,(T3,A,T60,F20.10))")&
        "Overlap energy of the core charge distribution:",ecore_overlap,&
        "Self energy of the core charge distribution:   ",ecore_self
    END IF

    IF (.NOT. scf_control%use_ot) THEN
     IF (scf_control%use_cholesky .OR. scf_control%use_arpack) THEN
       ndep=0
       CALL copy_real_matrix_to_full_matrix(s(1)%matrix,ortho)
       CALL cp_fm_cholesky_decompose(ortho)
     ELSE
       CALL calculate_ortho_matrix(ortho,s(1)%matrix,scf_work2,ndep,globenv,&
            scf_control%eps_eigval,scf_control%work_syevx)
     ENDIF
    ENDIF

    IF (ionode.AND.globenv%print%scf) THEN
      DO ispin=1,dft_control%nspins
        mo_set => c(ispin)%mo_set
        CALL get_mo_set(mo_set=mo_set,&
                        homo=homo, &
                        nelectron=nelectron_spin, &
                        nao=nao)
        WRITE (UNIT=output_unit,fmt=*) " For spin ",ispin
        WRITE (UNIT=output_unit,FMT="(/,(T3,A,I10))")&
           "Number of electrons:                    ",nelectron_spin,&
           "Number of occupied orbitals:            ",homo,&
           "Number of orbital functions:            ",nao,&
           "Number of independent orbital functions:",nao - ndep
      ENDDO
    END IF

    IF (globenv%print%ortho_matrix) THEN
      CALL write_blacs_matrix(ortho,4,6,qs_env,globenv)
    END IF

    ! we need a copy the core matrix for every spin..
    IF (dft_control%nspins.eq.2) THEN
       CALL copy_blacs_to_blacs_matrix(scf_work1(1)%matrix,scf_work1(2)%matrix)
    ENDIF
    CALL calculate_first_density_matrix(ortho,c,p_rmpv,s(1), &
             scf_work1,scf_work2,globenv,&
             scf_control=scf_control, id_nr=qs_env%id_nr)

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE init_scf_run

! *****************************************************************************

  SUBROUTINE calculate_ortho_matrix(ortho,s,work,ndep,globenv,eps_eigval,&
       work_syevx)

!   Purpose: Calculate the orthogonalization matrix (S**(-1/2))

!   History: - Creation (01.05.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(cp_full_matrix_type), POINTER        :: ortho,work
    TYPE(real_matrix_type), POINTER           :: s
    INTEGER, INTENT(OUT)                      :: ndep
    REAL(kind=wp), INTENT(in)                 :: eps_eigval,work_syevx

!   *** Local variables ***

    INTEGER :: handle,nao

!   ---------------------------------------------------------------------------

    CALL timeset("calculate_ortho_matrix","I","",handle)

    CALL copy_real_matrix_to_full_matrix(s,ortho)
    CALL power_blacs_matrix(ortho,work,-0.5_wp,eps_eigval,ndep,work_syevx)
    CALL symmetrise_blacs_matrix(ortho,work)

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE calculate_ortho_matrix

! *****************************************************************************

  SUBROUTINE eigensolver(ks,mo_set,ortho,work,do_level_shift,level_shift,&
                         globenv,use_cholesky, use_arpack,work_syevx,&
                         use_jacobi,eps_jacobi,jacobi_threshold)

!   Purpose: Diagonalise the Kohn-Sham matrix to get a new set of MO eigen-
!            vectors and MO eigenvalues.
!   ks will be modified

!   History: - Creation (01.05.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(cp_full_matrix_type), POINTER          :: ks,ortho,work
    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(mo_set_type), POINTER                :: mo_set
    LOGICAL, INTENT(INOUT)                    :: use_cholesky
    LOGICAL, INTENT(IN)                       :: do_level_shift,use_arpack,&
                                                 use_jacobi
    REAL(wp), INTENT(IN)                      :: level_shift,work_syevx,&
                                                 eps_jacobi,jacobi_threshold

!   *** Local variables ***
    INTEGER :: handle,imo,nao,nmo,homo,ncv,handle2
    TYPE(cp_full_matrix_type), POINTER :: mo_eigenvectors
    REAL(wp), DIMENSION(:), POINTER    :: mo_eigenvalues
    TYPE(op_symv_type)                 :: arpack_op

!   ---------------------------------------------------------------------------

    CALL timeset("eigensolver","I","",handle)

!   *** Diagonalise the Kohn-Sham matrix ***

    CALL get_mo_set(mo_set=mo_set,&
                    nao=nao,&
                    nmo=nmo,&
                    homo=homo, &
                    eigenvalues=mo_eigenvalues,&
                    eigenvectors=mo_eigenvectors)

    IF (do_level_shift) THEN
      IF (use_cholesky .OR. use_arpack) THEN
        CALL stop_program("eigensolver","level shift not implemented")
      ENDIF
    ENDIF

    IF (use_cholesky) THEN

       CALL cp_fm_cholesky_reduce(ks,ortho)
       IF (use_jacobi) THEN
         CALL blacs_syevx(ks,work,mo_eigenvalues,nao,work_syevx)
         use_cholesky = .FALSE.
         CALL cp_fm_cholesky_restore(work,nao,ortho,mo_eigenvectors,"SOLVE")
         CALL copy_blacs_to_blacs_matrix(mo_eigenvectors,ortho)
      ELSE
         CALL blacs_syevx(ks,work,mo_eigenvalues,homo,work_syevx)
         CALL cp_fm_cholesky_restore(work,homo,ortho,mo_eigenvectors,"SOLVE")
       END IF

    ELSEIF (use_arpack) THEN ! not really recommended, since rather slow, except for huge basis sets

!TC       ncv = nmo + 64
       ncv = nmo*2
       CALL cp_fm_cholesky_reduce(ks,ortho)
       
       CALL setup_op_symv(arpack_op, ks, work, .FALSE.)
       CALL arpack_diag_sy(work,mo_eigenvalues,nmo,ncv,.FALSE.,globenv, &
                           op_symv=arpack_op)

       CALL cp_fm_cholesky_restore(work,nmo,ortho,mo_eigenvectors,"SOLVE")

    ELSE

       CALL cp_fm_get_info(matrix=ks,nrow_global=nao)

!MK    CALL cp_fm_gemm("T","N",nao,nao,nao,1.0_wp,ks,ortho,0.0_wp,work)
       CALL cp_fm_symm("L","U",nao,nao,1.0_wp,ks,ortho,0.0_wp,work)
       CALL cp_fm_gemm("T","N",nao,nao,nao,1.0_wp,ortho,work,0.0_wp,ks)
 
       IF (do_level_shift) THEN
         DO imo=homo+1,nmo
            CALL blacs_add_to_element(ks,imo,imo,level_shift)
         END DO
       ENDIF

       IF (use_jacobi) THEN

!        epsi =1e-11
!        Anz_max=4
!        N=ks%matrix_struct%nrow_global

!        *** Scale the elements with the inverse ***
!        *** difference of the old eigenvalues   ***

!        CALL blacs_block_jacobi_scaled(ks,work,mo_eigenvalues,&
!                                       jacobi_threshold,homo+1)

!        *** Klassisches Jacobi-Block-Verfrahren ***
!        *** Rotiere Elemente nur, wenn sie groesser als thresh sind ***

         CALL blacs_block_jacobi_classic(ks,work,jacobi_threshold,homo+1)

       ELSE

         CALL blacs_syevx(ks,work,mo_eigenvalues,nmo,work_syevx)

       END IF

       CALL cp_fm_gemm("N","N",nao,nmo,nao,1.0_wp,ortho,work,0.0_wp,&
                        mo_eigenvectors)

       IF (do_level_shift) CALL correct_mo_eigenvalues(mo_set,level_shift)

!      *** Use last MO set as orthogonalization matrix ***

       IF ((do_level_shift).OR.(eps_jacobi /= 0.0_wp)) THEN
         CALL copy_blacs_to_blacs_matrix(mo_eigenvectors,ortho)
       END IF

    ENDIF
    ! write(6,*) "evals syevx",mo_eigenvalues(1:nmo)

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE eigensolver

! *****************************************************************************

  SUBROUTINE density_mixing(new_density,old_density,p_mix,delta,globenv)

!   Purpose: Perform a density mixing of the old (last SCF iteration) and the
!            new (current) density matrix.

!   History: - Creation (01.05.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(real_matrix_p_type), INTENT(INOUT)   :: new_density,old_density
    REAL(wp), INTENT(IN)                      :: p_mix
    REAL(wp), INTENT(OUT)                     :: delta

!   *** Local variables ***

    TYPE(real_block_node_type), POINTER :: block_node

    REAL(wp) :: r_mix
    INTEGER  :: handle,i,iblock_col,iblock_row,j,nblock_row

    REAL(wp), DIMENSION(:,:), POINTER :: p_new_block,p_old_block

!   ---------------------------------------------------------------------------

    CALL timeset("density_mixing","I","",handle)

    delta = 0.0_wp

    r_mix = 1.0_wp - p_mix

    CALL get_matrix_info(matrix=old_density%matrix,nblock_row=nblock_row)

    DO iblock_row=1,nblock_row

      block_node => first_block_node(matrix=new_density%matrix,&
                                     block_row=iblock_row)

      DO WHILE (ASSOCIATED(block_node))

        CALL get_block_node(block_node=block_node,&
                            block_col=iblock_col,&
                            BLOCK=p_new_block)

        CALL get_block_node(matrix=old_density%matrix,&
                            block_row=iblock_row,&
                            block_col=iblock_col,&
                            BLOCK=p_old_block)

        DO j=1,SIZE(p_new_block,2)
          DO i=1,SIZE(p_new_block,1)
            p_new_block(i,j) = p_new_block(i,j) - p_old_block(i,j)
            delta = MAX(delta,ABS(p_new_block(i,j)))
            p_old_block(i,j) = p_old_block(i,j) + p_mix*p_new_block(i,j)
            p_new_block(i,j) = p_old_block(i,j) + r_mix*p_new_block(i,j)
          END DO
        END DO

        block_node => next_block_node(block_node)

      END DO

    END DO

    CALL mp_max(delta,globenv%group)

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE density_mixing

! *****************************************************************************

  SUBROUTINE calculate_first_density_matrix(ortho,mo_array,p_rmpv,s_sparse,&
       work1,work2,globenv,scf_control,id_nr)

    TYPE(cp_full_matrix_type), POINTER        :: ortho,work2
    TYPE(cp_full_matrix_p_type), DIMENSION(:), POINTER        :: work1
    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(mo_set_p_type),  DIMENSION(:), POINTER               :: mo_array
    TYPE(real_matrix_p_type), INTENT(INOUT)   :: s_sparse
    TYPE(real_matrix_p_type), DIMENSION(:), INTENT(INOUT)   :: p_rmpv
    TYPE(scf_control_type), POINTER           :: scf_control
    INTEGER, INTENT(in)                       :: id_nr!FM ugly, should be removed

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name =&
      "calculate_first_density_matrix"
    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE "//routine_name//" (MODULE "//module_name//")"

!   *** Local variables ***

    TYPE(all_potential_type), POINTER :: all_potential
    TYPE(atomic_kind_type), POINTER   :: atomic_kind
    TYPE(cp_full_matrix_type), POINTER  :: mo_eigenvectors,sv
    TYPE(gth_potential_type), POINTER :: gth_potential
    TYPE(gto_basis_set_type), POINTER :: orb_basis_set

    REAL(wp) :: paa,saa,sum
    INTEGER  :: atom_a,first_sgfa,handle,iatom,ikind,iset,isgf,isgfa,ishell,&
                istat,la,last_sgfa,maxl,natom,nelec,nset,restart_unit,nmo

    INTEGER, DIMENSION(:), ALLOCATABLE :: econf

    INTEGER, DIMENSION(:), POINTER   :: atom_list,elec_conf,nshell
    INTEGER, DIMENSION(:,:), POINTER :: first_sgf,l,last_sgf
    INTEGER :: nspin,ispin,nelectron

!   ---------------------------------------------------------------------------

    nspin=size(mo_array)

    CALL timeset("calculate_first_density_matrix","I","",handle)

    IF (scf_control%use_ot .AND. .NOT.  &
           (scf_control%density_guess == "RANDOM" .OR.  &
              (scf_control%density_guess == "RESTART" .AND. &
              scf_control%level_shift == 0.0_wp))) THEN
       CALL stop_program("calculate_first_density_matrix","NYI")
    ENDIF

    IF (scf_control%density_guess == "RESTART") THEN

      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set)
      CALL read_mo_set(mo_array,atomic_kind_set,particle_set,globenv,&
           id_nr=id_nr)
      
      DO ispin=1,nspin
         IF (scf_control%level_shift /= 0.0_wp) THEN
             CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,&
                        eigenvectors=mo_eigenvectors)
             CALL copy_blacs_to_blacs_matrix(mo_eigenvectors,ortho)
         END IF

         CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,&
                           eigenvectors=mo_eigenvectors, nmo=nmo)
         ! ortho so that one can restart for different positions (basis sets?)
         CALL replicate_blacs_matrix(mo_eigenvectors,sv,"SV")
         CALL cp_sm_fm_multiply(s_sparse%matrix,mo_eigenvectors,sv,nmo,&
                                para_env=globenv%para_env)
         CALL blacs_make_basis(mo_eigenvectors,nmo, &
                                    ortho=sv,otype="SV")
         CALL cp_fm_release(sv)
         CALL calculate_density_matrix(mo_array(ispin)%mo_set,p_rmpv(ispin)%matrix,globenv)
      ENDDO

     ELSE IF (scf_control%density_guess == "RANDOM") THEN

      DO ispin=1,nspin
          CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,&
                            eigenvectors=mo_eigenvectors, nmo=nmo)
          CALL cp_fm_init_random(mo_eigenvectors,nmo)
          ! ortho so that one can restart for different positions (basis sets?)
          CALL replicate_blacs_matrix(mo_eigenvectors,sv,"SV")
          CALL cp_sm_fm_multiply(s_sparse%matrix,mo_eigenvectors,sv,nmo,&
               para_env=globenv%para_env)
          CALL blacs_make_basis(mo_eigenvectors,nmo, &
                                 ortho=sv,otype="SV")
          CALL cp_fm_release(sv)
          CALL calculate_density_matrix(mo_array(ispin)%mo_set,p_rmpv(ispin)%matrix,globenv)
      ENDDO

     ELSE IF (scf_control%density_guess == "CORE") THEN

      DO ispin=1,nspin
         CALL eigensolver(work1(ispin)%matrix,mo_array(ispin)%mo_set,ortho,work2, &
              .FALSE.,0.0_wp,globenv,&
              use_arpack=scf_control%use_arpack, &
              use_cholesky=scf_control%use_cholesky, &
              work_syevx=scf_control%work_syevx,&
              use_jacobi=.FALSE.,&
              eps_jacobi=scf_control%eps_jacobi,&
              jacobi_threshold=scf_control%jacobi_threshold)
         IF (scf_control%level_shift /= 0.0_wp) THEN
           CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,&
                           eigenvectors=mo_eigenvectors)
           CALL copy_blacs_to_blacs_matrix(mo_eigenvectors,ortho)
         END IF
         CALL calculate_density_matrix(mo_array(ispin)%mo_set,p_rmpv(ispin)%matrix,globenv)
      ENDDO

    ELSE IF (scf_control%density_guess == "ATOMIC") THEN

      DO ispin=1,nspin
          CALL get_mo_set(mo_array(ispin)%mo_set,nelectron=nelectron)
          CALL cp_fm_set_all(work1(ispin)%matrix,0.0_wp)

          sum = 0.0_wp

          CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,maxl=maxl)

          ALLOCATE (econf(0:maxl),STAT=istat)
          IF (istat /= 0) CALL stop_memory(routine,"econf",(maxl + 1)*int_size)

          CALL copy_real_matrix_to_full_matrix(s_sparse%matrix,work2)

          DO ikind=1,SIZE(atomic_kind_set)

            atomic_kind => atomic_kind_set(ikind)

            CALL get_atomic_kind(atomic_kind=atomic_kind,&
                                 natom=natom,&
                                 atom_list=atom_list,&
                                 all_potential=all_potential,&
                                 gth_potential=gth_potential,&
                                 orb_basis_set=orb_basis_set)

            IF (ASSOCIATED(all_potential)) THEN
              CALL get_potential(potential=all_potential,elec_conf=elec_conf)
            ELSE IF (ASSOCIATED(gth_potential)) THEN
              CALL get_potential(potential=gth_potential,elec_conf=elec_conf)
            ELSE
              CYCLE
            END IF

            maxl = SIZE(elec_conf) - 1
            econf(:) = 0
            econf(0:maxl) = elec_conf(0:maxl)

            CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                                   nset=nset,&
                                   nshell=nshell,&
                                   l=l,&
                                   first_sgf=first_sgf,&
                                   last_sgf=last_sgf)

            DO iset=1,nset
              DO ishell=1,nshell(iset)
                la = l(ishell,iset)
                nelec = 4*la + 2
                IF (econf(la) > 0) THEN
                  IF (econf(la) >= nelec) THEN
                    paa = 2.0_wp
                    econf(la) = econf(la) - nelec
                  ELSE
                    paa = REAL(2*econf(la),wp)/REAL(nelec,wp)
                    econf(la) = 0
                  END IF
                  first_sgfa = first_sgf(ishell,iset)
                  last_sgfa = last_sgf(ishell,iset)
                  DO isgfa=first_sgfa,last_sgfa
                    DO iatom=1,natom
                      atom_a = atom_list(iatom)
                      isgf = particle_set(atom_a)%first_sgf + isgfa - 1
                      CALL cp_fm_set_element(work1(ispin)%matrix,isgf,isgf,paa)
                      CALL blacs_get_element(work2,isgf,isgf,saa)
                      sum = sum + paa*saa
                    END DO
                  END DO
                END IF
              END DO
            END DO

          END DO

          DEALLOCATE (econf,STAT=istat)
          IF (istat /= 0) CALL stop_memory(routine,"econf")

          CALL mp_sum(sum,globenv%group)
    
          saa = REAL(nelectron,wp)/sum

          DO ikind=1,SIZE(atomic_kind_set)

            atomic_kind => atomic_kind_set(ikind)

            CALL get_atomic_kind(atomic_kind=atomic_kind,&
                                 natom=natom,&
                                 atom_list=atom_list,&
                                 orb_basis_set=orb_basis_set)
    
            CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                                   nset=nset,&
                                   nshell=nshell,&
                                   l=l,&
                                   first_sgf=first_sgf,&
                                   last_sgf=last_sgf)

            DO iset=1,nset
              DO ishell=1,nshell(iset)
                first_sgfa = first_sgf(ishell,iset)
                last_sgfa = last_sgf(ishell,iset)
                DO isgfa=first_sgfa,last_sgfa
                  DO iatom=1,natom
                    atom_a = atom_list(iatom)
                    isgf = particle_set(atom_a)%first_sgf + isgfa - 1
                    CALL blacs_get_element(work1(ispin)%matrix,isgf,isgf,paa)
                    paa = paa*saa
                    CALL cp_fm_set_element(work1(ispin)%matrix,isgf,isgf,paa)
                  END DO
                END DO
              END DO
            END DO
          END DO

          CALL copy_full_matrix_to_real_matrix(work1(ispin)%matrix,p_rmpv(ispin)%matrix)

      ENDDO

    ELSE

      CALL stop_program(routine,&
                        "An invalid keyword for the initial density "//&
                        "guess was specified")

    END IF

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE calculate_first_density_matrix

!!****f* qs_scf/pao_opt_create [1.0] *
!!
!!   NAME
!!     pao_opt_create
!!
!!   FUNCTION
!!    allocates and initializes the optimizer and the objective function
!!    
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE pao_opt_create(pao_optimizer,pao_obj_f,qs_env,global_env,&
     error)
  TYPE(cp_lbfgs_opt_gerd_f_type), POINTER :: pao_optimizer
  TYPE(pao_gerd_functional_type), POINTER :: pao_obj_f
  TYPE(qs_environment_type), INTENT(inout), TARGET :: qs_env
  TYPE(global_environment_type), INTENT(in), TARGET :: global_env
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

  CHARACTER(len=*), PARAMETER :: routine_name = "pao_opt_create"
  CHARACTER(len=*), PARAMETER :: routine = module_name//":"//routine_name

  LOGICAL :: failure,pao
  TYPE(dft_control_type), POINTER :: dft_control
  TYPE(pao_env_type), POINTER :: pao_env
  REAL(kind=wp), DIMENSION(:), POINTER :: angles
  TYPE(real_matrix_p_type), DIMENSION(:), POINTER :: p
  INTEGER :: ispin
  failure=.FALSE.

  ispin=1
  
  NULLIFY(pao_optimizer,pao_obj_f)

  ! checks if pao is active
  CALL get_qs_env(qs_env,dft_control=dft_control, pao_env=pao_env, p=p)
  pao=dft_control%qs_control%pao
  IF (pao) THEN
     CPPrecondition(ASSOCIATED(pao_env),cp_failure_level,routine,error,failure)
     IF (.NOT.failure) THEN
        ALLOCATE(pao_optimizer, pao_obj_f)
        CALL pao_g_ang_get(pao_env%angles_att,qs_env=qs_env,&
             global_env=global_env, angles=angles,error=error)
        IF (.NOT.ASSOCIATED(angles)) THEN
           CALL pao_g_ang_set(pao_env%angles_att,qs_env=qs_env,&
                global_env=global_env,angles_val=0.0_wp)
           CALL pao_g_ang_get(pao_env%angles_att,qs_env=qs_env,&
                global_env=global_env, angles=angles,error=error)
        END IF
        CALL pao_qs_env_did_change(qs_env,s_changed=.TRUE.,h_changed=.TRUE.,&
             full_reset=.FALSE., global_env=global_env)
        CALL pao_gerd_f_init(pao_obj_f,qs_env=qs_env, global_env=global_env,&
             scf_work1=scf_work1(ispin)%matrix, scf_work2=scf_work2, scf_work3=ortho,&
             error=error)
        CALL cp_opt_gerd_f_init(pao_optimizer, global_env=global_env,&
             obj_funct=pao_obj_f, x0=angles,wanted_relative_f_delta=0.0_wp,&
             wanted_projected_gradient=0.0_wp,max_f_per_iter=5,error=error)        
        ! map first_density to min_bas (should be calculated in the min_bas)
        CALL sm_matrix_transf_p_to_m(source_m=p(1)%matrix, &
             target_m=pao_env%min_density_m(1)%matrix%sm,&
             glob_angles=pao_env%angles_att,&
             qs_env=qs_env, global_env=global_env, error=error)
     END IF     
  END IF
END SUBROUTINE pao_opt_create
!***************************************************************************

!!****f* qs_scf/pao_opt_dealloc [1.0] *
!!
!!   NAME
!!     pao_opt_dealloc
!!
!!   FUNCTION
!!     deallocates the types needed for the optimization in the pao process:
!!     
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE pao_opt_dealloc(pao_optimizer,pao_obj_f,qs_env,global_env,&
     error)
  TYPE(cp_lbfgs_opt_gerd_f_type), POINTER :: pao_optimizer
  TYPE(pao_gerd_functional_type), POINTER :: pao_obj_f
  TYPE(qs_environment_type), INTENT(in), TARGET :: qs_env
  TYPE(global_environment_type), INTENT(in), TARGET :: global_env
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  CHARACTER(len=*), PARAMETER :: routine_name = "pao_opt_dealloc"
  CHARACTER(len=*), PARAMETER :: routine = module_name//":"//routine_name

  LOGICAL :: failure
  INTEGER :: stat
  failure=.FALSE.
  
  IF (ASSOCIATED(pao_optimizer).AND.ASSOCIATED(pao_obj_f)) THEN
     CALL cp_opt_gerd_f_dealloc_ref(pao_optimizer)
     CALL pao_gerd_f_dealloc_ref(pao_obj_f)
     DEALLOCATE(pao_optimizer, pao_obj_f, stat=stat)
     CPPostconditionNoFail(stat==0,cp_warning_level,routine,error)
  END IF
END SUBROUTINE pao_opt_dealloc
!***************************************************************************

END MODULE qs_scf
