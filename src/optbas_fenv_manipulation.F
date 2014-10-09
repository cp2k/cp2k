!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!
MODULE optbas_fenv_manipulation
  USE admm_types,                      ONLY: admm_env_create,&
                                             admm_type
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE basis_set_types,                 ONLY: gto_basis_set_type,&
                                             init_orb_basis_set
  USE cp_blacs_env,                    ONLY: cp_blacs_env_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_get_info,&
                                             cp_dbcsr_p_type,&
                                             cp_dbcsr_type
  USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm
  USE cp_fm_basic_linalg,              ONLY: cp_fm_upper_to_full
  USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose,&
                                             cp_fm_cholesky_invert
  USE cp_fm_pool_types,                ONLY: cp_fm_pool_p_type
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_release,&
                                             cp_fm_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE f77_blas
  USE input_constants,                 ONLY: debug_print_level,&
                                             do_admm_basis_projection,&
                                             do_admm_purify_mo_diag,&
                                             use_aux_fit_basis_set,&
                                             use_orb_basis_set
  USE input_section_types,             ONLY: section_vals_add_values,&
                                             section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get,&
                                             section_vals_val_set
  USE kinds,                           ONLY: default_string_length
  USE optimize_basis_types,            ONLY: basis_optimization_type,&
                                             flex_basis_type
  USE particle_types,                  ONLY: particle_type
  USE qs_energy_utils,                 ONLY: qs_energies_init
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_interactions,                 ONLY: init_interaction_radii
  USE qs_kind_types,                   ONLY: qs_kind_type
  USE qs_ks_methods,                   ONLY: qs_ks_did_change,&
                                             qs_ks_update_qs_env
  USE qs_ks_types,                     ONLY: qs_ks_env_type,&
                                             set_ks_env
  USE qs_matrix_pools,                 ONLY: mpools_get
  USE qs_mo_io,                        ONLY: read_mo_set
  USE qs_mo_methods,                   ONLY: calculate_density_matrix
  USE qs_mo_types,                     ONLY: init_mo_set,&
                                             mo_set_p_type
  USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
  USE qs_overlap,                      ONLY: build_overlap_matrix
  USE qs_rho_methods,                  ONLY: qs_rho_update_rho
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE scf_control_types,               ONLY: scf_control_type
  USE string_utilities,                ONLY: uppercase
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: modify_input_settings,&  
            allocate_mo_sets,&
            calculate_overlap,&
            create_opt_admm_env,&
            update_basis_set,&
            calculate_ks_matrix,&
            calculate_overlap_inverse
 
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'optbas_fenv_manipulation'  

CONTAINS

! *****************************************************************************
!> \brief change settings in the training input files to initialize
!>        all needed structures and adjust settings to basis optimization
!> \param basis_optimization ...
!> \param bas_id ...
!> \param input_file ...
!> \param error ...
!> \author Florian Schiffmann
! *****************************************************************************
  SUBROUTINE modify_input_settings(basis_optimization,bas_id,input_file,error)
    TYPE(basis_optimization_type)            :: basis_optimization
    INTEGER                                  :: bas_id
    TYPE(section_vals_type), POINTER         :: input_file
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'modify_input_settings', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: atom, orig_basis
    INTEGER                                  :: ibasis, ikind, jkind, nbasis, &
                                                nkind
    TYPE(section_vals_type), POINTER         :: admm_section, dft_section, &
                                                feval_section, kind_section, &
                                                subsys_section

     feval_section => section_vals_get_subs_vals(input_file,"FORCE_EVAL",error=error)
     dft_section => section_vals_get_subs_vals(feval_section,"DFT",error=error)
     admm_section => section_vals_get_subs_vals(dft_section,"AUXILIARY_DENSITY_MATRIX_METHOD",error=error)
     subsys_section => section_vals_get_subs_vals(feval_section,"SUBSYS",error=error)
     kind_section => section_vals_get_subs_vals(subsys_section,"KIND",error=error)


     CALL section_vals_val_set(feval_section,"PRINT%DISTRIBUTION%_SECTION_PARAMETERS_",&
                               i_val=debug_print_level,error=error)
     CALL section_vals_val_set(dft_section,"SCF%PRINT%TOTAL_DENSITIES%_SECTION_PARAMETERS_",&
                               i_val=debug_print_level,error=error)
     CALL section_vals_val_set(dft_section,"SCF%PRINT%DETAILED_ENERGY%_SECTION_PARAMETERS_",&
                               i_val=debug_print_level,error=error)
     CALL section_vals_add_values(admm_section,error)
     CALL section_vals_val_set(admm_section,"ADMM_PURIFICATION_METHOD",&
                               i_val=do_admm_purify_mo_diag,error=error)
     CALL section_vals_val_set(admm_section,"METHOD",i_val=do_admm_basis_projection,error=error)

     ! add the new basis file containing the templates to the basis file list 
     CALL section_vals_val_get(dft_section,"BASIS_SET_FILE_NAME",n_rep_val=nbasis,error=error)
     CALL section_vals_val_set(dft_section,"BASIS_SET_FILE_NAME",i_rep_val=nbasis+1,&
                               c_val=basis_optimization%work_basis_file,error=error)

     ! Set the auxilarry basis in the kind sections
     CALL section_vals_get(kind_section,n_repetition=nkind,error=error)
     DO ikind=1,nkind
       CALL section_vals_val_get(kind_section,"_SECTION_PARAMETERS_",&
            c_val=atom,i_rep_section=ikind,error=error)
       CALL uppercase(atom)
       CALL section_vals_val_get(kind_section,"BASIS_SET",&
             c_val=orig_basis,i_rep_section=ikind,error=error)
       CALL section_vals_val_set(kind_section,"AUX_FIT_BASIS_SET",&
             c_val=orig_basis,i_rep_section=ikind,error=error)

       DO jkind=1,basis_optimization%nkind
          IF(atom==basis_optimization%kind_basis(jkind)%element)THEN
         
             ibasis=basis_optimization%combination(bas_id,jkind)
             CALL section_vals_val_set(kind_section,"AUX_FIT_BASIS_SET",&
                 c_val=TRIM(ADJUSTL(basis_optimization%kind_basis(jkind)%flex_basis(ibasis)%basis_name)),&
                 i_rep_section=ikind,error=error)
           EXIT
           END IF
       END DO
     END DO

  END SUBROUTINE modify_input_settings

! *****************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE allocate_mo_sets(qs_env,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'allocate_mo_sets', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ispin
    INTEGER, DIMENSION(2)                    :: nelectron_spin
    LOGICAL                                  :: natom_mismatch
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_fm_pool_p_type), DIMENSION(:), &
      POINTER                                :: ao_mo_fm_pools, &
                                                ao_mo_fm_pools_aux_fit
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos, mos_aux_fit
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(section_vals_type), POINTER         :: dft_section

    CALL get_qs_env(qs_env=qs_env,dft_control=dft_control,&
                    mos=mos,nelectron_spin=nelectron_spin,&
                    mos_aux_fit=mos_aux_fit,&
                    atomic_kind_set=atomic_kind_set,&
                    qs_kind_set=qs_kind_set,&
                    scf_control=scf_control,&
                   particle_set=particle_set,&
!                   matrix_s=matrix_s,&
                    admm_env=admm_env,&
                    error=error)
    dft_section =>  section_vals_get_subs_vals(qs_env%input,"DFT",error=error)

    CALL mpools_get(qs_env%mpools, ao_mo_fm_pools=ao_mo_fm_pools,&
         error=error)
    DO ispin=1,dft_control%nspins
       IF (.NOT.ASSOCIATED(mos(ispin)%mo_set%mo_coeff)) THEN
          CALL init_mo_set(mos(ispin)%mo_set,&
               ao_mo_fm_pools(ispin)%pool,&
               name="qs_env"//TRIM(ADJUSTL(cp_to_string(qs_env%id_nr)))//&
               "%mo"//TRIM(ADJUSTL(cp_to_string(ispin))),&
               error=error)
       END IF
    END DO
 
    CALL read_mo_set(mos,atomic_kind_set,qs_kind_set,particle_set,qs_env%para_env,&
      id_nr=0,multiplicity=dft_control%multiplicity,dft_section=dft_section,&
      natom_mismatch=natom_mismatch,error=error)
   
    CALL mpools_get(qs_env%mpools_aux_fit, ao_mo_fm_pools=ao_mo_fm_pools_aux_fit,&
           error=error)
    DO ispin=1,dft_control%nspins
        IF (.NOT.ASSOCIATED(mos_aux_fit(ispin)%mo_set%mo_coeff)) THEN
           CALL init_mo_set(mos_aux_fit(ispin)%mo_set,&
                ao_mo_fm_pools_aux_fit(ispin)%pool,&
                name="qs_env"//TRIM(ADJUSTL(cp_to_string(qs_env%id_nr)))//&
                "%mo_aux_fit"//TRIM(ADJUSTL(cp_to_string(ispin))),&
                error=error)
        END IF
    END DO
    
   END SUBROUTINE allocate_mo_sets

! *****************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param error ...
! *****************************************************************************
   SUBROUTINE calculate_ks_matrix(qs_env,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calculate_ks_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ispin
    TYPE(qs_rho_type), POINTER               :: rho

    NULLIFY(rho)
    CALL qs_energies_init(qs_env,.FALSE.,error)
    CALL get_qs_env(qs_env, rho=rho, error=error)
    DO ispin=1,qs_env%dft_control%nspins
       CALL calculate_density_matrix(qs_env%mos(ispin)%mo_set,rho%rho_ao(ispin)%matrix,error=error)
    END DO
    CALL qs_rho_update_rho(rho, qs_env, error=error)
    CALL qs_ks_did_change(qs_env%ks_env,rho_changed=.TRUE.,error=error)
    qs_env%requires_mo_derivs=.FALSE.
    qs_env%dft_control%do_admm=.FALSE.
    CALL qs_ks_update_qs_env(qs_env%ks_env,qs_env=qs_env,calculate_forces=.FALSE.,&
            error=error)
    qs_env%dft_control%do_admm=.TRUE.

   END SUBROUTINE calculate_ks_matrix

! *****************************************************************************
!> \brief ...
!> \param matrix_s ...
!> \param matrix_s_inv ...
!> \param para_env ...
!> \param context ...
!> \param error ...
! *****************************************************************************
   SUBROUTINE calculate_overlap_inverse(matrix_s,matrix_s_inv,para_env,context,error)
    TYPE(cp_dbcsr_type), POINTER             :: matrix_s
    TYPE(cp_fm_type), POINTER                :: matrix_s_inv
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_blacs_env_type), POINTER         :: context
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calculate_overlap_inverse', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: nao
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: work1

    CALL cp_dbcsr_get_info(matrix_s,nfullrows_total=nao)
    CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=nao, ncol_global=nao,&
                             para_env=para_env, context=context,error=error)

    CALL cp_fm_create(matrix_s_inv,matrix_struct=fm_struct_tmp,error=error)
    CALL cp_fm_create(work1,matrix_struct=fm_struct_tmp,error=error)
    CALL copy_dbcsr_to_fm(matrix_s,matrix_s_inv,error=error)
    CALL cp_fm_upper_to_full(matrix_s_inv,work1,error=error)
    CALL cp_fm_cholesky_decompose(matrix_s_inv,error=error)
    CALL cp_fm_cholesky_invert(matrix_s_inv,error=error)
    CALL cp_fm_upper_to_full(matrix_s_inv,work1,error=error)
    CALL cp_fm_struct_release(fm_struct_tmp,error)
    CALL cp_fm_release(work1,error=error)

   END SUBROUTINE calculate_overlap_inverse

! *****************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param S_type ...
!> \param error ...
! *****************************************************************************
   SUBROUTINE calculate_overlap(qs_env,S_type,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    CHARACTER(LEN=*)                         :: S_type
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calculate_overlap', &
      routineP = moduleN//':'//routineN

    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s, matrix_s_aux_fit, &
                                                matrix_s_aux_fit_vs_orb
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_aux_fit, &
                                                sab_aux_fit_vs_orb, sab_orb
    TYPE(qs_ks_env_type), POINTER            :: ks_env

    NULLIFY(ks_env)

    IF(S_type=="S_AB")THEN
       CALL get_qs_env(qs_env,&
                       ks_env=ks_env,&
                       matrix_s=matrix_s,&
                       sab_orb=sab_orb,&
                       error=error)
       CALL build_overlap_matrix(qs_env,matrix_s=matrix_s,&
                          matrix_name="OVERLAP",&
                          basis_set_id_a=use_orb_basis_set,&
                          basis_set_id_b=use_orb_basis_set, &
                          sab_nl=sab_orb,&
                          error=error)
       CALL set_ks_env(ks_env,matrix_s=matrix_s,error=error)
    ELSE IF(S_type=="S_AB_AUX")THEN
       CALL get_qs_env(qs_env,&
                       ks_env=ks_env,&
                       matrix_s_aux_fit=matrix_s_aux_fit,&
                       matrix_s_aux_fit_vs_orb=matrix_s_aux_fit_vs_orb,&
                       sab_aux_fit=sab_aux_fit,&
                       sab_aux_fit_vs_orb=sab_aux_fit_vs_orb,&
                       error=error)
       CALL build_overlap_matrix(qs_env,matrix_s=matrix_s_aux_fit,&
                          matrix_name="AUX_FIT_OVERLAP",&
                          basis_set_id_a=use_aux_fit_basis_set,&
                          basis_set_id_b=use_aux_fit_basis_set, &
                          sab_nl=sab_aux_fit,&
                          error=error)
       CALL build_overlap_matrix(qs_env,matrix_s=matrix_s_aux_fit_vs_orb,&
                          matrix_name="MIXED_OVERLAP",&
                          basis_set_id_a=use_aux_fit_basis_set,&
                          basis_set_id_b=use_orb_basis_set, &
                          sab_nl=sab_aux_fit_vs_orb,&
                          error=error)
       CALL set_ks_env(ks_env,matrix_s_aux_fit_vs_orb=matrix_s_aux_fit_vs_orb,&
                       matrix_s_aux_fit=matrix_s_aux_fit,error=error)
    END IF

  END SUBROUTINE calculate_overlap

! *****************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE create_opt_admm_env(qs_env,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_opt_admm_env', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: natoms
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos, mos_aux_fit
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

    NULLIFY(admm_env, dft_control, mos, mos_aux_fit, particle_set)
    CALL get_qs_env(qs_env=qs_env, admm_env=admm_env, dft_control=dft_control,&
                    mos=mos,mos_aux_fit=mos_aux_fit, particle_set=particle_set, error=error)
    natoms = SIZE(particle_set,1)
    CALL admm_env_create(qs_env%admm_env, dft_control%admm_control, mos, mos_aux_fit,&
                         qs_env%para_env, natoms, error)

  END SUBROUTINE create_opt_admm_env

! *****************************************************************************
!> \brief ...
!> \param opt_bas ...
!> \param set_id ...
!> \param bas_id ...
!> \param qs_env ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE update_basis_set(opt_bas,set_id,bas_id,qs_env,error)
    TYPE(basis_optimization_type)            :: opt_bas
    INTEGER                                  :: set_id, bas_id
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'update_basis_set', &
      routineP = moduleN//':'//routineN

    CHARACTER(default_string_length)         :: elem
    INTEGER                                  :: ibasis, ikind, jkind
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set

    CALL get_qs_env(qs_env=qs_env,dft_control=dft_control,&
         atomic_kind_set=atomic_kind_set, qs_kind_set=qs_kind_set,error=error)
    DO ikind=1,SIZE(qs_kind_set)
       DO jkind=1,opt_bas%nkind
           CALL get_atomic_kind(atomic_kind_set(ikind), name=elem)
           CALL uppercase(elem)
           IF(elem==opt_bas%kind_basis(jkind)%element)THEN
              ibasis=opt_bas%combination(bas_id,jkind)
              CALL transfer_data_to_gto(qs_kind_set(ikind)%aux_fit_basis_set,&
                                        opt_bas%kind_basis(jkind)%flex_basis(ibasis),error)
              CALL init_orb_basis_set(qs_kind_set(ikind)%aux_fit_basis_set,error)
           END IF
       END DO
    END DO

    CALL init_interaction_radii(dft_control%qs_control,atomic_kind_set,qs_kind_set,error)

   END SUBROUTINE update_basis_set

! *****************************************************************************
!> \brief ...
!> \param gto_basis ...
!> \param basis ...
!> \param error ...
! *****************************************************************************
   SUBROUTINE transfer_data_to_gto(gto_basis,basis,error)
    TYPE(gto_basis_set_type), POINTER        :: gto_basis
    TYPE(flex_basis_type)                    :: basis
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'transfer_data_to_gto', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ipgf, iset, ishell

    DO iset=1,basis%nsets
       DO ishell=1,basis%subset(iset)%ncon_tot
          DO ipgf=1,basis%subset(iset)%nexp
             gto_basis%gcc(ipgf,ishell,iset) = basis%subset(iset)%coeff(ipgf,ishell)
          END DO
       END DO
       DO ipgf=1,basis%subset(iset)%nexp
          gto_basis%zet(ipgf,iset) = basis%subset(iset)%exps(ipgf)
       END DO
    END DO

  END SUBROUTINE transfer_data_to_gto

END MODULE optbas_fenv_manipulation
     
     

