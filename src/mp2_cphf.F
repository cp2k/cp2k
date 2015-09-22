!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Rountines to calculate CPHF like update and solve Z-vector equation 
!>        for MP2 gradients (only GPW)
!> \par History
!>      11.2013 created [Mauro Del Ben]
! *****************************************************************************
MODULE mp2_cphf
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE cp_blacs_env,                    ONLY: cp_blacs_env_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_add, cp_dbcsr_allocate_matrix_set, cp_dbcsr_col_block_sizes, &
       cp_dbcsr_copy, cp_dbcsr_create, cp_dbcsr_distribution, cp_dbcsr_init, &
       cp_dbcsr_p_type, cp_dbcsr_release, cp_dbcsr_row_block_sizes, &
       cp_dbcsr_scale, cp_dbcsr_set, dbcsr_type_symmetric
  USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                             copy_fm_to_dbcsr
  USE cp_fm_basic_linalg,              ONLY: cp_fm_upper_to_full
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_info,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_set_all,&
                                             cp_fm_to_fm_submat,&
                                             cp_fm_type
  USE cp_gemm_interface,               ONLY: cp_gemm
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE hfx_energy_potential,            ONLY: integrate_four_center
  USE hfx_types,                       ONLY: alloc_containers,&
                                             hfx_container_type,&
                                             hfx_init_container,&
                                             hfx_type
  USE input_constants,                 ONLY: hfx_do_eval_energy
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kahan_sum,                       ONLY: accurate_sum
  USE kinds,                           ONLY: dp
  USE linear_systems,                  ONLY: solve_system
  USE machine,                         ONLY: m_walltime
  USE mathconstants,                   ONLY: fourpi
  USE message_passing,                 ONLY: mp_sum
  USE mp2_types,                       ONLY: mp2_type
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_methods,                      ONLY: pw_axpy,&
                                             pw_copy,&
                                             pw_derive,&
                                             pw_integral_ab,&
                                             pw_scale,&
                                             pw_transfer
  USE pw_poisson_methods,              ONLY: pw_poisson_solve
  USE pw_poisson_types,                ONLY: pw_poisson_type
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_p_type,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_p_type
  USE qs_collocate_density,            ONLY: calculate_rho_elec
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_integrate_potential,          ONLY: integrate_v_core_rspace,&
                                             integrate_v_rspace
  USE qs_ks_types,                     ONLY: qs_ks_env_type,&
                                             set_ks_env
  USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
  USE qs_rho_methods,                  ONLY: qs_rho_rebuild
  USE qs_rho_types,                    ONLY: qs_rho_create,&
                                             qs_rho_get,&
                                             qs_rho_release,&
                                             qs_rho_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE virial_types,                    ONLY: virial_type

  !$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num, omp_get_num_threads

#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mp2_cphf'

  PUBLIC :: solve_z_vector_eq


  CONTAINS

! *****************************************************************************
!> \brief Solve Z-vector equations necessary for the calculation of the MP2
!>        gradients, in order to be consistent here the parameters for the 
!>        calculation of the CPHF like updats have to be exactly equal to the 
!>        SCF case
!> \param qs_env ...
!> \param mp2_env ...
!> \param para_env ...
!> \param dft_control ...
!> \param atomic_kind_set ...
!> \param mo_coeff ...
!> \param nmo ...
!> \param homo ...
!> \param Eigenval ...
!> \param unit_nr ...
!> \param Eigenval_beta ...
!> \param homo_beta ...
!> \param mo_coeff_beta ...
!> \author Mauro Del Ben, Vladimir Rybkin
! *****************************************************************************
  SUBROUTINE solve_z_vector_eq(qs_env,mp2_env,para_env,dft_control,&
                               atomic_kind_set,mo_coeff,nmo,homo,Eigenval,unit_nr,&
                               Eigenval_beta,homo_beta,mo_coeff_beta)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(mp2_type), POINTER                  :: mp2_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    INTEGER                                  :: nmo, homo
    REAL(KIND=dp), DIMENSION(:)              :: Eigenval
    INTEGER                                  :: unit_nr
    REAL(KIND=dp), DIMENSION(:), OPTIONAL    :: Eigenval_beta
    INTEGER, OPTIONAL                        :: homo_beta
    TYPE(cp_fm_type), OPTIONAL, POINTER      :: mo_coeff_beta

    CHARACTER(LEN=*), PARAMETER :: routineN = 'solve_z_vector_eq', &
      routineP = moduleN//':'//routineN

    INTEGER :: alpha, beta, bin, dimen, handle, handle2, i, i_global, &
      i_thread, iiB, ikind, irep, ispin, j_global, jjB, my_bin_size, &
      n_rep_hf, n_threads, ncol_local, nrow_local, transf_type_in, &
      transf_type_out, virtual, virtual_beta
    INTEGER, DIMENSION(3)                    :: comp
    INTEGER, DIMENSION(:), POINTER           :: col_indices, row_indices
    LOGICAL :: alpha_beta, do_dynamic_load_balancing, do_hfx, failure, &
      hfx_treat_lsd_in_core, use_virial
    REAL(KIND=dp)                            :: e_hartree, factor, out_alpha, &
                                                pair_energy, tot_rho_r
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: mat_deb
    REAL(KIND=dp), DIMENSION(3, 3)           :: h_stress
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    TYPE(cp_dbcsr_p_type)                    :: P_mu_nu
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: mat_mu_nu, matrix_ks, &
                                                matrix_p_mp2, matrix_s, &
                                                matrix_w_mp2, rho_ao
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER :: fm_back, fm_G_mu_nu, L_jb, L_jb_beta, &
      mo_coeff_o, mo_coeff_o_beta, mo_coeff_v, mo_coeff_v_beta, P_ia, &
      P_ia_beta, P_mo, P_mo_beta, W_mo, W_mo_beta
    TYPE(hfx_container_type), DIMENSION(:), &
      POINTER                                :: integral_containers
    TYPE(hfx_container_type), POINTER        :: maxval_container
    TYPE(hfx_type), POINTER                  :: actual_x_data
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type)                          :: dvg(3), pot_g, rho_g, rho_r, &
                                                temp_pw_g
    TYPE(pw_p_type), POINTER                 :: rho_core
    TYPE(pw_poisson_type), POINTER           :: poisson_env
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(qs_rho_type), POINTER               :: rho, rho_work
    TYPE(section_vals_type), POINTER         :: hfx_sections, input
    TYPE(virial_type), POINTER               :: virial

! Set alpha-beta 

    alpha_beta = .FALSE.
    IF (PRESENT(homo_beta) .AND. PRESENT(Eigenval_beta)  &
    .AND. PRESENT(mo_coeff_beta) ) alpha_beta = .TRUE.

    CALL timeset(routineN,handle)
    failure=.FALSE.

    ! start collecting stuff
    dimen=nmo
    virtual=dimen-homo
    IF (alpha_beta) virtual_beta = dimen - homo_beta
    NULLIFY(input,pw_env,matrix_s,blacs_env,rho,energy,force,virial,matrix_w_mp2,&
            matrix_p_mp2,matrix_ks,rho_core,sab_orb)
    CALL get_qs_env(qs_env,&
                    ks_env=ks_env,&
                    pw_env=pw_env,&
                    input=input,&
                    matrix_s=matrix_s,&
                    matrix_ks= matrix_ks,&
                    matrix_p_mp2=matrix_p_mp2,&
                    matrix_w_mp2=matrix_w_mp2,&
                    blacs_env=blacs_env,&
                    rho=rho,&
                    energy=energy,&
                    force=force,&
                    virial=virial,&
                    rho_core=rho_core,&
                    sab_orb=sab_orb)

    CALL qs_rho_get(rho, rho_ao=rho_ao)

    ! check if we have to calculate the virial
    use_virial = virial%pv_availability.AND.(.NOT.virial%pv_numer)

    ! mp2 matrices
    NULLIFY(P_mo, W_mo, L_jb)
    P_mo => mp2_env%ri_grad%P_mo
    W_mo => mp2_env%ri_grad%W_mo
    L_jb => mp2_env%ri_grad%L_jb
    IF (alpha_beta) THEN
       NULLIFY(P_mo_beta, W_mo_beta, L_jb_beta)
       P_mo_beta => mp2_env%ri_grad%P_mo_beta
       W_mo_beta => mp2_env%ri_grad%W_mo_beta
       L_jb_beta => mp2_env%ri_grad%L_jb_beta
    ENDIF

    ! pw stuff
    NULLIFY(poisson_env,pw_pools,auxbas_pw_pool)
    CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool,&
                    pw_pools=pw_pools, poisson_env=poisson_env)

    ! get some of the grids ready
    NULLIFY(rho_r%pw,rho_g%pw,pot_g%pw)
    CALL pw_pool_create_pw(auxbas_pw_pool,rho_r%pw,&
                           use_data=REALDATA3D,&
                           in_space=REALSPACE)
    CALL pw_pool_create_pw(auxbas_pw_pool,rho_g%pw,&
                           use_data=COMPLEXDATA1D,&
                           in_space=RECIPROCALSPACE)
    CALL pw_pool_create_pw(auxbas_pw_pool,pot_g%pw,&
                           use_data=COMPLEXDATA1D,&
                           in_space=RECIPROCALSPACE)

    ! hfx section
    NULLIFY(hfx_sections)
    hfx_sections => section_vals_get_subs_vals(input,"DFT%XC%HF")
    CALL section_vals_get(hfx_sections,explicit=do_hfx,n_repetition=n_rep_hf)
    IF( do_hfx ) THEN
      CALL section_vals_val_get(hfx_sections, "TREAT_LSD_IN_CORE", l_val=hfx_treat_lsd_in_core,&
                                i_rep_section=1)
    END IF

    ! create work array
    NULLIFY(mat_mu_nu)
    CALL cp_dbcsr_allocate_matrix_set(mat_mu_nu, dft_control%nspins)
    DO ispin=1,dft_control%nspins
      ALLOCATE(mat_mu_nu(ispin)%matrix)
      CALL cp_dbcsr_init(mat_mu_nu(ispin)%matrix)
      CALL cp_dbcsr_create(matrix=mat_mu_nu(ispin)%matrix,&
           name="T_mu_nu",&
           dist=cp_dbcsr_distribution(matrix_s(1)%matrix), matrix_type=dbcsr_type_symmetric,&
           row_blk_size=cp_dbcsr_row_block_sizes(matrix_s(1)%matrix),&
           col_blk_size=cp_dbcsr_col_block_sizes(matrix_s(1)%matrix),&
           nze=0)
      CALL cp_dbcsr_alloc_block_from_nbl(mat_mu_nu(ispin)%matrix,sab_orb)
      CALL cp_dbcsr_set(mat_mu_nu(ispin)%matrix,0.0_dp)
    END DO

    ALLOCATE(P_mu_nu%matrix)
    CALL cp_dbcsr_init(P_mu_nu%matrix)
    ! CALL cp_dbcsr_create(P_mu_nu%matrix,template=mat_mu_nu(1)%matrix)
    ! CALL cp_dbcsr_copy(P_mu_nu%matrix,mat_mu_nu(1)%matrix,name="P_mu_nu")
    CALL cp_dbcsr_copy(P_mu_nu%matrix,rho_ao(1)%matrix,name="P_mu_nu")
    CALL cp_dbcsr_set(P_mu_nu%matrix,0.0_dp)

    NULLIFY(fm_G_mu_nu, fm_struct_tmp)
    CALL cp_fm_struct_create(fm_struct_tmp,para_env=para_env,context=blacs_env, &
                             nrow_global=dimen,ncol_global=dimen)
    CALL cp_fm_create(fm_G_mu_nu, fm_struct_tmp, name="G_mu_nu")
    CALL cp_fm_create(fm_back, fm_struct_tmp, name="fm_back")
    CALL cp_fm_struct_release(fm_struct_tmp)
    CALL cp_fm_set_all(fm_G_mu_nu, 0.0_dp)
    CALL cp_fm_set_all(fm_back, 0.0_dp)

    NULLIFY(mo_coeff_o, fm_struct_tmp)
    CALL cp_fm_struct_create(fm_struct_tmp,para_env=para_env,context=blacs_env, &
                             nrow_global=dimen,ncol_global=homo)
    CALL cp_fm_create(mo_coeff_o, fm_struct_tmp, name="mo_coeff_o")
    CALL cp_fm_struct_release(fm_struct_tmp)
    CALL cp_fm_set_all(mo_coeff_o, 0.0_dp)
    CALL cp_fm_to_fm_submat(msource=mo_coeff, mtarget=mo_coeff_o, &
                            nrow=dimen, ncol=homo, &
                            s_firstrow=1, s_firstcol=1, &
                            t_firstrow=1, t_firstcol=1)
    
    NULLIFY(mo_coeff_v, fm_struct_tmp)
    CALL cp_fm_struct_create(fm_struct_tmp,para_env=para_env,context=blacs_env, &
                             nrow_global=dimen,ncol_global=virtual)
    CALL cp_fm_create(mo_coeff_v, fm_struct_tmp, name="mo_coeff_v")
    CALL cp_fm_struct_release(fm_struct_tmp)
    CALL cp_fm_set_all(mo_coeff_v, 0.0_dp)
    CALL cp_fm_to_fm_submat(msource=mo_coeff, mtarget=mo_coeff_v, &
                            nrow=dimen, ncol=virtual, & 
                            s_firstrow=1, s_firstcol=homo+1, &
                            t_firstrow=1, t_firstcol=1)
    
    IF (alpha_beta) THEN
       NULLIFY(mo_coeff_o_beta, fm_struct_tmp)
       CALL cp_fm_struct_create(fm_struct_tmp,para_env=para_env,context=blacs_env, &
                                nrow_global=dimen,ncol_global=homo_beta)
       CALL cp_fm_create(mo_coeff_o_beta, fm_struct_tmp, name="mo_coeff_o_beta")
       CALL cp_fm_struct_release(fm_struct_tmp)
       CALL cp_fm_set_all(mo_coeff_o_beta, 0.0_dp)
       CALL cp_fm_to_fm_submat(msource=mo_coeff_beta, mtarget=mo_coeff_o_beta, &
                               nrow=dimen, ncol=homo_beta, &
                               s_firstrow=1, s_firstcol=1, &
                               t_firstrow=1, t_firstcol=1)
       
       NULLIFY(mo_coeff_v_beta, fm_struct_tmp)
       CALL cp_fm_struct_create(fm_struct_tmp,para_env=para_env,context=blacs_env, &
                                nrow_global=dimen,ncol_global=virtual_beta)
       CALL cp_fm_create(mo_coeff_v_beta, fm_struct_tmp, name="mo_coeff_v_beta")
       CALL cp_fm_struct_release(fm_struct_tmp)
       CALL cp_fm_set_all(mo_coeff_v_beta, 0.0_dp)
       CALL cp_fm_to_fm_submat(msource=mo_coeff_beta, mtarget=mo_coeff_v_beta, &
                               nrow=dimen, ncol=virtual_beta, & 
                               s_firstrow=1, s_firstcol=homo_beta+1, &
                               t_firstrow=1, t_firstcol=1)
    ENDIF

    ! create a working rho environment
    NULLIFY(rho_work)
    CALL qs_rho_create(rho_work)
    CALL qs_rho_rebuild(rho=rho_work, qs_env=qs_env, rebuild_ao=.TRUE., rebuild_grids=.FALSE.)

    ! here we check if we have to reallocate the HFX container
    IF(mp2_env%ri_mp2%free_hfx_buffer) THEN
      CALL timeset(routineN//"_alloc_hfx",handle2)
      n_threads = 1
!$  n_threads = omp_get_max_threads()

      DO irep = 1, n_rep_hf
        DO i_thread = 0, n_threads-1
          actual_x_data => qs_env%x_data(irep, i_thread + 1)

          do_dynamic_load_balancing = .TRUE.
          IF( n_threads == 1 .OR. actual_x_data%memory_parameter%do_disk_storage ) do_dynamic_load_balancing = .FALSE.

          IF( do_dynamic_load_balancing ) THEN
            my_bin_size = SIZE(actual_x_data%distribution_energy)
          ELSE
            my_bin_size = 1
          END IF

          IF(.NOT. actual_x_data%memory_parameter%do_all_on_the_fly) THEN
           !  CALL dealloc_containers(actual_x_data, hfx_do_eval_energy)
            CALL alloc_containers(actual_x_data, my_bin_size, hfx_do_eval_energy)

            DO bin=1, my_bin_size
              maxval_container => actual_x_data%maxval_container(bin)
              integral_containers => actual_x_data%integral_containers(:,bin)
              CALL hfx_init_container(maxval_container, actual_x_data%memory_parameter%actual_memory_usage, .FALSE.)
              DO i=1,64
                CALL hfx_init_container(integral_containers(i), actual_x_data%memory_parameter%actual_memory_usage, .FALSE.)
              END DO
            END DO
          END IF
        END DO
      END DO
      CALL timestop(handle2)
    END IF

    ! not a good idea having screening on initial P
    IF(qs_env%x_data(1,1)%screening_parameter%do_initial_p_screening) THEN
      CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
                     "CPHF-like update of the Hartree Fock exchange part requested with the use "//&
                     "of SCREEN_ON_INITIAL_P. This may lead to unphysical results. "//&
                     "Set SCREEN_ON_INITIAL_P to false to avoid possible problems. "//&
                     CPSourceFileRef,&
                     only_ionode=.TRUE.)
    END IF

    ! update Lagrangian with the CPHF like update, occ-occ block, first call (recompute hfx integrals if needed)
    transf_type_in=1
    transf_type_out=1
    out_alpha=0.5_dp
    ! In alpha-beta case, L_bj_alpha has Coulomb and XC alpha-alpha part 
    ! and (only) Coulomb alpha-beta part and vice versa.

    ! Complete in closed shell case, alpha-alpha (Coulomb and XC) 
    ! part of L_bj(alpha) for open shell
    factor = 4.0_dp
    IF (alpha_beta) factor = 2.0_dp
    CALL cphf_like_update(qs_env,para_env,homo,virtual,dimen,&
                          mo_coeff,mo_coeff_o,mo_coeff_v,Eigenval,&
                          hfx_sections,energy,n_rep_hf,poisson_env,&
                          rho_work,pot_g,rho_g,rho_r,mat_mu_nu,P_mu_nu,&
                          P_mo,fm_G_mu_nu,fm_back,transf_type_in,out_alpha,&
                          L_jb,transf_type_out,&
                          recalc_hfx_integrals=mp2_env%ri_mp2%free_hfx_buffer,&
                          factor=factor)
    IF (alpha_beta) THEN
       ! Alpha-beta (Coulomb) part of L_jb(alpha)
       CALL cphf_like_update(qs_env,para_env,homo,virtual,dimen,&
                             mo_coeff,mo_coeff_o,mo_coeff_v,Eigenval,&
                             hfx_sections,energy,n_rep_hf,poisson_env,&
                             rho_work,pot_g,rho_g,rho_r,mat_mu_nu,P_mu_nu,&
                             P_mo_beta,fm_G_mu_nu,fm_back,transf_type_in,out_alpha,&
                             L_jb,transf_type_out,&
                             mo_coeff_beta=mo_coeff_beta, &
                             mo_coeff_o_beta=mo_coeff_o_beta, &
                             mo_coeff_v_beta=mo_coeff_v_beta, &
                             homo_beta=homo_beta, virtual_beta=virtual_beta,&
                             factor=factor)
       ! Beta-beta (Coulomb and XC) part of L_jb(beta)
       CALL cphf_like_update(qs_env,para_env,homo_beta,virtual_beta,dimen, &
                              mo_coeff_beta,mo_coeff_o_beta,mo_coeff_v_beta, &
                             Eigenval_beta,hfx_sections,energy,n_rep_hf,&
                             poisson_env,rho_work,pot_g,rho_g,rho_r,mat_mu_nu,&
                             P_mu_nu,P_mo_beta,fm_G_mu_nu,fm_back,transf_type_in,&
                             out_alpha,L_jb_beta,transf_type_out,factor=factor)
       ! Alpha-beta (Coulomb) part of L_jb(beta)
       CALL cphf_like_update(qs_env,para_env,homo_beta,virtual_beta,dimen,&
                             mo_coeff_beta,mo_coeff_o_beta,mo_coeff_v_beta, &
                             Eigenval_beta,&
                             hfx_sections,energy,n_rep_hf,poisson_env,&
                             rho_work,pot_g,rho_g,rho_r,mat_mu_nu,P_mu_nu,&
                             P_mo,fm_G_mu_nu,fm_back,transf_type_in,out_alpha,&
                             L_jb_beta,transf_type_out,&
                             mo_coeff_beta=mo_coeff, &
                             mo_coeff_o_beta=mo_coeff_o, &
                             mo_coeff_v_beta=mo_coeff_v, &
                             homo_beta=homo, virtual_beta=virtual,factor=factor)
    ENDIF

    ! update Lagrangian with the CPHF like update, virt-virt block
    transf_type_in=2
    transf_type_out=1
    out_alpha=0.5_dp
    ! Complete in closed shell case, alpha-alpha (Coulomb and XC) 
    ! part of L_bj(alpha) for open shell
    CALL cphf_like_update(qs_env,para_env,homo,virtual,dimen,&
                          mo_coeff,mo_coeff_o,mo_coeff_v,Eigenval,&
                          hfx_sections,energy,n_rep_hf,poisson_env,&
                          rho_work,pot_g,rho_g,rho_r,mat_mu_nu,P_mu_nu,&
                          P_mo,fm_G_mu_nu,fm_back,transf_type_in,out_alpha,&
                          L_jb,transf_type_out,factor=factor)
    IF (alpha_beta) THEN
       ! Alpha-beta (Coulomb) part of L_jb(alpha)
       CALL cphf_like_update(qs_env,para_env,homo,virtual,dimen,&
                             mo_coeff,mo_coeff_o,mo_coeff_v,Eigenval,&
                             hfx_sections,energy,n_rep_hf,poisson_env,&
                             rho_work,pot_g,rho_g,rho_r,mat_mu_nu,P_mu_nu,&
                             P_mo_beta,fm_G_mu_nu,fm_back,transf_type_in,out_alpha,&
                             L_jb,transf_type_out,&
                             mo_coeff_beta=mo_coeff_beta, &
                             mo_coeff_o_beta=mo_coeff_o_beta, &
                             mo_coeff_v_beta=mo_coeff_v_beta, &
                             homo_beta=homo_beta, virtual_beta=virtual_beta,&
                             factor=factor)
       ! Beta-beta (Coulomb and XC) part of L_jb(beta)
       CALL cphf_like_update(qs_env,para_env,homo_beta,virtual_beta,dimen, &
                              mo_coeff_beta,mo_coeff_o_beta,mo_coeff_v_beta, &
                             Eigenval_beta,hfx_sections,energy,n_rep_hf,&
                             poisson_env,rho_work,pot_g,rho_g,rho_r,mat_mu_nu,&
                             P_mu_nu,P_mo_beta,fm_G_mu_nu,fm_back,transf_type_in,&
                             out_alpha,L_jb_beta,transf_type_out,factor=factor)
       ! Alpha-beta (Coulomb) part of L_jb(beta)
       CALL cphf_like_update(qs_env,para_env,homo_beta,virtual_beta,dimen,&
                             mo_coeff_beta,mo_coeff_o_beta,mo_coeff_v_beta, &
                             Eigenval_beta,&
                             hfx_sections,energy,n_rep_hf,poisson_env,&
                             rho_work,pot_g,rho_g,rho_r,mat_mu_nu,P_mu_nu,&
                             P_mo,fm_G_mu_nu,fm_back,transf_type_in,out_alpha,&
                             L_jb_beta,transf_type_out,&
                             mo_coeff_beta=mo_coeff, &
                             mo_coeff_o_beta=mo_coeff_o, &
                             mo_coeff_v_beta=mo_coeff_v, &
                             homo_beta=homo, virtual_beta=virtual,factor=factor)
    ENDIF
    ! at this point Lagrnagian is completed ready to solve the Z-vector equations
    ! P_ia will contain the solution of these equations
    NULLIFY(P_ia, fm_struct_tmp)
    CALL cp_fm_struct_create(fm_struct_tmp,para_env=para_env,context=blacs_env, &
                             nrow_global=homo,ncol_global=virtual)
    CALL cp_fm_create(P_ia, fm_struct_tmp, name="P_ia")
    CALL cp_fm_struct_release(fm_struct_tmp)
    CALL cp_fm_set_all(P_ia, 0.0_dp)
    IF (alpha_beta) THEN
       NULLIFY(P_ia_beta, fm_struct_tmp)
       CALL cp_fm_struct_create(fm_struct_tmp,para_env=para_env,context=blacs_env, &
                         nrow_global=homo_beta,ncol_global=virtual_beta)
       CALL cp_fm_create(P_ia_beta, fm_struct_tmp, name="P_ia_beta")
       CALL cp_fm_struct_release(fm_struct_tmp)
       CALL cp_fm_set_all(P_ia_beta, 0.0_dp)
    ENDIF

    IF (.NOT. alpha_beta) THEN
       CALL solve_z_vector_eq_low(qs_env,mp2_env,para_env,homo,virtual,dimen,unit_nr,&
                               mo_coeff,mo_coeff_o,mo_coeff_v,Eigenval,blacs_env, &
                               hfx_sections,energy,n_rep_hf, &
                               poisson_env,rho_work,pot_g,rho_g,rho_r,mat_mu_nu, &
                               P_mu_nu,L_jb,fm_G_mu_nu,fm_back,P_ia)
    ELSE
       CALL solve_z_vector_eq_low(qs_env,mp2_env,para_env,homo,virtual,dimen,unit_nr,&
                                   mo_coeff,mo_coeff_o,mo_coeff_v,Eigenval,blacs_env,&
                                    hfx_sections,energy,n_rep_hf,&
                                   poisson_env,rho_work,pot_g,rho_g,rho_r,mat_mu_nu,&
                                   P_mu_nu,L_jb,fm_G_mu_nu,fm_back,P_ia,&
                                   homo_beta,Eigenval_beta,P_ia_beta, mo_coeff_beta,&
                                   mo_coeff_o_beta, mo_coeff_v_beta, L_jb_beta)
    ENDIF
    
    ! release Lagrangian
    CALL cp_fm_release(L_jb)
    IF (alpha_beta) CALL cp_fm_release(L_jb_beta)

    ! update the MP2-MO density matrix with the occ-virt block
    CALL cp_fm_to_fm_submat(msource=P_ia, mtarget=P_mo, &
                            nrow=homo, ncol=virtual, &
                            s_firstrow=1, s_firstcol=1, &
                            t_firstrow=1, t_firstcol=homo+1)
    CALL cp_fm_release(P_ia)
    ! transpose P_MO matrix (easy way to symmetrize)
    CALL cp_fm_set_all(fm_back, 0.0_dp)
    ! P_mo now is ready
    CALL cp_fm_upper_to_full(matrix=P_mo, work=fm_back)
    IF (alpha_beta) THEN
       ! update the MP2-MO density matrix with the occ-virt block
       CALL cp_fm_to_fm_submat(msource=P_ia_beta, mtarget=P_mo_beta, &
                               nrow=homo_beta, ncol=virtual_beta, &
                               s_firstrow=1, s_firstcol=1, &
                               t_firstrow=1, t_firstcol=homo_beta+1)
       CALL cp_fm_release(P_ia_beta)
       ! transpose P_MO matrix (easy way to symmetrize)
       CALL cp_fm_set_all(fm_back, 0.0_dp)
       ! P_mo now is ready
       CALL cp_fm_upper_to_full(matrix=P_mo_beta, work=fm_back)
    ENDIF

    ! do the final update to MP2 energy weighted matrix W_MO
    CALL cp_fm_get_info(matrix=W_mo,&
                        nrow_local=nrow_local,&
                        ncol_local=ncol_local,&
                        row_indices=row_indices,&
                        col_indices=col_indices)
    DO jjB=1, ncol_local
      j_global=col_indices(jjB)
      IF(j_global<=homo) THEN
        DO iiB=1, nrow_local
          i_global=row_indices(iiB)
          W_mo%local_data(iiB,jjB)=W_mo%local_data(iiB,jjB)-P_mo%local_data(iiB,jjB)*Eigenval(j_global)
        END DO
      ELSE
        DO iiB=1, nrow_local
          i_global=row_indices(iiB)
          IF(i_global<=homo) THEN
            ! virt-occ
            W_mo%local_data(iiB,jjB)=W_mo%local_data(iiB,jjB)-P_mo%local_data(iiB,jjB)*Eigenval(i_global)
          ELSE
            ! virt-virt
            W_mo%local_data(iiB,jjB)=W_mo%local_data(iiB,jjB)-P_mo%local_data(iiB,jjB)*Eigenval(j_global)
          END IF
        END DO
      END IF
    END DO
    IF (alpha_beta) THEN
       CALL cp_fm_get_info(matrix=W_mo_beta,&
                           nrow_local=nrow_local,&
                           ncol_local=ncol_local,&
                           row_indices=row_indices,&
                           col_indices=col_indices)
       DO jjB=1, ncol_local
         j_global=col_indices(jjB)
         IF(j_global<=homo_beta) THEN
           DO iiB=1, nrow_local
             i_global=row_indices(iiB)
             W_mo_beta%local_data(iiB,jjB)= &
             W_mo_beta%local_data(iiB,jjB)-P_mo_beta%local_data(iiB,jjB)*Eigenval_beta(j_global)
           END DO
         ELSE
           DO iiB=1, nrow_local
             i_global=row_indices(iiB)
             IF(i_global<=homo_beta) THEN
               ! virt-occ
               W_mo_beta%local_data(iiB,jjB)= &
               W_mo_beta%local_data(iiB,jjB)-P_mo_beta%local_data(iiB,jjB)*Eigenval_beta(i_global)
             ELSE
               ! virt-virt
               W_mo_beta%local_data(iiB,jjB)= &
               W_mo_beta%local_data(iiB,jjB)-P_mo_beta%local_data(iiB,jjB)*Eigenval_beta(j_global)
             END IF
           END DO
         END IF
       END DO
    ENDIF

    ! complete the occ-occ block of W_mo with a CPHF like update
    transf_type_in=4
    transf_type_out=2
    out_alpha=-0.5_dp
    ! Complete in closed shell case, alpha-alpha (Coulomb and XC) 
    ! part of W_mo(III)(alpha) for open shell
    CALL cphf_like_update(qs_env,para_env,homo,virtual,dimen,&
                          mo_coeff,mo_coeff_o,mo_coeff_v,Eigenval,&
                          hfx_sections,energy,n_rep_hf,poisson_env,&
                          rho_work,pot_g,rho_g,rho_r,mat_mu_nu,P_mu_nu,&
                          P_mo,fm_G_mu_nu,fm_back,transf_type_in,out_alpha,&
                          W_mo,transf_type_out,factor=factor)
    IF (alpha_beta) THEN
       ! Alpha-beta (Coulomb) part of W_mo(III)(alpha)
       CALL cphf_like_update(qs_env,para_env,homo,virtual,dimen,&
                             mo_coeff,mo_coeff_o,mo_coeff_v,Eigenval,&
                             hfx_sections,energy,n_rep_hf,poisson_env,&
                             rho_work,pot_g,rho_g,rho_r,mat_mu_nu,P_mu_nu,&
                             P_mo_beta,fm_G_mu_nu,fm_back,transf_type_in,out_alpha,&
                             W_mo,transf_type_out,&
                             mo_coeff_beta=mo_coeff_beta, &
                             mo_coeff_o_beta=mo_coeff_o_beta, &
                             mo_coeff_v_beta=mo_coeff_v_beta, &
                             homo_beta=homo_beta, virtual_beta=virtual_beta,&
                             factor=factor)
       ! Beta-beta (Coulomb and XC) part of W_mo(III)(beta)
       CALL cphf_like_update(qs_env,para_env,homo_beta,virtual_beta,dimen, &
                              mo_coeff_beta,mo_coeff_o_beta,mo_coeff_v_beta, &
                             Eigenval_beta,hfx_sections,energy,n_rep_hf,&
                             poisson_env,rho_work,pot_g,rho_g,rho_r,mat_mu_nu,&
                             P_mu_nu,P_mo_beta,fm_G_mu_nu,fm_back,transf_type_in,&
                             out_alpha,W_mo_beta,transf_type_out,factor=factor)
       ! Alpha-beta (Coulomb) part of W_mo(III)(beta)
       CALL cphf_like_update(qs_env,para_env,homo_beta,virtual_beta,dimen,&
                             mo_coeff_beta,mo_coeff_o_beta,mo_coeff_v_beta, &
                             Eigenval_beta,&
                             hfx_sections,energy,n_rep_hf,poisson_env,&
                             rho_work,pot_g,rho_g,rho_r,mat_mu_nu,P_mu_nu,&
                             P_mo,fm_G_mu_nu,fm_back,transf_type_in,out_alpha,&
                             W_mo_beta,transf_type_out,&
                             mo_coeff_beta=mo_coeff, &
                             mo_coeff_o_beta=mo_coeff_o, &
                             mo_coeff_v_beta=mo_coeff_v, &
                             homo_beta=homo, virtual_beta=virtual,factor=factor)
    ENDIF

    ! release DBCSR stuff
    DO ispin=1, dft_control%nspins
      CALL cp_dbcsr_release(mat_mu_nu(ispin)%matrix)
      DEALLOCATE(mat_mu_nu(ispin)%matrix)
    END DO
    DEALLOCATE(mat_mu_nu)
    CALL cp_dbcsr_release(P_mu_nu%matrix)
    DEALLOCATE(P_mu_nu%matrix)
    ! release fm stuff
    CALL cp_fm_release(fm_G_mu_nu)
    CALL cp_fm_release(mo_coeff_o)
    CALL cp_fm_release(mo_coeff_v)
    IF (alpha_beta) THEN
       CALL cp_fm_release(mo_coeff_o_beta)
       CALL cp_fm_release(mo_coeff_v_beta)
    ENDIF
    ! release rho stuff
    CALL qs_rho_release(rho_struct=rho_work)
 
    IF(.FALSE.) THEN
      ALLOCATE(mat_deb(dimen,dimen))
      mat_deb=0.0_dp
      CALL cp_fm_get_info(matrix=P_mo,&
                          nrow_local=nrow_local,&
                          ncol_local=ncol_local,&
                          row_indices=row_indices,&
                          col_indices=col_indices)
      DO jjB=1, ncol_local
        j_global=col_indices(jjB)
        DO iiB=1, nrow_local
          i_global=row_indices(iiB)
          mat_deb(i_global,j_global)=P_mo%local_data(iiB,jjB)
        END DO
      END DO
      CALL mp_sum(mat_deb,para_env%group)
      IF(para_env%mepos==0) CALL write_array(mat_deb(1:dimen,1:dimen))
      mat_deb=0.0_dp
      CALL cp_fm_get_info(matrix=W_mo,&
                          nrow_local=nrow_local,&
                          ncol_local=ncol_local,&
                          row_indices=row_indices,&
                          col_indices=col_indices)
      DO jjB=1, ncol_local
        j_global=col_indices(jjB)
        DO iiB=1, nrow_local
          i_global=row_indices(iiB)
          mat_deb(i_global,j_global)=W_mo%local_data(iiB,jjB)
        END DO
      END DO
      CALL mp_sum(mat_deb,para_env%group)
      IF(para_env%mepos==0) CALL write_array(mat_deb(1:dimen,1:dimen))
      DEALLOCATE(mat_deb)
      IF (alpha_beta) THEN
         ALLOCATE(mat_deb(dimen,dimen))
         mat_deb=0.0_dp
         CALL cp_fm_get_info(matrix=P_mo_beta,&
                             nrow_local=nrow_local,&
                             ncol_local=ncol_local,&
                             row_indices=row_indices,&
                             col_indices=col_indices)
         DO jjB=1, ncol_local
           j_global=col_indices(jjB)
           DO iiB=1, nrow_local
             i_global=row_indices(iiB)
             mat_deb(i_global,j_global)=P_mo_beta%local_data(iiB,jjB)
           END DO
         END DO
         CALL mp_sum(mat_deb,para_env%group)
         IF(para_env%mepos==0) CALL write_array(mat_deb(1:dimen,1:dimen))
         mat_deb=0.0_dp
         CALL cp_fm_get_info(matrix=W_mo_beta,&
                             nrow_local=nrow_local,&
                             ncol_local=ncol_local,&
                             row_indices=row_indices,&
                             col_indices=col_indices)
         DO jjB=1, ncol_local
           j_global=col_indices(jjB)
           DO iiB=1, nrow_local
             i_global=row_indices(iiB)
             mat_deb(i_global,j_global)=W_mo_beta%local_data(iiB,jjB)
           END DO
         END DO
         CALL mp_sum(mat_deb,para_env%group)
         IF(para_env%mepos==0) CALL write_array(mat_deb(1:dimen,1:dimen))
         DEALLOCATE(mat_deb)
      ENDIF
    END IF

    ! backtransform into AO basis, since P_mo and W_mo
    ! are symmetric (in principle), no need to symmetrize
    ! first W_mo
    CALL cp_fm_set_all(fm_back, 0.0_dp)
    CALL cp_gemm('N','N',dimen,dimen,dimen,1.0_dp,&
                    mo_coeff,W_mo,0.0_dp,fm_back,&
                    a_first_row=1,&
                    b_first_col=1,&
                    b_first_row=1,&
                    c_first_col=1,&
                    c_first_row=1)
    CALL cp_fm_set_all(W_mo, 0.0_dp)
    CALL cp_gemm('N','T',dimen,dimen,dimen,1.0_dp,&
                    fm_back,mo_coeff,0.0_dp,W_mo,&
                    a_first_row=1,&
                    b_first_col=1,&
                    b_first_row=1,&
                    c_first_col=1,&
                    c_first_row=1)
    ! and P_mo
    CALL cp_fm_set_all(fm_back, 0.0_dp)
    CALL cp_gemm('N','N',dimen,dimen,dimen,1.0_dp,&
                    mo_coeff,P_mo,0.0_dp,fm_back,&
                    a_first_row=1,&
                    b_first_col=1,&
                    b_first_row=1,&
                    c_first_col=1,&
                    c_first_row=1)
    CALL cp_fm_set_all(P_mo, 0.0_dp)
    CALL cp_gemm('N','T',dimen,dimen,dimen,1.0_dp,&
                    fm_back,mo_coeff,0.0_dp,P_mo,&
                    a_first_row=1,&
                    b_first_col=1,&
                    b_first_row=1,&
                    c_first_col=1,&
                    c_first_row=1)
    IF (alpha_beta) THEN
       CALL cp_fm_set_all(fm_back, 0.0_dp)
       CALL cp_gemm('N','N',dimen,dimen,dimen,1.0_dp,&
                       mo_coeff_beta,W_mo_beta,0.0_dp,fm_back,&
                       a_first_row=1,&
                       b_first_col=1,&
                       b_first_row=1,&
                       c_first_col=1,&
                       c_first_row=1)
       CALL cp_fm_set_all(W_mo_beta, 0.0_dp)
       CALL cp_gemm('N','T',dimen,dimen,dimen,1.0_dp,&
                       fm_back,mo_coeff_beta,0.0_dp,W_mo_beta,&
                       a_first_row=1,&
                       b_first_col=1,&
                       b_first_row=1,&
                       c_first_col=1,&
                       c_first_row=1)
       ! and P_mo
       CALL cp_fm_set_all(fm_back, 0.0_dp)
       CALL cp_gemm('N','N',dimen,dimen,dimen,1.0_dp,&
                       mo_coeff_beta,P_mo_beta,0.0_dp,fm_back,&
                       a_first_row=1,&
                       b_first_col=1,&
                       b_first_row=1,&
                       c_first_col=1,&
                       c_first_row=1)
       CALL cp_fm_set_all(P_mo_beta, 0.0_dp)
       CALL cp_gemm('N','T',dimen,dimen,dimen,1.0_dp,&
                       fm_back,mo_coeff_beta,0.0_dp,P_mo_beta,&
                       a_first_row=1,&
                       b_first_col=1,&
                       b_first_row=1,&
                       c_first_col=1,&
                       c_first_row=1)
    ENDIF
 
    ! copy W_mo into dbcsr
    CALL copy_fm_to_dbcsr(W_mo, matrix_w_mp2(1)%matrix, keep_sparsity=.TRUE.)
    IF (alpha_beta) THEN
        CALL copy_fm_to_dbcsr(W_mo_beta, matrix_w_mp2(2)%matrix, keep_sparsity=.TRUE.)
    ENDIF

    ! create mp2 DBCSR density
    CALL cp_dbcsr_allocate_matrix_set(matrix_p_mp2,dft_control%nspins)
    DO ispin=1, dft_control%nspins
       ALLOCATE(matrix_p_mp2(ispin)%matrix)
       CALL cp_dbcsr_init(matrix_p_mp2(ispin)%matrix)
       CALL cp_dbcsr_copy(matrix_p_mp2(ispin)%matrix,rho_ao(ispin)%matrix,&
            name="P MATRIX MP2")
       CALL cp_dbcsr_set(matrix_p_mp2(ispin)%matrix,0.0_dp)
       IF (ispin == 1) &
          CALL copy_fm_to_dbcsr(P_mo, matrix_p_mp2(ispin)%matrix,keep_sparsity=.TRUE.)
       IF (ispin == 2) &
          CALL copy_fm_to_dbcsr(P_mo_beta, matrix_p_mp2(ispin)%matrix,keep_sparsity=.TRUE.)
    END DO
    CALL set_ks_env(ks_env, matrix_p_mp2=matrix_p_mp2)

    CALL cp_fm_release(fm_back)

    ! release remaining fm stuff
    CALL cp_fm_release(W_mo)
    CALL cp_fm_release(P_mo)
    IF (alpha_beta) THEN
       CALL cp_fm_release(W_mo_beta)
       CALL cp_fm_release(P_mo_beta)
    ENDIF

    ! update the core-forces with the MP2-density contribution
    ! put MP2 density on the grid
    IF (alpha_beta) THEN ! In alpha_beta case, get the joint density
      CALL cp_dbcsr_add(matrix_p_mp2(1)%matrix,matrix_p_mp2(2)%matrix,1.0_dp,1.0_dp)
    ENDIF
    CALL calculate_rho_elec(matrix_p=matrix_p_mp2(1)%matrix,&
                            rho=rho_r,&
                            rho_gspace=rho_g,&
                            total_rho=tot_rho_r,&
                            ks_env=ks_env,&
                            soft_valid=.FALSE.)
    ! calculate the MP2 potential
    CALL pw_transfer(rho_r%pw, rho_g%pw)
    CALL pw_poisson_solve(poisson_env,rho_g%pw, pair_energy, pot_g%pw)
    CALL pw_transfer(pot_g%pw, rho_r%pw)
    CALL pw_scale(rho_r%pw,rho_r%pw%pw_grid%dvol)
 
    ! calculate core forces
    CALL integrate_v_core_rspace(rho_r, qs_env)
    DO ikind=1,SIZE(atomic_kind_set)
       force(ikind)%mp2_sep= force(ikind)%rho_core
       force(ikind)%rho_core=0.0_dp
    ENDDO
    ! right contribution
    IF (alpha_beta) THEN ! In alpha_beta case, get the joint density
      CALL cp_dbcsr_add(rho_ao(1)%matrix,rho_ao(2)%matrix,1.0_dp,1.0_dp)
      CALL cp_dbcsr_add(matrix_ks(1)%matrix,matrix_ks(2)%matrix,1.0_dp,1.0_dp)
    ENDIF
    CALL integrate_v_rspace(v_rspace=rho_r,pmat=rho_ao(1),hmat=matrix_ks(1),&
                            qs_env=qs_env,calculate_forces=.TRUE.)

    IF (alpha_beta) THEN ! In alpha_beta case, get the initial densities back
      CALL cp_dbcsr_add(rho_ao(1)%matrix,rho_ao(2)%matrix,1.0_dp,-1.0_dp)
      CALL cp_dbcsr_add(matrix_p_mp2(1)%matrix,matrix_p_mp2(2)%matrix,1.0_dp,-1.0_dp)
      CALL cp_dbcsr_add(matrix_ks(1)%matrix,matrix_ks(2)%matrix,1.0_dp,-1.0_dp)
    ENDIF

    IF(use_virial) THEN
      ! update virial if necessery with the volume term
      ! first create pw auxilliary stuff
      CALL timeset(routineN//"_Virial",handle2)
      NULLIFY(temp_pw_g%pw)
      CALL pw_pool_create_pw(auxbas_pw_pool,temp_pw_g%pw,&
                             use_data=COMPLEXDATA1D,&
                             in_space=RECIPROCALSPACE)
      DO i=1, 3
        NULLIFY(dvg(i)%pw)
        CALL pw_pool_create_pw(auxbas_pw_pool,dvg(i)%pw,&
                               use_data=COMPLEXDATA1D,&
                               in_space=RECIPROCALSPACE)
      END DO
     
      ! make a copy of the MP2 density in G space
      CALL pw_copy(rho_g%pw, temp_pw_g%pw)
      ! calculate MP2-like-hartree potential derivatives
      DO i=1, 3
        comp=0
        comp(i)=1
        CALL pw_copy(pot_g%pw, dvg(i)%pw)
        CALL pw_derive(dvg(i)%pw, comp)
      END DO
 
      ! calculate total SCF density and potential
      CALL calculate_rho_elec(matrix_p=rho_ao(1)%matrix,&
                              rho=rho_r,&
                              rho_gspace=rho_g,&
                              total_rho=tot_rho_r,&
                              ks_env=ks_env,&
                              soft_valid=.FALSE.)
      ! and associated potential
      CALL pw_transfer(rho_r%pw, rho_g%pw)
      ! don't forget the core density
      CALL pw_axpy(rho_core%pw, rho_g%pw)
      CALL pw_poisson_solve(poisson_env,rho_g%pw, pair_energy, pot_g%pw)

      ! finally update virial with the volume contribution
      e_hartree=pw_integral_ab(temp_pw_g%pw, pot_g%pw)
      h_stress=0.0_dp
      DO alpha=1, 3
        comp=0
        comp(alpha)=1
        CALL pw_copy(pot_g%pw, rho_g%pw)
        CALL pw_derive(rho_g%pw, comp)
        h_stress(alpha,alpha)=-e_hartree
        DO beta=alpha, 3
          h_stress(alpha,beta)=h_stress(alpha,beta) &
                     -2.0_dp*pw_integral_ab(rho_g%pw, dvg(beta)%pw)/fourpi
          h_stress (beta,alpha)=h_stress(alpha,beta)
        END DO
      END DO
      virial%pv_virial = virial%pv_virial + h_stress/REAL(para_env%num_pe,dp)

      ! free stuff
      CALL pw_pool_give_back_pw(auxbas_pw_pool,temp_pw_g%pw)
      DO i=1, 3
        CALL pw_pool_give_back_pw(auxbas_pw_pool,dvg(i)%pw)
      END DO
      CALL timestop(handle2)
    END IF

    DO ispin=1, dft_control%nspins
      CALL cp_dbcsr_add(rho_ao(ispin)%matrix, matrix_p_mp2(ispin)%matrix, 1.0_dp, 1.0_dp)
    END DO

    ! release stuff
    CALL pw_pool_give_back_pw(auxbas_pw_pool,rho_r%pw)
    CALL pw_pool_give_back_pw(auxbas_pw_pool,rho_g%pw)
    CALL pw_pool_give_back_pw(auxbas_pw_pool,pot_g%pw)

    CALL timestop(handle)

  END SUBROUTINE

! *****************************************************************************
!> \brief Here we performe the CPHF like update using GPW,
!>        transf_type_in  defines the type of transformation for the matrix in input
!>        transf_type_in = 1 -> occ-occ back transformation
!>        transf_type_in = 2 -> virt-virt back transformation
!>        transf_type_in = 3 -> occ-virt back transformation including the 
!>                              eigenvalues energy differences for the diagonal elements
!>        transf_type_in = 4 -> full range
!>        transf_type_out defines the type of transformation for the matrix in output
!>        transf_type_out = 1 -> occ-vit transformation 
!>        transf_type_out = 2 -> occ-occ transformation
!> \param qs_env ...
!> \param para_env ...
!> \param homo ...
!> \param virtual ...
!> \param dimen ...
!> \param mo_coeff ...
!> \param mo_coeff_o ...
!> \param mo_coeff_v ...
!> \param Eigenval ...
!> \param hfx_sections ...
!> \param energy ...
!> \param n_rep_hf ...
!> \param poisson_env ...
!> \param rho_work ...
!> \param pot_g ...
!> \param rho_g ...
!> \param rho_r ...
!> \param mat_mu_nu ...
!> \param P_mu_nu ...
!> \param fm_mo ...
!> \param fm_ao ...
!> \param fm_back ...
!> \param transf_type_in ...
!> \param out_alpha ...
!> \param fm_mo_out ...
!> \param transf_type_out ...
!> \param recalc_hfx_integrals ...
!> \param mo_coeff_beta ...
!> \param mo_coeff_o_beta ...
!> \param mo_coeff_v_beta ...
!> \param homo_beta ...
!> \param virtual_beta ...
!> \param factor ...
!> \author Mauro Del Ben, Vladimir Rybkin
! *****************************************************************************
  SUBROUTINE cphf_like_update(qs_env,para_env,homo,virtual,dimen,&
                              mo_coeff,mo_coeff_o,mo_coeff_v,Eigenval,&
                              hfx_sections,energy,n_rep_hf,poisson_env,&
                              rho_work,pot_g,rho_g,rho_r,mat_mu_nu,P_mu_nu,&
                              fm_mo,fm_ao,fm_back,transf_type_in,out_alpha,&
                              fm_mo_out,transf_type_out,recalc_hfx_integrals, &
                              mo_coeff_beta,mo_coeff_o_beta,mo_coeff_v_beta, &
                              homo_beta, virtual_beta, factor)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER                                  :: homo, virtual, dimen
    TYPE(cp_fm_type), POINTER                :: mo_coeff, mo_coeff_o, &
                                                mo_coeff_v
    REAL(KIND=dp), DIMENSION(:)              :: Eigenval
    TYPE(section_vals_type), POINTER         :: hfx_sections
    TYPE(qs_energy_type), POINTER            :: energy
    INTEGER                                  :: n_rep_hf
    TYPE(pw_poisson_type), POINTER           :: poisson_env
    TYPE(qs_rho_type), POINTER               :: rho_work
    TYPE(pw_p_type)                          :: pot_g, rho_g, rho_r
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: mat_mu_nu
    TYPE(cp_dbcsr_p_type)                    :: P_mu_nu
    TYPE(cp_fm_type), POINTER                :: fm_mo, fm_ao, fm_back
    INTEGER                                  :: transf_type_in
    REAL(KIND=dp)                            :: out_alpha
    TYPE(cp_fm_type), POINTER                :: fm_mo_out
    INTEGER                                  :: transf_type_out
    LOGICAL, OPTIONAL                        :: recalc_hfx_integrals
    TYPE(cp_fm_type), OPTIONAL, POINTER      :: mo_coeff_beta, &
                                                mo_coeff_o_beta, &
                                                mo_coeff_v_beta
    INTEGER, OPTIONAL                        :: homo_beta, virtual_beta
    REAL(KIND=dp), OPTIONAL                  :: factor

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cphf_like_update', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i_global, iiB, irep, &
                                                j_global, jjB, ncol_local, &
                                                nrow_local, ns
    INTEGER, DIMENSION(:), POINTER           :: col_indices, row_indices
    LOGICAL                                  :: alpha_beta, failure, &
                                                my_recalc_hfx_integrals
    REAL(KIND=dp)                            :: ex_energy, pair_energy, &
                                                total_rho
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: rho_work_ao
    TYPE(cp_dbcsr_p_type), DIMENSION(:, :), &
      POINTER                                :: mat_2d, rho_ao_2d
    TYPE(qs_ks_env_type), POINTER            :: ks_env

    alpha_beta = .FALSE.
    IF (PRESENT(mo_coeff_beta)) alpha_beta = .TRUE.

    NULLIFY(ks_env, rho_work_ao)
    CALL timeset(routineN,handle)
    failure=.FALSE.

    my_recalc_hfx_integrals=.FALSE.
    IF(PRESENT(recalc_hfx_integrals)) my_recalc_hfx_integrals=recalc_hfx_integrals

    CALL get_qs_env(qs_env, ks_env=ks_env)
    ! perform back transformation
    SELECT CASE(transf_type_in)
      CASE(1)
        IF (.NOT. alpha_beta) THEN
           ! occ-occ block
           CALL cp_fm_set_all(fm_back, 0.0_dp)
           CALL cp_gemm('N','N',dimen,homo,homo,1.0_dp,&
                           mo_coeff_o,fm_mo,0.0_dp,fm_back,&
                           a_first_row=1,&
                           b_first_col=1,&
                           b_first_row=1,&
                           c_first_col=1,&
                           c_first_row=1)
           CALL cp_fm_set_all(fm_ao, 0.0_dp)
           CALL cp_gemm('N','T',dimen,dimen,homo,1.0_dp,&
                           fm_back,mo_coeff_o,0.0_dp,fm_ao,&
                           a_first_row=1,&
                           b_first_col=1,&
                           b_first_row=1,&
                           c_first_col=1,&
                           c_first_row=1)
        ELSE
           ! occ-occ block
           CALL cp_fm_set_all(fm_back, 0.0_dp)
           CALL cp_gemm('N','N',dimen,homo_beta,homo_beta,1.0_dp,&
                           mo_coeff_o_beta,fm_mo,0.0_dp,fm_back,&
                           a_first_row=1,&
                           b_first_col=1,&
                           b_first_row=1,&
                           c_first_col=1,&
                           c_first_row=1)
           CALL cp_fm_set_all(fm_ao, 0.0_dp)
           CALL cp_gemm('N','T',dimen,dimen,homo_beta,1.0_dp,&
                           fm_back,mo_coeff_o_beta,0.0_dp,fm_ao,&
                           a_first_row=1,&
                           b_first_col=1,&
                           b_first_row=1,&
                           c_first_col=1,&
                           c_first_row=1)
        ENDIF

      CASE(2)
        IF (.NOT. alpha_beta) THEN
           ! virt-virt block
           CALL cp_fm_set_all(fm_back, 0.0_dp)
           CALL cp_gemm('N','N',dimen,virtual,virtual,1.0_dp,&
                           mo_coeff_v,fm_mo,0.0_dp,fm_back,&
                           a_first_row=1,&
                           b_first_col=homo+1,&
                           b_first_row=homo+1,&
                           c_first_col=1,&
                           c_first_row=1)
           CALL cp_fm_set_all(fm_ao, 0.0_dp)
           CALL cp_gemm('N','T',dimen,dimen,virtual,1.0_dp,&
                           fm_back,mo_coeff_v,0.0_dp,fm_ao,&
                           a_first_row=1,&
                           b_first_col=1,&
                           b_first_row=1,&
                           c_first_col=1,&
                           c_first_row=1)
        ELSE
           CALL cp_fm_set_all(fm_back, 0.0_dp)
           CALL cp_gemm('N','N',dimen,virtual_beta,virtual_beta,1.0_dp,&
                           mo_coeff_v_beta,fm_mo,0.0_dp,fm_back,&
                           a_first_row=1,&
                           b_first_col=homo_beta+1,&
                           b_first_row=homo_beta+1,&
                           c_first_col=1,&
                           c_first_row=1)
           CALL cp_fm_set_all(fm_ao, 0.0_dp)
           CALL cp_gemm('N','T',dimen,dimen,virtual_beta,1.0_dp,&
                           fm_back,mo_coeff_v_beta,0.0_dp,fm_ao,&
                           a_first_row=1,&
                           b_first_col=1,&
                           b_first_row=1,&
                           c_first_col=1,&
                           c_first_row=1)
        ENDIF

      CASE(3)
        IF (.NOT. alpha_beta) THEN
           CALL cp_fm_set_all(fm_back, 0.0_dp)
           CALL cp_gemm('N','N',dimen,virtual,homo,1.0_dp,&
                           mo_coeff_o,fm_mo,0.0_dp,fm_back,&
                           a_first_row=1,&
                           b_first_col=1,&
                           b_first_row=1,&
                           c_first_col=1,&
                           c_first_row=1)
           CALL cp_fm_set_all(fm_ao, 0.0_dp)
           CALL cp_gemm('N','T',dimen,dimen,virtual,1.0_dp,&
                           fm_back,mo_coeff_v,0.0_dp,fm_ao,&
                           a_first_row=1,&
                           b_first_col=1,&
                           b_first_row=1,&
                           c_first_col=1,&
                           c_first_row=1)
           ! and symmetrize (here again multiply instead of transposing)
           CALL cp_fm_set_all(fm_back, 0.0_dp)
           CALL cp_gemm('N','T',dimen,homo,virtual,1.0_dp,&
                           mo_coeff_v,fm_mo,0.0_dp,fm_back,&
                           a_first_row=1,&
                           b_first_col=1,&
                           b_first_row=1,&
                           c_first_col=1,&
                           c_first_row=1)
           CALL cp_gemm('N','T',dimen,dimen,homo,0.5_dp,&
                           fm_back,mo_coeff_o,0.5_dp,fm_ao,&
                           a_first_row=1,&
                           b_first_col=1,&
                           b_first_row=1,&
                           c_first_col=1,&
                           c_first_row=1)
           ! scale for the orbital energy differences for the diagonal elements
           fm_mo_out%local_data(:,:)=fm_mo%local_data(:,:)
           CALL cp_fm_get_info(matrix=fm_mo_out,&
                               nrow_local=nrow_local,&
                               ncol_local=ncol_local,&
                               row_indices=row_indices,&
                               col_indices=col_indices)
           DO jjB=1, ncol_local
             j_global=col_indices(jjB)
             DO iiB=1, nrow_local
               i_global=row_indices(iiB)
               fm_mo_out%local_data(iiB,jjB)=fm_mo_out%local_data(iiB,jjB)*&
                                             (Eigenval(j_global+homo)-Eigenval(i_global))
             END DO
           END DO
        ELSE
           CALL cp_fm_set_all(fm_back, 0.0_dp)
           CALL cp_gemm('N','N',dimen,virtual_beta,homo_beta,1.0_dp,&
                           mo_coeff_o_beta,fm_mo,0.0_dp,fm_back,&
                           a_first_row=1,&
                           b_first_col=1,&
                           b_first_row=1,&
                           c_first_col=1,&
                           c_first_row=1)
           CALL cp_fm_set_all(fm_ao, 0.0_dp)
           CALL cp_gemm('N','T',dimen,dimen,virtual_beta,1.0_dp,&
                           fm_back,mo_coeff_v_beta,0.0_dp,fm_ao,&
                           a_first_row=1,&
                           b_first_col=1,&
                           b_first_row=1,&
                           c_first_col=1,&
                           c_first_row=1)
           ! and symmetrize (here again multiply instead of transposing)
           CALL cp_fm_set_all(fm_back, 0.0_dp)
           CALL cp_gemm('N','T',dimen,homo_beta,virtual_beta,1.0_dp,&
                           mo_coeff_v_beta,fm_mo,0.0_dp,fm_back,&
                           a_first_row=1,&
                           b_first_col=1,&
                           b_first_row=1,&
                           c_first_col=1,&
                           c_first_row=1)
           CALL cp_gemm('N','T',dimen,dimen,homo_beta,0.5_dp,&
                           fm_back,mo_coeff_o_beta,0.5_dp,fm_ao,&
                           a_first_row=1,&
                           b_first_col=1,&
                           b_first_row=1,&
                           c_first_col=1,&
                           c_first_row=1)
        ENDIF

      CASE(4)
        ! all-all block
        IF (.NOT. alpha_beta) THEN
           CALL cp_fm_set_all(fm_back, 0.0_dp)
           CALL cp_gemm('N','N',dimen,dimen,dimen,1.0_dp,&
                           mo_coeff,fm_mo,0.0_dp,fm_back,&
                           a_first_row=1,&
                           b_first_col=1,&
                           b_first_row=1,&
                           c_first_col=1,&
                           c_first_row=1)
           CALL cp_fm_set_all(fm_ao, 0.0_dp)
           CALL cp_gemm('N','T',dimen,dimen,dimen,1.0_dp,&
                           fm_back,mo_coeff,0.0_dp,fm_ao,&
                           a_first_row=1,&
                           b_first_col=1,&
                           b_first_row=1,&
                           c_first_col=1,&
                           c_first_row=1)
        ELSE
           CALL cp_fm_set_all(fm_back, 0.0_dp)
           CALL cp_gemm('N','N',dimen,dimen,dimen,1.0_dp,&
                           mo_coeff_beta,fm_mo,0.0_dp,fm_back,&
                           a_first_row=1,&
                           b_first_col=1,&
                           b_first_row=1,&
                           c_first_col=1,&
                           c_first_row=1)
           CALL cp_fm_set_all(fm_ao, 0.0_dp)
           CALL cp_gemm('N','T',dimen,dimen,dimen,1.0_dp,&
                           fm_back,mo_coeff_beta,0.0_dp,fm_ao,&
                           a_first_row=1,&
                           b_first_col=1,&
                           b_first_row=1,&
                           c_first_col=1,&
                           c_first_row=1)
        ENDIF

      CASE DEFAULT
        ! nothing
    END SELECT

    ! copy fm into DBCSR
    CALL cp_dbcsr_set(P_mu_nu%matrix,0.0_dp)
    CALL copy_fm_to_dbcsr(fm_ao, P_mu_nu%matrix, keep_sparsity=.TRUE.)

    ! calculate associated density 
    CALL calculate_rho_elec(matrix_p=P_mu_nu%matrix,&
                            rho=rho_r,&
                            rho_gspace=rho_g,&
                            total_rho=total_rho,&
                            ks_env=ks_env)
    ! and calculate potential
    CALL pw_poisson_solve(poisson_env, rho_g%pw, pair_energy, pot_g%pw)
    CALL pw_transfer(pot_g%pw, rho_r%pw)
    CALL pw_scale(rho_r%pw,rho_r%pw%pw_grid%dvol)
    ! integrate the potential
    CALL cp_dbcsr_set(mat_mu_nu(1)%matrix,0.0_dp)
    CALL integrate_v_rspace(rho_r,hmat=mat_mu_nu(1), &
                            qs_env=qs_env,calculate_forces=.FALSE.,compute_tau=.FALSE.,gapw=.FALSE.)

    ! update with the exchange like contributions
    ! copy mat_mu_nu into rho_ao work

    ! Only for alpha-alpha and beta-beta
    IF (.NOT. alpha_beta) THEN
       CALL qs_rho_get(rho_work, rho_ao=rho_work_ao)
       CALL cp_dbcsr_set(rho_work_ao(1)%matrix,0.0_dp)
       CALL cp_dbcsr_copy(rho_work_ao(1)%matrix,P_mu_nu%matrix)
       ! save old EX energy
       ex_energy=energy%ex
       DO irep=1, n_rep_hf
         ns = SIZE(rho_work_ao)
         rho_ao_2d(1:ns,1:1) => rho_work_ao(1:ns)
         ns = SIZE(mat_mu_nu)
         mat_2d(1:ns,1:1) => mat_mu_nu(1:ns)
         CALL integrate_four_center(qs_env, mat_2d, energy, rho_ao_2d, hfx_sections,&
                                    para_env, my_recalc_hfx_integrals, irep, .TRUE.,&
                                    ispin=1)
       END DO
       ! restore original EX energy
       energy%ex=ex_energy
    ENDIF

    ! scale by a factor 4.0 (closed shell) or 2.0 (open-shell)
    CALL cp_dbcsr_scale(mat_mu_nu(1)%matrix,factor)
 
    ! copy back to fm
    CALL cp_fm_set_all(fm_ao, 0.0_dp)
    CALL copy_dbcsr_to_fm(matrix=mat_mu_nu(1)%matrix, fm=fm_ao)
    CALL cp_fm_set_all(fm_back, 0.0_dp)
    CALL cp_fm_upper_to_full(fm_ao, fm_back)

    ! transform to MO basis, here we always sum the result into the input matrix
    SELECT CASE(transf_type_out)
    ! In alpha-beta case the density is contracted with the orbitals of different
    ! spin. Thus, there's branching at forward trasform and no branching here,
    ! at back transform.
      CASE(1)
        ! occ-virt block
        CALL cp_fm_set_all(fm_back, 0.0_dp)
        CALL cp_gemm('T','N',homo,dimen,dimen,1.0_dp,&
                        mo_coeff_o,fm_ao,0.0_dp,fm_back,&
                        a_first_row=1,&
                        b_first_col=1,&
                        b_first_row=1,&
                        c_first_col=1,&
                        c_first_row=1)
        CALL cp_gemm('N','N',homo,virtual,dimen,out_alpha,&
                        fm_back,mo_coeff_v,1.0_dp,fm_mo_out,&
                        a_first_row=1,&
                        b_first_col=1,&
                        b_first_row=1,&
                        c_first_col=1,&
                        c_first_row=1)

      CASE(2)
        ! occ-occ block
        CALL cp_fm_set_all(fm_back, 0.0_dp)
        CALL cp_gemm('T','N',homo,dimen,dimen,1.0_dp,&
                        mo_coeff_o,fm_ao,0.0_dp,fm_back,&
                        a_first_row=1,&
                        b_first_col=1,&
                        b_first_row=1,&
                        c_first_col=1,&
                        c_first_row=1)
        CALL cp_gemm('N','N',homo,homo,dimen,out_alpha,&
                        fm_back,mo_coeff_o,1.0_dp,fm_mo_out,&
                        a_first_row=1,&
                        b_first_col=1,&
                        b_first_row=1,&
                        c_first_col=1,&
                        c_first_row=1)       

      CASE DEFAULT
        ! nothing
    END SELECT


    CALL timestop(handle)

  END SUBROUTINE cphf_like_update

! *****************************************************************************
!> \brief Low level subroutine for the iterative solution of a large 
!>        system of linear equation
!> \param qs_env ...
!> \param mp2_env ...
!> \param para_env ...
!> \param homo ...
!> \param virtual ...
!> \param dimen ...
!> \param unit_nr ...
!> \param mo_coeff ...
!> \param mo_coeff_o ...
!> \param mo_coeff_v ...
!> \param Eigenval ...
!> \param blacs_env ...
!> \param hfx_sections ...
!> \param energy ...
!> \param n_rep_hf ...
!> \param poisson_env ...
!> \param rho_work ...
!> \param pot_g ...
!> \param rho_g ...
!> \param rho_r ...
!> \param mat_mu_nu ...
!> \param P_mu_nu ...
!> \param L_jb ...
!> \param fm_G_mu_nu ...
!> \param fm_back ...
!> \param P_ia ...
!> \param homo_beta ...
!> \param Eigenval_beta ...
!> \param P_ia_beta ...
!> \param mo_coeff_beta ...
!> \param mo_coeff_o_beta ...
!> \param mo_coeff_v_beta ...
!> \param L_jb_beta ...
!> \author Mauro Del Ben, Vladimir Rybkin
! *****************************************************************************
  SUBROUTINE solve_z_vector_eq_low(qs_env,mp2_env,para_env,homo,virtual,dimen,unit_nr,&
                                   mo_coeff,mo_coeff_o,mo_coeff_v,Eigenval,blacs_env,&
                                   hfx_sections,energy,n_rep_hf,poisson_env,&
                                   rho_work,pot_g,rho_g,rho_r,mat_mu_nu,P_mu_nu,&
                                   L_jb,fm_G_mu_nu,fm_back,P_ia,homo_beta,Eigenval_beta,&
                                   P_ia_beta, mo_coeff_beta, mo_coeff_o_beta, mo_coeff_v_beta,&
                                   L_jb_beta)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(mp2_type), POINTER                  :: mp2_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER                                  :: homo, virtual, dimen, unit_nr
    TYPE(cp_fm_type), POINTER                :: mo_coeff, mo_coeff_o, &
                                                mo_coeff_v
    REAL(KIND=dp), DIMENSION(:)              :: Eigenval
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    TYPE(section_vals_type), POINTER         :: hfx_sections
    TYPE(qs_energy_type), POINTER            :: energy
    INTEGER                                  :: n_rep_hf
    TYPE(pw_poisson_type), POINTER           :: poisson_env
    TYPE(qs_rho_type), POINTER               :: rho_work
    TYPE(pw_p_type)                          :: pot_g, rho_g, rho_r
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: mat_mu_nu
    TYPE(cp_dbcsr_p_type)                    :: P_mu_nu
    TYPE(cp_fm_type), POINTER                :: L_jb, fm_G_mu_nu, fm_back, &
                                                P_ia
    INTEGER, OPTIONAL                        :: homo_beta
    REAL(KIND=dp), DIMENSION(:), OPTIONAL    :: Eigenval_beta
    TYPE(cp_fm_type), OPTIONAL, POINTER      :: P_ia_beta, mo_coeff_beta, &
                                                mo_coeff_o_beta, &
                                                mo_coeff_v_beta, L_jb_beta

    CHARACTER(LEN=*), PARAMETER :: routineN = 'solve_z_vector_eq_low', &
      routineP = moduleN//':'//routineN

    INTEGER :: cycle_counter, handle, i_global, iiB, iiter, j_global, jjB, &
      max_num_iter, ncol_local, nrow_local, transf_type_in, transf_type_out, &
      virtual_beta
    INTEGER, DIMENSION(:), POINTER           :: col_indices, row_indices
    LOGICAL                                  :: alpha_beta, converged, failure
    REAL(KIND=dp) :: alpha, beta, conv, eps_conv, factor, norm_b, &
      norm_b_beta, norms(3), norms_beta(3), out_alpha, rkrk, rkrk_beta, t1, t2
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: proj_bi_xj, temp_vals, &
                                                x_norm, xi_b
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: A_small, b_small, xi_Axi
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: Ax, xn
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: Ap, Ap_b, b_i, pk, pk_b, &
                                                precond, residual, rk, rk_b, &
                                                xk, xk_b

    CALL timeset(routineN,handle)
    failure=.FALSE.

    alpha_beta = .FALSE.
    factor = 4.0_dp
    IF (PRESENT(homo_beta) .AND. PRESENT(Eigenval_beta)) alpha_beta = .TRUE.
    IF (alpha_beta) THEN
       virtual_beta = dimen - homo_beta
       factor = 2.0_dp
    ENDIF

    max_num_iter=mp2_env%ri_grad%cphf_max_num_iter
    eps_conv=mp2_env%ri_grad%cphf_eps_conv

    IF (unit_nr>0) THEN
      WRITE(unit_nr,*)
      WRITE(unit_nr,'(T3,A)')           'MP2_CPHF| Iterative solution of Z-Vector equations'
      WRITE(unit_nr,'(T3,A,T45,ES8.1)') 'MP2_CPHF| Convergence threshold:', eps_conv
      WRITE(unit_nr,'(T3,A,T45,I8)')    'MP2_CPHF| Maximum number of iterations: ', max_num_iter
      WRITE(unit_nr,'(T4,A)') REPEAT("-",40)
      WRITE(unit_nr,'(T4,A,T15,A,T33,A)') 'Step','Time','Convergence'
      WRITE(unit_nr,'(T4,A)') REPEAT("-",40)
    END IF

    ! set the transformation type (equal for all methods all updates)
    transf_type_in=3
    transf_type_out=1
    out_alpha=1.0_dp

    ! set convergece flag
    converged=.FALSE.

    IF((.FALSE.) .OR. (alpha_beta)) THEN
      ! CG algorithm
      ! Used for open-shell systems
      ! create some work array
      NULLIFY(xk, pk, rk, Ap, fm_struct_tmp)
      CALL cp_fm_struct_create(fm_struct_tmp,para_env=para_env,context=blacs_env, &
                               nrow_global=homo,ncol_global=virtual)
      CALL cp_fm_create(xk, fm_struct_tmp, name="xk")
      CALL cp_fm_create(pk, fm_struct_tmp, name="pk")
      CALL cp_fm_create(rk, fm_struct_tmp, name="rk")
      CALL cp_fm_create(Ap, fm_struct_tmp, name="Ap")
      CALL cp_fm_struct_release(fm_struct_tmp)
      CALL cp_fm_set_all(xk, 0.0_dp)
      CALL cp_fm_set_all(pk, 0.0_dp)
      CALL cp_fm_set_all(rk, 0.0_dp)
      CALL cp_fm_set_all(Ap, 0.0_dp)

      ! copy -L_jb into pk and rk
      pk%local_data(:,:)=-L_jb%local_data(:,:)
      rk%local_data(:,:)=-L_jb%local_data(:,:)
      norm_b=0.0_dp
      norm_b=SUM(L_jb%local_data(:,:)*L_jb%local_data(:,:))
      CALL mp_sum(norm_b,para_env%group)

      IF (alpha_beta) THEN
         NULLIFY(xk_b, pk_b, rk_b, Ap_b, fm_struct_tmp)
         CALL cp_fm_struct_create(fm_struct_tmp,para_env=para_env,context=blacs_env, &
                            nrow_global=homo_beta,ncol_global=virtual_beta)
         CALL cp_fm_create(xk_b, fm_struct_tmp, name="xk")
         CALL cp_fm_create(pk_b, fm_struct_tmp, name="pk")
         CALL cp_fm_create(rk_b, fm_struct_tmp, name="rk")
         CALL cp_fm_create(Ap_b, fm_struct_tmp, name="Ap")
         CALL cp_fm_struct_release(fm_struct_tmp)
         CALL cp_fm_set_all(xk_b, 0.0_dp)
         CALL cp_fm_set_all(pk_b, 0.0_dp)
         CALL cp_fm_set_all(rk_b, 0.0_dp)
         CALL cp_fm_set_all(Ap_b, 0.0_dp)
 
         ! copy -L_jb_beta into pk_b and rk_b
         pk_b%local_data(:,:)=-L_jb_beta%local_data(:,:)
         rk_b%local_data(:,:)=-L_jb_beta%local_data(:,:)
         norm_b_beta=0.0_dp
         norm_b_beta=SUM(L_jb_beta%local_data(:,:)*L_jb_beta%local_data(:,:))
         CALL mp_sum(norm_b_beta,para_env%group)
         norm_b=norm_b + norm_b_beta
      ENDIF
      ! Set the norm
      norm_b=SQRT(norm_b)

      cycle_counter=0
      DO iiter=1, max_num_iter
        cycle_counter=cycle_counter+1
        t1 = m_walltime()

        ! calculate matrix-vector product
        CALL cp_fm_set_all(Ap, 0.0_dp)
        ! Full for closed shell. Alpha-alpha part of alpha for open shell.
        CALL cphf_like_update(qs_env,para_env,homo,virtual,dimen,&
                              mo_coeff,mo_coeff_o,mo_coeff_v,Eigenval,&
                              hfx_sections,energy,n_rep_hf,poisson_env,&
                              rho_work,pot_g,rho_g,rho_r,mat_mu_nu,P_mu_nu,&
                              pk,fm_G_mu_nu,fm_back,transf_type_in,out_alpha,&
                              Ap,transf_type_out,factor=factor)
        IF (alpha_beta) THEN
           ! Alpha-beta part of alpha.
           CALL cphf_like_update(qs_env,para_env,homo,virtual,dimen,&
                             mo_coeff,mo_coeff_o,mo_coeff_v,Eigenval,&
                             hfx_sections,energy,n_rep_hf,poisson_env,&
                             rho_work,pot_g,rho_g,rho_r,mat_mu_nu,P_mu_nu,&
                             pk_b,fm_G_mu_nu,fm_back,transf_type_in,out_alpha,&
                             Ap,transf_type_out,factor=factor, &
                             mo_coeff_beta=mo_coeff_beta, &
                             mo_coeff_o_beta=mo_coeff_o_beta, &
                             mo_coeff_v_beta=mo_coeff_v_beta, &
                             homo_beta=homo_beta, virtual_beta=virtual_beta)
           ! Beta-beta part (Coulomb and XC) of beta.                     
           CALL cp_fm_set_all(Ap_b, 0.0_dp)
           CALL cphf_like_update(qs_env,para_env,homo_beta,virtual_beta,dimen, &
                              mo_coeff_beta,mo_coeff_o_beta,mo_coeff_v_beta, &
                             Eigenval_beta,hfx_sections,energy,n_rep_hf,&
                             poisson_env,rho_work,pot_g,rho_g,rho_r,mat_mu_nu,&
                             P_mu_nu,Pk_b,fm_G_mu_nu,fm_back,transf_type_in,&
                             out_alpha,Ap_b,transf_type_out,factor=factor)
           ! Beta-alpha part of beta.
            CALL cphf_like_update(qs_env,para_env,homo_beta,virtual_beta,dimen,&
                             mo_coeff_beta,mo_coeff_o_beta,mo_coeff_v_beta, &
                             Eigenval_beta,&
                             hfx_sections,energy,n_rep_hf,poisson_env,&
                             rho_work,pot_g,rho_g,rho_r,mat_mu_nu,P_mu_nu,&
                             pk,fm_G_mu_nu,fm_back,transf_type_in,out_alpha,&
                             Ap_b,transf_type_out,factor=factor,&
                             mo_coeff_beta=mo_coeff, &
                             mo_coeff_o_beta=mo_coeff_o, &
                             mo_coeff_v_beta=mo_coeff_v, &
                             homo_beta=homo, virtual_beta=virtual)
        ENDIF

        norms=0.0_dp
        norms(1)=SUM(rk%local_data(:,:)*rk%local_data(:,:))
        norms(2)=SUM(rk%local_data(:,:)*pk%local_data(:,:))
        norms(3)=SUM(pk%local_data(:,:)*Ap%local_data(:,:))
        CALL mp_sum(norms,para_env%group)
        IF (alpha_beta) THEN
           norms_beta=0.0_dp
           norms_beta(1)=SUM(rk_b%local_data(:,:)*rk_b%local_data(:,:))
           norms_beta(2)=SUM(rk_b%local_data(:,:)*pk_b%local_data(:,:))
           norms_beta(3)=SUM(pk_b%local_data(:,:)*Ap_b%local_data(:,:))
           CALL mp_sum(norms_beta,para_env%group)
           norms = norms + norms_beta
        ENDIF
        alpha=norms(1)/norms(3)

        xk%local_data(:,:)=xk%local_data(:,:)+alpha*pk%local_data(:,:)
        rk%local_data(:,:)=rk%local_data(:,:)-alpha*Ap%local_data(:,:)
        IF (alpha_beta) THEN
           xk_b%local_data(:,:)=xk_b%local_data(:,:)+alpha*pk_b%local_data(:,:)
           rk_b%local_data(:,:)=rk_b%local_data(:,:)-alpha*Ap_b%local_data(:,:)
        ENDIF

        rkrk=0.0_dp
        rkrk=SUM(rk%local_data(:,:)*rk%local_data(:,:))
        CALL mp_sum(rkrk,para_env%group) 

        IF (alpha_beta) THEN
           rkrk_beta=0.0_dp
           rkrk_beta = SUM(rk_b%local_data(:,:)*rk_b%local_data(:,:))
           CALL mp_sum(rkrk_beta,para_env%group) 
           rkrk = rkrk + rkrk_beta
        ENDIF
        beta=rkrk/norms(2)

        pk%local_data(:,:)=rk%local_data(:,:)+beta*pk%local_data(:,:)
        IF (alpha_beta) THEN
           pk_b%local_data(:,:)=rk_b%local_data(:,:)+beta*pk_b%local_data(:,:)
        ENDIF

        conv=SQRT(rkrk)/norm_b
 
        t2 = m_walltime()

        IF (unit_nr>0) THEN
          WRITE(unit_nr,'(T3,I5,T13,F6.1,11X,F14.8)') iiter, t2-t1, conv
        END IF             

        IF(conv<=eps_conv) THEN
          converged=.TRUE.
          EXIT
        END IF

      END DO

      P_ia%local_data(:,:)=xk%local_data(:,:)
      IF (alpha_beta) THEN
         P_ia_beta%local_data(:,:)=xk_b%local_data(:,:)
      ENDIF

      CALL cp_fm_release(xk)
      CALL cp_fm_release(pk)
      CALL cp_fm_release(rk)
      CALL cp_fm_release(Ap)
      IF (alpha_beta) THEN
         CALL cp_fm_release(xk_b)
         CALL cp_fm_release(pk_b)
         CALL cp_fm_release(rk_b)
         CALL cp_fm_release(Ap_b)
      ENDIF

    ELSE
      ! Pople method
      ! change sign to L_jb
      L_jb%local_data(:,:)=-L_jb%local_data(:,:)

      ! allocate stuff
      ALLOCATE(xn(1:max_num_iter))
      ALLOCATE(Ax(1:max_num_iter))

      ! create fm structure
      NULLIFY(fm_struct_tmp)
      CALL cp_fm_struct_create(fm_struct_tmp,para_env=para_env,context=blacs_env, &
                               nrow_global=homo,ncol_global=virtual)


      ! create preconditioner (for now only orbital energy differences)
      NULLIFY(precond)
      CALL cp_fm_create(precond, fm_struct_tmp, name="precond")
      CALL cp_fm_set_all(precond,1.0_dp)
      CALL cp_fm_get_info(matrix=precond,&
                          nrow_local=nrow_local,&
                          ncol_local=ncol_local,&
                          row_indices=row_indices,&
                          col_indices=col_indices)
      DO jjB=1, ncol_local
        j_global=col_indices(jjB)
        DO iiB=1, nrow_local
          i_global=row_indices(iiB)
          precond%local_data(iiB,jjB)=precond%local_data(iiB,jjB)/&
                                      (Eigenval(j_global+homo)-Eigenval(i_global))
        END DO
      END DO

      ! x(iiter) vector
      NULLIFY(b_i)
      CALL cp_fm_create(b_i, fm_struct_tmp, name="b_i")
      CALL cp_fm_set_all(b_i, 0.0_dp)
      b_i%local_data(:,:)=precond%local_data(:,:)*L_jb%local_data(:,:)

      ! create the residual vector (r), we check convergence on the norm of 
      ! this vector r=(Ax-b) 
      NULLIFY(residual)
      CALL cp_fm_create(residual, fm_struct_tmp, name="residual")
      CALL cp_fm_set_all(residual, 0.0_dp)

      ! allocate array containing the various scalar products
      ALLOCATE(x_norm(1:max_num_iter))
      ALLOCATE(xi_b(1:max_num_iter))
      ALLOCATE(xi_Axi(1:max_num_iter,0:max_num_iter))
      x_norm=0.0_dp
      xi_b=0.0_dp
      xi_Axi=0.0_dp

      cycle_counter=0
      DO iiter=1, max_num_iter
        cycle_counter=cycle_counter+1
 
        t1 = m_walltime()

        ! create and update x_i (orthogonalization with previous vectors)
        NULLIFY(xn(iiter)%matrix)
        CALL cp_fm_create(xn(iiter)%matrix, fm_struct_tmp, name="xi")
        CALL cp_fm_set_all(xn(iiter)%matrix, 0.0_dp)

        ! first compute the projection of the actual b_i into all previous x_i
        ! already scaled with the norm of each x_i
        ALLOCATE(proj_bi_xj(iiter-1))
        DO iiB=1, iiter-1
          proj_bi_xj(iiB)=0.0_dp
          proj_bi_xj(iiB)=accurate_sum(b_i%local_data(1:nrow_local,1:ncol_local)*&
                                       xn(iiB)%matrix%local_data(1:nrow_local,1:ncol_local))
          proj_bi_xj(iiB)=proj_bi_xj(iiB)/x_norm(iiB)
        END DO
        CALL mp_sum(proj_bi_xj,para_env%group)

        ! update actual x_i
        xn(iiter)%matrix%local_data(:,:)=b_i%local_data(:,:)
        DO iiB=1, iiter-1
          xn(iiter)%matrix%local_data(:,:)=xn(iiter)%matrix%local_data(:,:)-&
                                           xn(iiB)%matrix%local_data(:,:)*proj_bi_xj(iiB)
        END DO
        DEALLOCATE(proj_bi_xj)
         
        ! create Ax(iiter) that will store the matrix vector product for this cycle
        NULLIFY(Ax(iiter)%matrix)
        CALL cp_fm_create(Ax(iiter)%matrix, fm_struct_tmp, name="Ai")
        CALL cp_fm_set_all(Ax(iiter)%matrix, 0.0_dp)
        ! perform the matrix-vector product (CPHF like update) 
        CALL cphf_like_update(qs_env,para_env,homo,virtual,dimen,&
                              mo_coeff,mo_coeff_o,mo_coeff_v,Eigenval,&
                              hfx_sections,energy,n_rep_hf,poisson_env,&
                              rho_work,pot_g,rho_g,rho_r,mat_mu_nu,P_mu_nu,&
                              xn(iiter)%matrix,fm_G_mu_nu,fm_back,transf_type_in,out_alpha,&
                              Ax(iiter)%matrix,transf_type_out,factor=factor)

        ! in order to reduce the number of calls to mp_sum here we
        ! cluster all necessary scalar products into a sigle vector
        ! temp_vals contains:
        ! 1:iiter -> <Ax_i|x_j>
        ! iiter+1 -> <x_i|b>
        ! iiter+2 -> <x_i|x_i>
        ALLOCATE(temp_vals(iiter+2))
        temp_vals=0.0_dp
        ! <Ax_i|x_j>
        DO iiB=1, iiter
          temp_vals(iiB)=accurate_sum(Ax(iiter)%matrix%local_data(1:nrow_local,1:ncol_local)*&
                                      xn(iiB)%matrix%local_data(1:nrow_local,1:ncol_local))
        END DO        
        ! <x_i|b>
        temp_vals(iiter+1)=accurate_sum(xn(iiter)%matrix%local_data(1:nrow_local,1:ncol_local)*&
                                                    L_jb%local_data(1:nrow_local,1:ncol_local))
        ! norm
        temp_vals(iiter+2)=accurate_sum(xn(iiter)%matrix%local_data(1:nrow_local,1:ncol_local)*&
                                        xn(iiter)%matrix%local_data(1:nrow_local,1:ncol_local))
        CALL mp_sum(temp_vals,para_env%group)
        ! update <Ax_i|x_j>,  <x_i|b> and norm <x_i|x_i>
        xi_Axi(iiter,1:iiter)=temp_vals(1:iiter)
        xi_Axi(1:iiter,iiter)=temp_vals(1:iiter) 
        xi_b(iiter)   = temp_vals(iiter+1)
        x_norm(iiter) = temp_vals(iiter+2)
        ! deallocate temp_vals
        DEALLOCATE(temp_vals)

        ! solve reduced system 
        IF(ALLOCATED(A_small)) DEALLOCATE(A_small)
        IF(ALLOCATED(b_small)) DEALLOCATE(b_small)
        ALLOCATE(A_small(iiter,iiter))
        ALLOCATE(b_small(iiter,1))
        A_small(1:iiter,1:iiter)=xi_Axi(1:iiter,1:iiter)
        b_small(1:iiter,1)=xi_b(1:iiter)

        CALL solve_system(matrix=A_small, mysize=iiter, eigenvectors=b_small)

        ! check for convergence
        CALL cp_fm_set_all(residual, 0.0_dp)
        DO iiB=1, iiter
          residual%local_data(1:nrow_local,1:ncol_local)=residual%local_data(1:nrow_local,1:ncol_local)+&
                                    b_small(iiB,1)*Ax(iiB)%matrix%local_data(1:nrow_local,1:ncol_local)
        END DO
        residual%local_data(1:nrow_local,1:ncol_local)=residual%local_data(1:nrow_local,1:ncol_local)-&
                                                           L_jb%local_data(1:nrow_local,1:ncol_local)
        conv=0.0_dp
        conv=accurate_sum(residual%local_data(1:nrow_local,1:ncol_local)*&
                          residual%local_data(1:nrow_local,1:ncol_local))
        CALL mp_sum(conv,para_env%group)
        conv=SQRT(conv)

        t2 = m_walltime()

        IF (unit_nr>0) THEN
          WRITE(unit_nr,'(T3,I5,T13,F6.1,11X,F14.8)') iiter, t2-t1, conv
        END IF

        IF(conv<=eps_conv) THEN
          converged=.TRUE.
          EXIT
        END IF

        ! update b_i for the next round
        b_i%local_data(:,:)=precond%local_data(:,:)*Ax(iiter)%matrix%local_data(:,:)

      END DO
 
      ! store solution into P_ia
      DO iiter=1, cycle_counter
        P_ia%local_data(1:nrow_local,1:ncol_local)=P_ia%local_data(1:nrow_local,1:ncol_local)+&
                        b_small(iiter,1)*xn(iiter)%matrix%local_data(1:nrow_local,1:ncol_local)
      END DO

      DEALLOCATE(x_norm)
      DEALLOCATE(xi_b)
      DEALLOCATE(xi_Axi)

      CALL cp_fm_release(precond)
      CALL cp_fm_release(b_i)
      CALL cp_fm_release(residual)
      CALL cp_fm_struct_release(fm_struct_tmp)

      ! release Ax, xn
      DO iiter=1, cycle_counter
        CALL cp_fm_release(Ax(iiter)%matrix)
        CALL cp_fm_release(xn(iiter)%matrix)
      END DO
      DEALLOCATE(xn)
      DEALLOCATE(Ax)

    END IF

    IF (unit_nr>0) THEN
      WRITE(unit_nr,'(T4,A)') REPEAT("-",40)
      IF(converged) THEN
        WRITE(unit_nr,'(T3,A,I5,A)')'Z-Vector equations converged in',cycle_counter,' steps'
      ELSE
        WRITE(unit_nr,'(T3,A,I5,A)')'Z-Vector equations NOT converged in',cycle_counter,' steps'
      END IF
    END IF

    CALL timestop(handle)

  END SUBROUTINE solve_z_vector_eq_low

! *****************************************************************************
!> \brief ...
!> \param mat ...
!> \param unitout ...
! *****************************************************************************
  SUBROUTINE write_array(mat,unitout)
    REAL(KIND=dp), DIMENSION(:, :)           :: mat
    INTEGER, OPTIONAL                        :: unitout

    INTEGER                                  :: iii, jjj

     WRITE(*,*)
     DO iii=1, SIZE(mat,1)
       WRITE(*,*) iii
       DO jjj=1, SIZE(mat,2), 20
         IF(PRESENT(unitout)) THEN
           WRITE(1000+unitout,'(1000F10.5)') mat(iii,jjj:MIN(SIZE(mat,2),jjj+19))
         ELSE
           WRITE(*,'(1000F10.5)') mat(iii,jjj:MIN(SIZE(mat,2),jjj+19))
         END IF
       END DO
       WRITE(*,*)
     END DO
     WRITE(*,*)
  END SUBROUTINE

END MODULE mp2_cphf
