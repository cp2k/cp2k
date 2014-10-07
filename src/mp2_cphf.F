!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Rountines to calculate CPHF like update and solve Z-vector equation 
!>        for MP2 gradients (only GPW)
!> \par History
!>      11.2013 created [Mauro Del Ben]
! *****************************************************************************
MODULE mp2_cphf
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE cell_types,                      ONLY: cell_type
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
! USE f77_blas
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
  USE input_constants,                 ONLY: hfx_do_eval_energy,&
                                             use_orb_basis_set
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
  USE particle_types,                  ONLY: particle_type
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
  USE qs_rho_methods,                  ONLY: qs_rho_rebuild
  USE qs_rho_types,                    ONLY: qs_rho_release,&
                                             qs_rho_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE virial_types,                    ONLY: virial_type
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
!> \param cell ...
!> \param particle_set ...
!> \param atomic_kind_set ...
!> \param mo_coeff ...
!> \param nmo ...
!> \param homo ...
!> \param Eigenval ...
!> \param unit_nr ...
!> \param error ...
!> \author Mauro Del Ben 
! *****************************************************************************
  SUBROUTINE solve_z_vector_eq(qs_env,mp2_env,para_env,dft_control,cell,particle_set,&
                               atomic_kind_set,mo_coeff,nmo,homo,Eigenval,unit_nr,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(mp2_type), POINTER                  :: mp2_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(cell_type), POINTER                 :: cell
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    INTEGER                                  :: nmo, homo
    REAL(KIND=dp), DIMENSION(:)              :: Eigenval
    INTEGER                                  :: unit_nr
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'solve_z_vector_eq', &
      routineP = moduleN//':'//routineN

    INTEGER :: alpha, beta, bin, dimen, handle, handle2, i, i_global, &
      i_thread, iiB, ikind, irep, ispin, j_global, jjB, my_bin_size, &
      n_rep_hf, n_threads, ncol_local, nrow_local, stat, transf_type_in, &
      transf_type_out, virtual
    INTEGER, DIMENSION(3)                    :: comp
    INTEGER, DIMENSION(:), POINTER           :: col_indices, row_indices
    LOGICAL                                  :: do_dynamic_load_balancing, &
                                                do_hfx, failure, &
                                                hfx_treat_lsd_in_core, &
                                                use_virial
    REAL(KIND=dp)                            :: e_hartree, out_alpha, &
                                                pair_energy, tot_rho_r
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: mat_deb
    REAL(KIND=dp), DIMENSION(3, 3)           :: h_stress
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    TYPE(cp_dbcsr_p_type)                    :: P_mu_nu
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: mat_mu_nu, matrix_ks, &
                                                matrix_p_mp2, matrix_s, &
                                                matrix_w_mp2
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: fm_back, fm_G_mu_nu, L_jb, &
                                                mo_coeff_o, mo_coeff_v, P_ia, &
                                                P_mo, W_mo
    TYPE(hfx_container_type), DIMENSION(:), &
      POINTER                                :: integral_containers
    TYPE(hfx_container_type), POINTER        :: maxval_container
    TYPE(hfx_type), POINTER                  :: actual_x_data
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type)                          :: dvg(3), pot_g, rho_g, rho_r, &
                                                temp_pw_g
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

!$  INTEGER :: omp_get_max_threads

    CALL timeset(routineN,handle)
    failure=.FALSE.

    ! start collecting stuff
    dimen=nmo
    virtual=dimen-homo
    NULLIFY(input,pw_env,matrix_s,blacs_env,rho,energy,force,virial,matrix_w_mp2,&
            matrix_p_mp2,matrix_ks)
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
                    error=error)

    ! check if we have to calculate the virial
    use_virial = virial%pv_availability.AND.(.NOT.virial%pv_numer)

    ! mp2 matrices
    NULLIFY(P_mo, W_mo, L_jb)
    P_mo => mp2_env%ri_grad%P_mo
    W_mo => mp2_env%ri_grad%W_mo
    L_jb => mp2_env%ri_grad%L_jb

    ! pw stuff
    NULLIFY(poisson_env,pw_pools,auxbas_pw_pool)
    CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool,&
                    pw_pools=pw_pools, poisson_env=poisson_env, error=error)

    ! get some of the grids ready
    NULLIFY(rho_r%pw,rho_g%pw,pot_g%pw)
    CALL pw_pool_create_pw(auxbas_pw_pool,rho_r%pw,&
                           use_data=REALDATA3D,&
                           in_space=REALSPACE,error=error)
    CALL pw_pool_create_pw(auxbas_pw_pool,rho_g%pw,&
                           use_data=COMPLEXDATA1D,&
                           in_space=RECIPROCALSPACE,error=error)
    CALL pw_pool_create_pw(auxbas_pw_pool,pot_g%pw,&
                           use_data=COMPLEXDATA1D,&
                           in_space=RECIPROCALSPACE,error=error)

    ! hfx section
    NULLIFY(hfx_sections)
    hfx_sections => section_vals_get_subs_vals(input,"DFT%XC%HF",error=error)
    CALL section_vals_get(hfx_sections,explicit=do_hfx,n_repetition=n_rep_hf,error=error)
    IF( do_hfx ) THEN
      CALL section_vals_val_get(hfx_sections, "TREAT_LSD_IN_CORE", l_val=hfx_treat_lsd_in_core,&
                                i_rep_section=1,error=error)
    END IF

    ! create work array
    NULLIFY(mat_mu_nu)
    CALL cp_dbcsr_allocate_matrix_set(mat_mu_nu, dft_control%nspins, error)
    DO ispin=1,dft_control%nspins
      ALLOCATE(mat_mu_nu(ispin)%matrix)
      CALL cp_dbcsr_init(mat_mu_nu(ispin)%matrix,error=error)
      CALL cp_dbcsr_create(matrix=mat_mu_nu(ispin)%matrix,&
           name="T_mu_nu",&
           dist=cp_dbcsr_distribution(matrix_s(1)%matrix), matrix_type=dbcsr_type_symmetric,&
           row_blk_size=cp_dbcsr_row_block_sizes(matrix_s(1)%matrix),&
           col_blk_size=cp_dbcsr_col_block_sizes(matrix_s(1)%matrix),&
           nblks=0, nze=0, error=error)
      CALL cp_dbcsr_alloc_block_from_nbl(mat_mu_nu(ispin)%matrix,qs_env%sab_orb,error=error)
      CALL cp_dbcsr_set(mat_mu_nu(ispin)%matrix,0.0_dp,error=error)
    END DO

    ALLOCATE(P_mu_nu%matrix)
    CALL cp_dbcsr_init(P_mu_nu%matrix,error=error)
    ! CALL cp_dbcsr_create(P_mu_nu%matrix,template=mat_mu_nu(1)%matrix,error=error)
    ! CALL cp_dbcsr_copy(P_mu_nu%matrix,mat_mu_nu(1)%matrix,name="P_mu_nu",error=error)
    CALL cp_dbcsr_copy(P_mu_nu%matrix,rho%rho_ao(1)%matrix,name="P_mu_nu",error=error)
    CALL cp_dbcsr_set(P_mu_nu%matrix,0.0_dp,error=error)

    NULLIFY(fm_G_mu_nu, fm_struct_tmp)
    CALL cp_fm_struct_create(fm_struct_tmp,para_env=para_env,context=blacs_env, &
                             nrow_global=dimen,ncol_global=dimen,error=error)
    CALL cp_fm_create(fm_G_mu_nu, fm_struct_tmp, name="G_mu_nu",error=error)
    CALL cp_fm_create(fm_back, fm_struct_tmp, name="fm_back",error=error)
    CALL cp_fm_struct_release(fm_struct_tmp,error=error)
    CALL cp_fm_set_all(fm_G_mu_nu, 0.0_dp,error=error)
    CALL cp_fm_set_all(fm_back, 0.0_dp,error=error)

    NULLIFY(mo_coeff_o, fm_struct_tmp)
    CALL cp_fm_struct_create(fm_struct_tmp,para_env=para_env,context=blacs_env, &
                             nrow_global=dimen,ncol_global=homo,error=error)
    CALL cp_fm_create(mo_coeff_o, fm_struct_tmp, name="mo_coeff_o",error=error)
    CALL cp_fm_struct_release(fm_struct_tmp,error=error)
    CALL cp_fm_set_all(mo_coeff_o, 0.0_dp,error=error)
    CALL cp_fm_to_fm_submat(msource=mo_coeff, mtarget=mo_coeff_o, &
                            nrow=dimen, ncol=homo, &
                            s_firstrow=1, s_firstcol=1, &
                            t_firstrow=1, t_firstcol=1, error=error)
    
    NULLIFY(mo_coeff_v, fm_struct_tmp)
    CALL cp_fm_struct_create(fm_struct_tmp,para_env=para_env,context=blacs_env, &
                             nrow_global=dimen,ncol_global=virtual,error=error)
    CALL cp_fm_create(mo_coeff_v, fm_struct_tmp, name="mo_coeff_v",error=error)
    CALL cp_fm_struct_release(fm_struct_tmp,error=error)
    CALL cp_fm_set_all(mo_coeff_v, 0.0_dp,error=error)
    CALL cp_fm_to_fm_submat(msource=mo_coeff, mtarget=mo_coeff_v, &
                            nrow=dimen, ncol=virtual, & 
                            s_firstrow=1, s_firstcol=homo+1, &
                            t_firstrow=1, t_firstcol=1, error=error)

    ! create a working rho environment
    NULLIFY(rho_work)
    CALL qs_rho_rebuild(rho=rho_work, qs_env=qs_env, rebuild_ao=.TRUE., rebuild_grids=.FALSE., error=error)

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
           !  CALL dealloc_containers(actual_x_data, hfx_do_eval_energy, error)
            CALL alloc_containers(actual_x_data, my_bin_size, hfx_do_eval_energy, error)
           
            DO bin=1, my_bin_size
              maxval_container => actual_x_data%maxval_container(bin)
              integral_containers => actual_x_data%integral_containers(:,bin)
              CALL hfx_init_container(maxval_container, actual_x_data%memory_parameter%actual_memory_usage, .FALSE., error)
              DO i=1,64
                CALL hfx_init_container(integral_containers(i), actual_x_data%memory_parameter%actual_memory_usage, .FALSE., error)
              END DO
            END DO
          END IF
        END DO
      END DO
      CALL timestop(handle2)
    END IF

    ! update lagrangian with the CPHF like update, occ-occ block, first call (recompute hfx integrals if needed)
    transf_type_in=1
    transf_type_out=1
    out_alpha=0.5_dp
    CALL cphf_like_update(qs_env,mp2_env,para_env,homo,virtual,dimen,unit_nr,&
                          mo_coeff,mo_coeff_o,mo_coeff_v,Eigenval,dft_control,&
                          hfx_sections,energy,n_rep_hf,pw_env,poisson_env,&
                          rho_work,pot_g,rho_g,rho_r,mat_mu_nu,P_mu_nu,&
                          P_mo,fm_G_mu_nu,fm_back,transf_type_in,out_alpha,&
                          L_jb,transf_type_out,error,&
                          recalc_hfx_integrals=mp2_env%ri_mp2%free_hfx_buffer)

    ! update lagrangian with the CPHF like update, virt-virt block
    transf_type_in=2
    transf_type_out=1
    out_alpha=0.5_dp
    CALL cphf_like_update(qs_env,mp2_env,para_env,homo,virtual,dimen,unit_nr,&
                          mo_coeff,mo_coeff_o,mo_coeff_v,Eigenval,dft_control,&
                          hfx_sections,energy,n_rep_hf,pw_env,poisson_env,&
                          rho_work,pot_g,rho_g,rho_r,mat_mu_nu,P_mu_nu,&
                          P_mo,fm_G_mu_nu,fm_back,transf_type_in,out_alpha,&
                          L_jb,transf_type_out,error)
    ! at this point Lagrnagian is completed ready to solve the Z-vector equations
    ! P_ia will contain the solution of these equations
    NULLIFY(P_ia, fm_struct_tmp)
    CALL cp_fm_struct_create(fm_struct_tmp,para_env=para_env,context=blacs_env, &
                             nrow_global=homo,ncol_global=virtual,error=error)
    CALL cp_fm_create(P_ia, fm_struct_tmp, name="P_ia",error=error)
    CALL cp_fm_struct_release(fm_struct_tmp,error=error)
    CALL cp_fm_set_all(P_ia, 0.0_dp,error=error)

    CALL solve_z_vector_eq_low(qs_env,mp2_env,para_env,homo,virtual,dimen,unit_nr,&
                               mo_coeff,mo_coeff_o,mo_coeff_v,Eigenval,blacs_env,dft_control,&
                               hfx_sections,energy,n_rep_hf,pw_env,poisson_env,&
                               rho_work,pot_g,rho_g,rho_r,mat_mu_nu,P_mu_nu,&
                               L_jb,fm_G_mu_nu,fm_back,P_ia,error)
    
    ! release Lagrangian
    CALL cp_fm_release(L_jb,error=error)

    ! update the MP2-MO density matrix with the occ-virt block
    CALL cp_fm_to_fm_submat(msource=P_ia, mtarget=P_mo, &
                            nrow=homo, ncol=virtual, &
                            s_firstrow=1, s_firstcol=1, &
                            t_firstrow=1, t_firstcol=homo+1, error=error)
    CALL cp_fm_release(P_ia,error=error)
    ! transpose P_MO matrix (easy way to symmetrize)
    CALL cp_fm_set_all(fm_back, 0.0_dp,error=error)
    ! P_mo now is ready
    CALL cp_fm_upper_to_full(matrix=P_mo, work=fm_back, error=error)

    ! do the final update to MP2 energy weighted matrix W_MO
    CALL cp_fm_get_info(matrix=W_mo,&
                        nrow_local=nrow_local,&
                        ncol_local=ncol_local,&
                        row_indices=row_indices,&
                        col_indices=col_indices,&
                        error=error)
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

    ! complete the occ-occ block of W_mo with a CPHF like update
    transf_type_in=4
    transf_type_out=2
    out_alpha=-0.5_dp
    CALL cphf_like_update(qs_env,mp2_env,para_env,homo,virtual,dimen,unit_nr,&
                          mo_coeff,mo_coeff_o,mo_coeff_v,Eigenval,dft_control,&
                          hfx_sections,energy,n_rep_hf,pw_env,poisson_env,&
                          rho_work,pot_g,rho_g,rho_r,mat_mu_nu,P_mu_nu,&
                          P_mo,fm_G_mu_nu,fm_back,transf_type_in,out_alpha,&
                          W_mo,transf_type_out,error)

    ! release DBCSR stuff
    DO ispin=1, dft_control%nspins
      CALL cp_dbcsr_release(mat_mu_nu(ispin)%matrix,error=error)
      DEALLOCATE(mat_mu_nu(ispin)%matrix,STAT=stat)
    END DO
    DEALLOCATE(mat_mu_nu,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL cp_dbcsr_release(P_mu_nu%matrix,error=error)
    DEALLOCATE(P_mu_nu%matrix,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ! release fm stuff
    CALL cp_fm_release(fm_G_mu_nu,error=error)
    CALL cp_fm_release(mo_coeff_o,error=error)
    CALL cp_fm_release(mo_coeff_v,error=error)
    ! release rho stuff
    CALL qs_rho_release(rho_struct=rho_work,error=error)
 
    IF(.FALSE.) THEN
      ALLOCATE(mat_deb(dimen,dimen))
      mat_deb=0.0_dp
      CALL cp_fm_get_info(matrix=P_mo,&
                          nrow_local=nrow_local,&
                          ncol_local=ncol_local,&
                          row_indices=row_indices,&
                          col_indices=col_indices,&
                          error=error)
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
                          col_indices=col_indices,&
                          error=error)
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
    END IF

    ! backtransform into AO basis, since P_mo and W_mo
    ! are symmetric (in principle), not need to symmetrize
    ! first W_mo
    CALL cp_fm_set_all(fm_back, 0.0_dp,error=error)
    CALL cp_gemm('N','N',dimen,dimen,dimen,1.0_dp,&
                    mo_coeff,W_mo,0.0_dp,fm_back,&
                    error=error,&
                    a_first_row=1,&
                    b_first_col=1,&
                    b_first_row=1,&
                    c_first_col=1,&
                    c_first_row=1)
    CALL cp_fm_set_all(W_mo, 0.0_dp,error=error)
    CALL cp_gemm('N','T',dimen,dimen,dimen,1.0_dp,&
                    fm_back,mo_coeff,0.0_dp,W_mo,&
                    error=error,&
                    a_first_row=1,&
                    b_first_col=1,&
                    b_first_row=1,&
                    c_first_col=1,&
                    c_first_row=1)
    ! and P_mo
    CALL cp_fm_set_all(fm_back, 0.0_dp,error=error)
    CALL cp_gemm('N','N',dimen,dimen,dimen,1.0_dp,&
                    mo_coeff,P_mo,0.0_dp,fm_back,&
                    error=error,&
                    a_first_row=1,&
                    b_first_col=1,&
                    b_first_row=1,&
                    c_first_col=1,&
                    c_first_row=1)
    CALL cp_fm_set_all(P_mo, 0.0_dp,error=error)
    CALL cp_gemm('N','T',dimen,dimen,dimen,1.0_dp,&
                    fm_back,mo_coeff,0.0_dp,P_mo,&
                    error=error,&
                    a_first_row=1,&
                    b_first_col=1,&
                    b_first_row=1,&
                    c_first_col=1,&
                    c_first_row=1)
 
    ! copy W_mo into dbcsr
    CALL copy_fm_to_dbcsr(W_mo, matrix_w_mp2(1)%matrix, keep_sparsity=.TRUE.,error=error)

    ! create mp2 DBCSR density
    CALL cp_dbcsr_allocate_matrix_set(matrix_p_mp2,qs_env%dft_control%nspins,error=error)
    DO ispin=1, dft_control%nspins
       ALLOCATE(matrix_p_mp2(ispin)%matrix)
       CALL cp_dbcsr_init(matrix_p_mp2(ispin)%matrix,error=error)
       CALL cp_dbcsr_copy(matrix_p_mp2(ispin)%matrix,rho%rho_ao(ispin)%matrix,&
            name="P MATRIX MP2",error=error)
       CALL cp_dbcsr_set(matrix_p_mp2(ispin)%matrix,0.0_dp,error=error)
       CALL copy_fm_to_dbcsr(P_mo, matrix_p_mp2(ispin)%matrix,keep_sparsity=.TRUE.,error=error)
    END DO
    CALL set_ks_env(ks_env, matrix_p_mp2=matrix_p_mp2, error=error)

    ! release remaining fm stuff
    CALL cp_fm_release(W_mo,error=error)
    CALL cp_fm_release(P_mo,error=error)
    CALL cp_fm_release(fm_back,error=error)

    ! update the core-forces with the MP2-density contribution
    ! put MP2 density on the grid
    CALL calculate_rho_elec(matrix_p=matrix_p_mp2(1)%matrix,&
                            rho=rho_r,&
                            rho_gspace=rho_g,&
                            total_rho=tot_rho_r,&
                            qs_env=qs_env,&
                            soft_valid=.FALSE.,&
                            error=error)
    ! calculate the MP2 potential
    CALL pw_transfer(rho_r%pw, rho_g%pw, error=error)
    CALL pw_poisson_solve(poisson_env,rho_g%pw, pair_energy, pot_g%pw,error=error)
    CALL pw_transfer(pot_g%pw, rho_r%pw, error=error)
    CALL pw_scale(rho_r%pw,rho_r%pw%pw_grid%dvol, error=error)
 
    ! calculate core forces
    CALL integrate_v_core_rspace(rho_r, qs_env, error=error)
    DO ikind=1,SIZE(atomic_kind_set)
       force(ikind)%mp2_sep=force(ikind)%rho_core
       force(ikind)%rho_core=0.0_dp
    ENDDO
    ! right contribution
    CALL integrate_v_rspace(v_rspace=rho_r,p=qs_env%rho%rho_ao(1),h=matrix_ks(1),&
                            qs_env=qs_env,calculate_forces=.TRUE.,error=error)

    IF(use_virial) THEN
      ! update virial if necessery with the volume term
      ! first create pw auxilliary stuff
      CALL timeset(routineN//"_Virial",handle2)
      NULLIFY(temp_pw_g%pw)
      CALL pw_pool_create_pw(auxbas_pw_pool,temp_pw_g%pw,&
                             use_data=COMPLEXDATA1D,&
                             in_space=RECIPROCALSPACE,error=error)
      DO i=1, 3
        NULLIFY(dvg(i)%pw)
        CALL pw_pool_create_pw(auxbas_pw_pool,dvg(i)%pw,&
                               use_data=COMPLEXDATA1D,&
                               in_space=RECIPROCALSPACE,error=error)
      END DO
     
      ! make a copy of the MP2 density in G space
      CALL pw_copy(rho_g%pw, temp_pw_g%pw, error=error)
      ! calculate MP2-like-hartree potential derivatives
      DO i=1, 3
        comp=0
        comp(i)=1
        CALL pw_copy(pot_g%pw, dvg(i)%pw, error=error)
        CALL pw_derive(dvg(i)%pw, comp, error=error)
      END DO
 
      ! calculate total SCF density and potential
      CALL calculate_rho_elec(matrix_p=qs_env%rho%rho_ao(1)%matrix,&
                              rho=rho_r,&
                              rho_gspace=rho_g,&
                              total_rho=tot_rho_r,&
                              qs_env=qs_env,&
                              soft_valid=.FALSE.,&
                              error=error)
      ! and associated potential
      CALL pw_transfer(rho_r%pw, rho_g%pw, error=error)
      ! don't forget the core density
      CALL pw_axpy(qs_env%rho_core%pw, rho_g%pw, error=error)
      CALL pw_poisson_solve(poisson_env,rho_g%pw, pair_energy, pot_g%pw,error=error)

      ! finally update virial with the volume contribution
      e_hartree=pw_integral_ab(temp_pw_g%pw, pot_g%pw, error=error)
      h_stress=0.0_dp
      DO alpha=1, 3
        comp=0
        comp(alpha)=1
        CALL pw_copy(pot_g%pw, rho_g%pw, error=error)
        CALL pw_derive(rho_g%pw, comp, error=error)
        h_stress(alpha,alpha)=-e_hartree
        DO beta=alpha, 3
          h_stress(alpha,beta)=h_stress(alpha,beta) &
                     -2.0_dp*pw_integral_ab(rho_g%pw, dvg(beta)%pw, error=error)/fourpi
          h_stress (beta,alpha)=h_stress(alpha,beta)
        END DO
      END DO
      virial%pv_virial = virial%pv_virial + h_stress/REAL(para_env%num_pe,dp)

      ! free stuff
      CALL pw_pool_give_back_pw(auxbas_pw_pool,temp_pw_g%pw,error=error)
      DO i=1, 3
        CALL pw_pool_give_back_pw(auxbas_pw_pool,dvg(i)%pw,error=error)
      END DO
      CALL timestop(handle2)
    END IF

    DO ispin=1, qs_env%dft_control%nspins
      CALL cp_dbcsr_add(qs_env%rho%rho_ao(ispin)%matrix, matrix_p_mp2(ispin)%matrix, 1.0_dp,  1.0_dp, error)
    END DO

    ! release stuff
    CALL pw_pool_give_back_pw(auxbas_pw_pool,rho_r%pw,error=error)
    CALL pw_pool_give_back_pw(auxbas_pw_pool,rho_g%pw,error=error)
    CALL pw_pool_give_back_pw(auxbas_pw_pool,pot_g%pw,error=error)

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
!> \param dft_control ...
!> \param hfx_sections ...
!> \param energy ...
!> \param n_rep_hf ...
!> \param pw_env ...
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
!> \param error ...
!> \param recalc_hfx_integrals ...
!> \author Mauro Del Ben 
! *****************************************************************************
  SUBROUTINE cphf_like_update(qs_env,mp2_env,para_env,homo,virtual,dimen,unit_nr,&
                              mo_coeff,mo_coeff_o,mo_coeff_v,Eigenval,dft_control,&
                              hfx_sections,energy,n_rep_hf,pw_env,poisson_env,&
                              rho_work,pot_g,rho_g,rho_r,mat_mu_nu,P_mu_nu,&
                              fm_mo,fm_ao,fm_back,transf_type_in,out_alpha,&
                              fm_mo_out,transf_type_out,error,recalc_hfx_integrals)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(mp2_type), POINTER                  :: mp2_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER                                  :: homo, virtual, dimen, unit_nr
    TYPE(cp_fm_type), POINTER                :: mo_coeff, mo_coeff_o, &
                                                mo_coeff_v
    REAL(KIND=dp), DIMENSION(:)              :: Eigenval
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(section_vals_type), POINTER         :: hfx_sections
    TYPE(qs_energy_type), POINTER            :: energy
    INTEGER                                  :: n_rep_hf
    TYPE(pw_env_type), POINTER               :: pw_env
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
    TYPE(cp_error_type), INTENT(inout)       :: error
    LOGICAL, OPTIONAL                        :: recalc_hfx_integrals

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cphf_like_update', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i_global, iiB, irep, &
                                                j_global, jjB, ncol_local, &
                                                nrow_local
    INTEGER, DIMENSION(:), POINTER           :: col_indices, row_indices
    LOGICAL                                  :: failure, &
                                                my_recalc_hfx_integrals
    REAL(KIND=dp)                            :: ex_energy, pair_energy, &
                                                total_rho

    CALL timeset(routineN,handle)
    failure=.FALSE.

    my_recalc_hfx_integrals=.FALSE.
    IF(PRESENT(recalc_hfx_integrals)) my_recalc_hfx_integrals=recalc_hfx_integrals

    ! perform back transformation
    SELECT CASE(transf_type_in)
      CASE(1)
        ! occ-occ block
        CALL cp_fm_set_all(fm_back, 0.0_dp,error=error)
        CALL cp_gemm('N','N',dimen,homo,homo,1.0_dp,&
                        mo_coeff_o,fm_mo,0.0_dp,fm_back,&
                        error=error,&
                        a_first_row=1,&
                        b_first_col=1,&
                        b_first_row=1,&
                        c_first_col=1,&
                        c_first_row=1)
        CALL cp_fm_set_all(fm_ao, 0.0_dp,error=error)
        CALL cp_gemm('N','T',dimen,dimen,homo,1.0_dp,&
                        fm_back,mo_coeff_o,0.0_dp,fm_ao,&
                        error=error,&
                        a_first_row=1,&
                        b_first_col=1,&
                        b_first_row=1,&
                        c_first_col=1,&
                        c_first_row=1)

      CASE(2)
        ! virt-virt block
        CALL cp_fm_set_all(fm_back, 0.0_dp,error=error)
        CALL cp_gemm('N','N',dimen,virtual,virtual,1.0_dp,&
                        mo_coeff_v,fm_mo,0.0_dp,fm_back,&
                        error=error,&
                        a_first_row=1,&
                        b_first_col=homo+1,&
                        b_first_row=homo+1,&
                        c_first_col=1,&
                        c_first_row=1)
        CALL cp_fm_set_all(fm_ao, 0.0_dp,error=error)
        CALL cp_gemm('N','T',dimen,dimen,virtual,1.0_dp,&
                        fm_back,mo_coeff_v,0.0_dp,fm_ao,&
                        error=error,&
                        a_first_row=1,&
                        b_first_col=1,&
                        b_first_row=1,&
                        c_first_col=1,&
                        c_first_row=1)

      CASE(3)
        CALL cp_fm_set_all(fm_back, 0.0_dp,error=error)
        CALL cp_gemm('N','N',dimen,virtual,homo,1.0_dp,&
                        mo_coeff_o,fm_mo,0.0_dp,fm_back,&
                        error=error,&
                        a_first_row=1,&
                        b_first_col=1,&
                        b_first_row=1,&
                        c_first_col=1,&
                        c_first_row=1)
        CALL cp_fm_set_all(fm_ao, 0.0_dp,error=error)
        CALL cp_gemm('N','T',dimen,dimen,virtual,1.0_dp,&
                        fm_back,mo_coeff_v,0.0_dp,fm_ao,&
                        error=error,&
                        a_first_row=1,&
                        b_first_col=1,&
                        b_first_row=1,&
                        c_first_col=1,&
                        c_first_row=1)
        ! and symmetrize (here again multiply instead of transposing)
        CALL cp_fm_set_all(fm_back, 0.0_dp,error=error)
        CALL cp_gemm('N','T',dimen,homo,virtual,1.0_dp,&
                        mo_coeff_v,fm_mo,0.0_dp,fm_back,&
                        error=error,&
                        a_first_row=1,&
                        b_first_col=1,&
                        b_first_row=1,&
                        c_first_col=1,&
                        c_first_row=1)
        CALL cp_gemm('N','T',dimen,dimen,homo,0.5_dp,&
                        fm_back,mo_coeff_o,0.5_dp,fm_ao,&
                        error=error,&
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
                            col_indices=col_indices,&
                            error=error)
        DO jjB=1, ncol_local
          j_global=col_indices(jjB)
          DO iiB=1, nrow_local
            i_global=row_indices(iiB)
            fm_mo_out%local_data(iiB,jjB)=fm_mo_out%local_data(iiB,jjB)*&
                                          (Eigenval(j_global+homo)-Eigenval(i_global))
          END DO
        END DO

      CASE(4)
        ! all-all block
        CALL cp_fm_set_all(fm_back, 0.0_dp,error=error)
        CALL cp_gemm('N','N',dimen,dimen,dimen,1.0_dp,&
                        mo_coeff,fm_mo,0.0_dp,fm_back,&
                        error=error,&
                        a_first_row=1,&
                        b_first_col=1,&
                        b_first_row=1,&
                        c_first_col=1,&
                        c_first_row=1)
        CALL cp_fm_set_all(fm_ao, 0.0_dp,error=error)
        CALL cp_gemm('N','T',dimen,dimen,dimen,1.0_dp,&
                        fm_back,mo_coeff,0.0_dp,fm_ao,&
                        error=error,&
                        a_first_row=1,&
                        b_first_col=1,&
                        b_first_row=1,&
                        c_first_col=1,&
                        c_first_row=1)

      CASE DEFAULT
        ! nothing
    END SELECT

    ! copy fm into DBCSR
    CALL cp_dbcsr_set(P_mu_nu%matrix,0.0_dp,error=error)
    CALL copy_fm_to_dbcsr(fm_ao, P_mu_nu%matrix, keep_sparsity=.TRUE., error=error)

    ! calculate associated density 
    CALL calculate_rho_elec(matrix_p=P_mu_nu%matrix,&
                            rho=rho_r,&
                            rho_gspace=rho_g,&
                            total_rho=total_rho,&
                            qs_env=qs_env,error=error)
    ! and calculate potential
    CALL pw_poisson_solve(poisson_env, rho_g%pw, pair_energy, pot_g%pw, error=error)
    CALL pw_transfer(pot_g%pw, rho_r%pw, error=error)
    CALL pw_scale(rho_r%pw,rho_r%pw%pw_grid%dvol, error=error)
    ! integrate the potential
    CALL cp_dbcsr_set(mat_mu_nu(1)%matrix,0.0_dp,error=error)
    CALL integrate_v_rspace(rho_r, h=mat_mu_nu(1), &
                            qs_env=qs_env,calculate_forces=.FALSE.,compute_tau=.FALSE.,gapw=.FALSE.,&
                            basis_set_id=use_orb_basis_set,&
                            error=error)
    ! update with the exchange like contributions
    ! copy mat_mu_nu into rho_ao work
    CALL cp_dbcsr_set(rho_work%rho_ao(1)%matrix,0.0_dp,error=error)
    CALL cp_dbcsr_copy(rho_work%rho_ao(1)%matrix,P_mu_nu%matrix,error=error)
    ! save old EX energy
    ex_energy=energy%ex
    DO irep=1, n_rep_hf
      CALL integrate_four_center(qs_env, mat_mu_nu, energy, rho_work%rho_ao, hfx_sections,&
                                 para_env, my_recalc_hfx_integrals, irep, .TRUE.,&
                                 ispin=1, error=error)
    END DO
    ! restore original EX energy
    energy%ex=ex_energy

    ! scale by factor 4.0
    CALL cp_dbcsr_scale(mat_mu_nu(1)%matrix,4.0_dp,error=error)

    ! copy back to fm
    CALL cp_fm_set_all(fm_ao, 0.0_dp,error=error)
    CALL copy_dbcsr_to_fm(matrix=mat_mu_nu(1)%matrix, fm=fm_ao, error=error)
    CALL cp_fm_set_all(fm_back, 0.0_dp,error=error)
    CALL cp_fm_upper_to_full(fm_ao, fm_back, error)

    ! transform to MO basis, here we always sum the result into the input matrix
    SELECT CASE(transf_type_out)
      CASE(1)
        ! occ-virt block
        CALL cp_fm_set_all(fm_back, 0.0_dp,error=error)
        CALL cp_gemm('T','N',homo,dimen,dimen,1.0_dp,&
                        mo_coeff_o,fm_ao,0.0_dp,fm_back,&
                        error=error,&
                        a_first_row=1,&
                        b_first_col=1,&
                        b_first_row=1,&
                        c_first_col=1,&
                        c_first_row=1)
        CALL cp_gemm('N','N',homo,virtual,dimen,out_alpha,&
                        fm_back,mo_coeff_v,1.0_dp,fm_mo_out,&
                        error=error,&
                        a_first_row=1,&
                        b_first_col=1,&
                        b_first_row=1,&
                        c_first_col=1,&
                        c_first_row=1)

      CASE(2)
        ! occ-occ block
        CALL cp_fm_set_all(fm_back, 0.0_dp,error=error)
        CALL cp_gemm('T','N',homo,dimen,dimen,1.0_dp,&
                        mo_coeff_o,fm_ao,0.0_dp,fm_back,&
                        error=error,&
                        a_first_row=1,&
                        b_first_col=1,&
                        b_first_row=1,&
                        c_first_col=1,&
                        c_first_row=1)
        CALL cp_gemm('N','N',homo,homo,dimen,out_alpha,&
                        fm_back,mo_coeff_o,1.0_dp,fm_mo_out,&
                        error=error,&
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
!> \param dft_control ...
!> \param hfx_sections ...
!> \param energy ...
!> \param n_rep_hf ...
!> \param pw_env ...
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
!> \param error ...
!> \author Mauro Del Ben 
! *****************************************************************************
  SUBROUTINE solve_z_vector_eq_low(qs_env,mp2_env,para_env,homo,virtual,dimen,unit_nr,&
                                   mo_coeff,mo_coeff_o,mo_coeff_v,Eigenval,blacs_env,dft_control,&
                                   hfx_sections,energy,n_rep_hf,pw_env,poisson_env,&
                                   rho_work,pot_g,rho_g,rho_r,mat_mu_nu,P_mu_nu,&
                                   L_jb,fm_G_mu_nu,fm_back,P_ia,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(mp2_type), POINTER                  :: mp2_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER                                  :: homo, virtual, dimen, unit_nr
    TYPE(cp_fm_type), POINTER                :: mo_coeff, mo_coeff_o, &
                                                mo_coeff_v
    REAL(KIND=dp), DIMENSION(:)              :: Eigenval
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(section_vals_type), POINTER         :: hfx_sections
    TYPE(qs_energy_type), POINTER            :: energy
    INTEGER                                  :: n_rep_hf
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_poisson_type), POINTER           :: poisson_env
    TYPE(qs_rho_type), POINTER               :: rho_work
    TYPE(pw_p_type)                          :: pot_g, rho_g, rho_r
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: mat_mu_nu
    TYPE(cp_dbcsr_p_type)                    :: P_mu_nu
    TYPE(cp_fm_type), POINTER                :: L_jb, fm_G_mu_nu, fm_back, &
                                                P_ia
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'solve_z_vector_eq_low', &
      routineP = moduleN//':'//routineN

    INTEGER :: cycle_counter, handle, i_global, iiB, iiter, j_global, jjB, &
      max_num_iter, ncol_local, nrow_local, stat, transf_type_in, &
      transf_type_out
    INTEGER, DIMENSION(:), POINTER           :: col_indices, row_indices
    LOGICAL                                  :: converged, failure
    REAL(KIND=dp)                            :: alpha, beta, conv, eps_conv, &
                                                norm_b, norms(3), out_alpha, &
                                                rkrk, t1, t2
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: proj_bi_xj, temp_vals, &
                                                x_norm, xi_b
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: A_small, b_small, xi_Axi
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: Ax, xn
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: Ap, b_i, pk, precond, &
                                                residual, rk, xk

    CALL timeset(routineN,handle)
    failure=.FALSE.

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

    IF(.FALSE.) THEN
      ! CG algorithm
      ! create some work array
      NULLIFY(xk, pk, rk, Ap, fm_struct_tmp)
      CALL cp_fm_struct_create(fm_struct_tmp,para_env=para_env,context=blacs_env, &
                               nrow_global=homo,ncol_global=virtual,error=error)
      CALL cp_fm_create(xk, fm_struct_tmp, name="xk",error=error)
      CALL cp_fm_create(pk, fm_struct_tmp, name="pk",error=error)
      CALL cp_fm_create(rk, fm_struct_tmp, name="rk",error=error)
      CALL cp_fm_create(Ap, fm_struct_tmp, name="Ap",error=error)
      CALL cp_fm_struct_release(fm_struct_tmp,error=error)
      CALL cp_fm_set_all(xk, 0.0_dp,error=error)
      CALL cp_fm_set_all(pk, 0.0_dp,error=error)
      CALL cp_fm_set_all(rk, 0.0_dp,error=error)
      CALL cp_fm_set_all(Ap, 0.0_dp,error=error)

      ! copy -L_jb into pk and rk
      pk%local_data(:,:)=-L_jb%local_data(:,:)
      rk%local_data(:,:)=-L_jb%local_data(:,:)
      norm_b=0.0_dp
      norm_b=SUM(L_jb%local_data(:,:)*L_jb%local_data(:,:))
      CALL mp_sum(norm_b,para_env%group)
      norm_b=SQRT(norm_b)

      cycle_counter=0
      DO iiter=1, max_num_iter
        cycle_counter=cycle_counter+1
        t1 = m_walltime()

        ! calculate matrix-vector product
        CALL cp_fm_set_all(Ap, 0.0_dp,error=error)
        CALL cphf_like_update(qs_env,mp2_env,para_env,homo,virtual,dimen,unit_nr,&
                              mo_coeff,mo_coeff_o,mo_coeff_v,Eigenval,dft_control,&
                              hfx_sections,energy,n_rep_hf,pw_env,poisson_env,&
                              rho_work,pot_g,rho_g,rho_r,mat_mu_nu,P_mu_nu,&
                              pk,fm_G_mu_nu,fm_back,transf_type_in,out_alpha,&
                              Ap,transf_type_out,error)

        norms=0.0_dp
        norms(1)=SUM(rk%local_data(:,:)*rk%local_data(:,:))
        norms(2)=SUM(rk%local_data(:,:)*pk%local_data(:,:))
        norms(3)=SUM(pk%local_data(:,:)*Ap%local_data(:,:))
        CALL mp_sum(norms,para_env%group)
        alpha=norms(1)/norms(3)

        xk%local_data(:,:)=xk%local_data(:,:)+alpha*pk%local_data(:,:)
        rk%local_data(:,:)=rk%local_data(:,:)-alpha*Ap%local_data(:,:)

        rkrk=0.0_dp
        rkrk=SUM(rk%local_data(:,:)*rk%local_data(:,:))
        CALL mp_sum(rkrk,para_env%group) 
        beta=rkrk/norms(2)

        pk%local_data(:,:)=rk%local_data(:,:)+beta*pk%local_data(:,:)

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

      CALL cp_fm_release(xk,error=error)
      CALL cp_fm_release(pk,error=error)
      CALL cp_fm_release(rk,error=error)
      CALL cp_fm_release(Ap,error=error)

    ELSE
      ! Pople method
      ! change sign to L_jb
      L_jb%local_data(:,:)=-L_jb%local_data(:,:)

      ! allocate stuff
      ALLOCATE(xn(1:max_num_iter),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(Ax(1:max_num_iter),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

      ! create fm structure
      NULLIFY(fm_struct_tmp)
      CALL cp_fm_struct_create(fm_struct_tmp,para_env=para_env,context=blacs_env, &
                               nrow_global=homo,ncol_global=virtual,error=error)

      ! create preconditioner (for now only orbital energy differences)
      NULLIFY(precond)
      CALL cp_fm_create(precond, fm_struct_tmp, name="precond",error=error)
      CALL cp_fm_set_all(precond,1.0_dp,error=error)
      CALL cp_fm_get_info(matrix=precond,&
                          nrow_local=nrow_local,&
                          ncol_local=ncol_local,&
                          row_indices=row_indices,&
                          col_indices=col_indices,&
                          error=error)
      DO jjB=1, ncol_local
        j_global=col_indices(jjB)
        DO iiB=1, nrow_local
          i_global=row_indices(iiB)
          precond%local_data(iiB,jjB)=precond%local_data(iiB,jjB)/&
                                      (Eigenval(j_global+homo)-Eigenval(i_global))
        END DO
      END DO

      ! create b_i, work array needed for the orthogonalization of the 
      ! x(iiter) vector
      NULLIFY(b_i)
      CALL cp_fm_create(b_i, fm_struct_tmp, name="b_i",error=error)
      CALL cp_fm_set_all(b_i, 0.0_dp, error=error)
      b_i%local_data(:,:)=precond%local_data(:,:)*L_jb%local_data(:,:)

      ! create the residual vector (r), we check convergence on the norm of 
      ! this vector r=(Ax-b) 
      NULLIFY(residual)
      CALL cp_fm_create(residual, fm_struct_tmp, name="residual",error=error)
      CALL cp_fm_set_all(residual, 0.0_dp, error=error)

      ! allocate array containing the various scalar products
      ALLOCATE(x_norm(1:max_num_iter),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(xi_b(1:max_num_iter),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(xi_Axi(1:max_num_iter,0:max_num_iter),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      x_norm=0.0_dp
      xi_b=0.0_dp
      xi_Axi=0.0_dp

      cycle_counter=0
      DO iiter=1, max_num_iter
        cycle_counter=cycle_counter+1
 
        t1 = m_walltime()

        ! create and update x_i (orthogonalization with previous vectors)
        NULLIFY(xn(iiter)%matrix)
        CALL cp_fm_create(xn(iiter)%matrix, fm_struct_tmp, name="xi",error=error)
        CALL cp_fm_set_all(xn(iiter)%matrix, 0.0_dp, error=error)

        ! first compute the projection of the actual b_i into all previous x_i
        ! already scaled with the norm of each x_i
        ALLOCATE(proj_bi_xj(iiter-1),STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
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
        DEALLOCATE(proj_bi_xj,STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

        ! create Ax(iiter) that will store the matrix vector product for this cycle
        NULLIFY(Ax(iiter)%matrix)
        CALL cp_fm_create(Ax(iiter)%matrix, fm_struct_tmp, name="Ai",error=error)
        CALL cp_fm_set_all(Ax(iiter)%matrix, 0.0_dp, error=error)
        ! performe the matrix-vector product (CPHF like update) 
        CALL cphf_like_update(qs_env,mp2_env,para_env,homo,virtual,dimen,unit_nr,&
                              mo_coeff,mo_coeff_o,mo_coeff_v,Eigenval,dft_control,&
                              hfx_sections,energy,n_rep_hf,pw_env,poisson_env,&
                              rho_work,pot_g,rho_g,rho_r,mat_mu_nu,P_mu_nu,&
                              xn(iiter)%matrix,fm_G_mu_nu,fm_back,transf_type_in,out_alpha,&
                              Ax(iiter)%matrix,transf_type_out,error)

        ! in order to reduce the number of calls to mp_sum here we
        ! cluster all necessary scalar products into a sigle vector
        ! temp_vals contains:
        ! 1:iiter -> <Ax_i|x_j>
        ! iiter+1 -> <x_i|b>
        ! iiter+2 -> <x_i|x_i>
        ALLOCATE(temp_vals(iiter+2),STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
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
        DEALLOCATE(temp_vals,STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

        ! solve reduced system 
        IF(ALLOCATED(A_small)) DEALLOCATE(A_small)
        IF(ALLOCATED(b_small)) DEALLOCATE(b_small)
        ALLOCATE(A_small(iiter,iiter),STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        ALLOCATE(b_small(iiter,1),STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        A_small(1:iiter,1:iiter)=xi_Axi(1:iiter,1:iiter)
        b_small(1:iiter,1)=xi_b(1:iiter)

        CALL solve_system(matrix=A_small, mysize=iiter, eigenvectors=b_small)

        ! check for convergence
        CALL cp_fm_set_all(residual, 0.0_dp, error=error)
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

      DEALLOCATE(x_norm,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      DEALLOCATE(xi_b,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      DEALLOCATE(xi_Axi,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

      CALL cp_fm_release(precond,error=error)
      CALL cp_fm_release(b_i,error=error)
      CALL cp_fm_release(residual,error=error)
      CALL cp_fm_struct_release(fm_struct_tmp,error=error)

      ! release Ax, xn
      DO iiter=1, cycle_counter
        CALL cp_fm_release(Ax(iiter)%matrix,error=error)
        CALL cp_fm_release(xn(iiter)%matrix,error=error)
      END DO
      DEALLOCATE(xn,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      DEALLOCATE(Ax,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

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
