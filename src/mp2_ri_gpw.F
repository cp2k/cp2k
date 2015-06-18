!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Rountines to calculate RI-GPW-MP2 energy using pw
!> \par History
!>      06.2012 created [Mauro Del Ben]
! *****************************************************************************
MODULE mp2_ri_gpw
  USE ai_coulomb,                      ONLY: coulomb2_new
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: gto_basis_set_type
  USE bibliography,                    ONLY: DelBen2013,&
                                             cite_reference
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                             cp_blacs_env_release,&
                                             cp_blacs_env_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_create, cp_dbcsr_get_info, cp_dbcsr_init, cp_dbcsr_multiply, &
       cp_dbcsr_p_type, cp_dbcsr_release, cp_dbcsr_set, cp_dbcsr_type, &
       dbcsr_type_no_symmetry
  USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                             cp_dbcsr_m_by_n_from_template
  USE cp_fm_basic_linalg,              ONLY: cp_fm_triangular_invert
  USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose
  USE cp_fm_diag,                      ONLY: cp_fm_syevx
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_info,&
                                             cp_fm_release,&
                                             cp_fm_set_all,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_para_env,                     ONLY: cp_para_env_create,&
                                             cp_para_env_release
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE gaussian_gridlevels,             ONLY: gaussian_gridlevel
  USE input_constants,                 ONLY: use_orb_basis_set,&
                                             use_ri_aux_basis_set
  USE kinds,                           ONLY: dp,&
                                             int_8
  USE machine,                         ONLY: m_flush,&
                                             m_memory,&
                                             m_walltime
  USE mathconstants,                   ONLY: fourpi,&
                                             pi
  USE message_passing,                 ONLY: mp_allgather,&
                                             mp_comm_split_direct,&
                                             mp_max,&
                                             mp_min,&
                                             mp_sendrecv,&
                                             mp_sum,&
                                             mp_sync
  USE molecule_kind_types,             ONLY: molecule_kind_type
  USE molecule_types_new,              ONLY: molecule_type
  USE mp2_ri_grad_util,                ONLY: complete_gamma
  USE mp2_types,                       ONLY: mp2_type
  USE orbital_pointers,                ONLY: ncoset
  USE particle_types,                  ONLY: particle_type
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_grid_types,                   ONLY: pw_grid_type
  USE pw_methods,                      ONLY: pw_copy,&
                                             pw_scale,&
                                             pw_transfer
  USE pw_poisson_methods,              ONLY: pw_poisson_solve
  USE pw_poisson_types,                ONLY: pw_poisson_type
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_type
  USE pw_types,                        ONLY: REALDATA3D,&
                                             REALSPACE,&
                                             pw_p_type,&
                                             pw_type
  USE qs_collocate_density,            ONLY: calculate_wavefunction
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_integrate_potential,          ONLY: integrate_pgf_product_rspace,&
                                             integrate_v_rspace
  USE qs_kind_types,                   ONLY: qs_kind_type
  USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                             neighbor_list_iterate,&
                                             neighbor_list_iterator_create,&
                                             neighbor_list_iterator_p_type,&
                                             neighbor_list_iterator_release,&
                                             neighbor_list_set_p_type
  USE realspace_grid_types,            ONLY: realspace_grid_desc_p_type,&
                                             realspace_grid_p_type,&
                                             rs_grid_release,&
                                             rs_grid_retain
  USE rs_pw_interface,                 ONLY: potential_pw2rs
  USE task_list_types,                 ONLY: task_list_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: get_limit
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mp2_ri_gpw'

  PUBLIC :: mp2_ri_gpw_compute_in,&
            mp2_ri_gpw_compute_en


  CONTAINS

! *****************************************************************************
!> \brief with ri mp2 gpw
!> \param BIb_C ...
!> \param BIb_C_gw ...
!> \param ends_array ...
!> \param ends_B_virtual ...
!> \param sizes_array ...
!> \param sizes_B_virtual ...
!> \param starts_array ...
!> \param starts_B_virtual ...
!> \param dimen_RI ...
!> \param qs_env ...
!> \param para_env ...
!> \param para_env_sub ...
!> \param color_sub ...
!> \param dft_control ...
!> \param cell ...
!> \param particle_set ...
!> \param atomic_kind_set ...
!> \param qs_kind_set ...
!> \param mo_coeff ...
!> \param nmo ...
!> \param homo ...
!> \param rho_r ...
!> \param rho_g ...
!> \param pot_g ...
!> \param mat_munu ...
!> \param sab_orb_sub ...
!> \param pw_env_sub ...
!> \param poisson_env ...
!> \param auxbas_pw_pool ...
!> \param task_list_sub ...
!> \param mo_coeff_o ...
!> \param mo_coeff_v ...
!> \param mo_coeff_all ...
!> \param mo_coeff_gw ...
!> \param eps_filter ...
!> \param unit_nr ...
!> \param mp2_memory ...
!> \param calc_PQ_cond_num ...
!> \param calc_forces ...
!> \param error ...
!> \param blacs_env_sub ...
!> \param my_do_gw ...
!> \param starts_B_all ...
!> \param sizes_B_all ...
!> \param ends_B_all ...
!> \param gw_corr_lev_occ ...
!> \param gw_corr_lev_virt ...
!> \param rlength ...
!> \param do_truncation_gw ...
!> \param error_sub ...
!> \param BIb_C_beta ...
!> \param BIb_C_gw_beta ...
!> \param ends_B_virtual_beta ...
!> \param sizes_B_virtual_beta ...
!> \param starts_B_virtual_beta ...
!> \param homo_beta ...
!> \param mo_coeff_o_beta ...
!> \param mo_coeff_v_beta ...
!> \param mo_coeff_all_beta ...
!> \param mo_coeff_gw_beta ...
!> \author Mauro Del Ben 
! *****************************************************************************
  SUBROUTINE mp2_ri_gpw_compute_in(BIb_C,BIb_C_gw,ends_array,ends_B_virtual,&
                                   sizes_array,sizes_B_virtual,starts_array,starts_B_virtual,&
                                   dimen_RI,qs_env,para_env,para_env_sub,color_sub,dft_control,&
                                   cell,particle_set,atomic_kind_set,qs_kind_set,mo_coeff,nmo,homo,&
                                   rho_r,rho_g,pot_g,mat_munu,sab_orb_sub,pw_env_sub,&
                                   poisson_env,auxbas_pw_pool,task_list_sub,mo_coeff_o,mo_coeff_v,mo_coeff_all,&
                                   mo_coeff_gw,eps_filter,unit_nr,&
                                   mp2_memory,calc_PQ_cond_num,calc_forces,error,blacs_env_sub,my_do_gw,&
                                   starts_B_all,sizes_B_all,ends_B_all,gw_corr_lev_occ,gw_corr_lev_virt,&
                                   rlength,do_truncation_gw,error_sub,&
                                   BIb_C_beta,BIb_C_gw_beta,ends_B_virtual_beta,sizes_B_virtual_beta,starts_B_virtual_beta,&
                                   homo_beta,mo_coeff_o_beta,mo_coeff_v_beta,mo_coeff_all_beta,mo_coeff_gw_beta)
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: BIb_C, BIb_C_gw
    INTEGER, ALLOCATABLE, DIMENSION(:) :: ends_array, ends_B_virtual, &
      sizes_array, sizes_B_virtual, starts_array, starts_B_virtual
    INTEGER                                  :: dimen_RI
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_para_env_type), POINTER          :: para_env, para_env_sub
    INTEGER                                  :: color_sub
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(cell_type), POINTER                 :: cell
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    INTEGER                                  :: nmo, homo
    TYPE(pw_p_type)                          :: rho_r, rho_g, pot_g
    TYPE(cp_dbcsr_p_type)                    :: mat_munu
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb_sub
    TYPE(pw_env_type), POINTER               :: pw_env_sub
    TYPE(pw_poisson_type), POINTER           :: poisson_env
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(task_list_type), POINTER            :: task_list_sub
    TYPE(cp_dbcsr_type), POINTER             :: mo_coeff_o, mo_coeff_v, &
                                                mo_coeff_all, mo_coeff_gw
    REAL(KIND=dp)                            :: eps_filter
    INTEGER                                  :: unit_nr
    REAL(KIND=dp)                            :: mp2_memory
    LOGICAL                                  :: calc_PQ_cond_num, calc_forces
    TYPE(cp_error_type), INTENT(inout)       :: error
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env_sub
    LOGICAL                                  :: my_do_gw
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: starts_B_all, sizes_B_all, &
                                                ends_B_all
    INTEGER                                  :: gw_corr_lev_occ, &
                                                gw_corr_lev_virt
    REAL(KIND=dp)                            :: rlength
    LOGICAL                                  :: do_truncation_gw
    TYPE(cp_error_type), INTENT(inout)       :: error_sub
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :), OPTIONAL           :: BIb_C_beta, BIb_C_gw_beta
    INTEGER, ALLOCATABLE, DIMENSION(:), &
      OPTIONAL                               :: ends_B_virtual_beta, &
                                                sizes_B_virtual_beta, &
                                                starts_B_virtual_beta
    INTEGER, OPTIONAL                        :: homo_beta
    TYPE(cp_dbcsr_type), OPTIONAL, POINTER   :: mo_coeff_o_beta, &
                                                mo_coeff_v_beta, &
                                                mo_coeff_all_beta, &
                                                mo_coeff_gw_beta

    CHARACTER(LEN=*), PARAMETER :: routineN = 'mp2_ri_gpw_compute_in', &
      routineP = moduleN//':'//routineN

    INTEGER :: gw_corr_lev_total, handle, handle2, handle3, i, i_counter, ig, &
      iiB, itmp(2), jjB, LLL, max_row_col_local, max_row_col_local_beta, &
      max_row_col_local_gw, my_B_all_end, my_B_all_size, my_B_all_start, &
      my_B_size, my_B_size_beta, my_B_virtual_end, my_B_virtual_end_beta, &
      my_B_virtual_start, my_B_virtual_start_beta, my_group_L_end, &
      my_group_L_size, my_group_L_start, ncol_local, nfullcols_total, &
      nfullrows_total, ngroup, nrow_local, num_small_eigen, stat, virtual, &
      virtual_beta
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: sub_proc_map
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: local_col_row_info, &
                                                local_col_row_info_beta, &
                                                local_col_row_info_gw
    LOGICAL                                  :: do_alpha_beta, failure
    REAL(KIND=dp)                            :: cond_num, g2, g3d, gg, &
                                                mem_for_iaK, pair_energy, &
                                                wfn_size
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: my_Lrows
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: BIb_C_gw_trunc
    TYPE(cp_dbcsr_type) :: matrix_ia_jb, matrix_ia_jb_beta, matrix_ia_jnu, &
      matrix_ia_jnu_beta, matrix_in_jm, matrix_in_jm_beta, matrix_in_jnu, &
      matrix_in_jnu_beta
    TYPE(cp_fm_type), POINTER                :: fm_BIb_gw, fm_BIb_gw_beta, &
                                                fm_BIb_jb, fm_BIb_jb_beta, &
                                                fm_matrix_L
    TYPE(cp_para_env_type), POINTER          :: para_env_L
    TYPE(pw_grid_type), POINTER              :: grid
    TYPE(pw_p_type)                          :: psi_L
    TYPE(pw_type), POINTER                   :: gf

    CALL timeset(routineN,handle)

    failure=.FALSE.


    CALL cite_reference(DelBen2013)


    do_alpha_beta=.FALSE.
    IF(PRESENT(BIb_C_beta).AND.&
       PRESENT(ends_B_virtual_beta).AND.&
       PRESENT(sizes_B_virtual_beta).AND.& 
       PRESENT(starts_B_virtual_beta).AND.& 
       PRESENT(homo_beta).AND.&
       PRESENT(mo_coeff_o_beta).AND.&
       PRESENT(mo_coeff_v_beta)) do_alpha_beta=.TRUE.

    virtual=nmo-homo

    CALL create_intermediate_matrices(matrix_ia_jnu,matrix_ia_jb,mo_coeff_o,virtual,homo,&
                                      fm_BIb_jb,"fm_BIb_jb",max_row_col_local,&
                                      nfullrows_total,nfullcols_total,nrow_local, &
                                      ncol_local,blacs_env_sub,&
                                      para_env_sub,local_col_row_info,error_sub,error)

    CALL create_parallelization_arrays(para_env_sub,ends_B_virtual,sizes_B_virtual,&
                                       starts_B_virtual,my_B_virtual_start,&
                                       my_B_virtual_end,my_B_size,&
                                       virtual,error)


    IF(do_alpha_beta) THEN

      virtual_beta=nmo-homo_beta

      CALL create_intermediate_matrices(matrix_ia_jnu_beta,matrix_ia_jb_beta,mo_coeff_o_beta,&
                                        virtual_beta,homo_beta,&
                                        fm_BIb_jb_beta,"fm_BIb_jb_beta",&
                                        max_row_col_local_beta,&
                                        nfullrows_total,nfullcols_total,nrow_local, &
                                        ncol_local,blacs_env_sub,&
                                        para_env_sub,local_col_row_info_beta,error_sub,error)

      CALL create_parallelization_arrays(para_env_sub,ends_B_virtual_beta,sizes_B_virtual_beta,&
                                         starts_B_virtual_beta,my_B_virtual_start_beta,&
                                         my_B_virtual_end_beta,my_B_size_beta,&
                                         virtual_beta,error)

    END IF

    ! in the case of G0W0, we need (K|nm), n,m may be occ or virt (m restricted to corrected levels)
    IF(my_do_gw) THEN

      gw_corr_lev_total=gw_corr_lev_virt+gw_corr_lev_occ    

      CALL create_intermediate_matrices(matrix_in_jnu,matrix_in_jm,mo_coeff_gw,&
                                        nmo,gw_corr_lev_total,&
                                        fm_BIb_gw,"fm_BIb_gw",&
                                        max_row_col_local_gw,&
                                        nfullrows_total,nfullcols_total,nrow_local, &
                                        ncol_local,blacs_env_sub,&
                                        para_env_sub,local_col_row_info_gw,error_sub,error)

      CALL create_parallelization_arrays(para_env_sub,ends_B_all,sizes_B_all,&
                                         starts_B_all,my_B_all_start,&
                                         my_B_all_end,my_B_all_size,&
                                         nmo,error)

      IF(do_alpha_beta) THEN
        ! deallocate local_col_row_info_gw, otherwise it gets twice allocated in create_intermediate_m
        DEALLOCATE(local_col_row_info_gw)
        CALL create_intermediate_matrices(matrix_in_jnu_beta,matrix_in_jm_beta,mo_coeff_gw_beta,&
                                          nmo,gw_corr_lev_total,&
                                          fm_BIb_gw_beta,"fm_BIb_gw_beta",&
                                          max_row_col_local_gw,&
                                          nfullrows_total,nfullcols_total,nrow_local, &
                                          ncol_local,blacs_env_sub,&
                                          para_env_sub,local_col_row_info_gw,error_sub,error)

        ! we don"t need parallelization arrays for beta since the matrix sizes of B_nm^P is the same
        ! for the beta case and therefore the parallelization of beta is the same than for alpha

      END IF
    END IF


    wfn_size=REAL(SIZE(rho_r%pw%cr3d),KIND=dp)
    CALL mp_max(wfn_size,para_env%group)

    ngroup=para_env%num_pe/para_env_sub%num_pe

    ALLOCATE(sub_proc_map(-para_env_sub%num_pe:2*para_env_sub%num_pe-1))
    DO i=0,para_env_sub%num_pe-1
      sub_proc_map(i)=i
      sub_proc_map(-i-1)=para_env_sub%num_pe-i-1
      sub_proc_map(para_env_sub%num_pe+i)=i
    END DO

    ! start real calculation
    NULLIFY(psi_L%pw)
    CALL pw_pool_create_pw(auxbas_pw_pool,psi_L%pw,&
                            use_data=REALDATA3D,&
                            in_space=REALSPACE,error=error_sub)

    ! calculate L^{-1}
    CALL calculate_Lmin1(qs_env,para_env,para_env_sub,para_env_L,mp2_memory,&
                         eps_filter,fm_matrix_L,ngroup,color_sub,dimen_RI,&
                         mo_coeff,dft_control,psi_L,rho_r,rho_g,pot_g,pw_env_sub,poisson_env,&
                         my_group_L_size,my_group_L_start,my_group_L_end,sab_orb_sub,&
                         sizes_array,starts_array,ends_array,calc_PQ_cond_num,cond_num,&
                         num_small_eigen,auxbas_pw_pool,error_sub,error)

    IF (unit_nr>0) WRITE (UNIT=unit_nr,FMT="(T3,A,T75,i6)")&
                      "RI_INFO| Cholesky decomposition group size:", para_env_L%num_pe
    IF(calc_PQ_cond_num) THEN
      IF (unit_nr>0) WRITE (UNIT=unit_nr,FMT="(T3,A,T67,ES14.5)")&
                        "RI_INFO| Condition number of the (P|Q):", cond_num
      IF (unit_nr>0) WRITE (UNIT=unit_nr,FMT="(T3,A,T75,i6)")&
                        "RI_INFO| Number of Eigenvalue of (P|Q) smaller than 10^(-3):",num_small_eigen
    END IF

    ! replicate the necessary row of the L^{-1} matrix on each proc
    CALL grep_Lcols(para_env_L,dimen_RI,fm_matrix_L,&
                    my_group_L_start,my_group_L_end,my_group_L_size,my_Lrows,&
                    error)
    ! clean the L^{-1} matrix
    CALL cp_fm_release(fm_matrix_L, error=error)
    CALL cp_para_env_release(para_env_L,error=error_sub)

    IF(calc_forces) THEN
      ! we need (P|Q)^(-1/2) for future use, just save it
      ! in a fully (home made) distributed way 
      itmp=get_limit(dimen_RI,para_env_sub%num_pe,para_env_sub%mepos)
      lll=itmp(2)-itmp(1)+1
      ALLOCATE(qs_env%mp2_env%ri_grad%PQ_half(lll,my_group_L_size),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      qs_env%mp2_env%ri_grad%PQ_half(:,:)=my_Lrows(itmp(1):itmp(2),1:my_group_L_size)
    END IF

    IF (unit_nr>0) THEN
      WRITE (UNIT=unit_nr,FMT="(T3,A,T75,i6)")&
               "RI_INFO| Occupied  basis set size:", homo,&
               "RI_INFO| Virtual   basis set size:", virtual,&
               "RI_INFO| Auxiliary basis set size:", dimen_RI 

      mem_for_iaK=dimen_RI*REAL(homo,KIND=dp)*virtual*8.0_dp/(1024_dp**2)
      IF(do_alpha_beta) mem_for_iaK=mem_for_iaK+&
                        dimen_RI*REAL(homo_beta,KIND=dp)*(nmo-homo_beta)*8.0_dp/(1024_dp**2)

      WRITE(unit_nr,'(T3,A,T67,F11.2,A3)') 'RI_INFO| Total memory for (ia|K) integrals:',&
                                           mem_for_iaK, ' MB'
      IF(my_do_gw) THEN
        mem_for_iaK=dimen_RI*REAL(nmo,KIND=dp)*gw_corr_lev_total*8.0_dp/(1024_dp**2)

        WRITE(unit_nr,'(T3,A,T67,F11.2,A3)') 'RI_INFO| Total memory for G0W0-(nm|K) integrals:',&
                                           mem_for_iaK, ' MB'
      END IF
      CALL m_flush(unit_nr)
    ENDIF

    CALL timeset(routineN//"_alloc",handle2)
    CALL mp_sync(para_env%group) ! sync to see memory output

    ! array that will store the (ia|K) integrals
    ALLOCATE(BIb_C(my_group_L_size,my_B_size,homo),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
    BIb_C=0.0_dp

    IF(do_alpha_beta) THEN
      ALLOCATE(BIb_C_beta(my_group_L_size,my_B_size_beta,homo_beta),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
      BIb_C_beta=0.0_dp
    END IF

    ! in the case of GW, we also need (nm|K) 
    IF(my_do_gw) THEN
      ALLOCATE(BIb_C_gw(my_group_L_size,my_B_all_size,gw_corr_lev_total),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
      BIb_C_gw=0.0_dp
      IF(do_alpha_beta) THEN
        ALLOCATE(BIb_C_gw_beta(my_group_L_size,my_B_all_size,gw_corr_lev_total),STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
        BIb_C_gw_beta=0.0_dp
      END IF
    END IF


    CALL timestop(handle2)

    CALL timeset(routineN//"_loop",handle2)
    i_counter=0
    DO LLL=my_group_L_start, my_group_L_end
      i_counter=i_counter+1

      ! pseudo psi_L
      CALL calculate_wavefunction(mo_coeff,1,psi_L,rho_g, atomic_kind_set, &
                                  qs_kind_set,cell,dft_control,particle_set, &
                                  pw_env_sub,basis_set_id=use_ri_aux_basis_set,&
                                  external_vector=my_Lrows(:,LLL-my_group_L_start+1),&
                                  error=error_sub)

      CALL timeset(routineN//"_pot",handle3)
      rho_r%pw%cr3d = psi_L%pw%cr3d
      CALL pw_transfer(rho_r%pw, rho_g%pw, error=error_sub)
      CALL pw_poisson_solve(poisson_env,rho_g%pw, pair_energy, pot_g%pw,error=error_sub)
      CALL pw_transfer(pot_g%pw, rho_r%pw, error=error_sub)
      CALL pw_scale(rho_r%pw,rho_r%pw%pw_grid%dvol, error=error_sub)
      CALL timestop(handle3)

      ! and finally (K|mu nu)
      CALL timeset(routineN//"_int",handle3)
      CALL cp_dbcsr_set(mat_munu%matrix,0.0_dp,error=error_sub)
      CALL integrate_v_rspace(rho_r, h=mat_munu,qs_env=qs_env,calculate_forces=.FALSE.,compute_tau=.FALSE.,gapw=.FALSE.,&
         basis_set_id=use_orb_basis_set, pw_env_external=pw_env_sub, task_list_external=task_list_sub, error=error_sub)
      CALL timestop(handle3)

      ! multiply and goooooooo ...
      CALL timeset(routineN//"_mult_o",handle3)
      CALL cp_dbcsr_multiply("N", "N", 1.0_dp, mat_munu%matrix, mo_coeff_o, &
                              0.0_dp, matrix_ia_jnu, filter_eps=eps_filter, error=error_sub) 
      IF(do_alpha_beta) THEN
        ! transform orbitals using the beta coeff matrix
        CALL cp_dbcsr_multiply("N", "N", 1.0_dp, mat_munu%matrix, mo_coeff_o_beta, &
                                0.0_dp, matrix_ia_jnu_beta, filter_eps=eps_filter, error=error_sub)
      END IF
      CALL timestop(handle3)
      CALL timeset(routineN//"_mult_v",handle3)
      CALL cp_dbcsr_multiply("T", "N", 1.0_dp, matrix_ia_jnu, mo_coeff_v, &
                              0.0_dp, matrix_ia_jb, filter_eps=eps_filter, error=error_sub)
      IF(do_alpha_beta) THEN
        ! transform orbitals using the beta coeff matrix
        CALL cp_dbcsr_multiply("T", "N", 1.0_dp, matrix_ia_jnu_beta, mo_coeff_v_beta, &
                                0.0_dp, matrix_ia_jb_beta, filter_eps=eps_filter, error=error_sub)
      END IF

      CALL timestop(handle3)

      ! now fill the matrix
      CALL timeset(routineN//"_E_Ex_1",handle3)
      CALL copy_dbcsr_to_fm(matrix_ia_jb, fm_BIb_jb, error=error_sub)

      CALL grep_my_integrals(para_env_sub,fm_BIb_jb,BIb_C(i_counter,1:my_B_size,1:homo),max_row_col_local,&
                             homo,virtual,sub_proc_map,local_col_row_info,&
                             my_B_virtual_end,my_B_virtual_start,my_B_size,&
                             error_sub)
      IF(do_alpha_beta) THEN
        CALL copy_dbcsr_to_fm(matrix_ia_jb_beta, fm_BIb_jb_beta, error=error_sub)
        CALL grep_my_integrals(para_env_sub,fm_BIb_jb_beta,&
                               BIb_C_beta(i_counter,1:my_B_size_beta,1:homo_beta),max_row_col_local_beta,&
                               homo_beta,virtual_beta,sub_proc_map,local_col_row_info_beta,&
                               my_B_virtual_end_beta,my_B_virtual_start_beta,my_B_size_beta,&
                               error_sub)
      END IF

      CALL timestop(handle3)


      IF(my_do_gw) THEN
        ! transform (K|mu nu) to (K|nm), n corresponds to corrected GW levels, m is in nmo
        CALL timeset(routineN//"_mult_gw",handle3)
        CALL cp_dbcsr_multiply("N", "N", 1.0_dp, mat_munu%matrix, mo_coeff_gw, &
                                0.0_dp, matrix_in_jnu, filter_eps=eps_filter, error=error_sub)
        CALL cp_dbcsr_multiply("T", "N", 1.0_dp, matrix_in_jnu, mo_coeff_all, &
                                0.0_dp,matrix_in_jm,filter_eps=eps_filter, error=error_sub)

        CALL timestop(handle3)

        CALL timeset(routineN//"_E_Ex_2_gw",handle3)

        CALL copy_dbcsr_to_fm(matrix_in_jm, fm_BIb_gw, error=error_sub)
        CALL grep_my_integrals(para_env_sub,fm_BIb_gw,BIb_C_gw(i_counter,1:my_B_all_size,1:gw_corr_lev_total),&
                               max_row_col_local_gw,&
                               gw_corr_lev_total,nmo,sub_proc_map,local_col_row_info_gw,&
                               my_B_all_end,my_B_all_start,my_B_all_size,&
                               error_sub)

        ! the same for beta
        IF(do_alpha_beta) THEN
          CALL cp_dbcsr_multiply("N", "N", 1.0_dp, mat_munu%matrix, mo_coeff_gw_beta, &
                                  0.0_dp, matrix_in_jnu_beta, filter_eps=eps_filter, error=error_sub)
          CALL cp_dbcsr_multiply("T", "N", 1.0_dp, matrix_in_jnu_beta, mo_coeff_all_beta, &
                                  0.0_dp,matrix_in_jm_beta,filter_eps=eps_filter, error=error_sub)
  
          CALL timestop(handle3)
  
          CALL timeset(routineN//"_E_Ex_2_gw",handle3)
  
          CALL copy_dbcsr_to_fm(matrix_in_jm_beta, fm_BIb_gw_beta, error=error_sub)
          CALL grep_my_integrals(para_env_sub,fm_BIb_gw_beta,BIb_C_gw_beta(i_counter,1:my_B_all_size,1:gw_corr_lev_total),&
                                 max_row_col_local_gw,&
                                 gw_corr_lev_total,nmo,sub_proc_map,local_col_row_info_gw,&
                                 my_B_all_end,my_B_all_start,my_B_all_size,&
                                 error_sub)
        END IF

        CALL timestop(handle3)
  
      END IF

    END DO


    ! use truncated operator for diagonal elements B_nn^P in case of periodic calculations 
    IF(my_do_gw .AND. do_truncation_gw) THEN

      ! truncation
      IF(my_group_L_start<=my_group_L_end)  THEN
        grid => poisson_env%green_fft%influence_fn%pw_grid
        gf   => poisson_env%green_fft%influence_fn
        DO ig = grid%first_gne0, grid%ngpts_cut_local
          g2 = grid%gsq(ig)
          gg = SQRT(g2)
          g3d = fourpi/g2
          gf%cc(ig) = g3d*(1.0_dp-COS(rlength*gg))
        END DO
        IF(grid%have_g0) &
          gf%cc(1)=0.5_dp*fourpi*rlength*rlength
      END IF


      i_counter=0
      DO LLL=my_group_L_start, my_group_L_end
        i_counter=i_counter+1
  
        ! pseudo psi_L
        CALL calculate_wavefunction(mo_coeff,1,psi_L,rho_g, atomic_kind_set, &
                                    qs_kind_set,cell,dft_control,particle_set, &
                                    pw_env_sub,basis_set_id=use_ri_aux_basis_set,&
                                    external_vector=my_Lrows(:,LLL-my_group_L_start+1),&
                                    error=error_sub)
  
        ALLOCATE(BIb_C_gw_trunc(my_group_L_size,my_B_all_size,gw_corr_lev_total),STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
        BIb_C_gw_trunc=0.0_dp

        CALL timeset(routineN//"_pot_trunc",handle3)
        rho_r%pw%cr3d = psi_L%pw%cr3d
        CALL pw_transfer(rho_r%pw, rho_g%pw, error=error_sub)

        CALL pw_poisson_solve(poisson_env,rho_g%pw, pair_energy, pot_g%pw,error=error_sub)
        CALL pw_transfer(pot_g%pw, rho_r%pw, error=error_sub)
        CALL pw_scale(rho_r%pw,rho_r%pw%pw_grid%dvol, error=error_sub)
        CALL timestop(handle3)
 
        ! and finally (K|mu nu)
        CALL timeset(routineN//"_int_trunc",handle3)
        CALL cp_dbcsr_set(mat_munu%matrix,0.0_dp,error=error_sub)
        CALL integrate_v_rspace(rho_r, h=mat_munu,qs_env=qs_env,calculate_forces=.FALSE.,compute_tau=.FALSE.,gapw=.FALSE.,&
           basis_set_id=use_orb_basis_set, pw_env_external=pw_env_sub, task_list_external=task_list_sub, error=error_sub)
        CALL timestop(handle3)
 
        CALL timeset(routineN//"_mult_gw_trunc",handle3)
        CALL cp_dbcsr_multiply("N", "N", 1.0_dp, mat_munu%matrix, mo_coeff_gw, &
                                0.0_dp, matrix_in_jnu, filter_eps=eps_filter, error=error_sub)
        CALL cp_dbcsr_multiply("T", "N", 1.0_dp, matrix_in_jnu, mo_coeff_all, &
                                0.0_dp,matrix_in_jm,filter_eps=eps_filter, error=error_sub)

        CALL timestop(handle3)

        CALL timeset(routineN//"_E_Ex_2_trunc",handle3)

        CALL copy_dbcsr_to_fm(matrix_in_jm, fm_BIb_gw, error=error_sub)
        CALL grep_my_integrals(para_env_sub,fm_BIb_gw,BIb_C_gw_trunc(i_counter,1:my_B_all_size,1:gw_corr_lev_total),&
                               max_row_col_local_gw,&
                               gw_corr_lev_total,nmo,sub_proc_map,local_col_row_info_gw,&
                               my_B_all_end,my_B_all_start,my_B_all_size,&
                               error_sub)

        DO iiB=1,my_B_all_size
          DO jjB=1,gw_corr_lev_total
!          DO jjB=homo-gw_corr_lev_occ+1,homo+gw_corr_lev_virt
            ! only use the truncation for the diagonal elements B^nn_P
            IF(iiB + my_B_all_start - 1 == homo - gw_corr_lev_occ + jjB) THEN
              BIb_C_gw(i_counter,iiB,jjB)=BIb_C_gw_trunc(i_counter,iiB,jjB)
            END IF
          END DO
        END DO

        DEALLOCATE(BIb_C_gw_trunc)

        CALL timestop(handle3)

      END DO

      ! set this flag so that the truncation of above will be overwritten the next time the poisson_env is used
      IF(my_group_L_start<=my_group_L_end)  &
        poisson_env%rebuild=.TRUE. 

    END IF ! use truncated operator for GW
    CALL timestop(handle2)

    DEALLOCATE(my_Lrows)

    CALL cp_fm_release(fm_BIb_jb, error=error_sub)
    DEALLOCATE(local_col_row_info)

    CALL cp_dbcsr_release(matrix_ia_jnu,error=error_sub)
    CALL cp_dbcsr_release(matrix_ia_jb,error=error_sub)
    IF(do_alpha_beta) THEN
      CALL cp_dbcsr_release(matrix_ia_jnu_beta,error=error_sub)
      CALL cp_dbcsr_release(matrix_ia_jb_beta,error=error_sub)
      CALL cp_fm_release(fm_BIb_jb_beta, error=error_sub)
      DEALLOCATE(local_col_row_info_beta)
    END IF

    IF(my_do_gw) THEN
      CALL cp_dbcsr_release(matrix_in_jnu,error=error_sub)
      CALL cp_dbcsr_release(matrix_in_jm,error=error_sub)
      CALL cp_fm_release(fm_BIb_gw, error=error_sub)
      DEALLOCATE(local_col_row_info_gw)
      IF(do_alpha_beta) THEN
        CALL cp_dbcsr_release(matrix_in_jnu_beta,error=error_sub)
        CALL cp_dbcsr_release(matrix_in_jm_beta,error=error_sub)
        CALL cp_fm_release(fm_BIb_gw_beta, error=error_sub)
      END IF
    END IF

    CALL pw_pool_give_back_pw(auxbas_pw_pool,psi_L%pw,error=error_sub)

    DEALLOCATE(sub_proc_map)

    CALL timestop(handle)

  END SUBROUTINE mp2_ri_gpw_compute_in


! *****************************************************************************
!> \brief Encapsulate building of intermediate matrices matrix_ia_jnu(_beta
!>         matrix_ia_jb(_beta),fm_BIb_jb(_beta),matrix_in_jnu(for G0W0) and
!>         fm_BIb_all(for G0W0)
!> \param matrix_ia_jnu ...
!> \param matrix_ia_jb ...
!> \param mo_coeff_templ ...
!> \param size_1 ...
!> \param size_2 ...
!> \param fm_BIb_jb ...
!> \param matrix_name_2 ...
!> \param max_row_col_local ...
!> \param nfullrows_total ...
!> \param nfullcols_total ...
!> \param nrow_local ...
!> \param ncol_local ...
!> \param blacs_env_sub ...
!> \param para_env_sub ...
!> \param local_col_row_info ...
!> \param error_sub ...
!> \param error ...
!> \author Jan Wilhelm
! *****************************************************************************
  SUBROUTINE create_intermediate_matrices(matrix_ia_jnu,matrix_ia_jb,mo_coeff_templ,size_1,size_2,&
                                          fm_BIb_jb,matrix_name_2,max_row_col_local,&
                                          nfullrows_total,nfullcols_total,nrow_local, &
                                          ncol_local,blacs_env_sub,&
                                          para_env_sub,local_col_row_info,error_sub,error)

    TYPE(cp_dbcsr_type)                      :: matrix_ia_jnu, matrix_ia_jb, &
                                                mo_coeff_templ
    INTEGER                                  :: size_1, size_2
    TYPE(cp_fm_type), POINTER                :: fm_BIb_jb
    CHARACTER(LEN=*)                         :: matrix_name_2
    INTEGER                                  :: max_row_col_local, &
                                                nfullrows_total, &
                                                nfullcols_total, nrow_local, &
                                                ncol_local
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env_sub
    TYPE(cp_para_env_type), POINTER          :: para_env_sub
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: local_col_row_info
    TYPE(cp_error_type), INTENT(inout)       :: error_sub, error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'create_intermediate_matrices', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, stat
    INTEGER, DIMENSION(:), POINTER           :: col_indices, row_indices
    LOGICAL                                  :: failure
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct

    CALL timeset(routineN,handle)
    failure=.FALSE.

    ! initialize and create the matrix (K|jnu)
    CALL cp_dbcsr_init(matrix_ia_jnu,error=error_sub)
    CALL cp_dbcsr_create(matrix_ia_jnu,template=mo_coeff_templ,error=error_sub)

    ! Allocate Sparse matrices: (K|jb)
    CALL cp_dbcsr_init(matrix_ia_jb,error=error_sub)
    CALL cp_dbcsr_m_by_n_from_template(matrix_ia_jb,template=mo_coeff_templ,m=size_2,n=size_1,&
                                       sym=dbcsr_type_no_symmetry,error=error)

    ! set all to zero in such a way that the memory is actually allocated
    CALL cp_dbcsr_set(matrix_ia_jnu,0.0_dp,error=error_sub)
    CALL cp_dbcsr_set(matrix_ia_jb,0.0_dp,error=error_sub)

    ! create the analogous of matrix_ia_jb in fm type
    NULLIFY(fm_BIb_jb)
    NULLIFY(fm_struct)
    CALL cp_dbcsr_get_info(matrix_ia_jb,nfullrows_total=nfullrows_total,nfullcols_total=nfullcols_total)
    CALL cp_fm_struct_create(fm_struct,context=blacs_env_sub,nrow_global=nfullrows_total,&
                             ncol_global=nfullcols_total,para_env=para_env_sub,error=error_sub)
    CALL cp_fm_create(fm_BIb_jb,fm_struct,name=matrix_name_2,error=error_sub)

    CALL copy_dbcsr_to_fm(matrix_ia_jb, fm_BIb_jb, error=error_sub)
    CALL cp_fm_struct_release(fm_struct,error=error_sub)

    CALL cp_fm_get_info(matrix=fm_BIb_jb,&
                        nrow_local=nrow_local,&
                        ncol_local=ncol_local,&
                        row_indices=row_indices,&
                        col_indices=col_indices,&
                        error=error_sub)

    max_row_col_local=MAX(nrow_local,ncol_local)
    CALL mp_max(max_row_col_local,para_env_sub%group)

    ALLOCATE(local_col_row_info(0:max_row_col_local,2),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
    local_col_row_info=0
    ! 0,1 nrows
    local_col_row_info(0,1)=nrow_local
    local_col_row_info(1:nrow_local,1)=row_indices(1:nrow_local)
    ! 0,2 ncols
    local_col_row_info(0,2)=ncol_local
    local_col_row_info(1:ncol_local,2)=col_indices(1:ncol_local)


    CALL timestop(handle)

  END SUBROUTINE create_intermediate_matrices 


! *****************************************************************************
!> \brief Encapsulate building of parallelization arrays
!> \param para_env_sub ...
!> \param ends_B ...
!> \param sizes_B ...
!> \param starts_B ...
!> \param my_B_start ...
!> \param my_B_end ...
!> \param my_B_size ...
!> \param size_1 ...
!> \param error ...
!> \author Jan Wilhelm
! *****************************************************************************
  SUBROUTINE create_parallelization_arrays(para_env_sub,ends_B,sizes_B,starts_B,&
                                           my_B_start,my_B_end,my_B_size,size_1,error)

    TYPE(cp_para_env_type), POINTER          :: para_env_sub
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: ends_B, sizes_B, starts_B
    INTEGER                                  :: my_B_start, my_B_end, &
                                                my_B_size, size_1
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'create_parallelization_arrays', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iproc, itmp(2), stat
    LOGICAL                                  :: failure

    CALL timeset(routineN,handle)
    failure=.FALSE.

    ! divide the b states in the sub_group in such a way to create
    ! b_start and b_end for each proc inside the sub_group
    ALLOCATE(sizes_B(0:para_env_sub%num_pe-1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    sizes_B=0
    ALLOCATE(starts_B(0:para_env_sub%num_pe-1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    starts_B=0
    ALLOCATE(ends_B(0:para_env_sub%num_pe-1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ends_B=0

    DO iproc=0, para_env_sub%num_pe-1
      itmp=get_limit(size_1,para_env_sub%num_pe,iproc)
      starts_B(iproc)=itmp(1)
      ends_B(iproc)=itmp(2)
      sizes_B(iproc)=itmp(2)-itmp(1)+1
    END DO

    my_B_start=starts_B(para_env_sub%mepos)
    my_B_end=ends_B(para_env_sub%mepos)
    my_B_size=sizes_B(para_env_sub%mepos)


    CALL timestop(handle)

  END SUBROUTINE create_parallelization_arrays


! *****************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param para_env ...
!> \param para_env_sub ...
!> \param para_env_L ...
!> \param mp2_memory ...
!> \param eps_filter ...
!> \param fm_matrix_L ...
!> \param ngroup ...
!> \param color_sub ...
!> \param dimen_RI ...
!> \param mo_coeff ...
!> \param dft_control ...
!> \param psi_L ...
!> \param rho_r ...
!> \param rho_g ...
!> \param pot_g ...
!> \param pw_env_sub ...
!> \param poisson_env ...
!> \param my_group_L_size ...
!> \param my_group_L_start ...
!> \param my_group_L_end ...
!> \param sab_orb_sub ...
!> \param sizes_array ...
!> \param starts_array ...
!> \param ends_array ...
!> \param calc_PQ_cond_num ...
!> \param cond_num ...
!> \param num_small_eigen ...
!> \param auxbas_pw_pool ...
!> \param error_sub ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE calculate_Lmin1(qs_env,para_env,para_env_sub,para_env_L,mp2_memory,&
                             eps_filter,fm_matrix_L,ngroup,color_sub,dimen_RI,&
                             mo_coeff,dft_control,psi_L,rho_r,rho_g,pot_g,pw_env_sub,poisson_env,&
                             my_group_L_size,my_group_L_start,my_group_L_end,sab_orb_sub,&
                             sizes_array,starts_array,ends_array,calc_PQ_cond_num,cond_num,&
                             num_small_eigen,auxbas_pw_pool,error_sub,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_para_env_type), POINTER          :: para_env, para_env_sub, &
                                                para_env_L
    REAL(KIND=dp)                            :: mp2_memory, eps_filter
    TYPE(cp_fm_type), POINTER                :: fm_matrix_L
    INTEGER                                  :: ngroup, color_sub, dimen_RI
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(pw_p_type)                          :: psi_L, rho_r, rho_g, pot_g
    TYPE(pw_env_type), POINTER               :: pw_env_sub
    TYPE(pw_poisson_type), POINTER           :: poisson_env
    INTEGER                                  :: my_group_L_size, &
                                                my_group_L_start, &
                                                my_group_L_end
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb_sub
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: sizes_array, starts_array, &
                                                ends_array
    LOGICAL                                  :: calc_PQ_cond_num
    REAL(KIND=dp)                            :: cond_num
    INTEGER                                  :: num_small_eigen
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(cp_error_type), INTENT(inout)       :: error_sub, error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_Lmin1', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: debug = .FALSE.

    INTEGER :: best_group_size, color_exc, color_L, comm_exc, comm_exchange, &
      comm_L, dir, group_size, handle, handle2, handle3, i, i_bf, i_counter, &
      i_counter_diff_bf, i_counter_diff_bf_at_any_atom, i_counter_diff_sets, &
      i_global, iatom, igrid_level, igroup, iiB, ikind, info_chol, inode, &
      ipgf, iproc, iset, istat, itmp(2), j, j_global, j_ref_atom, j_ref_bf, &
      j_ref_kind, j_ref_set, jatom, jjB, jkind, jset, lb(3), LLL, &
      location(3), m, my_group_L_q_end, my_group_L_q_size, &
      my_group_L_q_start, my_group_pot_end, my_group_pot_size, &
      my_group_pot_start, na1, na2, natom, ncoa, ncoa_a, ncoa_b, ncol_local
    INTEGER :: nkind, nrow_local, nseta, nsetb, offset, offsetb, &
      proc_receive, proc_receive_static, proc_send, proc_send_static, &
      proc_shift, rec_L_end, rec_L_size, rec_L_start, sgfa, sgfb, stat, &
      strat_group_size, sub_sub_color, tp(3), ub(3)
    INTEGER, ALLOCATABLE, DIMENSION(:) :: diffbf_atom_arbitr_atom, &
      diffbf_noinbas, diffbfatanyatom_noinbas, ends_array_pot, kind_of, &
      proc_map, sizes_array_pot, starts_array_pot, sub_ends_array, &
      sub_sizes_array, sub_starts_array
    INTEGER, ALLOCATABLE, DIMENSION(:, :) :: atom_set_afo_diffbf, &
      atom_set_afo_diffbf_atanyatom, diffsets_atom_set_bf, ref_matrix
    INTEGER, ALLOCATABLE, DIMENSION(:, :, :) :: naolimits, ref_at_se_sgf_ofs
    INTEGER, DIMENSION(:), POINTER           :: col_indices, la_max, la_min, &
                                                lb_max, lb_min, npgfa, npgfb, &
                                                nsgfa, nsgfb, row_indices
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb
    LOGICAL                                  :: check_i, check_j, do_mult_2c, &
                                                failure, map_it_here
    REAL(KIND=dp)                            :: min_mem_for_QK, pair_energy, &
                                                rab2, rac2, scaling, zet_tmp
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: egen_L, f, wf_vector
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: I_ab, I_tmp, L_external_col, &
                                                L_local_col, L_local_col_q, &
                                                vac
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: v
    REAL(KIND=dp), DIMENSION(1)              :: exp_q
    REAL(KIND=dp), DIMENSION(3)              :: ra, rab, rac
    REAL(KIND=dp), DIMENSION(:), POINTER     :: set_radius_a, set_radius_b
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: I_tmp2, rpgfa, rpgfb, sphi_a, &
                                                sphi_b, zeta, zetb
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env_L
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct
    TYPE(cp_fm_type), POINTER                :: fm_matrix_L_diag
    TYPE(cp_para_env_type), POINTER          :: para_env_exc, &
                                                para_env_exchange
    TYPE(gto_basis_set_type), POINTER        :: basis_set_a, basis_set_b
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: nl_iterator
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_p_type), ALLOCATABLE, &
      DIMENSION(:)                           :: rho_r_array
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(realspace_grid_desc_p_type), &
      DIMENSION(:), POINTER                  :: rs_descs
    TYPE(realspace_grid_p_type), &
      DIMENSION(:), POINTER                  :: rs_v

    CALL timeset(routineN,handle)
    failure=.FALSE.


    ! get stuff
    CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, qs_kind_set=qs_kind_set,&
                    cell=cell, molecule_set=molecule_set,particle_set=particle_set,&
                    molecule_kind_set=molecule_kind_set,error=error)

    ! blacs_env => qs_env%blacs_env

    nkind = SIZE(qs_kind_set)
    natom = SIZE(particle_set)

    ALLOCATE (kind_of(natom),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL get_atomic_kind_set(atomic_kind_set,kind_of=kind_of)

    DO ikind=1, nkind
      CPPrecondition(ASSOCIATED(qs_kind_set(ikind)%ri_aux_basis_set),cp_failure_level,routineP,error,failure)
      ! qs_kind_set(ikind)%ri_aux_basis_set%norm_type=2
      ! CALL init_orb_basis_set(qs_kind_set(ikind)%ri_aux_basis_set,error)
    END DO

    dimen_RI=0
    DO iatom=1, natom
      ikind=kind_of(iatom)
      dimen_RI=dimen_RI+qs_kind_set(ikind)%ri_aux_basis_set%nsgf
    END DO

    ! calculate wich rows of L^{-1} to have
    ALLOCATE(sizes_array(0:ngroup-1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    sizes_array=0
    ALLOCATE(starts_array(0:ngroup-1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    starts_array=0
    ALLOCATE(ends_array(0:ngroup-1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ends_array=0

    DO igroup=0,ngroup-1
       itmp=get_limit(dimen_RI,ngroup,igroup)
       starts_array(igroup)=itmp(1)
       ends_array(igroup)=itmp(2)
       sizes_array(igroup)=itmp(2)-itmp(1)+1
    ENDDO

    my_group_L_size=sizes_array(color_sub)
    my_group_L_start=starts_array(color_sub)
    my_group_L_end=ends_array(color_sub)

    CALL timeset(routineN//"_loop_lm",handle2)

    do_mult_2c =  qs_env%mp2_env%mp2_gpw%do_mult_2c

    IF(do_mult_2c) THEN

      ! write arrays with number(=label) of every basis function, its quantum number (l,m)
      ! and the reference to the first basis function with (l,m) 
      ALLOCATE(ref_matrix(dimen_RI,10),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

      ! naolimits(iatom,iset,1/2)=begin/end reference entries in L_local_q 
      ALLOCATE(naolimits(natom,100,2),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
 
      ! array with reference (atom,set) for arbitrary (atom,set);
      ALLOCATE(ref_at_se_sgf_ofs(natom,100,4),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
 
 
      i_counter=0
      i_counter_diff_bf=0
      i_counter_diff_sets=0
      i_counter_diff_bf_at_any_atom=0
 
      ! 1) build up reference matrices to get the reference q_alpha for every basis function alpha
 
      offset = 0
 
      DO iatom=1, natom
        ikind=kind_of(iatom)
        basis_set_a => qs_kind_set(ikind)%ri_aux_basis_set
  
        la_max       =>  basis_set_a%lmax
        la_min       =>  basis_set_a%lmin
        npgfa        =>  basis_set_a%npgf
        nseta        =   basis_set_a%nset
        nsgfa        =>  basis_set_a%nsgf_set
        zeta         =>  basis_set_a%zet
        first_sgfa   =>  basis_set_a%first_sgf
 
 
        DO iset=1, nseta
 
          DO ipgf=1, npgfa(iset)
             
            !pseudo-m quantum number runs from 1 to 2l+1
            DO m=1,nsgfa(iset)
              i_counter=i_counter+1
              !first entry: l; only in case of la_max(iset)==la_min(iset) the shortcut will work, see IF below
              ref_matrix(i_counter,1) = la_max(iset)   
              !second entry: pseudo-m quantum number runs from 1 to 2l+1
              ref_matrix(i_counter,2) = m
              ref_matrix(i_counter,3) = iatom
              ref_matrix(i_counter,4) = iset
              ref_matrix(i_counter,5) = ipgf
              ref_matrix(i_counter,8) = first_sgfa(1,iset)
 
 
              !Check, whether (l,m) combination already exists
              IF(i_counter>=2) THEN
                DO i_bf=1, i_counter-1
                  ! 6th row of ref_matrix contains number of first basis function with identical (l,m)
                  ! for set with la_min=/la_max at the same atom 
                  IF( ref_matrix(i_bf,1)==la_max(iset) .AND. ref_matrix(i_bf,2)==m .AND. &
                      la_max(iset)==la_min(iset) .AND. ref_matrix(i_bf,3)==iatom) THEN
                    ref_matrix(i_counter,6) = i_bf
    
                    ! we are done with checking the basis functions in DO when a (l,m) match was found.
                    EXIT 
  
                  END IF
                    ! 0 if no basis function has this (l,m) combination
                    ref_matrix(i_counter,6) = 0
  
                END DO !loop already existing basis functions
  
                DO i_bf=1, i_counter-1
                  ! 6th row of ref_matrix contains number of first basis function with identical (l,m)
                  ! for set with la_min=/la_max AT ANY ATOM (!!DIFFERENCE!!) 
                  IF( ref_matrix(i_bf,1)==la_max(iset) .AND. ref_matrix(i_bf,2)==m .AND. &
                      la_max(iset)==la_min(iset)) THEN
                    ref_matrix(i_counter,9) = i_bf
  
                    ! we are done with checking the basis functions in DO when a (l,m) match was found.
                    EXIT
  
                  END IF
                    ! 0 if no basis function has this (l,m) combination
                    ref_matrix(i_counter,9) = 0
  
                END DO !loop already existing basis functions
  
  
  
              ELSE 
                ref_matrix(1,6) = 0
                ref_matrix(1,9) = 0
  
              END IF 
              ! for new (l,m) combination increase i_counter_diff_bf by 1 and set reference to itself 
              IF(ref_matrix(i_counter,6) == 0) THEN
                i_counter_diff_bf = i_counter_diff_bf + 1
                ref_matrix(i_counter,7) = i_counter_diff_bf
              ! for already existing (l,m) combination, 7th entry refers to number/label of first (l,m) combination
              ELSE
                ref_matrix(i_counter,7) = &
                         ref_matrix(ref_matrix(i_counter,6),7)
              END IF
  
              ! for new (l,m) at any atom combination increase i_counter_diff_bf_at_any_atom by 1
              !  and set reference to itself 
              IF(ref_matrix(i_counter,9) == 0) THEN
                i_counter_diff_bf_at_any_atom = i_counter_diff_bf_at_any_atom + 1
                ref_matrix(i_counter,10) = i_counter_diff_bf_at_any_atom
              ! for already existing (l,m) combination, 7th entry refers to number/label of first (l,m) combination
              ELSE
                ref_matrix(i_counter,10) = &
                         ref_matrix(ref_matrix(i_counter,9),10)
              END IF
  
  
            END DO !m
  
          END DO !prim. gf.
  
          IF(ref_matrix(i_counter,6) == 0) THEN
            i_counter_diff_sets = i_counter_diff_sets + 1
            ref_at_se_sgf_ofs(iatom,iset,1) = iatom
            ref_at_se_sgf_ofs(iatom,iset,2) = iset
            ref_at_se_sgf_ofs(iatom,iset,3) = first_sgfa(1,iset)
  
          ELSE
            ref_at_se_sgf_ofs(iatom,iset,1) = ref_matrix(ref_matrix(i_counter,6),3)
            ref_at_se_sgf_ofs(iatom,iset,2) = ref_matrix(ref_matrix(i_counter,6),4)
            ref_at_se_sgf_ofs(iatom,iset,3) = ref_matrix(ref_matrix(i_counter,6),8)
          END IF
  
          naolimits(iatom,iset,1) = ref_matrix(i_counter,7)-nsgfa(iset)+1
          naolimits(iatom,iset,2) = ref_matrix(i_counter,7)
  
          ref_at_se_sgf_ofs(iatom,iset,4) = offset
          offset=offset+nsgfa(iset)
  
  
       END DO !sets
     
      END DO !atoms
   
      ! array with number of q_a basis function at an atom
      ALLOCATE(diffbf_noinbas(i_counter_diff_bf),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
   
      ! array with number of a basis function with first unique (l,m) combination
      ALLOCATE(diffbfatanyatom_noinbas(i_counter_diff_bf_at_any_atom),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
   
      ! array with number of q_a (at arbitrary atom) with q_a at specific atom
      ALLOCATE(diffbf_atom_arbitr_atom(i_counter_diff_bf),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
   
      ! array with atom and set number for every q_a bf
      ALLOCATE(atom_set_afo_diffbf(i_counter_diff_bf,2),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
   
      ! array with atom and set number for every q_a bf
      ALLOCATE(atom_set_afo_diffbf_atanyatom(i_counter_diff_bf,2),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
   
      j=1
      DO i=1,dimen_RI
        IF(ref_matrix(i,9)==0) THEN
          diffbfatanyatom_noinbas(j)=i
          atom_set_afo_diffbf_atanyatom(j,1)=ref_matrix(i,3)
          atom_set_afo_diffbf_atanyatom(j,2)=ref_matrix(i,4)
          j= j+1
        END  IF
      END DO 
   
     j=1
     DO i=1,dimen_RI
       IF(ref_matrix(i,6)==0) THEN
       diffbf_noinbas(j)=i
       atom_set_afo_diffbf(j,1)=ref_matrix(i,3)
       atom_set_afo_diffbf(j,2)=ref_matrix(i,4)
       diffbf_atom_arbitr_atom(j)=ref_matrix(i,10)
       j=j+1
       END IF
     END DO
  
     ! array with number of first set with quantum number l
     ALLOCATE(diffsets_atom_set_bf(i_counter_diff_sets,3),STAT=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  
     ! first set with new quantum number l is on atom 1, set 1
     diffsets_atom_set_bf(1,1)=1
     diffsets_atom_set_bf(1,2)=1
     diffsets_atom_set_bf(1,3)=1
  
     j=2
     DO i=2,dimen_RI
       IF(ref_matrix(i,6)==0 .AND. ( diffsets_atom_set_bf(j-1,1)/=ref_matrix(i,3) .OR. &
          diffsets_atom_set_bf(j-1,2)/=ref_matrix(i,4)) ) THEN
          diffsets_atom_set_bf(j,1)=ref_matrix(i,3)
          diffsets_atom_set_bf(j,2)=ref_matrix(i,4)
          diffsets_atom_set_bf(j,3)=ref_matrix(i,7)
          j=j+1
       END IF
     END DO
  
     IF(debug) THEN  
       PRINT *, 'i_counter_diff_bf_at_any_atom = ', i_counter_diff_bf_at_any_atom
   
       DO i=1,dimen_RI
         WRITE(*,*)(ref_matrix(i,j),j=1,10)
       END DO
   
       PRINT *, 'diffbf_noinbas'
       DO i=1,i_counter_diff_bf
         WRITE(*,*)(diffbf_noinbas(i))
       END DO
   
       PRINT *, 'diffbfatanyatom_noinbas'
       DO i=1,i_counter_diff_bf_at_any_atom
         WRITE(*,*)(diffbfatanyatom_noinbas(i))
       END DO
   
       PRINT *, 'diffbf_atom_arbitr_atom'
       DO i=1,i_counter_diff_bf
         WRITE(*,*)(diffbf_atom_arbitr_atom(i))
       END DO
    
       PRINT *, 'atom_set_afo_diffbf'
       DO i=1,i_counter_diff_bf
         PRINT *, atom_set_afo_diffbf(i,1), atom_set_afo_diffbf(i,2)
       END DO
   
       PRINT *, 'atom_set_afo_diffbf_atanyatom'
       DO i=1,i_counter_diff_bf_at_any_atom
         PRINT *, atom_set_afo_diffbf_atanyatom(i,1), atom_set_afo_diffbf_atanyatom(i,2)
       END DO
   
       PRINT *, 'diffsets_atom_set_bf'
       DO i=1,i_counter_diff_sets
         WRITE(*,*)(diffsets_atom_set_bf(i,j),j=1,3)
       END DO
   
       PRINT *, 'naolimits'
       DO iatom=1,natom
          ikind=kind_of(iatom)
          basis_set_a => qs_kind_set(ikind)%ri_aux_basis_set
          nseta        =   basis_set_a%nset
    
          DO iset=1, nseta
            WRITE(*,*)(naolimits(iatom,iset,j),j=1,2)
          END DO
       END DO
   
       PRINT *, 'ref_at_se_sgf_ofs'
       DO iatom=1,natom
          ikind=kind_of(iatom)
          basis_set_a => qs_kind_set(ikind)%ri_aux_basis_set
          nseta        =   basis_set_a%nset
    
          DO iset=1, nseta
            WRITE(*,*)(ref_at_se_sgf_ofs(iatom,iset,j),j=1,4)
          END DO
       END DO
   
       PRINT *, 'i_counter_diff_bf(end) = ', i_counter_diff_bf
       PRINT *, 'i_counter_diff_sets(end) = ', i_counter_diff_sets
       PRINT *, 'i_counter(end) = ', i_counter
       PRINT *, 'dimen_RI = ', dimen_RI
       PRINT *, 'my_group_L_end = ' ,my_group_L_end
       PRINT *, 'my_group_L_start = ' ,my_group_L_start 

     END IF !debug
 
     ! For parallelization: Check, which elements of L_local_col_q    
     ! are needed for setting up L_local_q
  
     my_group_L_q_start = naolimits(ref_matrix(my_group_L_start,3),1,1)
     IF(ref_matrix(my_group_L_end,3)/=natom) THEN
       my_group_L_q_end = naolimits(ref_matrix(my_group_L_end,3)+1,1,1)-1
     ELSE
       nkind=kind_of(natom)
       basis_set_a => qs_kind_set(nkind)%ri_aux_basis_set
       nseta       =  basis_set_a%nset
       my_group_L_q_end = naolimits(ref_matrix(my_group_L_end,3),nseta,2)
     END IF
     my_group_L_q_size = my_group_L_q_end - my_group_L_q_start + 1
  
     ! calculate which potentials to calculate
     ALLOCATE(sizes_array_pot(0:ngroup-1),STAT=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     sizes_array_pot=0
     ALLOCATE(starts_array_pot(0:ngroup-1),STAT=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     starts_array_pot=0
     ALLOCATE(ends_array_pot(0:ngroup-1),STAT=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     ends_array_pot=0
  
     DO igroup=0,ngroup-1
        itmp=get_limit(i_counter_diff_bf_at_any_atom,ngroup,igroup)
        starts_array_pot(igroup)=itmp(1)
        ends_array_pot(igroup)=itmp(2)
        sizes_array_pot(igroup)=itmp(2)-itmp(1)+1
     ENDDO
  
     my_group_pot_size=sizes_array_pot(color_sub)
     my_group_pot_start=starts_array_pot(color_sub)
     my_group_pot_end=ends_array_pot(color_sub)
  
  
     DEALLOCATE(ref_matrix)
  
     ! 2.0) calculate potential (q_alpha|_p for first (l,m) combination at an atom
     i_counter = 0
     exp_q = 0.5_dp
  
     ALLOCATE(rho_r_array(i_counter_diff_bf_at_any_atom),STAT=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
 
     ALLOCATE(wf_vector(dimen_RI),STAT=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
 
  
     DO j=my_group_pot_start,my_group_pot_end
        wf_vector=0.0_dp
        ! only take one basis function unique l,m-combination at any atom
        wf_vector(diffbfatanyatom_noinbas(j))=1.0_dp
  
        ! set exponent of aux basis set temporarily to exp_q since reference wavefunction
        ! q_a/b should have exponent exp_q. Then integration (q_a|q_b) is numerically
        ! feasible in reciprocal space      
        jatom = atom_set_afo_diffbf_atanyatom(j,1)
        jset  = atom_set_afo_diffbf_atanyatom(j,2)
        jkind = kind_of(jatom)
        basis_set_b  =>  qs_kind_set(jkind)%ri_aux_basis_set
        zetb         =>  basis_set_b%zet
        lb_max       =>  basis_set_b%lmax
  
        zet_tmp = zetb(1,jset)
        zetb(1,jset) = exp_q(1)
  
        ! pseudo psi_L
        CALL calculate_wavefunction(mo_coeff,1,psi_L,rho_g, atomic_kind_set,&
                                    qs_kind_set,cell,dft_control,particle_set, &
                                    pw_env_sub,basis_set_id=use_ri_aux_basis_set,&
                                    external_vector=wf_vector,&
                                    error=error_sub)
  
        ! set back the exponent
        zetb(1,jset) = zet_tmp
  
        rho_r%pw%cr3d = psi_L%pw%cr3d
        CALL pw_transfer(rho_r%pw, rho_g%pw, error=error_sub)
        CALL pw_poisson_solve(poisson_env,rho_g%pw, pair_energy, pot_g%pw,error=error_sub)
        CALL pw_transfer(pot_g%pw, rho_r%pw, error=error_sub)
        CALL pw_scale(rho_r%pw,rho_r%pw%pw_grid%dvol, error=error_sub)
  
        ! Writing array of potential copied from xc/xc_rho_set_types.F
  
        NULLIFY(rho_r_array(j)%pw)
  
        rho_r_array(j)%pw => rho_r%pw
  
        CALL pw_pool_create_pw(auxbas_pw_pool,rho_r_array(j)%pw,&
             use_data=REALDATA3D, in_space=REALSPACE, &
             error=error)
        CALL pw_copy(rho_r%pw,rho_r_array(j)%pw,error=error)
  
     END DO
  
     ! Fill pw_grid and zero-components cr3d, aso. into rho_r_array for parallization
     ! If no potential was calculated because of my_group_pot_start>my_group_pot_end, 
     ! calculate potential of my_group_pot_start-component to get a pw_grid  
     IF(my_group_pot_start>my_group_pot_end) THEN
        wf_vector=0.0_dp
        ! pseudo psi_L 
        CALL calculate_wavefunction(mo_coeff,1,psi_L,rho_g, atomic_kind_set,&
                                    qs_kind_set,cell,dft_control,particle_set, &
                                    pw_env_sub,basis_set_id=use_ri_aux_basis_set,&
                                    external_vector=wf_vector,&
                                    error=error_sub)
  
        rho_r%pw%cr3d = psi_L%pw%cr3d
        CALL pw_transfer(rho_r%pw, rho_g%pw, error=error_sub)
        CALL pw_poisson_solve(poisson_env,rho_g%pw, pair_energy, pot_g%pw,error=error_sub)
        CALL pw_transfer(pot_g%pw, rho_r%pw, error=error_sub)
        CALL pw_scale(rho_r%pw,rho_r%pw%pw_grid%dvol, error=error_sub)
  
     END IF
     
     DEALLOCATE(wf_vector)
  
     ! Fill pw_grid and 0.0_dp potential in every component of rho_r_array
     DO j=1,i_counter_diff_bf_at_any_atom
       IF(j<my_group_pot_start .OR. j>my_group_pot_end) THEN
  
         NULLIFY(rho_r_array(j)%pw)
         rho_r%pw%cr3d = 0.0_dp
         rho_r_array(j)%pw => rho_r%pw
  
         CALL pw_pool_create_pw(auxbas_pw_pool,rho_r_array(j)%pw,&
              use_data=REALDATA3D, in_space=REALSPACE,&
              error=error)
         CALL pw_copy(rho_r%pw,rho_r_array(j)%pw,error=error)
  
       END IF
     END DO
   
     color_exc=para_env_sub%mepos
     CALL mp_comm_split_direct(para_env%group,comm_exc,color_exc)
     NULLIFY(para_env_exc)
     CALL cp_para_env_create(para_env_exc,comm_exc,error=error)
  
     ! collect all calculated potentials
     DO j=1,i_counter_diff_bf_at_any_atom
       CALL mp_sum(rho_r_array(j)%pw%cr3d,para_env_exc%group)
     END DO
  
     CALL cp_para_env_release(para_env_exc,error=error)
  
     ALLOCATE(L_local_col_q(i_counter_diff_bf,my_group_L_q_size),STAT=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  
     L_local_col_q=0.0_dp
  
     i_counter = 0
  
  
     ! 2.1) calculate (q_alpha|q_beta)_p for every (l_alpha,m_alpha), (l_beta,m_beta)
     !      at every atom. Use the potential (q_alpha|_p calculated in 2.0)
  
     DO j=my_group_L_q_start, my_group_L_q_end
        i_counter=i_counter+1
  
        jatom = atom_set_afo_diffbf(j,1)
        jset  = atom_set_afo_diffbf(j,2)
        
        j_ref_bf   = diffbf_atom_arbitr_atom(j)
        j_ref_atom = atom_set_afo_diffbf_atanyatom(j_ref_bf,1)
        j_ref_set  = atom_set_afo_diffbf_atanyatom(j_ref_bf,2)
  
        j_ref_kind = kind_of(j_ref_atom)
        basis_set_b  =>  qs_kind_set(j_ref_kind)%ri_aux_basis_set
        zetb         =>  basis_set_b%zet
        lb_max       =>  basis_set_b%lmax
  
  
        NULLIFY(rs_v)
        NULLIFY(rs_descs)
        CALL pw_env_get(pw_env_sub, rs_descs=rs_descs, rs_grids=rs_v, error=error_sub)
        DO i=1,SIZE(rs_v)
          ! allocation
          CALL rs_grid_retain(rs_v(i)%rs_grid,error=error_sub)
        END DO
        CALL potential_pw2rs(rs_v,rho_r_array(j_ref_bf),pw_env_sub,error_sub)
  
        ! integrate the little bastards
        offset=0
        ! only take one single set once for every (l,m) combination
        DO i=1,i_counter_diff_sets
           iatom = diffsets_atom_set_bf(i,1)
           iset  = diffsets_atom_set_bf(i,2)
  
           ikind=kind_of(iatom)
           basis_set_a => qs_kind_set(ikind)%ri_aux_basis_set
  
           first_sgfa   =>  basis_set_a%first_sgf
           la_max       =>  basis_set_a%lmax
           la_min       =>  basis_set_a%lmin
           npgfa        =>  basis_set_a%npgf
           nseta        =   basis_set_a%nset
           nsgfa        =>  basis_set_a%nsgf_set
           rpgfa        =>  basis_set_a%pgf_radius
           set_radius_a =>  basis_set_a%set_radius
           sphi_a       =>  basis_set_a%sphi
           zeta         =>  basis_set_a%zet
  
           ra(:) = pbc(particle_set(iatom)%r,cell)-&
                   pbc(particle_set(jatom)%r,cell)+&
                   pbc(particle_set(j_ref_atom)%r,cell)
           rab=0.0_dp
           rab2=0.0_dp
  
           ncoa = npgfa(iset)*ncoset(la_max(iset))
           sgfa = first_sgfa(1,iset)
  
           ALLOCATE(I_tmp2(ncoa,1))
           I_tmp2=0.0_dp
           ALLOCATE(I_ab(nsgfa(iset),1))
           I_ab=0.0_dp
  
           igrid_level = gaussian_gridlevel(pw_env_sub%gridlevel_info,MINVAL(exp_q(:)))
           map_it_here=.FALSE.
           IF (.NOT. ALL (rs_v(igrid_level)%rs_grid%desc%perd == 1)) THEN
              DO dir = 1,3
                    ! bounds of local grid (i.e. removing the 'wings'), if periodic
                    tp(dir) = FLOOR(DOT_PRODUCT(cell%h_inv(dir,:),ra)*rs_v(igrid_level)%rs_grid%desc%npts(dir))
                    tp(dir) = MODULO ( tp(dir), rs_v(igrid_level)%rs_grid%desc%npts(dir) )
                    IF (rs_v(igrid_level)%rs_grid%desc%perd(dir) .NE. 1) THEN
                       lb(dir) = rs_v(igrid_level)%rs_grid%lb_local ( dir ) + rs_v(igrid_level)%rs_grid%desc%border
                       ub(dir) = rs_v(igrid_level)%rs_grid%ub_local ( dir ) - rs_v(igrid_level)%rs_grid%desc%border
                    ELSE
                       lb(dir) = rs_v(igrid_level)%rs_grid%lb_local ( dir )
                       ub(dir) = rs_v(igrid_level)%rs_grid%ub_local ( dir )
                    ENDIF
                    ! distributed grid, only map if it is local to the grid
                    location(dir)=tp(dir)+rs_v(igrid_level)%rs_grid%desc%lb(dir)
              ENDDO
              IF  (lb(1)<=location(1) .AND. location(1)<=ub(1) .AND. &
                   lb(2)<=location(2) .AND. location(2)<=ub(2) .AND. &
                   lb(3)<=location(3) .AND. location(3)<=ub(3)) THEN
                 map_it_here=.TRUE.
              ENDIF
           ELSE
              ! not distributed, just a round-robin distribution over the full set of CPUs
              IF (MODULO(offset,para_env_sub%num_pe)==para_env_sub%mepos) map_it_here=.TRUE.
           ENDIF
  
           offset=offset+nsgfa(iset)
  
           IF (map_it_here) THEN
             DO ipgf=1, npgfa(iset)
               na1=(ipgf - 1)*ncoset(la_max(iset)) + 1
               na2=ipgf*ncoset(la_max(iset))
               igrid_level = gaussian_gridlevel(pw_env_sub%gridlevel_info,exp_q(1))
  
               CALL integrate_pgf_product_rspace(la_max=la_max(iset),zeta=exp_q(1)/2.0_dp,la_min=la_min(iset),&
                                                 lb_max=0,zetb=exp_q(1)/2.0_dp,lb_min=0,&
                                                 ra=ra,rab=rab,rab2=rab2,&
                                                 rsgrid=rs_v(igrid_level)%rs_grid,&
                                                 cell=cell,&
                                                 cube_info=pw_env_sub%cube_info(igrid_level),&
                                                 hab=I_tmp2,&
                                                 o1=na1-1,&
                                                 o2=0,&
                                                 map_consistent=.TRUE.,&
                                                 eps_gvg_rspace=dft_control%qs_control%eps_gvg_rspace,&
                                                 calculate_forces=.FALSE.,&
                                                 error=error_sub)
             END DO
  
            CALL dgemm("T","N",nsgfa(iset),1,ncoa,&
                         1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                         I_tmp2(1,1),SIZE(I_tmp2,1),&
                         1.0_dp,I_ab(1,1),SIZE(I_ab,1))
  
             ! integration result has to be scaled since from basis set differing exponent exp_q
             ! causes different normalization
             L_local_col_q(offset-nsgfa(iset)+1:offset,i_counter)=I_ab(1:nsgfa(iset),1)*&
                                   (exp_q(1)/zeta(1,iset))**(la_max(iset)*0.5_dp+3.0_dp/4.0_dp)*&
                                   (exp_q(1)/zetb(1,j_ref_set))**(lb_max(j_ref_set)*0.5_dp+3.0_dp/4.0_dp)
           END IF
  
           DEALLOCATE(I_tmp2)
           DEALLOCATE(I_ab)
  
  
        END DO ! q_beta)-set loop      
  
        ! deallocation
        DO i=1,SIZE(rs_v)
          CALL rs_grid_release(rs_v(i)%rs_grid, error=error_sub)
        END DO
  
     END DO !(q_alpha| loop
  
     ! post cleanup
     DO j=1,i_counter_diff_bf_at_any_atom
       CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_r_array(j)%pw,&
                 accept_non_compatible=.TRUE., error=error)
     END DO
  
     DEALLOCATE(rho_r_array)
  
     DEALLOCATE(diffbf_noinbas)
     DEALLOCATE(diffbfatanyatom_noinbas)
  
  
  !  3) Distribute (q_alpha|q_beta)_p to L_local_col
  
     ALLOCATE(L_local_col(dimen_RI,my_group_L_size),STAT=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     L_local_col=0.0_dp
  
     DO  iatom=1, natom
       ikind=kind_of(iatom)
       basis_set_a  =>  qs_kind_set(ikind)%ri_aux_basis_set
       first_sgfa   =>  basis_set_a%first_sgf
       la_max       =>  basis_set_a%lmax
       la_min       =>  basis_set_a%lmin
       npgfa        =>  basis_set_a%npgf
       nseta        =   basis_set_a%nset
       nsgfa        =>  basis_set_a%nsgf_set
       rpgfa        =>  basis_set_a%pgf_radius
       set_radius_a =>  basis_set_a%set_radius
       sphi_a       =>  basis_set_a%sphi
       zeta         =>  basis_set_a%zet
  
       DO iset=1, nseta
  
         offsetb = 0
  
         DO jatom=1, natom
           jkind=kind_of(jatom)
           basis_set_b  => qs_kind_set(jkind)%ri_aux_basis_set
           first_sgfb   =>  basis_set_b%first_sgf
           lb_max       =>  basis_set_b%lmax
           lb_min       =>  basis_set_b%lmin
           npgfb        =>  basis_set_b%npgf
           nsetb        =   basis_set_b%nset
           nsgfb        =>  basis_set_b%nsgf_set
           rpgfb        =>  basis_set_b%pgf_radius
           set_radius_b =>  basis_set_b%set_radius
           sphi_b       =>  basis_set_b%sphi
           zetb         =>  basis_set_b%zet
  
           DO jset=1, nsetb
  
             ! for parallelization: check, if jset is in paralleliz. range
             IF(ref_at_se_sgf_ofs(jatom,jset,4) + nsgfb(jset) >= my_group_L_start .AND. &
                ref_at_se_sgf_ofs(jatom,jset,4) + 1 <= my_group_L_end ) THEN
  
               ! add contribution from (q_a|q_b)_p
               scaling = 1.0_dp/(zeta(1,iset))**(la_max(iset)*0.5_dp+3.0_dp/4.0_dp)/&
                       (zetb(1,jset))**(lb_max(jset)*0.5_dp+3.0_dp/4.0_dp)*&
                       exp_q(1)**(la_max(iset)*0.5_dp+3.0_dp/4.0_dp)*&
                       exp_q(1)**(lb_max(jset)*0.5_dp+3.0_dp/4.0_dp)
  
               DO i=1,nsgfa(iset)
                 DO j=1,nsgfb(jset)
                   ! for parallelization: check, if sgf_j is in paralleliz. range
                   IF(ref_at_se_sgf_ofs(jatom,jset,4) + j >= my_group_L_start .AND. &
                      ref_at_se_sgf_ofs(jatom,jset,4) + j <= my_group_L_end ) THEN
  
                     ! Scale (q_a|q_b) by c_a*c_b, only implemented for uncontracted Gaussians and lmin=lmax
                     L_local_col(ref_at_se_sgf_ofs(iatom,iset,4)+i,&
                                     ref_at_se_sgf_ofs(jatom,jset,4)+j-my_group_L_start+1) =&
                         L_local_col_q(naolimits(iatom,iset,1)+i-1,&
                                       naolimits(jatom,jset,1)+j-1-my_group_L_q_start+1)*&
                         scaling
                     
                   END IF
                 END DO
               END DO
   
               ! for s-functions alpha and beta, we need additional correction V_a(G=0)-V_q(G=0))*beta(G=0)
               ! and V_b(G=0)-V_q(G=0))*alpha(G=0), only add once for every subgroup (for mepos==0)
               IF (la_max(iset)==0 .AND. lb_max(jset)==0 .AND. para_env_sub%mepos == 0) THEN
                 L_local_col(ref_at_se_sgf_ofs(iatom,iset,4)+1,&
                                 ref_at_se_sgf_ofs(jatom,jset,4)-my_group_L_start+2) = &
                   L_local_col(ref_at_se_sgf_ofs(iatom,iset,4)+1,&
                                   ref_at_se_sgf_ofs(jatom,jset,4)-my_group_L_start+2) - &
                   Pi/(rho_r%pw%pw_grid%vol)*(4.0_dp*Pi**2.0_dp/(zeta(1,iset)*zetb(1,jset)))**0.75_dp*&
                   (2.0_dp/exp_q(1)-1.0_dp/zeta(1,iset)-1.0_dp/zetb(1,jset))
               END IF
  
             END IF
  
           END DO !jset
         END DO !jatom
       END DO !iset
     END DO !iatom
  
     DEALLOCATE(naolimits)
     DEALLOCATE(L_local_col_q) 
  
  
     ! End of calculation (q_alpha|q_beta) for every (l_alpha,m_alpha), (l_beta,m_beta) combination
 
 
     ! 4) Calculate (a|b)_s - (q_alpha|q_beta)_s and add them to (q_a|q_b)_p
  
     ! Loop over neighbor list
     CALL neighbor_list_iterator_create(nl_iterator,sab_orb_sub)
     DO WHILE (neighbor_list_iterate(nl_iterator)==0)
       CALL get_iterator_info(nl_iterator,ikind=ikind,jkind=jkind,inode=inode,&
                              iatom=iatom,jatom=jatom,r=rac)
  
       ikind=kind_of(iatom)
       basis_set_a  =>  qs_kind_set(ikind)%ri_aux_basis_set
       first_sgfa   =>  basis_set_a%first_sgf
       la_max       =>  basis_set_a%lmax
       la_min       =>  basis_set_a%lmin
       npgfa        =>  basis_set_a%npgf
       nseta        =   basis_set_a%nset
       nsgfa        =>  basis_set_a%nsgf_set
       rpgfa        =>  basis_set_a%pgf_radius
       set_radius_a =>  basis_set_a%set_radius
       sphi_a       =>  basis_set_a%sphi
       zeta         =>  basis_set_a%zet
  
       DO iset=1, nseta
  
           jkind=kind_of(jatom)
           basis_set_b  => qs_kind_set(jkind)%ri_aux_basis_set
           first_sgfb   =>  basis_set_b%first_sgf
           lb_max       =>  basis_set_b%lmax
           lb_min       =>  basis_set_b%lmin
           npgfb        =>  basis_set_b%npgf
           nsetb        =   basis_set_b%nset
           nsgfb        =>  basis_set_b%nsgf_set
           rpgfb        =>  basis_set_b%pgf_radius
           set_radius_b =>  basis_set_b%set_radius
           sphi_b       =>  basis_set_b%sphi
           zetb         =>  basis_set_b%zet
  
           DO jset=1, nsetb
  
             ! for parallelization: check, if (jatom,jset) is in paralleliz. range
             check_j = ref_at_se_sgf_ofs(jatom,jset,4) + nsgfb(jset) >= my_group_L_start .AND. &
                       ref_at_se_sgf_ofs(jatom,jset,4) + 1 <= my_group_L_end  
  
             check_i = ref_at_se_sgf_ofs(iatom,iset,4) + nsgfa(iset) >= my_group_L_start .AND. &
                       ref_at_se_sgf_ofs(iatom,iset,4) + 1 <= my_group_L_end .AND. &
                       iatom/=jatom 
  
  
             IF(check_j .OR. check_i) THEN
  
               ncoa_a = npgfa(iset)*ncoset(la_max(iset))
               ncoa_b = npgfb(jset)*ncoset(lb_max(jset))
               sgfa = first_sgfa(1,iset)
               sgfb = first_sgfb(1,jset)
  
  
               ALLOCATE(I_tmp(ncoa_a,nsgfb(jset)),STAT=istat)
               CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
               ALLOCATE(I_tmp2(nsgfa(iset),nsgfb(jset)),STAT=istat)
               CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  
               ALLOCATE(f(0:la_max(iset)+lb_max(jset)+2),&
                        v(ncoa_a,ncoa_b,la_max(iset)+lb_max(jset)+1),&
                        vac(ncoa_a,ncoa_b),STAT=istat)
               CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
               f=0.0_dp
               v=0.0_dp
               vac=0.0_dp
               I_tmp=0.0_dp
               I_tmp2=0.0_dp
  
               rac2 = rac(1)*rac(1) + rac(2)*rac(2) + rac(3)*rac(3)
  
  
               ! calculate spherical matrix elements (a|b)_s 
               CALL coulomb2_new(la_max=la_max(iset),npgfa=npgfa(iset),zeta=zeta(:,iset),&
                                 la_min=la_min(iset), &
                                 lc_max=lb_max(jset),npgfc=npgfb(jset),zetc=zetb(:,jset),&
                                 lc_min=lb_min(jset),rac=rac,&
                                 rac2=rac2,vac=vac,v=v,f=f)
  
  
               ! trafo to spherical and add (a|b)_s to L_local_col
               I_tmp(1:ncoa_a,1:nsgfb(jset)) = MATMUL(vac(1:ncoa_a,1:ncoa_b),sphi_b(1:ncoa_b,sgfb:(sgfb+nsgfb(jset)-1)))
               I_tmp2(1:nsgfa(iset),1:nsgfb(jset)) = MATMUL(TRANSPOSE(sphi_a(1:ncoa_a,sgfa:(sgfa+nsgfa(iset)-1))),&
                                                     I_tmp(1:ncoa_a,1:nsgfb(jset)))
  
               f=0.0_dp
               v=0.0_dp
               vac=0.0_dp
  
               ! calculate -(q_a|q_b)_s
               CALL coulomb2_new(la_max=la_max(iset),npgfa=npgfa(iset),zeta=exp_q(:),&
                                 la_min=la_min(iset), &
                                 lc_max=lb_max(jset),npgfc=npgfb(jset),zetc=exp_q(:),&
                                 lc_min=lb_min(jset),rac=rac,&
                                 rac2=rac2,vac=vac,v=v,f=f)
  
  
               ! trafo to spherical and scale by c_a*c_b, only implemented for uncontracted Gaussians and lmin=lmax
               I_tmp(1:ncoa_a,1:nsgfb(jset)) = MATMUL(vac(1:ncoa_a,1:ncoa_b),&
                                   sphi_b(1:ncoa_b,ref_at_se_sgf_ofs(jatom,jset,3):&
                                   (ref_at_se_sgf_ofs(jatom,jset,3)+nsgfb(jset)-1)))
               I_tmp2(1:nsgfa(iset),1:nsgfb(jset)) = I_tmp2(1:nsgfa(iset),1:nsgfb(jset))- &
                   MATMUL(TRANSPOSE(sphi_a(1:ncoa_a,ref_at_se_sgf_ofs(iatom,iset,3):&
                                           (ref_at_se_sgf_ofs(iatom,iset,3)+nsgfa(iset)-1))),&
                          I_tmp(1:ncoa_a,1:nsgfb(jset)))/&
                          (zeta(1,iset))**(la_max(iset)*0.5_dp+3.0_dp/4.0_dp)/&
                          (zetb(1,jset))**(lb_max(jset)*0.5_dp+3.0_dp/4.0_dp)*&
                          exp_q(1)**(la_max(iset)*1.0_dp+3.0_dp/2.0_dp)*&
                          exp_q(1)**(lb_max(jset)*1.0_dp+3.0_dp/2.0_dp)/&
                          (zeta(1,ref_at_se_sgf_ofs(iatom,iset,2)))**(la_max(iset)*0.5_dp+3.0_dp/4.0_dp)/&
                          (zetb(1,ref_at_se_sgf_ofs(jatom,jset,2)))**(lb_max(jset)*0.5_dp+3.0_dp/4.0_dp)
  
  
               ! add contribution from (a|b)_s-(q_a|q_b)_s
               DO i=1,nsgfa(iset)
                 DO j=1,nsgfb(jset)
                   ! for parallelization: check, if sgf_j is in paralleliz. range
                   IF(ref_at_se_sgf_ofs(jatom,jset,4) + j >= my_group_L_start .AND. &
                      ref_at_se_sgf_ofs(jatom,jset,4) + j <= my_group_L_end   .AND. &
                      check_j) THEN
  
                       L_local_col(ref_at_se_sgf_ofs(iatom,iset,4)+i,&
                                       ref_at_se_sgf_ofs(jatom,jset,4)+j-my_group_L_start+1) = &
                                    L_local_col(ref_at_se_sgf_ofs(iatom,iset,4)+i,&
                                                    ref_at_se_sgf_ofs(jatom,jset,4)+j-my_group_L_start+1)+&
                                    I_tmp2(i,j) 
                   END IF
  
                   IF(ref_at_se_sgf_ofs(iatom,iset,4) + i >= my_group_L_start .AND. &
                      ref_at_se_sgf_ofs(iatom,iset,4) + i <= my_group_L_end   .AND. &
                      check_i) THEN
                     L_local_col(ref_at_se_sgf_ofs(jatom,jset,4)+j,&
                                     ref_at_se_sgf_ofs(iatom,iset,4)+i-my_group_L_start+1) = &
                                   L_local_col(ref_at_se_sgf_ofs(jatom,jset,4)+j,&
                                                   ref_at_se_sgf_ofs(iatom,iset,4)+i-my_group_L_start+1)+&
                                   I_tmp2(i,j)
                   END IF
                 END DO
               END DO  
  
  
               DEALLOCATE(f,v,vac,I_tmp,I_tmp2,STAT=istat)
            
             END IF !parallelization
  
           END DO !jset
         END DO !iset
  
     END DO
  
     CALL neighbor_list_iterator_release(nl_iterator)
 
     DEALLOCATE(diffsets_atom_set_bf,ref_at_se_sgf_ofs) 
     DEALLOCATE(atom_set_afo_diffbf) 
     DEALLOCATE(sizes_array_pot,starts_array_pot,ends_array_pot)
 
     CALL mp_sum(L_local_col,para_env_sub%group)
  
     ! end if of new algorithm
   ELSE  
     ! old algorithm

     ALLOCATE(L_local_col(dimen_RI,my_group_L_size),STAT=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     L_local_col=0.0_dp
  
     ALLOCATE(wf_vector(dimen_RI),STAT=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  
     i_counter=0
     DO LLL=my_group_L_start, my_group_L_end
        i_counter=i_counter+1
  
        wf_vector=0.0_dp
        wf_vector(LLL)=1.0_dp
  
        ! pseudo psi_L
        CALL calculate_wavefunction(mo_coeff,1,psi_L,rho_g, atomic_kind_set,&
                                    qs_kind_set,cell,dft_control,particle_set, &
                                    pw_env_sub,basis_set_id=use_ri_aux_basis_set,&
                                    external_vector=wf_vector,&
                                    error=error_sub)
  
        CALL timeset(routineN//"_pot_lm",handle3)
        rho_r%pw%cr3d = psi_L%pw%cr3d
        CALL pw_transfer(rho_r%pw, rho_g%pw, error=error_sub)
        CALL pw_poisson_solve(poisson_env,rho_g%pw, pair_energy, pot_g%pw,error=error_sub)
        CALL pw_transfer(pot_g%pw, rho_r%pw, error=error_sub)
        CALL pw_scale(rho_r%pw,rho_r%pw%pw_grid%dvol, error=error_sub)
  
        NULLIFY(rs_v)
        NULLIFY(rs_descs)
        CALL pw_env_get(pw_env_sub, rs_descs=rs_descs, rs_grids=rs_v, error=error_sub)
        DO i=1,SIZE(rs_v)
          CALL rs_grid_retain(rs_v(i)%rs_grid,error=error_sub)
        END DO
        CALL potential_pw2rs(rs_v,rho_r,pw_env_sub,error_sub)
  
        CALL timestop(handle3)
   
        ! integrate the little bastards
        offset=0
        DO iatom=1, natom
          ikind=kind_of(iatom)
          basis_set_a => qs_kind_set(ikind)%ri_aux_basis_set
  
          first_sgfa   =>  basis_set_a%first_sgf
          la_max       =>  basis_set_a%lmax
          la_min       =>  basis_set_a%lmin
          npgfa        =>  basis_set_a%npgf
          nseta        =   basis_set_a%nset
          nsgfa        =>  basis_set_a%nsgf_set
          rpgfa        =>  basis_set_a%pgf_radius
          set_radius_a =>  basis_set_a%set_radius
          sphi_a       =>  basis_set_a%sphi
          zeta         =>  basis_set_a%zet
  
          ra(:) = pbc(particle_set(iatom)%r,cell)
          rab=0.0_dp
          rab2=0.0_dp
  
          DO iset=1, nseta
           ncoa = npgfa(iset)*ncoset(la_max(iset))
           sgfa = first_sgfa(1,iset)
  
           ALLOCATE(I_tmp2(ncoa,1))
           I_tmp2=0.0_dp
           ALLOCATE(I_ab(nsgfa(iset),1))
           I_ab=0.0_dp
  
  
           igrid_level = gaussian_gridlevel(pw_env_sub%gridlevel_info,MINVAL(zeta(:,iset)))
           map_it_here=.FALSE.
           IF (.NOT. ALL (rs_v(igrid_level)%rs_grid%desc%perd == 1)) THEN
              DO dir = 1,3
                    ! bounds of local grid (i.e. removing the 'wings'), if periodic
                    tp(dir) = FLOOR(DOT_PRODUCT(cell%h_inv(dir,:),ra)*rs_v(igrid_level)%rs_grid%desc%npts(dir))
                    tp(dir) = MODULO ( tp(dir), rs_v(igrid_level)%rs_grid%desc%npts(dir) )
                    IF (rs_v(igrid_level)%rs_grid%desc%perd(dir) .NE. 1) THEN
                       lb(dir) = rs_v(igrid_level)%rs_grid%lb_local ( dir ) + rs_v(igrid_level)%rs_grid%desc%border
                       ub(dir) = rs_v(igrid_level)%rs_grid%ub_local ( dir ) - rs_v(igrid_level)%rs_grid%desc%border
                    ELSE
                       lb(dir) = rs_v(igrid_level)%rs_grid%lb_local ( dir )
                       ub(dir) = rs_v(igrid_level)%rs_grid%ub_local ( dir )
                    ENDIF
                    ! distributed grid, only map if it is local to the grid
                    location(dir)=tp(dir)+rs_v(igrid_level)%rs_grid%desc%lb(dir)
              ENDDO
              IF  (lb(1)<=location(1) .AND. location(1)<=ub(1) .AND. &
                   lb(2)<=location(2) .AND. location(2)<=ub(2) .AND. &
                   lb(3)<=location(3) .AND. location(3)<=ub(3)) THEN
                 map_it_here=.TRUE.
              ENDIF
           ELSE
              ! not distributed, just a round-robin distribution over the full set of CPUs
              IF (MODULO(offset,para_env_sub%num_pe)==para_env_sub%mepos) map_it_here=.TRUE.
           ENDIF
  
           offset=offset+nsgfa(iset)
  
           IF (map_it_here) THEN
             DO ipgf=1, npgfa(iset)
               sgfa = first_sgfa(1,iset)
  
               na1=(ipgf - 1)*ncoset(la_max(iset)) + 1
               na2=ipgf*ncoset(la_max(iset))
               igrid_level = gaussian_gridlevel(pw_env_sub%gridlevel_info,zeta(ipgf,iset))
  
               CALL integrate_pgf_product_rspace(la_max=la_max(iset),zeta=zeta(ipgf,iset)/2.0_dp,la_min=la_min(iset),&
                                                 lb_max=0,zetb=zeta(ipgf,iset)/2.0_dp,lb_min=0,&
                                                 ra=ra,rab=rab,rab2=rab2,&
                                                 rsgrid=rs_v(igrid_level)%rs_grid,&
                                                 cell=cell,&
                                                 cube_info=pw_env_sub%cube_info(igrid_level),&
                                                 hab=I_tmp2,&
                                                 o1=na1-1,&
                                                 o2=0,&
                                                 map_consistent=.TRUE.,&
                                                 eps_gvg_rspace=dft_control%qs_control%eps_gvg_rspace,&
                                                 calculate_forces=.FALSE.,&
                                                 error=error_sub)
             END DO
           
             CALL dgemm("T","N",nsgfa(iset),1,ncoa,&
                         1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                         I_tmp2(1,1),SIZE(I_tmp2,1),&
                         1.0_dp,I_ab(1,1),SIZE(I_ab,1))
   
             L_local_col(offset-nsgfa(iset)+1:offset,i_counter)=I_ab(1:nsgfa(iset),1)
           END IF
  
           DEALLOCATE(I_tmp2)
           DEALLOCATE(I_ab)
  
          END DO
        END DO
  
        DO i=1,SIZE(rs_v)
          CALL rs_grid_release(rs_v(i)%rs_grid, error=error_sub)
        END DO
  
      END DO
 
      DEALLOCATE(wf_vector)
 
      CALL mp_sum(L_local_col,para_env_sub%group)
  
    END IF  !end old code

    CALL timestop(handle2)

    ! split the total number of proc in a subgroup of the size of ~1/10 of the 
    ! total num of proc
    best_group_size=para_env%num_pe

    strat_group_size=MAX(1,para_env%num_pe/10)

    min_mem_for_QK=REAL(dimen_RI,KIND=dp)*dimen_RI*3.0_dp*8.0_dp/1024_dp/1024_dp

    group_size=strat_group_size-1
    DO iproc=strat_group_size, para_env%num_pe
      group_size=group_size+1
      ! check that group_size is a multiple of sub_group_size and a divisor of
      ! the total num of proc
      IF(MOD(para_env%num_pe,group_size)/=0.OR.MOD(group_size,para_env_sub%num_pe)/=0) CYCLE

      ! check for memory
      IF(REAL(group_size,KIND=dp)*mp2_memory<min_mem_for_QK) CYCLE

      best_group_size=group_size
      EXIT
    END DO

    ! create the L group
    color_L=para_env%mepos/best_group_size
    CALL mp_comm_split_direct(para_env%group,comm_L,color_L)
    NULLIFY(para_env_L)
    CALL cp_para_env_create(para_env_L,comm_L,error=error)

    ! create the blacs_L
    NULLIFY(blacs_env_L)
    CALL cp_blacs_env_create(blacs_env=blacs_env_L, para_env=para_env_L, error=error)

    ! now create the exchange group (for communication only between members not belonging to the
    ! same group
    sub_sub_color=para_env_sub%mepos
    CALL mp_comm_split_direct(para_env_L%group,comm_exchange,sub_sub_color)
    NULLIFY(para_env_exchange)
    CALL cp_para_env_create(para_env_exchange,comm_exchange,error=error_sub)

    ! crate the proc maps
    ALLOCATE(proc_map(-para_env_exchange%num_pe:2*para_env_exchange%num_pe-1))
    DO i=0,para_env_exchange%num_pe-1
      proc_map(i)=i
      proc_map(-i-1)=para_env_exchange%num_pe-i-1
      proc_map(para_env_exchange%num_pe+i)=i
    END DO

    ! create the information array
    ALLOCATE(sub_sizes_array(0:para_env_exchange%num_pe-1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    sub_sizes_array=0
    ALLOCATE(sub_starts_array(0:para_env_exchange%num_pe-1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    sub_starts_array=0
    ALLOCATE(sub_ends_array(0:para_env_exchange%num_pe-1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    sub_ends_array=0

    sub_sizes_array(para_env_exchange%mepos)=my_group_L_size
    sub_starts_array(para_env_exchange%mepos)=my_group_L_start
    sub_ends_array(para_env_exchange%mepos)=my_group_L_end

    CALL mp_sum(sub_sizes_array,para_env_exchange%group)
    CALL mp_sum(sub_starts_array,para_env_exchange%group)
    CALL mp_sum(sub_ends_array,para_env_exchange%group)

    ! create the full matrix L defined in the L group
    NULLIFY(fm_matrix_L)
    NULLIFY(fm_struct)
    CALL cp_fm_struct_create(fm_struct,context=blacs_env_L,nrow_global=dimen_RI,&
                             ncol_global=dimen_RI,para_env=para_env_L,error=error)
    CALL cp_fm_create(fm_matrix_L,fm_struct,name="fm_matrix_L",error=error)
    CALL cp_fm_struct_release(fm_struct,error=error)

    CALL cp_fm_set_all(matrix=fm_matrix_L,alpha=0.0_dp,error=error)

    CALL cp_fm_get_info(matrix=fm_matrix_L,&
                        nrow_local=nrow_local,&
                        ncol_local=ncol_local,&
                        row_indices=row_indices,&
                        col_indices=col_indices,&
                        error=error)
   
    DO jjB=1, ncol_local
      j_global=col_indices(jjB)
      IF(j_global>=my_group_L_start.AND.j_global<=my_group_L_end) THEN
        DO iiB=1, nrow_local
          i_global=row_indices(iiB)
          fm_matrix_L%local_data(iiB,jjB)=L_local_col(i_global,j_global-my_group_L_start+1)
        END DO
      END IF
    END DO
 
    proc_send_static=proc_map(para_env_exchange%mepos+1)
    proc_receive_static=proc_map(para_env_exchange%mepos-1)
 
    DO proc_shift=1, para_env_exchange%num_pe-1
      proc_send=proc_map(para_env_exchange%mepos+proc_shift)
      proc_receive=proc_map(para_env_exchange%mepos-proc_shift)
     
      rec_L_size=sub_sizes_array(proc_receive)
      rec_L_start=sub_starts_array(proc_receive)
      rec_L_end=sub_ends_array(proc_receive)

      ALLOCATE(L_external_col(dimen_RI,rec_L_size))
      L_external_col=0.0_dp

      CALL  mp_sendrecv(L_local_col,proc_send_static,L_external_col,proc_receive_static,para_env_exchange%group)

      DO jjB=1, ncol_local
        j_global=col_indices(jjB)
        IF(j_global>=rec_L_start.AND.j_global<=rec_L_end) THEN
          DO iiB=1, nrow_local
            i_global=row_indices(iiB)
            fm_matrix_L%local_data(iiB,jjB)=L_external_col(i_global,j_global-rec_L_start+1)
          END DO
        END IF
      END DO

      DEALLOCATE(L_local_col)
      ALLOCATE(L_local_col(dimen_RI,rec_L_size),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      L_local_col(:,:)=L_external_col

      DEALLOCATE(L_external_col)

    END DO

    DEALLOCATE(L_local_col)

    ! free the old exchange group stuff
    DEALLOCATE(proc_map)
    CALL cp_para_env_release(para_env_exchange,error=error_sub)

    DEALLOCATE(sub_sizes_array)
    DEALLOCATE(sub_starts_array)
    DEALLOCATE(sub_ends_array)

    ! create the new group for the mp_sum of the local data
    sub_sub_color=para_env_L%mepos
    CALL mp_comm_split_direct(para_env%group,comm_exchange,sub_sub_color)
    NULLIFY(para_env_exchange)
    CALL cp_para_env_create(para_env_exchange,comm_exchange,error=error_sub)

    CALL mp_sum(fm_matrix_L%local_data,para_env_exchange%group)

    CALL cp_para_env_release(para_env_exchange,error=error_sub)

    cond_num=1.0_dp
    num_small_eigen=0
    IF(calc_PQ_cond_num) THEN
      ! calculate the condition number of the (P|Q) matrix
      ! create a copy of the matrix
      NULLIFY(fm_matrix_L_diag)
      NULLIFY(fm_struct)
      CALL cp_fm_struct_create(fm_struct,context=blacs_env_L,nrow_global=dimen_RI,&
                               ncol_global=dimen_RI,para_env=para_env_L,error=error)
      CALL cp_fm_create(fm_matrix_L_diag,fm_struct,name="fm_matrix_L_diag",error=error)
      CALL cp_fm_struct_release(fm_struct,error=error)

      CALL cp_fm_set_all(matrix=fm_matrix_L_diag,alpha=0.0_dp,error=error) 

      CALL cp_fm_to_fm(source=fm_matrix_L,destination=fm_matrix_L_diag,error=error)

      ALLOCATE(egen_L(dimen_RI),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

      egen_L=0.0_dp
      CALL cp_fm_syevx(matrix=fm_matrix_L_diag,eigenvalues=egen_L,error=error) 
   
      num_small_eigen=0
      DO iiB=1, dimen_RI
        IF(ABS(egen_L(iiB))<0.001_dp) num_small_eigen=num_small_eigen+1
      END DO

      cond_num=MAXVAL(ABS(egen_L))/MINVAL(ABS(egen_L))

      CALL cp_fm_release(fm_matrix_L_diag, error=error)

      DEALLOCATE(egen_L)
    END IF

    ! do cholesky decomposition
    CALL cp_fm_cholesky_decompose(matrix=fm_matrix_L, n=dimen_RI, info_out=info_chol, error=error)
    CPPostcondition(info_chol==0,cp_failure_level,routineP,error,failure)

    CALL cp_fm_triangular_invert(matrix_a=fm_matrix_L,uplo_tr='U',error=error)

    ! clean the lower part of the L^{-1} matrix (just to not have surprises afterwards)
    CALL cp_fm_get_info(matrix=fm_matrix_L,&
                        nrow_local=nrow_local,&
                        ncol_local=ncol_local,&
                        row_indices=row_indices,&
                        col_indices=col_indices,&
                        error=error)
    DO iiB=1, nrow_local
      i_global=row_indices(iiB)
      DO jjB=1, ncol_local
        j_global=col_indices(jjB)
        IF(j_global<i_global) fm_matrix_L%local_data(iiB,jjB)=0.0_dp
      END DO
    END DO

    ! release blacs_env
    CALL cp_blacs_env_release(blacs_env_L,error=error)

    CALL timestop(handle)

  END SUBROUTINE calculate_Lmin1

! *****************************************************************************
!> \brief ...
!> \param para_env ...
!> \param dimen_RI ...
!> \param fm_matrix_L ...
!> \param my_group_L_start ...
!> \param my_group_L_end ...
!> \param my_group_L_size ...
!> \param my_Lrows ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE grep_Lcols(para_env,dimen_RI,fm_matrix_L,&
                        my_group_L_start,my_group_L_end,my_group_L_size,my_Lrows,&
                        error)
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER                                  :: dimen_RI
    TYPE(cp_fm_type), POINTER                :: fm_matrix_L
    INTEGER                                  :: my_group_L_start, &
                                                my_group_L_end, &
                                                my_group_L_size
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: my_Lrows
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'grep_Lcols', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, handle2, i_global, iiB, j_global, jjB, &
      max_row_col_local, ncol_local, ncol_rec, nrow_local, nrow_rec, &
      proc_receive, proc_receive_static, proc_send, proc_send_static, &
      proc_shift, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: proc_map
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: local_col_row_info, &
                                                rec_col_row_info
    INTEGER, DIMENSION(:), POINTER           :: col_indices, col_indices_rec, &
                                                row_indices, row_indices_rec
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: local_L, local_L_internal, &
                                                rec_L

    CALL timeset(routineN,handle)

    ALLOCATE(my_Lrows(dimen_RI,my_group_L_size),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    my_Lrows=0.0_dp

    ! proc_map, vector that replicate the processor numbers also
    ! for negative and positive number > num_pe
    ! needed to know which is the processor, to respect to another one,
    ! for a given shift
    ALLOCATE(proc_map(-para_env%num_pe:2*para_env%num_pe-1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DO iiB=0,para_env%num_pe-1
      proc_map(iiB)=iiB
      proc_map(-iiB-1)=para_env%num_pe-iiB-1
      proc_map(para_env%num_pe+iiB)=iiB
    END DO

    CALL cp_fm_get_info(matrix=fm_matrix_L,&
                        nrow_local=nrow_local,&
                        ncol_local=ncol_local,&
                        row_indices=row_indices,&
                        col_indices=col_indices,&
                        local_data=local_L_internal,&
                        error=error)

    ALLOCATE(local_L(nrow_local,ncol_local),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    local_L=local_L_internal(1:nrow_local,1:ncol_local)

    max_row_col_local=MAX(nrow_local,ncol_local)
    CALL mp_max(max_row_col_local,para_env%group)

    ALLOCATE(local_col_row_info(0:max_row_col_local,2),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    local_col_row_info=0
    ! 0,1 nrows
    local_col_row_info(0,1)=nrow_local
    local_col_row_info(1:nrow_local,1)=row_indices(1:nrow_local)
    ! 0,2 ncols
    local_col_row_info(0,2)=ncol_local
    local_col_row_info(1:ncol_local,2)=col_indices(1:ncol_local)

    ALLOCATE(rec_col_row_info(0:max_row_col_local,2),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ! accumulate data on my_Lrows starting from myself
    DO jjB=1, ncol_local
      j_global=col_indices(jjB)
      IF(j_global>=my_group_L_start.AND.j_global<=my_group_L_end) THEN
        DO iiB=1, nrow_local
          i_global=row_indices(iiB)
          my_Lrows(i_global,j_global-my_group_L_start+1)=local_L(iiB,jjB)
        END DO
      END IF
    END DO    

    proc_send_static=proc_map(para_env%mepos+1)
    proc_receive_static=proc_map(para_env%mepos-1)

    CALL timeset(routineN//"_comm",handle2)

    DO proc_shift=1, para_env%num_pe-1
      proc_send=proc_map(para_env%mepos+proc_shift)
      proc_receive=proc_map(para_env%mepos-proc_shift)

      ! first exchange information on the local data
      rec_col_row_info=0
      CALL  mp_sendrecv(local_col_row_info,proc_send_static,rec_col_row_info,proc_receive_static,para_env%group)
      nrow_rec=rec_col_row_info(0,1)
      ncol_rec=rec_col_row_info(0,2)

      ALLOCATE(row_indices_rec(nrow_rec),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      row_indices_rec=rec_col_row_info(1:nrow_rec,1)

      ALLOCATE(col_indices_rec(ncol_rec),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      col_indices_rec=rec_col_row_info(1:ncol_rec,2)

      ALLOCATE(rec_L(nrow_rec,ncol_rec),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      rec_L=0.0_dp

      ! then send and receive the real data
      CALL  mp_sendrecv(local_L,proc_send_static,rec_L,proc_receive_static,para_env%group)

      ! accumulate the received data on my_Lrows
      DO jjB=1, ncol_rec
        j_global=col_indices_rec(jjB)
        IF(j_global>=my_group_L_start.AND.j_global<=my_group_L_end) THEN
          DO iiB=1, nrow_rec
            i_global=row_indices_rec(iiB)
            my_Lrows(i_global,j_global-my_group_L_start+1)=rec_L(iiB,jjB)
          END DO
        END IF
      END DO

      local_col_row_info(:,:)=rec_col_row_info
      DEALLOCATE(local_L)
      ALLOCATE(local_L(nrow_rec,ncol_rec),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      local_L=rec_L

      DEALLOCATE(col_indices_rec)
      DEALLOCATE(row_indices_rec)
      DEALLOCATE(rec_L)
    END DO
    CALL timestop(handle2)

    DEALLOCATE(local_col_row_info)
    DEALLOCATE(rec_col_row_info)
    DEALLOCATE(proc_map)
    DEALLOCATE(local_L)

    CALL timestop(handle)

  END SUBROUTINE

! *****************************************************************************
!> \brief ...
!> \param para_env_sub ...
!> \param fm_BIb_jb ...
!> \param BIb_jb ...
!> \param max_row_col_local ...
!> \param homo ...
!> \param virtual ...
!> \param proc_map ...
!> \param local_col_row_info ...
!> \param my_B_virtual_end ...
!> \param my_B_virtual_start ...
!> \param my_B_size ...
!> \param error_sub ...
! *****************************************************************************
  SUBROUTINE grep_my_integrals(para_env_sub,fm_BIb_jb,BIb_jb,max_row_col_local,&
                               homo,virtual,proc_map,local_col_row_info,&
                               my_B_virtual_end,my_B_virtual_start,my_B_size,&
                               error_sub)
    TYPE(cp_para_env_type), POINTER          :: para_env_sub
    TYPE(cp_fm_type), POINTER                :: fm_BIb_jb
    REAL(KIND=dp), DIMENSION(:, :)           :: BIb_jb
    INTEGER                                  :: max_row_col_local, homo, &
                                                virtual
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: proc_map
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: local_col_row_info
    INTEGER                                  :: my_B_virtual_end, &
                                                my_B_virtual_start, my_B_size
    TYPE(cp_error_type), INTENT(inout)       :: error_sub

    CHARACTER(LEN=*), PARAMETER :: routineN = 'grep_my_integrals', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i_global, iiB, j_global, jjB, &
                                                ncol_rec, nrow_rec, &
                                                proc_receive, proc_send, &
                                                proc_shift, stat
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: rec_col_row_info
    INTEGER, DIMENSION(:), POINTER           :: col_indices_rec, &
                                                row_indices_rec
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: local_BI, rec_BI

    failure=.FALSE.

    ALLOCATE(rec_col_row_info(0:max_row_col_local,2),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)

    rec_col_row_info(:,:)=local_col_row_info

    nrow_rec=rec_col_row_info(0,1)
    ncol_rec=rec_col_row_info(0,2)

    ALLOCATE(row_indices_rec(nrow_rec),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
    row_indices_rec=rec_col_row_info(1:nrow_rec,1)

    ALLOCATE(col_indices_rec(ncol_rec),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
    col_indices_rec=rec_col_row_info(1:ncol_rec,2)

    ! accumulate data on BIb_jb buffer starting from myself
    DO jjB=1, ncol_rec
      j_global=col_indices_rec(jjB)
      IF(j_global>=my_B_virtual_start.AND.j_global<=my_B_virtual_end) THEN
        DO iiB=1, nrow_rec
          i_global=row_indices_rec(iiB)
          BIb_jb(j_global-my_B_virtual_start+1,i_global)=fm_BIb_jb%local_data(iiB,jjB)
        END DO
      END IF
    END DO

    DEALLOCATE(row_indices_rec)
    DEALLOCATE(col_indices_rec)

    IF(para_env_sub%num_pe>1) THEN
      ALLOCATE(local_BI(nrow_rec,ncol_rec))
      local_BI(1:nrow_rec,1:ncol_rec)=fm_BIb_jb%local_data(1:nrow_rec,1:ncol_rec)

      DO proc_shift=1, para_env_sub%num_pe-1
        proc_send=proc_map(para_env_sub%mepos+proc_shift)
        proc_receive=proc_map(para_env_sub%mepos-proc_shift)

        ! first exchange information on the local data
        rec_col_row_info=0
        CALL  mp_sendrecv(local_col_row_info,proc_send,rec_col_row_info,proc_receive,para_env_sub%group)
        nrow_rec=rec_col_row_info(0,1)
        ncol_rec=rec_col_row_info(0,2)

        ALLOCATE(row_indices_rec(nrow_rec),STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
        row_indices_rec=rec_col_row_info(1:nrow_rec,1)

        ALLOCATE(col_indices_rec(ncol_rec),STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
        col_indices_rec=rec_col_row_info(1:ncol_rec,2)

        ALLOCATE(rec_BI(nrow_rec,ncol_rec),STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
        rec_BI=0.0_dp

        ! then send and receive the real data
        CALL  mp_sendrecv(local_BI,proc_send,rec_BI,proc_receive,para_env_sub%group)

        ! accumulate the received data on BIb_jb buffer 
        DO jjB=1, ncol_rec
          j_global=col_indices_rec(jjB)
          IF(j_global>=my_B_virtual_start.AND.j_global<=my_B_virtual_end) THEN
            DO iiB=1, nrow_rec
              i_global=row_indices_rec(iiB)
              BIb_jb(j_global-my_B_virtual_start+1,i_global)=rec_BI(iiB,jjB)
            END DO
          END IF
        END DO

        DEALLOCATE(col_indices_rec)
        DEALLOCATE(row_indices_rec)
        DEALLOCATE(rec_BI)
      END DO

      DEALLOCATE(local_BI)
    END IF

    DEALLOCATE(rec_col_row_info)

  END SUBROUTINE grep_my_integrals

! *****************************************************************************
!> \brief ...
!> \param Emp2 ...
!> \param Emp2_Cou ...
!> \param Emp2_EX ...
!> \param BIb_C ...
!> \param mp2_env ...
!> \param para_env ...
!> \param para_env_sub ...
!> \param color_sub ...
!> \param ends_array ...
!> \param ends_B_virtual ...
!> \param sizes_array ...
!> \param sizes_B_virtual ...
!> \param starts_array ...
!> \param starts_B_virtual ...
!> \param Eigenval ...
!> \param nmo ...
!> \param homo ...
!> \param dimen_RI ...
!> \param unit_nr ...
!> \param calc_forces ...
!> \param calc_ex ...
!> \param error ...
!> \param error_sub ...
!> \param open_shell_SS ...
!> \param BIb_C_beta ...
!> \param homo_beta ...
!> \param Eigenval_beta ...
!> \param ends_B_virtual_beta ...
!> \param sizes_B_virtual_beta ...
!> \param starts_B_virtual_beta ...
! *****************************************************************************
  SUBROUTINE mp2_ri_gpw_compute_en(Emp2,Emp2_Cou,Emp2_EX,BIb_C,mp2_env,para_env,para_env_sub,color_sub,&
                                   ends_array,ends_B_virtual,sizes_array,sizes_B_virtual,starts_array,starts_B_virtual,&
                                   Eigenval,nmo,homo,dimen_RI,unit_nr,calc_forces,calc_ex,error,error_sub,&
                                   open_shell_SS,BIb_C_beta,homo_beta,Eigenval_beta,&
                                   ends_B_virtual_beta,sizes_B_virtual_beta, starts_B_virtual_beta)
    REAL(KIND=dp)                            :: Emp2, Emp2_Cou, Emp2_EX
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: BIb_C
    TYPE(mp2_type), POINTER                  :: mp2_env
    TYPE(cp_para_env_type), POINTER          :: para_env, para_env_sub
    INTEGER                                  :: color_sub
    INTEGER, ALLOCATABLE, DIMENSION(:) :: ends_array, ends_B_virtual, &
      sizes_array, sizes_B_virtual, starts_array, starts_B_virtual
    REAL(KIND=dp), DIMENSION(:)              :: Eigenval
    INTEGER                                  :: nmo, homo, dimen_RI, unit_nr
    LOGICAL, INTENT(IN)                      :: calc_forces
    LOGICAL                                  :: calc_ex
    TYPE(cp_error_type), INTENT(inout)       :: error, error_sub
    LOGICAL, OPTIONAL                        :: open_shell_SS
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :), OPTIONAL           :: BIb_C_beta
    INTEGER, OPTIONAL                        :: homo_beta
    REAL(KIND=dp), DIMENSION(:), OPTIONAL    :: Eigenval_beta
    INTEGER, ALLOCATABLE, DIMENSION(:), &
      OPTIONAL                               :: ends_B_virtual_beta, &
                                                sizes_B_virtual_beta, &
                                                starts_B_virtual_beta

    CHARACTER(LEN=*), PARAMETER :: routineN = 'mp2_ri_gpw_compute_en', &
      routineP = moduleN//':'//routineN

    INTEGER :: a, a_global, b, b_global, best_block_size, &
      best_integ_group_size, block_size, comm_exchange, comm_P, end_point, &
      handle, handle2, handle3, iiB, ij_counter, ij_counter_rec, &
      ij_counter_send, ij_index, ijk_counter, ijk_counter_send, ijk_index, &
      integ_group_size, irep, jjB, kkB, kkk, Lend_pos, lll, Lstart_pos, &
      max_ij_pairs, max_ijk, min_integ_group_size, my_B_size, my_B_size_beta, &
      my_B_virtual_end, my_B_virtual_end_beta, my_B_virtual_start, &
      my_B_virtual_start_beta, my_block_size, my_group_L_end, &
      my_group_L_size, my_group_L_size_orig, my_group_L_start, my_homo_beta, &
      my_i, my_ij_pairs
    INTEGER :: my_ijk, my_j, my_k, my_new_group_L_size, my_num_dgemm_call, &
      ngroup, num_IJ_blocks, num_integ_group, num_sing_ij, pos_integ_group, &
      proc_receive, proc_send, proc_shift, rec_B_size, rec_B_virtual_end, &
      rec_B_virtual_start, rec_block_size, rec_i, rec_ij_index, rec_j, &
      rec_L_size, send_B_size, send_B_virtual_end, send_B_virtual_start, &
      send_block_size, send_i, send_ij_index, send_ijk_index, send_j, send_k, &
      send_L_size, start_point, stat, sub_P_color, sub_sub_color, &
      total_ij_pairs, total_ijk, virtual, virtual_beta
    INTEGER, ALLOCATABLE, DIMENSION(:) :: integ_group_pos2color_sub, &
      new_sizes_array, num_ij_pairs, num_ijk, proc_map, proc_map_rep, &
      sizes_array_orig, sub_proc_map
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: ij_map, ijk_map
    INTEGER, ALLOCATABLE, DIMENSION(:, :, :) :: ranges_info_array
    LOGICAL                                  :: failure, my_alpha_beta_case, &
                                                my_open_shell_SS
    REAL(KIND=dp) :: actual_flop_rate, mem_for_aK, mem_for_comm, mem_for_iaK, &
      mem_for_rep, mem_min, mem_per_group, mem_real, my_flop_rate, &
      null_mat_rec(2,2,2), null_mat_send(2,2,2), P_ij_elem, sym_fac, t_end, &
      t_start
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: external_ab, external_i_aL, &
                                                local_ab, send_ab, t_ab
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: B_ia_Q, BI_C_rec, local_i_aL, &
                                                local_j_aL, local_k_aL, &
                                                Y_i_aP, Y_j_aP
    TYPE(cp_para_env_type), POINTER          :: para_env_exchange, &
                                                para_env_P, para_env_rep

    CALL timeset(routineN,handle)
    failure=.FALSE.

    my_open_shell_SS=.FALSE.
    IF(PRESENT(open_shell_SS)) my_open_shell_SS=open_shell_SS

    my_alpha_beta_case=.FALSE.
    IF(PRESENT(BIb_C_beta).AND.&
       PRESENT(ends_B_virtual_beta).AND.&
       PRESENT(sizes_B_virtual_beta).AND.&
       PRESENT(starts_B_virtual_beta).AND.&
       PRESENT(homo_beta).AND.&
       PRESENT(Eigenval_beta)) my_alpha_beta_case=.TRUE.

    virtual=nmo-homo
    IF(my_alpha_beta_case) virtual_beta=nmo-homo_beta

    CALL mp2_ri_get_sizes(mp2_env,para_env, para_env_sub,ends_array, ends_B_virtual,sizes_array, sizes_B_virtual,&
    starts_array, starts_B_virtual,homo, dimen_RI, unit_nr,color_sub,best_block_size,best_integ_group_size, block_size, &
    integ_group_size,min_integ_group_size, my_B_size,my_B_virtual_end, my_B_virtual_start, my_group_L_size, &
    my_group_L_start, my_group_L_end, ngroup, num_IJ_blocks, num_integ_group,pos_integ_group, virtual,my_alpha_beta_case, &
    my_open_shell_SS, actual_flop_rate, mem_for_aK, mem_for_comm, mem_for_iaK,mem_for_rep, mem_min, mem_per_group, mem_real)

    IF(my_alpha_beta_case) THEN
      my_B_virtual_start_beta=starts_B_virtual_beta(para_env_sub%mepos)
      my_B_virtual_end_beta=ends_B_virtual_beta(para_env_sub%mepos)
      my_B_size_beta=sizes_B_virtual_beta(para_env_sub%mepos)
      my_homo_beta=homo_beta
    ELSE
      my_B_virtual_start_beta=my_B_virtual_start
      my_B_virtual_end_beta=my_B_virtual_end
      my_B_size_beta=my_B_size
      my_homo_beta=homo
    END IF

    ! now create a group that contains all the proc that have the same virtual starting point
    ! in the integ group
    ! sub_sub_color=para_env_sub%mepos
    CALL mp2_ri_create_group(BIb_C,para_env,para_env_sub,homo,color_sub, ends_array, ends_B_virtual, &
                 sizes_array,sizes_B_virtual,starts_array,starts_B_virtual,calc_forces,error, error_sub,&
                 comm_exchange,integ_group_size, my_B_size,iiB, my_B_virtual_end,my_B_virtual_start, my_group_L_end,&
                 my_group_L_size, my_group_L_size_orig, my_group_L_start,my_new_group_L_size, ngroup, &
                 sub_sub_color,integ_group_pos2color_sub,new_sizes_array, proc_map, proc_map_rep, sizes_array_orig,&
                 sub_proc_map,ranges_info_array,para_env_exchange,para_env_rep,num_integ_group)

    ! *****************************************************************
    ! **********  REPLICATION-BLOCKED COMMUNICATION SCHEME  ***********
    ! *****************************************************************
    ! introduce block size, the number of occupied orbitals has to be a
    ! multiple of the block size

    ! Calculate the maximum number of ij pairs that have to be computed
    ! among groups
    CALL mp2_ri_communication(my_alpha_beta_case,total_ij_pairs,homo,homo_beta,num_IJ_blocks,&
               block_size,ngroup,ij_map,color_sub,my_ij_pairs,error_sub,my_open_shell_SS,unit_nr)
    
    ALLOCATE(num_ij_pairs(0:para_env_exchange%num_pe-1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
    num_ij_pairs=0
    num_ij_pairs(para_env_exchange%mepos)=my_ij_pairs
    CALL mp_sum(num_ij_pairs,para_env_exchange%group)

    max_ij_pairs=MAXVAL(num_ij_pairs)

    ! start real stuff
    CALL mp2_ri_allocate(local_ab, t_ab,mp2_env,homo,virtual,dimen_RI,my_B_size,&
                                block_size,my_B_size_beta,my_group_L_size,Y_i_aP,&
                                Y_j_aP,local_i_aL, local_j_aL,calc_forces,error_sub)

    CALL timeset(routineN//"_RI_loop",handle2)
    null_mat_rec=0.0_dp
    null_mat_send=0.0_dp
    Emp2=0.0_dp
    Emp2_Cou=0.0_dp
    Emp2_EX=0.0_dp
    my_num_dgemm_call=0
    my_flop_rate=0.0_dp
    DO ij_index=1, max_ij_pairs

      IF(ij_index<=my_ij_pairs) THEN
        ! We have work to do
        ij_counter=(ij_index-MIN(1,color_sub))*ngroup+color_sub
        my_i=ij_map(ij_counter,1)      
        my_j=ij_map(ij_counter,2) 
        my_block_size=ij_map(ij_counter,3)

        local_i_aL=0.0_dp
        ! local_i_aL(my_group_L_start:my_group_L_end,1:my_B_size)=BIb_C(1:my_group_L_size,1:my_B_size,my_i)
        DO irep=0, num_integ_group-1
          Lstart_pos=ranges_info_array(1,irep,para_env_exchange%mepos)
          Lend_pos=ranges_info_array(2,irep,para_env_exchange%mepos)
          start_point=ranges_info_array(3,irep,para_env_exchange%mepos)
          end_point=ranges_info_array(4,irep,para_env_exchange%mepos)

          local_i_aL(Lstart_pos:Lend_pos,1:my_B_size,1:my_block_size)=&
                                                  BIb_C(start_point:end_point,1:my_B_size,my_i:my_i+my_block_size-1)
        END DO

        local_j_aL=0.0_dp
        ! local_j_aL(my_group_L_start:my_group_L_end,1:my_B_size)=BIb_C(1:my_group_L_size,1:my_B_size,my_j)
        DO irep=0, num_integ_group-1
          Lstart_pos=ranges_info_array(1,irep,para_env_exchange%mepos)
          Lend_pos=ranges_info_array(2,irep,para_env_exchange%mepos)
          start_point=ranges_info_array(3,irep,para_env_exchange%mepos)
          end_point=ranges_info_array(4,irep,para_env_exchange%mepos)

          IF(.NOT.my_alpha_beta_case) THEN
            local_j_aL(Lstart_pos:Lend_pos,1:my_B_size,1:my_block_size)=&
                                                  BIb_C(start_point:end_point,1:my_B_size,my_j:my_j+my_block_size-1)
          ELSE
            local_j_aL(Lstart_pos:Lend_pos,1:my_B_size_beta,1:my_block_size)=&
                                                  BIb_C_beta(start_point:end_point,1:my_B_size_beta,my_j:my_j+my_block_size-1)
          END IF
        END DO

        ! collect data from other proc
        CALL timeset(routineN//"_comm",handle3)
        DO proc_shift=1, para_env_exchange%num_pe-1
          proc_send=proc_map(para_env_exchange%mepos+proc_shift)
          proc_receive=proc_map(para_env_exchange%mepos-proc_shift)

          send_ij_index=num_ij_pairs(proc_send)

          rec_L_size=sizes_array(proc_receive)
          ALLOCATE(BI_C_rec(rec_L_size,MAX(my_B_size,my_B_size_beta),my_block_size),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
 
          IF(ij_index<=send_ij_index) THEN
            ! ij_counter_send=(ij_index-MIN(1,proc_send))*ngroup+proc_send
            ij_counter_send=(ij_index-MIN(1,integ_group_pos2color_sub(proc_send)))*ngroup+integ_group_pos2color_sub(proc_send)
            send_i=ij_map(ij_counter_send,1)
            send_j=ij_map(ij_counter_send,2)
            send_block_size=ij_map(ij_counter_send,3)

            ! occupied i
            BI_C_rec=0.0_dp
            CALL  mp_sendrecv(BIb_C(1:my_group_L_size,1:my_B_size,send_i:send_i+send_block_size-1),proc_send,&
                              BI_C_rec(1:rec_L_size,1:my_B_size,1:my_block_size),proc_receive,&
                              para_env_exchange%group)
            DO irep=0, num_integ_group-1
              Lstart_pos=ranges_info_array(1,irep,proc_receive)
              Lend_pos=ranges_info_array(2,irep,proc_receive)
              start_point=ranges_info_array(3,irep,proc_receive)
              end_point=ranges_info_array(4,irep,proc_receive)
 
              local_i_aL(Lstart_pos:Lend_pos,1:my_B_size,1:my_block_size)=&
                                                  BI_C_rec(start_point:end_point,1:my_B_size,1:my_block_size)

            END DO

            ! occupied j
            BI_C_rec=0.0_dp
            IF(.NOT.my_alpha_beta_case) THEN
              CALL  mp_sendrecv(BIb_C(1:my_group_L_size,1:my_B_size,send_j:send_j+send_block_size-1),proc_send,&
                                BI_C_rec(1:rec_L_size,1:my_B_size,1:my_block_size),proc_receive,&
                                para_env_exchange%group)
            ELSE
              CALL  mp_sendrecv(BIb_C_beta(1:my_group_L_size,1:my_B_size_beta,send_j:send_j+send_block_size-1),proc_send,&
                                BI_C_rec(1:rec_L_size,1:my_B_size_beta,1:my_block_size),proc_receive,&
                                para_env_exchange%group)
            END IF

            DO irep=0, num_integ_group-1
              Lstart_pos=ranges_info_array(1,irep,proc_receive)
              Lend_pos=ranges_info_array(2,irep,proc_receive)
              start_point=ranges_info_array(3,irep,proc_receive)
              end_point=ranges_info_array(4,irep,proc_receive)
 
              IF(.NOT.my_alpha_beta_case) THEN
                local_j_aL(Lstart_pos:Lend_pos,1:my_B_size,1:my_block_size)=&
                                                  BI_C_rec(start_point:end_point,1:my_B_size,1:my_block_size)
              ELSE
                local_j_aL(Lstart_pos:Lend_pos,1:my_B_size_beta,1:my_block_size)=&
                                                  BI_C_rec(start_point:end_point,1:my_B_size_beta,1:my_block_size)
              END IF

            END DO

          ELSE
            ! we send the null matrix while we know that we have to receive something

            ! occupied i
            BI_C_rec=0.0_dp
            CALL  mp_sendrecv(null_mat_send,proc_send,&
                              BI_C_rec(1:rec_L_size,1:my_B_size,1:my_block_size),proc_receive,&
                              para_env_exchange%group)

            DO irep=0, num_integ_group-1
              Lstart_pos=ranges_info_array(1,irep,proc_receive)
              Lend_pos=ranges_info_array(2,irep,proc_receive)
              start_point=ranges_info_array(3,irep,proc_receive)
              end_point=ranges_info_array(4,irep,proc_receive)

              local_i_aL(Lstart_pos:Lend_pos,1:my_B_size,1:my_block_size)=&
                                                  BI_C_rec(start_point:end_point,1:my_B_size,1:my_block_size)

            END DO

            ! occupied j
            BI_C_rec=0.0_dp
            IF(.NOT.my_alpha_beta_case) THEN
              CALL  mp_sendrecv(null_mat_send,proc_send,&
                                BI_C_rec(1:rec_L_size,1:my_B_size,1:my_block_size),proc_receive,&
                                para_env_exchange%group)
            ELSE
              CALL  mp_sendrecv(null_mat_send,proc_send,&
                                BI_C_rec(1:rec_L_size,1:my_B_size_beta,1:my_block_size),proc_receive,&
                                para_env_exchange%group)
            END IF
            DO irep=0, num_integ_group-1
              Lstart_pos=ranges_info_array(1,irep,proc_receive)
              Lend_pos=ranges_info_array(2,irep,proc_receive)
              start_point=ranges_info_array(3,irep,proc_receive)
              end_point=ranges_info_array(4,irep,proc_receive)
        
              IF(.NOT.my_alpha_beta_case) THEN
                local_j_aL(Lstart_pos:Lend_pos,1:my_B_size,1:my_block_size)=&
                                                  BI_C_rec(start_point:end_point,1:my_B_size,1:my_block_size)
              ELSE
                local_j_aL(Lstart_pos:Lend_pos,1:my_B_size_beta,1:my_block_size)=&
                                                  BI_C_rec(start_point:end_point,1:my_B_size_beta,1:my_block_size)
              END IF

            END DO

          END IF

          DEALLOCATE(BI_C_rec)

        END DO
        CALL timestop(handle3)
 
        IF(calc_forces) THEN
          Y_i_aP=0.0_dp
          Y_j_aP=0.0_dp
        END IF

        ! loop over the block elements
        DO iiB=1, my_block_size
          DO jjB=1, my_block_size
            CALL timeset(routineN//"_expansion",handle3)
            ! calculate the integrals (ia|jb) strating from my local data ...
            local_ab=0.0_dp
            t_start=m_walltime()
            CALL dgemm('T','N',my_B_size,my_B_size_beta,dimen_RI,1.0_dp,&
                       local_i_aL(:,:,iiB),dimen_RI,local_j_aL(:,:,jjB),dimen_RI,&
                       0.0_dp,local_ab(my_B_virtual_start:my_B_virtual_end,1:my_B_size_beta),my_B_size)
            t_end=m_walltime()
            actual_flop_rate=2.0_dp*my_B_size*my_B_size_beta*REAL(dimen_RI,KIND=dp)/(MAX(0.01_dp,t_end-t_start))
            my_flop_rate=my_flop_rate+actual_flop_rate
            my_num_dgemm_call=my_num_dgemm_call+1
            ! ... and from the other of my subgroup
            DO proc_shift=1, para_env_sub%num_pe-1
              proc_send=sub_proc_map(para_env_sub%mepos+proc_shift)
              proc_receive=sub_proc_map(para_env_sub%mepos-proc_shift)

              rec_B_size=sizes_B_virtual(proc_receive)
              rec_B_virtual_end=ends_B_virtual(proc_receive)
              rec_B_virtual_start=starts_B_virtual(proc_receive)
 
              ALLOCATE(external_i_aL(dimen_RI,rec_B_size),STAT=stat)
              CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
              external_i_aL=0.0_dp

              CALL  mp_sendrecv(local_i_aL(:,:,iiB),proc_send,&
                                external_i_aL,proc_receive,&
                                para_env_sub%group)

              ! local_ab(rec_B_virtual_start:rec_B_virtual_end,1:my_B_size)=MATMUL(TRANSPOSE(external_i_aL),local_j_aL)
              t_start=m_walltime()
              CALL dgemm('T','N',rec_B_size,my_B_size_beta,dimen_RI,1.0_dp,&
                         external_i_aL,dimen_RI,local_j_aL(:,:,jjB),dimen_RI,&
                         0.0_dp,local_ab(rec_B_virtual_start:rec_B_virtual_end,1:my_B_size_beta),rec_B_size)
              t_end=m_walltime()
              actual_flop_rate=2.0_dp*rec_B_size*my_B_size_beta*REAL(dimen_RI,KIND=dp)/(MAX(0.01_dp,t_end-t_start))
              my_flop_rate=my_flop_rate+actual_flop_rate
              my_num_dgemm_call=my_num_dgemm_call+1

              DEALLOCATE(external_i_aL)
            END DO
            CALL timestop(handle3)

            CALL timeset(routineN//"_ener",handle3)
            ! calculate coulomb only MP2
            sym_fac=2.0_dp
            IF(my_i==my_j) sym_fac=1.0_dp
            IF(.NOT.my_alpha_beta_case) THEN
              ! IF(my_open_shell_SS) sym_fac=sym_fac/2.0_dp
              DO b=1, my_B_size
                b_global=b+my_B_virtual_start-1
                DO a=1, virtual
                  Emp2_Cou=Emp2_Cou-sym_fac*2.0_dp*local_ab(a,b)**2/&
                            (Eigenval(homo+a)+Eigenval(homo+b_global)-Eigenval(my_i+iiB-1)-Eigenval(my_j+jjB-1))
                END DO
              END DO
            ELSE
              DO b=1, my_B_size_beta
                b_global=b+my_B_virtual_start_beta-1
                DO a=1, virtual
                  Emp2_Cou=Emp2_Cou-local_ab(a,b)**2/&
                            (Eigenval(homo+a)+Eigenval_beta(homo_beta+b_global)-Eigenval(my_i+iiB-1)-Eigenval_beta(my_j+jjB-1))
                END DO
              END DO 
            END IF
            
            IF(calc_ex) THEN
              ! contract integrals with orbital energies for exchange MP2 energy
              ! starting with local ...
              ! IF(my_open_shell_SS) sym_fac=sym_fac*2.0_dp
              IF(calc_forces) t_ab=0.0_dp
              DO b=1, my_B_size
                b_global=b+my_B_virtual_start-1
                DO a=1, my_B_size
                  a_global=a+my_B_virtual_start-1
                  Emp2_Ex=Emp2_Ex+sym_fac*local_ab(a_global,b)*local_ab(b_global,a)/&
                            (Eigenval(homo+a_global)+Eigenval(homo+b_global)-Eigenval(my_i+iiB-1)-Eigenval(my_j+jjB-1))
                  IF(calc_forces) t_ab(a_global,b)=-(2.0_dp*local_ab(a_global,b)-local_ab(b_global,a))/&
                                  (Eigenval(homo+a_global)+Eigenval(homo+b_global)-Eigenval(my_i+iiB-1)-Eigenval(my_j+jjB-1))
                END DO
              END DO
              ! ... and then with external data
              DO proc_shift=1, para_env_sub%num_pe-1
                proc_send=sub_proc_map(para_env_sub%mepos+proc_shift)
                proc_receive=sub_proc_map(para_env_sub%mepos-proc_shift)

                rec_B_size=sizes_B_virtual(proc_receive)
                rec_B_virtual_end=ends_B_virtual(proc_receive)
                rec_B_virtual_start=starts_B_virtual(proc_receive)

                send_B_size=sizes_B_virtual(proc_send)
                send_B_virtual_end=ends_B_virtual(proc_send)
                send_B_virtual_start=starts_B_virtual(proc_send)

                ALLOCATE(external_ab(my_B_size,rec_B_size),STAT=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
                external_ab=0.0_dp

                CALL  mp_sendrecv(local_ab(send_B_virtual_start:send_B_virtual_end,1:my_B_size),proc_send,&
                                  external_ab(1:my_B_size,1:rec_B_size),proc_receive,&
                                  para_env_sub%group)

                DO b=1, my_B_size
                  b_global=b+my_B_virtual_start-1
                  DO a=1, rec_B_size
                    a_global=a+rec_B_virtual_start-1
                    Emp2_Ex=Emp2_Ex+sym_fac*local_ab(a_global,b)*external_ab(b,a)/&
                            (Eigenval(homo+a_global)+Eigenval(homo+b_global)-Eigenval(my_i+iiB-1)-Eigenval(my_j+jjB-1))
                    IF(calc_forces) t_ab(a_global,b)=-(2.0_dp*local_ab(a_global,b)-external_ab(b,a))/&
                                    (Eigenval(homo+a_global)+Eigenval(homo+b_global)-Eigenval(my_i+iiB-1)-Eigenval(my_j+jjB-1))
                  END DO
                END DO

                DEALLOCATE(external_ab)
              END DO
            END IF
            CALL timestop(handle3)

            IF(calc_forces) THEN
              ! update P_ab, Gamma_P_ia
              CALL mp2_update_P_gamma(mp2_env,para_env_sub,ends_B_virtual,sizes_B_virtual, starts_B_virtual,&
              Eigenval, homo, dimen_RI, error_sub, b_global, iiB, jjB, my_B_size, &
              my_B_virtual_end, my_B_virtual_start, my_i, my_j, proc_receive, proc_send,proc_shift, &
              rec_B_size, rec_B_virtual_end, rec_B_virtual_start, send_B_size, send_B_virtual_end, & 
              send_B_virtual_start, virtual, proc_map, sub_proc_map,external_ab, local_ab,&
              send_ab, t_ab,local_i_aL, local_j_aL, Y_i_aP, Y_j_aP)

            END IF

          END DO ! jjB
        END DO ! iiB

      ELSE
        ! No work to do and we know that we have to receive nothing, but send something
        ! send data to other proc
        DO proc_shift=1, para_env_exchange%num_pe-1
          proc_send=proc_map(para_env_exchange%mepos+proc_shift)
          proc_receive=proc_map(para_env_exchange%mepos-proc_shift)

          send_ij_index=num_ij_pairs(proc_send)

          IF(ij_index<=send_ij_index) THEN
            ! something to send
            ! ij_counter_send=(ij_index-MIN(1,proc_send))*ngroup+proc_send
            ij_counter_send=(ij_index-MIN(1,integ_group_pos2color_sub(proc_send)))*ngroup+integ_group_pos2color_sub(proc_send)
            send_i=ij_map(ij_counter_send,1)
            send_j=ij_map(ij_counter_send,2)
            send_block_size=ij_map(ij_counter_send,3)
         
            ! occupied i
            CALL  mp_sendrecv(BIb_C(1:my_group_L_size,1:my_B_size,send_i:send_i+send_block_size-1),proc_send,&
                              null_mat_rec,proc_receive,&
                              para_env_exchange%group)
            ! occupied j
            IF(.NOT.my_alpha_beta_case) THEN
              CALL  mp_sendrecv(BIb_C(1:my_group_L_size,1:my_B_size,send_j:send_j+send_block_size-1),proc_send,&
                                null_mat_rec,proc_receive,&
                                para_env_exchange%group)            
            ELSE
              CALL  mp_sendrecv(BIb_C_beta(1:my_group_L_size,1:my_B_size_beta,send_j:send_j+send_block_size-1),proc_send,&
                                null_mat_rec,proc_receive,&
                                para_env_exchange%group)
            END IF

          ELSE
            ! nothing to send 
            ! occupied i
            CALL  mp_sendrecv(null_mat_send,proc_send,&
                              null_mat_rec,proc_receive,&
                              para_env_exchange%group)
            ! occupied j
            CALL  mp_sendrecv(null_mat_send,proc_send,&
                              null_mat_rec,proc_receive,&
                              para_env_exchange%group)

          END IF
        END DO
      END IF

      ! redistribute gamma
      IF(calc_forces) THEN
        CALL   mp2_redistribute_gamma(mp2_env,error_sub,end_point, iiB, jjb, ij_counter_rec, &
               ij_index, irep, kkk, Lend_pos, lll, Lstart_pos, my_B_size, my_block_size, my_group_L_size, my_i,&
               my_ij_pairs, my_j, ngroup, num_integ_group, proc_receive, proc_send, proc_shift, rec_block_size,&
               rec_i, rec_ij_index, rec_j, send_L_size, start_point, integ_group_pos2color_sub,num_ij_pairs, &
               proc_map, ij_map, ranges_info_array, Y_i_aP, Y_j_aP, para_env_exchange, &
               null_mat_rec, null_mat_send, sizes_array )
      END IF

    END DO
    CALL timestop(handle2)   

    DEALLOCATE(local_i_aL)
    DEALLOCATE(local_j_aL)
    DEALLOCATE(ij_map)
    DEALLOCATE(num_ij_pairs)

    IF(calc_forces) THEN
      DEALLOCATE(Y_i_aP, Y_j_aP, t_ab, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)

      CALL timeset(routineN//"_ij_sing",handle2)
      ! here we check if there are almost degenerate ij
      ! pairs and we update P_ij with these contribution.
      ! If all pairs are degenerate with each other this step will scale O(N^6),
      ! if the number of degenerate pairs scales linearly with the system size
      ! this step will scale O(N^5).
      ! Start counting the number of almost degenerate ij pairs according
      ! to eps_canonical
      num_sing_ij=0
      DO iiB=1, homo
        ! diagonal elements already updated
        DO jjB=iiB+1, homo
          IF(ABS(Eigenval(jjB)-Eigenval(iiB))<mp2_env%ri_mp2%eps_canonical) num_sing_ij=num_sing_ij+1
        END DO
      END DO

      IF (unit_nr>0) THEN
         WRITE (UNIT=unit_nr,FMT="(T3,A,T75,i6)")&
                                 "MO_INFO| Number of ij pairs below EPS_CANONICAL:",num_sing_ij
      END IF

      ! total number of elemets that have to be computed
      total_ijk=homo*num_sing_ij
      ALLOCATE(ijk_map(total_ijk,3),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
      ijk_map=0
      my_ijk=0
      ijk_counter=0
      DO iiB=1, homo
        ! diagonal elements already updated
        DO jjB=iiB+1, homo
          IF(ABS(Eigenval(jjB)-Eigenval(iiB))>=mp2_env%ri_mp2%eps_canonical) CYCLE
          DO kkB=1, homo
            ijk_counter=ijk_counter+1
            ijk_map(ijk_counter,1)=iiB
            ijk_map(ijk_counter,2)=jjB
            ijk_map(ijk_counter,3)=kkB
            IF (MOD(ijk_counter,ngroup)==color_sub) my_ijk=my_ijk+1
          END DO
        END DO
      END DO

      ALLOCATE(num_ijk(0:para_env_exchange%num_pe-1),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
      num_ijk=0
      num_ijk(para_env_exchange%mepos)=my_ijk
      CALL mp_sum(num_ijk,para_env_exchange%group)
      max_ijk=MAXVAL(num_ijk)

      ALLOCATE(local_i_aL(dimen_RI,my_B_size,1),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
      ALLOCATE(local_j_aL(dimen_RI,my_B_size,1),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
      ALLOCATE(local_k_aL(dimen_RI,my_B_size,1),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
      ALLOCATE(t_ab(virtual,my_B_size),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)

      DO ijk_index=1, max_ijk

        IF(ijk_index<=my_ijk) THEN
          ! work to be done
          ijk_counter=(ijk_index-MIN(1,color_sub))*ngroup+color_sub
          my_i=ijk_map(ijk_counter,1)
          my_j=ijk_map(ijk_counter,2)
          my_k=ijk_map(ijk_counter,3)

          local_i_aL=0.0_dp
          local_j_al=0.0_dp
          local_k_al=0.0_dp
          DO irep=0, num_integ_group-1
            Lstart_pos=ranges_info_array(1,irep,para_env_exchange%mepos)
            Lend_pos=ranges_info_array(2,irep,para_env_exchange%mepos)
            start_point=ranges_info_array(3,irep,para_env_exchange%mepos)
            end_point=ranges_info_array(4,irep,para_env_exchange%mepos)

            local_i_aL(Lstart_pos:Lend_pos,1:my_B_size,1)=BIb_C(start_point:end_point,1:my_B_size,my_i)
            local_j_aL(Lstart_pos:Lend_pos,1:my_B_size,1)=BIb_C(start_point:end_point,1:my_B_size,my_j)
            local_k_aL(Lstart_pos:Lend_pos,1:my_B_size,1)=BIb_C(start_point:end_point,1:my_B_size,my_k)
          END DO

          DO proc_shift=1, para_env_exchange%num_pe-1
            proc_send=proc_map(para_env_exchange%mepos+proc_shift)
            proc_receive=proc_map(para_env_exchange%mepos-proc_shift)

            send_ijk_index=num_ijk(proc_send)

            rec_L_size=sizes_array(proc_receive)
            ALLOCATE(BI_C_rec(rec_L_size,my_B_size,1),STAT=stat)
            CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)

            IF(ijk_index<=send_ijk_index) THEN
              ! something to send
              ijk_counter_send=(ijk_index-MIN(1,integ_group_pos2color_sub(proc_send)))*ngroup+integ_group_pos2color_sub(proc_send)
              send_i=ijk_map(ijk_counter_send,1)
              send_j=ijk_map(ijk_counter_send,2)
              send_k=ijk_map(ijk_counter_send,3)
            END IF

            ! occupied i
            BI_C_rec=0.0_dp
            IF(ijk_index<=send_ijk_index) THEN
              CALL  mp_sendrecv(BIb_C(1:my_group_L_size,1:my_B_size,send_i),proc_send,&
                                BI_C_rec(1:rec_L_size,1:my_B_size,1),proc_receive,&
                                para_env_exchange%group)
            ELSE
              ! nothing to send
              CALL  mp_sendrecv(null_mat_send,proc_send,&
                                BI_C_rec(1:rec_L_size,1:my_B_size,1:1),proc_receive,&
                                para_env_exchange%group)
            END IF
            DO irep=0, num_integ_group-1
              Lstart_pos=ranges_info_array(1,irep,proc_receive)
              Lend_pos=ranges_info_array(2,irep,proc_receive)
              start_point=ranges_info_array(3,irep,proc_receive)
              end_point=ranges_info_array(4,irep,proc_receive)

              local_i_aL(Lstart_pos:Lend_pos,1:my_B_size,1)=BI_C_rec(start_point:end_point,1:my_B_size,1)
            END DO

            ! occupied j
            BI_C_rec=0.0_dp
            IF(ijk_index<=send_ijk_index) THEN
              CALL  mp_sendrecv(BIb_C(1:my_group_L_size,1:my_B_size,send_j),proc_send,&
                                BI_C_rec(1:rec_L_size,1:my_B_size,1),proc_receive,&
                                para_env_exchange%group)
            ELSE
              ! nothing to send
              CALL  mp_sendrecv(null_mat_send,proc_send,&
                                BI_C_rec(1:rec_L_size,1:my_B_size,1:1),proc_receive,&
                                para_env_exchange%group)
            END IF
            DO irep=0, num_integ_group-1
              Lstart_pos=ranges_info_array(1,irep,proc_receive)
              Lend_pos=ranges_info_array(2,irep,proc_receive)
              start_point=ranges_info_array(3,irep,proc_receive)
              end_point=ranges_info_array(4,irep,proc_receive)

              local_j_aL(Lstart_pos:Lend_pos,1:my_B_size,1)=BI_C_rec(start_point:end_point,1:my_B_size,1)
            END DO

            ! occupied k
            BI_C_rec=0.0_dp
            IF(ijk_index<=send_ijk_index) THEN
              CALL  mp_sendrecv(BIb_C(1:my_group_L_size,1:my_B_size,send_k),proc_send,&
                                BI_C_rec(1:rec_L_size,1:my_B_size,1),proc_receive,&
                                para_env_exchange%group)
            ELSE
              ! nothing to send
              CALL  mp_sendrecv(null_mat_send,proc_send,&
                                BI_C_rec(1:rec_L_size,1:my_B_size,1:1),proc_receive,&
                                para_env_exchange%group)
            END IF
            DO irep=0, num_integ_group-1
              Lstart_pos=ranges_info_array(1,irep,proc_receive)
              Lend_pos=ranges_info_array(2,irep,proc_receive)
              start_point=ranges_info_array(3,irep,proc_receive)
              end_point=ranges_info_array(4,irep,proc_receive)

              local_k_aL(Lstart_pos:Lend_pos,1:my_B_size,1)=BI_C_rec(start_point:end_point,1:my_B_size,1)
            END DO

            DEALLOCATE(BI_C_rec)
          END DO

          ! expand integrals
          CALL timeset(routineN//"_exp_ik",handle3)
          local_ab=0.0_dp
          CALL dgemm('T','N',my_B_size,my_B_size,dimen_RI,1.0_dp,&
                         local_i_aL(:,:,1),dimen_RI,local_k_aL(:,:,1),dimen_RI,&
                         0.0_dp,local_ab(my_B_virtual_start:my_B_virtual_end,1:my_B_size),my_B_size)
          DO proc_shift=1, para_env_sub%num_pe-1
            proc_send=sub_proc_map(para_env_sub%mepos+proc_shift)
            proc_receive=sub_proc_map(para_env_sub%mepos-proc_shift)

            rec_B_size=sizes_B_virtual(proc_receive)
            rec_B_virtual_end=ends_B_virtual(proc_receive)
            rec_B_virtual_start=starts_B_virtual(proc_receive)

            ALLOCATE(external_i_aL(dimen_RI,rec_B_size),STAT=stat)
            CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
            external_i_aL=0.0_dp

            CALL  mp_sendrecv(local_i_aL(:,:,1),proc_send,&
                              external_i_aL,proc_receive,&
                              para_env_sub%group)

            CALL dgemm('T','N',rec_B_size,my_B_size,dimen_RI,1.0_dp,&
                           external_i_aL,dimen_RI,local_k_aL(:,:,1),dimen_RI,&
                           0.0_dp,local_ab(rec_B_virtual_start:rec_B_virtual_end,1:my_B_size),rec_B_size)

            DEALLOCATE(external_i_aL)
          END DO
          CALL timestop(handle3)

          CALL timeset(routineN//"_tab",handle3)
          t_ab=0.0_dp
          DO b=1, my_B_size
            b_global=b+my_B_virtual_start-1
            DO a=1, my_B_size
              a_global=a+my_B_virtual_start-1
              t_ab(a_global,b)=(2.0_dp*local_ab(a_global,b)-local_ab(b_global,a))/&
                               (Eigenval(my_i)+Eigenval(my_k)-Eigenval(homo+a_global)-Eigenval(homo+b_global))
            END DO
          END DO
          DO proc_shift=1, para_env_sub%num_pe-1
            proc_send=sub_proc_map(para_env_sub%mepos+proc_shift)
            proc_receive=sub_proc_map(para_env_sub%mepos-proc_shift)

            rec_B_size=sizes_B_virtual(proc_receive)
            rec_B_virtual_end=ends_B_virtual(proc_receive)
            rec_B_virtual_start=starts_B_virtual(proc_receive)

            send_B_size=sizes_B_virtual(proc_send)
            send_B_virtual_end=ends_B_virtual(proc_send)
            send_B_virtual_start=starts_B_virtual(proc_send)

            ALLOCATE(external_ab(my_B_size,rec_B_size),STAT=stat)
            CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
            external_ab=0.0_dp

            CALL  mp_sendrecv(local_ab(send_B_virtual_start:send_B_virtual_end,1:my_B_size),proc_send,&
                              external_ab(1:my_B_size,1:rec_B_size),proc_receive,&
                              para_env_sub%group)

            DO b=1, my_B_size
              b_global=b+my_B_virtual_start-1
              DO a=1, rec_B_size
                a_global=a+rec_B_virtual_start-1
                t_ab(a_global,b)=(2.0_dp*local_ab(a_global,b)-external_ab(b,a))/&
                                 (Eigenval(my_i)+Eigenval(my_k)-Eigenval(homo+a_global)-Eigenval(homo+b_global))
              END DO
            END DO

            DEALLOCATE(external_ab)
          END DO
          CALL timestop(handle3)

          CALL timeset(routineN//"_exp_jk",handle3)
          local_ab=0.0_dp
          CALL dgemm('T','N',my_B_size,my_B_size,dimen_RI,1.0_dp,&
                         local_j_aL(:,:,1),dimen_RI,local_k_aL(:,:,1),dimen_RI,&
                         0.0_dp,local_ab(my_B_virtual_start:my_B_virtual_end,1:my_B_size),my_B_size)
          DO proc_shift=1, para_env_sub%num_pe-1
            proc_send=sub_proc_map(para_env_sub%mepos+proc_shift)
            proc_receive=sub_proc_map(para_env_sub%mepos-proc_shift)

            rec_B_size=sizes_B_virtual(proc_receive)
            rec_B_virtual_end=ends_B_virtual(proc_receive)
            rec_B_virtual_start=starts_B_virtual(proc_receive)

            ALLOCATE(external_i_aL(dimen_RI,rec_B_size),STAT=stat)
            CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
            external_i_aL=0.0_dp

            CALL  mp_sendrecv(local_j_aL(:,:,1),proc_send,&
                              external_i_aL,proc_receive,&
                              para_env_sub%group)

            CALL dgemm('T','N',rec_B_size,my_B_size,dimen_RI,1.0_dp,&
                           external_i_aL,dimen_RI,local_k_aL(:,:,1),dimen_RI,&
                           0.0_dp,local_ab(rec_B_virtual_start:rec_B_virtual_end,1:my_B_size),rec_B_size)

            DEALLOCATE(external_i_aL)
          END DO
          CALL timestop(handle3)

          CALL timeset(routineN//"_Pij",handle3)
          DO b=1, my_B_size
            b_global=b+my_B_virtual_start-1
            DO a=1, virtual
              local_ab(a,b)=local_ab(a,b)/&
                            (Eigenval(my_j)+Eigenval(my_k)-Eigenval(homo+a)-Eigenval(homo+b_global))
            END DO
          END DO
          P_ij_elem=SUM(local_ab*t_ab)*2.0_dp
          mp2_env%ri_grad%P_ij(my_i,my_j)=mp2_env%ri_grad%P_ij(my_i,my_j)-P_ij_elem
          mp2_env%ri_grad%P_ij(my_j,my_i)=mp2_env%ri_grad%P_ij(my_j,my_i)-P_ij_elem
          CALL timestop(handle3)

        ELSE
          ! no work to be done, possible messeges to be exchanged
          DO proc_shift=1, para_env_exchange%num_pe-1
            proc_send=proc_map(para_env_exchange%mepos+proc_shift)
            proc_receive=proc_map(para_env_exchange%mepos-proc_shift)

            send_ijk_index=num_ijk(proc_send)

            IF(ijk_index<=send_ijk_index) THEN
              ! somethig to send 
              ijk_counter_send=(ijk_index-MIN(1,integ_group_pos2color_sub(proc_send)))*ngroup+integ_group_pos2color_sub(proc_send)
              send_i=ijk_map(ijk_counter_send,1)
              send_j=ijk_map(ijk_counter_send,2)
              send_k=ijk_map(ijk_counter_send,3)
              ! occupied i
              CALL  mp_sendrecv(BIb_C(1:my_group_L_size,1:my_B_size,send_i:send_i),proc_send,&
                                null_mat_rec,proc_receive,&
                                para_env_exchange%group)
              ! occupied j
              CALL  mp_sendrecv(BIb_C(1:my_group_L_size,1:my_B_size,send_j:send_j),proc_send,&
                                null_mat_rec,proc_receive,&
                                para_env_exchange%group)
              ! occupied k
              CALL  mp_sendrecv(BIb_C(1:my_group_L_size,1:my_B_size,send_k:send_k),proc_send,&
                                null_mat_rec,proc_receive,&
                                para_env_exchange%group)
            ELSE
              ! nothing to send 
              ! occupied i
              CALL  mp_sendrecv(null_mat_send,proc_send,&
                                null_mat_rec,proc_receive,&
                                para_env_exchange%group)
              ! occupied j
              CALL  mp_sendrecv(null_mat_send,proc_send,&
                                null_mat_rec,proc_receive,&
                                para_env_exchange%group)
              ! occupied k
              CALL  mp_sendrecv(null_mat_send,proc_send,&
                                null_mat_rec,proc_receive,&
                                para_env_exchange%group)
            END IF

          END DO ! proc loop
        END IF
      END DO ! ijk_index loop
      CALL timestop(handle2)

      DEALLOCATE(local_i_aL, local_j_aL, local_k_aL)
      DEALLOCATE(ijk_map)
      DEALLOCATE(num_ijk)
      DEALLOCATE(t_ab)

    END IF

    DEALLOCATE(integ_group_pos2color_sub)
    DEALLOCATE(local_ab)

    CALL mp_sum(Emp2_Cou,para_env%group)
    CALL mp_sum(Emp2_Ex,para_env%group)
  
    IF(calc_forces) THEN
      ! sum P_ab
      mp2_env%ri_grad%P_ab(:,:) = mp2_env%ri_grad%P_ab(:,:)*2.0_dp
      sub_P_color=para_env_sub%mepos
      CALL mp_comm_split_direct(para_env%group,comm_P,sub_P_color)
      NULLIFY(para_env_P)
      CALL cp_para_env_create(para_env_P,comm_P,error=error_sub)
      CALL mp_sum(mp2_env%ri_grad%P_ab,para_env_P%group)
      ! release para_env_P
      CALL cp_para_env_release(para_env_P,error=error_sub)
     
      ! sum P_ij (later)
      ! mp2_env%ri_grad%P_ij=mp2_env%ri_grad%P_ij*2.0_dp
      ! CALL mp_sum(mp2_env%ri_grad%P_ij,para_env%group)
      ! WRITE(*,*) mp2_env%ri_mp2%eps_canonical

      ! XXXXXXXXXXXXXX
      ! Write P_ab
      ! ALLOCATE(external_ab(virtual,virtual))
      ! external_ab=0.0_dp
      ! external_ab(my_B_virtual_start:my_B_virtual_end,1:virtual)=mp2_env%ri_grad%P_ab
      ! CALL mp_sum(external_ab,para_env_sub%group)
      ! IF(para_env%mepos==0) CALL write_array(external_ab)
      ! DEALLOCATE(external_ab) 
      ! ! Write P_ij
      ! IF(para_env%mepos==0) CALL write_array(mp2_env%ri_grad%P_ij)
      ! XXXXXXXXXXXXXX

      ! recover original information (before replication)
      DEALLOCATE(sizes_array)
      iiB=SIZE(sizes_array_orig)
      ALLOCATE(sizes_array(0:iiB-1),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
      sizes_array(:)=sizes_array_orig
      DEALLOCATE(sizes_array_orig)

      ! make a copy of the original integrals (ia|Q)
      my_group_L_size=my_group_L_size_orig
      ALLOCATE(B_ia_Q(homo,my_B_size,my_group_L_size))
      B_ia_Q=0.0_dp
      DO jjB=1, homo
        DO iiB=1, my_B_size
          B_ia_Q(jjB,iiB,1:my_group_L_size)=BIb_C(1:my_group_L_size,iiB,jjB)
        END DO
      END DO
      DEALLOCATE(BIb_C)

      ! sum Gamma and dereplicate
      ALLOCATE(BIb_C(homo,my_B_size,my_group_L_size))
      DO proc_shift=1, para_env_rep%num_pe-1
        ! invert order
        proc_send=proc_map_rep(para_env_rep%mepos-proc_shift)
        proc_receive=proc_map_rep(para_env_rep%mepos+proc_shift)

        start_point=ranges_info_array(3,proc_shift,para_env_exchange%mepos)
        end_point=ranges_info_array(4,proc_shift,para_env_exchange%mepos)

        BIb_C=0.0_dp

        CALL  mp_sendrecv(mp2_env%ri_grad%Gamma_P_ia(1:homo,1:my_B_size,start_point:end_point),proc_send,&
                          BIb_C,proc_receive,&
                          para_env_rep%group)

        mp2_env%ri_grad%Gamma_P_ia(1:homo,1:my_B_size,1:my_group_L_size)=&
                                    mp2_env%ri_grad%Gamma_P_ia(1:homo,1:my_B_size,1:my_group_L_size)+BIb_C

      END DO
      BIb_C(:,:,:)=mp2_env%ri_grad%Gamma_P_ia(1:homo,1:my_B_size,1:my_group_L_size)
      DEALLOCATE(mp2_env%ri_grad%Gamma_P_ia)
      ALLOCATE(mp2_env%ri_grad%Gamma_P_ia(homo,my_B_size,my_group_L_size),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
      mp2_env%ri_grad%Gamma_P_ia(:,:,:)=BIb_C
      DEALLOCATE(BIb_C)

      CALL complete_gamma(mp2_env,B_ia_Q,dimen_RI,homo,virtual,para_env,para_env_sub,ngroup,color_sub,&
                          my_group_L_size,my_group_L_start,my_group_L_end,&
                          my_B_size,my_B_virtual_start,my_B_virtual_end,&
                          ends_array,ends_B_virtual,sizes_array,sizes_B_virtual,starts_array,starts_B_virtual,&     
                          sub_proc_map,error,error_sub)                   

      ! ALLOCATE(Y_i_aP(homo,virtual,dimen_RI))
      ! Y_i_aP=0.0_dp
      ! Y_i_aP(1:homo,my_B_virtual_start:my_B_virtual_end,my_group_L_start:my_group_L_end)=mp2_env%ri_grad%Gamma_P_ia
      ! CALL mp_sum(Y_i_aP,para_env%group)
      ! IF(para_env%mepos==0) THEN
      !   DO iiB=1, dimen_RI
      !     WRITE(*,*) iiB
      !     CALL write_array(Y_i_aP(1:homo,1:virtual,iiB))
      !     WRITE(*,*) 
      !   END DO
      ! END IF
      ! DEALLOCATE(Y_i_aP)

    END IF

    Emp2=Emp2_Cou+Emp2_EX

    DEALLOCATE(proc_map)
    DEALLOCATE(sub_proc_map) 
    DEALLOCATE(proc_map_rep)
    DEALLOCATE(ranges_info_array)

    IF(.NOT.my_open_shell_SS) THEN
      ! keep the array for the next calculations
      IF(ALLOCATED(BIb_C)) DEALLOCATE(BIb_C)
      DEALLOCATE(sizes_array)
      DEALLOCATE(starts_array)
      DEALLOCATE(ends_array)
      DEALLOCATE(starts_B_virtual)
      DEALLOCATE(ends_B_virtual)
      DEALLOCATE(sizes_B_virtual)
      IF(my_alpha_beta_case) THEN
        DEALLOCATE(BIb_C_beta)
        DEALLOCATE(starts_B_virtual_beta)
        DEALLOCATE(ends_B_virtual_beta)
        DEALLOCATE(sizes_B_virtual_beta)
      END IF
    END IF

    CALL cp_para_env_release(para_env_exchange,error=error_sub)
    CALL cp_para_env_release(para_env_rep,error=error_sub)

    my_flop_rate=my_flop_rate/REAL(MAX(my_num_dgemm_call,1),KIND=dp)/1.0E9_dp
    CALL mp_sum(my_flop_rate,para_env%group)
    my_flop_rate=my_flop_rate/para_env%num_pe
    IF (unit_nr>0) WRITE (UNIT=unit_nr,FMT="(T3,A,T66,F15.2)")&
                          "PERFORMANCE| DGEMM flop rate (Gflops / MPI rank):", my_flop_rate

    CALL timestop(handle)

    END SUBROUTINE mp2_ri_gpw_compute_en

! *****************************************************************************
!> \brief ...
!> \param BIb_C ...
!> \param para_env ...
!> \param para_env_sub ...
!> \param para_env_exchange ...
!> \param para_env_rep ...
!> \param ngroup ...
!> \param integ_group_size ...
!> \param num_integ_group ...
!> \param homo ...
!> \param proc_map ...
!> \param proc_map_rep ...
!> \param ends_array ...
!> \param ends_B_virtual ...
!> \param sizes_array ...
!> \param sizes_B_virtual ...
!> \param starts_array ...
!> \param starts_B_virtual ...
!> \param my_B_size ...
!> \param my_B_virtual_start ...
!> \param my_B_virtual_end ...
!> \param my_group_L_size ...
!> \param my_group_L_start ...
!> \param my_group_L_end ...
!> \param my_new_group_L_size ...
!> \param new_sizes_array ...
!> \param ranges_info_array ...
!> \param error_sub ...
! *****************************************************************************
    SUBROUTINE replicate_iaK_2intgroup(BIb_C,para_env,para_env_sub,para_env_exchange,para_env_rep,&
                                       ngroup,integ_group_size,num_integ_group,homo,proc_map,proc_map_rep,&
                                       ends_array,ends_B_virtual,sizes_array,sizes_B_virtual,starts_array,starts_B_virtual,&
                                       my_B_size,my_B_virtual_start,my_B_virtual_end,&
                                       my_group_L_size,my_group_L_start,my_group_L_end,&
                                       my_new_group_L_size,new_sizes_array,ranges_info_array,&
                                       error_sub)
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: BIb_C
    TYPE(cp_para_env_type), POINTER          :: para_env, para_env_sub, &
                                                para_env_exchange, &
                                                para_env_rep
    INTEGER                                  :: ngroup, integ_group_size, &
                                                num_integ_group, homo
    INTEGER, ALLOCATABLE, DIMENSION(:) :: proc_map, proc_map_rep, ends_array, &
      ends_B_virtual, sizes_array, sizes_B_virtual, starts_array, &
      starts_B_virtual
    INTEGER :: my_B_size, my_B_virtual_start, my_B_virtual_end, &
      my_group_L_size, my_group_L_start, my_group_L_end, my_new_group_L_size
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: new_sizes_array
    INTEGER, ALLOCATABLE, DIMENSION(:, :, :) :: ranges_info_array
    TYPE(cp_error_type), INTENT(inout)       :: error_sub

    CHARACTER(LEN=*), PARAMETER :: routineN = 'replicate_iaK_2intgroup', &
      routineP = moduleN//':'//routineN

    INTEGER :: comm_rep, end_point, handle, i, max_L_size, proc_receive, &
      proc_send, proc_shift, start_point, stat, sub_sub_color
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: rep_ends_array, &
                                                rep_sizes_array, &
                                                rep_starts_array
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: BIb_C_copy
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :)                  :: BIb_C_gather

    CALL timeset(routineN,handle)
    failure=.FALSE.                                      

    ! create the replication group
    sub_sub_color=para_env_sub%mepos*para_env_exchange%num_pe+para_env_exchange%mepos
    CALL mp_comm_split_direct(para_env%group,comm_rep,sub_sub_color)
    NULLIFY(para_env_rep)
    CALL cp_para_env_create(para_env_rep,comm_rep,error=error_sub)

    ! crate the proc maps
    ALLOCATE(proc_map_rep(-para_env_rep%num_pe:2*para_env_rep%num_pe-1))
    DO i=0,para_env_rep%num_pe-1
      proc_map_rep(i)=i
      proc_map_rep(-i-1)=para_env_rep%num_pe-i-1
      proc_map_rep(para_env_rep%num_pe+i)=i
    END DO

    ! create the new limits for K according to the size
    ! of the integral group
    ALLOCATE(new_sizes_array(0:para_env_exchange%num_pe-1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
    new_sizes_array=0
    ALLOCATE(ranges_info_array(4,0:para_env_rep%num_pe-1,0:para_env_exchange%num_pe-1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
    ranges_info_array=0

    ! info array for replication
    ALLOCATE(rep_ends_array(0:para_env_rep%num_pe-1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
    rep_ends_array=0
    ALLOCATE(rep_starts_array(0:para_env_rep%num_pe-1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
    rep_starts_array=0
    ALLOCATE(rep_sizes_array(0:para_env_rep%num_pe-1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
    rep_sizes_array=0

    rep_sizes_array(para_env_rep%mepos)=my_group_L_size
    rep_starts_array(para_env_rep%mepos)=my_group_L_start
    rep_ends_array(para_env_rep%mepos)=my_group_L_end

    CALL mp_sum(rep_sizes_array,para_env_rep%group)
    CALL mp_sum(rep_starts_array,para_env_rep%group)
    CALL mp_sum(rep_ends_array,para_env_rep%group)

    ! calculate my_new_group_L_size according to sizes_array
    my_new_group_L_size=my_group_L_size
    ranges_info_array(1,0,para_env_exchange%mepos)=my_group_L_start
    ranges_info_array(2,0,para_env_exchange%mepos)=my_group_L_end
    ranges_info_array(3,0,para_env_exchange%mepos)=1
    ranges_info_array(4,0,para_env_exchange%mepos)=my_group_L_size

    DO proc_shift=1, para_env_rep%num_pe-1
      proc_send=proc_map_rep(para_env_rep%mepos+proc_shift)
      proc_receive=proc_map_rep(para_env_rep%mepos-proc_shift)

      my_new_group_L_size=my_new_group_L_size+rep_sizes_array(proc_receive)

      ranges_info_array(1,proc_shift,para_env_exchange%mepos)=rep_starts_array(proc_receive)
      ranges_info_array(2,proc_shift,para_env_exchange%mepos)=rep_ends_array(proc_receive)
      ranges_info_array(3,proc_shift,para_env_exchange%mepos)=ranges_info_array(4,proc_shift-1,para_env_exchange%mepos)+1
      ranges_info_array(4,proc_shift,para_env_exchange%mepos)=my_new_group_L_size

    END DO
    new_sizes_array(para_env_exchange%mepos)=my_new_group_L_size

    CALL mp_sum(new_sizes_array,para_env_exchange%group)
    CALL mp_sum(ranges_info_array,para_env_exchange%group)

    IF(.FALSE.) THEN
      ! replication scheme using mp_sendrecv
      ALLOCATE(BIb_C_copy(my_group_L_size,my_B_size,homo),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
      BIb_C_copy(:,:,:)=BIb_C

      DEALLOCATE(BIb_C)

      ALLOCATE(BIb_C(my_new_group_L_size,my_B_size,homo),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
      BIb_C=0.0_dp

      start_point=ranges_info_array(3,0,para_env_exchange%mepos)
      end_point=ranges_info_array(4,0,para_env_exchange%mepos)

      BIb_C(start_point:end_point,1:my_B_size,1:homo)=BIb_C_copy

      DO proc_shift=1, para_env_rep%num_pe-1
        proc_send=proc_map_rep(para_env_rep%mepos+proc_shift)
        proc_receive=proc_map_rep(para_env_rep%mepos-proc_shift)

        start_point=ranges_info_array(3,proc_shift,para_env_exchange%mepos)
        end_point=ranges_info_array(4,proc_shift,para_env_exchange%mepos)

        CALL  mp_sendrecv(BIb_C_copy(1:my_group_L_size,1:my_B_size,1:homo),proc_send,&
                          BIb_C(start_point:end_point,1:my_B_size,1:homo),proc_receive,&
                          para_env_rep%group)

      END DO

      DEALLOCATE(BIb_C_copy)

    ELSE
      ! replication scheme using mp_allgather
      ! get the max L size of the 
      max_L_size=MAXVAL(sizes_array)

      ALLOCATE(BIb_C_copy(max_L_size,my_B_size,homo),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
      BIb_C_copy=0.0_dp
      BIb_C_copy(1:my_group_L_size,1:my_B_size,1:homo)=BIb_C

      DEALLOCATE(BIb_C)

      ALLOCATE(BIb_C_gather(max_L_size,my_B_size,homo,0:para_env_rep%num_pe-1),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
      BIb_C_gather=0.0_dp

      CALL mp_allgather(BIb_C_copy,BIb_C_gather,para_env_rep%group)

      DEALLOCATE(BIb_C_copy)

      ALLOCATE(BIb_C(my_new_group_L_size,my_B_size,homo),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
      BIb_C=0.0_dp

      ! reorder data
      DO proc_shift=0, para_env_rep%num_pe-1
        proc_send=proc_map_rep(para_env_rep%mepos+proc_shift)
        proc_receive=proc_map_rep(para_env_rep%mepos-proc_shift)

        start_point=ranges_info_array(3,proc_shift,para_env_exchange%mepos)
        end_point=ranges_info_array(4,proc_shift,para_env_exchange%mepos)

        BIb_C(start_point:end_point,1:my_B_size,1:homo)=&
                     BIb_C_gather(1:end_point-start_point+1,1:my_B_size,1:homo,proc_receive)

      END DO

      DEALLOCATE(BIb_C_gather)

    END IF

    ! DEALLOCATE(proc_map_rep)
    DEALLOCATE(rep_sizes_array)
    DEALLOCATE(rep_starts_array)
    DEALLOCATE(rep_ends_array)

    ! CALL  cp_para_env_release(para_env_rep,error=error_sub)

    CALL timestop(handle)

    END SUBROUTINE replicate_iaK_2intgroup

! *****************************************************************************
!> \brief ...
!> \param mat ...
! *****************************************************************************
      SUBROUTINE write_array(mat)
    REAL(KIND=dp), DIMENSION(:, :)           :: mat

    INTEGER                                  :: iii, jjj

      WRITE(*,*)
      DO iii=1, SIZE(mat,1)
        DO jjj=1, SIZE(mat,2), 20
          WRITE(*,'(1000F10.5)') mat(iii,jjj:MIN(SIZE(mat,2),jjj+19))
        END DO
        WRITE(*,*)
      END DO
      WRITE(*,*)
     END SUBROUTINE

! *****************************************************************************
!> \brief ...
!> \param local_ab ...
!> \param t_ab ...
!> \param mp2_env ...
!> \param homo ...
!> \param virtual ...
!> \param dimen_RI ...
!> \param my_B_size ...
!> \param block_size ...
!> \param my_B_size_beta ...
!> \param my_group_L_size ...
!> \param Y_i_aP ...
!> \param Y_j_aP ...
!> \param local_i_aL ...
!> \param local_j_aL ...
!> \param calc_forces ...
!> \param error_sub ...
! *****************************************************************************
    SUBROUTINE mp2_ri_allocate(local_ab, t_ab,mp2_env,homo,virtual,dimen_RI,my_B_size,&
                                block_size,my_B_size_beta,my_group_L_size,Y_i_aP,&
                                Y_j_aP,local_i_aL, local_j_aL,calc_forces,error_sub)
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: local_ab, t_ab
    TYPE(mp2_type), POINTER                  :: mp2_env
    INTEGER                                  :: homo, virtual, dimen_RI, &
                                                my_B_size, block_size, &
                                                my_B_size_beta, &
                                                my_group_L_size
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: Y_i_aP, Y_j_aP, local_i_aL, &
                                                local_j_aL
    LOGICAL                                  :: calc_forces
    TYPE(cp_error_type), INTENT(inout)       :: error_sub

    CHARACTER(LEN=*), PARAMETER :: routineN = 'mp2_ri_allocate', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, stat
    LOGICAL                                  :: failure

    CALL timeset(routineN,handle)
    failure=.FALSE.

    ALLOCATE(local_i_aL(dimen_RI,my_B_size,block_size),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
    ALLOCATE(local_j_aL(dimen_RI,my_B_size_beta,block_size),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
    ALLOCATE(local_ab(virtual,my_B_size_beta),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)

    IF(calc_forces) THEN
      ALLOCATE(mp2_env%ri_grad%P_ij(homo,homo),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
      ALLOCATE(mp2_env%ri_grad%P_ab(my_B_size,virtual),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
      ALLOCATE(mp2_env%ri_grad%Gamma_P_ia(homo,my_B_size,my_group_L_size),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
      mp2_env%ri_grad%P_ij=0.0_dp
      mp2_env%ri_grad%P_ab=0.0_dp
      mp2_env%ri_grad%Gamma_P_ia=0.0_dp
      ALLOCATE(Y_i_aP(my_B_size,dimen_RI,block_size),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
      ALLOCATE(Y_j_aP(my_B_size,dimen_RI,block_size),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
      ALLOCATE(t_ab(virtual,my_B_size_beta),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
    END IF

    CALL timestop(handle)

    END SUBROUTINE mp2_ri_allocate

! *****************************************************************************
!> \brief ...
!> \param my_alpha_beta_case ...
!> \param total_ij_pairs ...
!> \param homo ...
!> \param homo_beta ...
!> \param num_IJ_blocks ...
!> \param block_size ...
!> \param ngroup ...
!> \param ij_map ...
!> \param color_sub ...
!> \param my_ij_pairs ...
!> \param error_sub ...
!> \param my_open_shell_SS ...
!> \param unit_nr ...
! *****************************************************************************
    SUBROUTINE mp2_ri_communication(my_alpha_beta_case,total_ij_pairs,homo,homo_beta,num_IJ_blocks,&
               block_size,ngroup,ij_map,color_sub,my_ij_pairs,error_sub,my_open_shell_SS,unit_nr)
    LOGICAL                                  :: my_alpha_beta_case
    INTEGER                                  :: total_ij_pairs, homo, &
                                                homo_beta, num_IJ_blocks, &
                                                block_size, ngroup
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: ij_map
    INTEGER                                  :: color_sub, my_ij_pairs
    TYPE(cp_error_type), INTENT(inout)       :: error_sub
    LOGICAL                                  :: my_open_shell_SS
    INTEGER                                  :: unit_nr

    CHARACTER(LEN=*), PARAMETER :: routineN = 'mp2_ri_communication', &
      routineP = moduleN//':'//routineN

    INTEGER :: assigned_blocks, first_I_block, first_J_block, handle, iiB, &
      ij_block_counter, ij_counter, jjB, last_i_block, last_J_block, &
      num_block_per_group, stat, total_ij_block, total_ij_pairs_blocks
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: ij_marker
    LOGICAL                                  :: failure

! Calculate the maximum number of ij pairs that have to be computed
! among groups

    CALL timeset(routineN,handle)
    failure=.FALSE.

    IF(.NOT.my_alpha_beta_case) THEN
      total_ij_pairs=homo*(1+homo)/2
      num_IJ_blocks=homo/block_size-1

      first_I_block=1
      last_i_block=block_size*(num_IJ_blocks-1)

      first_J_block=block_size+1
      last_J_block=block_size*(num_IJ_blocks+1)

      ij_block_counter=0
      DO iiB=first_I_block, last_i_block, block_size
        DO jjB=iiB+block_size, last_J_block, block_size
          ij_block_counter=ij_block_counter+1
        END DO
      END DO

      total_ij_block=ij_block_counter
      num_block_per_group=total_ij_block/ngroup
      assigned_blocks=num_block_per_group*ngroup

      total_ij_pairs_blocks=assigned_blocks+(total_ij_pairs-assigned_blocks*(block_size**2))

      ALLOCATE(ij_marker(homo,homo),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
      ij_marker=0
      ALLOCATE(ij_map(total_ij_pairs_blocks,3),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
      ij_map=0
      ij_counter=0
      my_ij_pairs=0
      DO iiB=first_I_block, last_i_block, block_size
        DO jjB=iiB+block_size, last_J_block, block_size
          IF(ij_counter+1>assigned_blocks) EXIT
          ij_counter=ij_counter+1
          ij_marker(iiB:iiB+block_size-1,jjB:jjB+block_size-1)=1
          ij_map(ij_counter,1)=iiB
          ij_map(ij_counter,2)=jjB
          ij_map(ij_counter,3)=block_size
          IF (MOD(ij_counter,ngroup)==color_sub) my_ij_pairs=my_ij_pairs+1
        END DO
      END DO
      DO iiB=1, homo
        DO jjB=iiB, homo
          IF(ij_marker(iiB,jjB)==0) THEN
            ij_counter=ij_counter+1
            ij_map(ij_counter,1)=iiB
            ij_map(ij_counter,2)=jjB
            ij_map(ij_counter,3)=1
            IF (MOD(ij_counter,ngroup)==color_sub) my_ij_pairs=my_ij_pairs+1
          END IF
        END DO
      END DO
      DEALLOCATE(ij_marker)

      IF((.NOT.my_open_shell_SS)) THEN
        IF (unit_nr>0) THEN
          IF(block_size==1) THEN
            WRITE (UNIT=unit_nr,FMT="(T3,A,T66,F15.1)")&
                        "RI_INFO| Percentage of ij pairs communicated with block size 1:", 100.0_dp
          ELSE
            WRITE (UNIT=unit_nr,FMT="(T3,A,T66,F15.1)")&
                        "RI_INFO| Percentage of ij pairs communicated with block size 1:", &
                        100.0_dp*REAL((total_ij_pairs-assigned_blocks*(block_size**2)),KIND=dp)/REAL(total_ij_pairs,KIND=dp)
          END IF
          CALL m_flush(unit_nr)
        END IF
      END IF

    ELSE
      ! alpha-beta case no index symmetry
      total_ij_pairs=homo*homo_beta
      ALLOCATE(ij_map(total_ij_pairs,3),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
      ij_map=0
      ij_counter=0
      my_ij_pairs=0
      DO iiB=1, homo 
        DO jjB=1, homo_beta
          ij_counter=ij_counter+1
          ij_map(ij_counter,1)=iiB
          ij_map(ij_counter,2)=jjB
          ij_map(ij_counter,3)=1
          IF (MOD(ij_counter,ngroup)==color_sub) my_ij_pairs=my_ij_pairs+1
        END DO
      END DO
    END IF

    CALL timestop(handle)
    
    END SUBROUTINE mp2_ri_communication
    
! *****************************************************************************
!> \brief ...
!> \param BIb_C ...
!> \param para_env ...
!> \param para_env_sub ...
!> \param homo ...
!> \param color_sub ...
!> \param ends_array ...
!> \param ends_B_virtual ...
!> \param sizes_array ...
!> \param sizes_B_virtual ...
!> \param starts_array ...
!> \param starts_B_virtual ...
!> \param calc_forces ...
!> \param error ...
!> \param error_sub ...
!> \param comm_exchange ...
!> \param integ_group_size ...
!> \param my_B_size ...
!> \param iiB ...
!> \param my_B_virtual_end ...
!> \param my_B_virtual_start ...
!> \param my_group_L_end ...
!> \param my_group_L_size ...
!> \param my_group_L_size_orig ...
!> \param my_group_L_start ...
!> \param my_new_group_L_size ...
!> \param ngroup ...
!> \param sub_sub_color ...
!> \param integ_group_pos2color_sub ...
!> \param new_sizes_array ...
!> \param proc_map ...
!> \param proc_map_rep ...
!> \param sizes_array_orig ...
!> \param sub_proc_map ...
!> \param ranges_info_array ...
!> \param para_env_exchange ...
!> \param para_env_rep ...
!> \param num_integ_group ...
! *****************************************************************************
    SUBROUTINE mp2_ri_create_group(BIb_C,para_env,para_env_sub,homo,color_sub, ends_array, ends_B_virtual, &
                 sizes_array,sizes_B_virtual,starts_array,starts_B_virtual,calc_forces,error, error_sub,&
                 comm_exchange,integ_group_size, my_B_size,iiB, my_B_virtual_end,my_B_virtual_start, my_group_L_end,&
                 my_group_L_size, my_group_L_size_orig, my_group_L_start,my_new_group_L_size, ngroup, &
                 sub_sub_color,integ_group_pos2color_sub,new_sizes_array, proc_map, proc_map_rep, sizes_array_orig,&
                 sub_proc_map,ranges_info_array,para_env_exchange,para_env_rep, num_integ_group)
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: BIb_C
    TYPE(cp_para_env_type), POINTER          :: para_env, para_env_sub
    INTEGER                                  :: homo, color_sub
    INTEGER, ALLOCATABLE, DIMENSION(:) :: ends_array, ends_B_virtual, &
      sizes_array, sizes_B_virtual, starts_array, starts_B_virtual
    LOGICAL, INTENT(IN)                      :: calc_forces
    TYPE(cp_error_type), INTENT(inout)       :: error, error_sub
    INTEGER :: comm_exchange, integ_group_size, my_B_size, iiB, &
      my_B_virtual_end, my_B_virtual_start, my_group_L_end, my_group_L_size, &
      my_group_L_size_orig, my_group_L_start, my_new_group_L_size, ngroup, &
      sub_sub_color
    INTEGER, ALLOCATABLE, DIMENSION(:) :: integ_group_pos2color_sub, &
      new_sizes_array, proc_map, proc_map_rep, sizes_array_orig, sub_proc_map
    INTEGER, ALLOCATABLE, DIMENSION(:, :, :) :: ranges_info_array
    TYPE(cp_para_env_type), POINTER          :: para_env_exchange, &
                                                para_env_rep
    INTEGER                                  :: num_integ_group

    CHARACTER(LEN=*), PARAMETER :: routineN = 'mp2_ri_create_group', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, stat
    LOGICAL                                  :: failure

!

    CALL timeset(routineN,handle)
    failure=.FALSE.
    !
    sub_sub_color=para_env_sub%mepos*num_integ_group+color_sub/integ_group_size
    CALL mp_comm_split_direct(para_env%group,comm_exchange,sub_sub_color)
    NULLIFY(para_env_exchange)
    CALL cp_para_env_create(para_env_exchange,comm_exchange,error=error_sub)

    ! create the proc maps
    ALLOCATE(proc_map(-para_env_exchange%num_pe:2*para_env_exchange%num_pe-1))
    DO i=0,para_env_exchange%num_pe-1
      proc_map(i)=i
      proc_map(-i-1)=para_env_exchange%num_pe-i-1
      proc_map(para_env_exchange%num_pe+i)=i
    END DO

    ALLOCATE(sub_proc_map(-para_env_sub%num_pe:2*para_env_sub%num_pe-1))
    DO i=0,para_env_sub%num_pe-1
      sub_proc_map(i)=i
      sub_proc_map(-i-1)=para_env_sub%num_pe-i-1
      sub_proc_map(para_env_sub%num_pe+i)=i
    END DO

    CALL replicate_iaK_2intgroup(BIb_C,para_env,para_env_sub,para_env_exchange,para_env_rep,&
                                 ngroup,integ_group_size,num_integ_group,homo,proc_map,proc_map_rep,&
                                 ends_array,ends_B_virtual,sizes_array,sizes_B_virtual,starts_array,starts_B_virtual,&
                                 my_B_size,my_B_virtual_start,my_B_virtual_end,&
                                 my_group_L_size,my_group_L_start,my_group_L_end,&
                                 my_new_group_L_size,new_sizes_array,ranges_info_array,&
                                 error_sub)

    ALLOCATE(integ_group_pos2color_sub(0:para_env_exchange%num_pe-1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
    integ_group_pos2color_sub=0
    integ_group_pos2color_sub(para_env_exchange%mepos)=color_sub
    CALL mp_sum(integ_group_pos2color_sub,para_env_exchange%group)

    IF(calc_forces) THEN
      iiB=SIZE(sizes_array)
      ALLOCATE(sizes_array_orig(0:iiB-1),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
      sizes_array_orig(:)=sizes_array
    END IF  
 
    my_group_L_size_orig=my_group_L_size
    my_group_L_size=my_new_group_L_size
    DEALLOCATE(sizes_array)

    ALLOCATE(sizes_array(0:integ_group_size-1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
    sizes_array(:)=new_sizes_array
      
    DEALLOCATE(new_sizes_array)
    !
    CALL timestop(handle)

    END SUBROUTINE mp2_ri_create_group

! *****************************************************************************
!> \brief ...
!> \param mp2_env ...
!> \param para_env ...
!> \param para_env_sub ...
!> \param ends_array ...
!> \param ends_B_virtual ...
!> \param sizes_array ...
!> \param sizes_B_virtual ...
!> \param starts_array ...
!> \param starts_B_virtual ...
!> \param homo ...
!> \param dimen_RI ...
!> \param unit_nr ...
!> \param color_sub ...
!> \param best_block_size ...
!> \param best_integ_group_size ...
!> \param block_size ...
!> \param integ_group_size ...
!> \param min_integ_group_size ...
!> \param my_B_size ...
!> \param my_B_virtual_end ...
!> \param my_B_virtual_start ...
!> \param my_group_L_size ...
!> \param my_group_L_start ...
!> \param my_group_L_end ...
!> \param ngroup ...
!> \param num_IJ_blocks ...
!> \param num_integ_group ...
!> \param pos_integ_group ...
!> \param virtual ...
!> \param my_alpha_beta_case ...
!> \param my_open_shell_SS ...
!> \param actual_flop_rate ...
!> \param mem_for_aK ...
!> \param mem_for_comm ...
!> \param mem_for_iaK ...
!> \param mem_for_rep ...
!> \param mem_min ...
!> \param mem_per_group ...
!> \param mem_real ...
! *****************************************************************************
    SUBROUTINE mp2_ri_get_sizes(mp2_env,para_env, para_env_sub,ends_array, ends_B_virtual,sizes_array, sizes_B_virtual,&
    starts_array, starts_B_virtual,homo, dimen_RI, unit_nr,color_sub,best_block_size,best_integ_group_size, block_size, &
    integ_group_size,min_integ_group_size, my_B_size,my_B_virtual_end, my_B_virtual_start, my_group_L_size, &
    my_group_L_start, my_group_L_end, ngroup, num_IJ_blocks, num_integ_group,pos_integ_group, virtual,my_alpha_beta_case,&
    my_open_shell_SS,actual_flop_rate, mem_for_aK, mem_for_comm, mem_for_iaK,mem_for_rep, mem_min, mem_per_group, mem_real)
    TYPE(mp2_type), POINTER                  :: mp2_env
    TYPE(cp_para_env_type), POINTER          :: para_env, para_env_sub
    INTEGER, ALLOCATABLE, DIMENSION(:) :: ends_array, ends_B_virtual, &
      sizes_array, sizes_B_virtual, starts_array, starts_B_virtual
    INTEGER :: homo, dimen_RI, unit_nr, color_sub, best_block_size, &
      best_integ_group_size, block_size, integ_group_size, &
      min_integ_group_size, my_B_size, my_B_virtual_end, my_B_virtual_start, &
      my_group_L_size, my_group_L_start, my_group_L_end, ngroup, &
      num_IJ_blocks, num_integ_group, pos_integ_group, virtual
    LOGICAL                                  :: my_alpha_beta_case, &
                                                my_open_shell_SS
    REAL(KIND=dp) :: actual_flop_rate, mem_for_aK, mem_for_comm, mem_for_iaK, &
      mem_for_rep, mem_min, mem_per_group, mem_real

    CHARACTER(LEN=*), PARAMETER :: routineN = 'mp2_ri_get_sizes', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iiB
    INTEGER(KIND=int_8)                      :: mem
    LOGICAL                                  :: failure

    CALL timeset(routineN,handle)
    failure=.FALSE.

    ngroup=para_env%num_pe/para_env_sub%num_pe

    ! Calculate available memory and create integral group according to that
    ! mem_for_iaK is the memory needed for storing the 3 centre integrals
    mem_for_iaK=REAL(homo,KIND=dp)*virtual*dimen_RI*8.0_dp/(1024_dp**2)
    mem_for_aK=REAL(virtual,KIND=dp)*dimen_RI*8.0_dp/(1024_dp**2)

    CALL m_memory(mem)
    mem_real=(mem+1024*1024-1)/(1024*1024)
    ! mp_min .... a hack.. it should be mp_max, but as it turns out, on some processes the previously freed memory (hfx)
    ! has not been given back to the OS yet. 
    CALL mp_min(mem_real,para_env%group)

    mem_min=2.0_dp*REAL(homo,KIND=dp)*MAXVAL(sizes_B_virtual)*MAXVAL(sizes_array)*8.0_dp/(1024**2)
    mem_min=mem_min+3.0_dp*MAXVAL(sizes_B_virtual)*REAL(dimen_RI,KIND=dp)*8.0_dp/(1024**2)
 
    IF((.NOT.my_open_shell_SS).AND.(.NOT.my_alpha_beta_case)) THEN
      IF (unit_nr>0) WRITE(unit_nr,'(T3,A,T69,F9.2,A3)') 'RI_INFO| Minimum required memory per MPI process:',&
                                                          mem_min, ' MB'
    END IF
   
    mem_real=mp2_env%mp2_memory

    mem_per_group=mem_real*para_env_sub%num_pe
    
    ! here we try to find the best block_size and integ_group_size
    best_integ_group_size=ngroup
    best_block_size=1

    ! in the open shell case no replication and no block communication is done 
    IF((.NOT.my_open_shell_SS).AND.(.NOT.my_alpha_beta_case)) THEN
      ! Here we split the memory half for the communication, half for replication
      IF(mp2_env%ri_mp2%block_size>0) THEN
        best_block_size=mp2_env%ri_mp2%block_size
        mem_for_rep=MAX(mem_min,mem_per_group-2.0_dp*mem_for_aK*best_block_size)
      ELSE
        mem_for_rep=mem_per_group/2.0_dp
      END IF
      ! calculate the minimum replication group size according to the available memory
      min_integ_group_size=CEILING(2.0_dp*mem_for_iaK/mem_for_rep)

      integ_group_size=MIN(min_integ_group_size,ngroup)-1
      DO iiB=min_integ_group_size+1, ngroup
        integ_group_size=integ_group_size+1
        ! check that the ngroup is a multiple of  integ_group_size
        IF(MOD(ngroup,integ_group_size)/=0) CYCLE
        ! check that the integ group size is not too small (10% is empirical for now)
        IF(REAL(integ_group_size,KIND=dp)/REAL(ngroup,KIND=dp)<0.1_dp) CYCLE

        best_integ_group_size=integ_group_size
        EXIT
      END DO

      IF(.NOT.(mp2_env%ri_mp2%block_size>0)) THEN
        mem_for_comm=mem_per_group-2.0_dp*mem_for_iaK/best_integ_group_size
        ! calculate the maximum block size and make sure it is not smaller than 1
        best_block_size=MAX(INT(mem_for_comm/(2.0_dp*mem_for_aK)),1)
        DO 
          num_IJ_blocks=(homo/best_block_size)
          num_IJ_blocks=(num_IJ_blocks*num_IJ_blocks-num_IJ_blocks)/2
          IF(num_IJ_blocks>ngroup.OR.best_block_size==1) THEN
            EXIT
          ELSE
            best_block_size=best_block_size-1
          END IF
        END DO
      END IF
 
      ! check that best_block_size is not bigger than homo/2-1
      best_block_size=MIN(MAX(homo/2-1+MOD(homo,2),1),best_block_size)
    END IF

    integ_group_size=best_integ_group_size
    block_size=best_block_size

    IF((.NOT.my_open_shell_SS).AND.(.NOT.my_alpha_beta_case)) THEN
      IF (unit_nr>0) THEN
         WRITE (UNIT=unit_nr,FMT="(T3,A,T75,i6)")&
                                 "RI_INFO| Group size for integral replication:", integ_group_size*para_env_sub%num_pe
         WRITE (UNIT=unit_nr,FMT="(T3,A,T75,i6)")&
                                 "RI_INFO| Block size:", block_size
         CALL m_flush(unit_nr)
      END IF
    END IF
    
    num_integ_group=ngroup/integ_group_size

    pos_integ_group=MOD(color_sub,integ_group_size)

    my_group_L_size=sizes_array(color_sub)
    my_group_L_start=starts_array(color_sub)
    my_group_L_end=ends_array(color_sub)

    my_B_virtual_start=starts_B_virtual(para_env_sub%mepos)
    my_B_virtual_end=ends_B_virtual(para_env_sub%mepos)
    my_B_size=sizes_B_virtual(para_env_sub%mepos)

    CALL timestop(handle)

    END SUBROUTINE mp2_ri_get_sizes

! *****************************************************************************
!> \brief ...
!> \param mp2_env ...
!> \param para_env_sub ...
!> \param ends_B_virtual ...
!> \param sizes_B_virtual ...
!> \param starts_B_virtual ...
!> \param Eigenval ...
!> \param homo ...
!> \param dimen_RI ...
!> \param error_sub ...
!> \param b_global ...
!> \param iiB ...
!> \param jjB ...
!> \param my_B_size ...
!> \param my_B_virtual_end ...
!> \param my_B_virtual_start ...
!> \param my_i ...
!> \param my_j ...
!> \param proc_receive ...
!> \param proc_send ...
!> \param proc_shift ...
!> \param rec_B_size ...
!> \param rec_B_virtual_end ...
!> \param rec_B_virtual_start ...
!> \param send_B_size ...
!> \param send_B_virtual_end ...
!> \param send_B_virtual_start ...
!> \param virtual ...
!> \param proc_map ...
!> \param sub_proc_map ...
!> \param external_ab ...
!> \param local_ab ...
!> \param send_ab ...
!> \param t_ab ...
!> \param local_i_aL ...
!> \param local_j_aL ...
!> \param Y_i_aP ...
!> \param Y_j_aP ...
! *****************************************************************************
    SUBROUTINE mp2_update_P_gamma(mp2_env,para_env_sub,ends_B_virtual,sizes_B_virtual, starts_B_virtual,&
               Eigenval, homo, dimen_RI, error_sub, b_global, iiB, jjB, my_B_size, &
               my_B_virtual_end, my_B_virtual_start, my_i, my_j, proc_receive, proc_send,proc_shift, &
               rec_B_size, rec_B_virtual_end, rec_B_virtual_start, send_B_size, send_B_virtual_end, & 
               send_B_virtual_start, virtual, proc_map, sub_proc_map,external_ab, local_ab,&
               send_ab, t_ab,local_i_aL, local_j_aL, Y_i_aP, Y_j_aP)
    TYPE(mp2_type), POINTER                  :: mp2_env
    TYPE(cp_para_env_type), POINTER          :: para_env_sub
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: ends_B_virtual, &
                                                sizes_B_virtual, &
                                                starts_B_virtual
    REAL(KIND=dp), DIMENSION(:)              :: Eigenval
    INTEGER                                  :: homo, dimen_RI
    TYPE(cp_error_type), INTENT(inout)       :: error_sub
    INTEGER :: b_global, iiB, jjB, my_B_size, my_B_virtual_end, &
      my_B_virtual_start, my_i, my_j, proc_receive, proc_send, proc_shift, &
      rec_B_size, rec_B_virtual_end, rec_B_virtual_start, send_B_size, &
      send_B_virtual_end, send_B_virtual_start, virtual
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: proc_map, sub_proc_map
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: external_ab, local_ab, &
                                                send_ab, t_ab
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: local_i_aL, local_j_aL, &
                                                Y_i_aP, Y_j_aP

    CHARACTER(LEN=*), PARAMETER :: routineN = 'mp2_update_P_gamma', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: a, b, handle, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: P_ij_diag

! update P_ab, Gamma_P_ia
! First, P_ab

    CALL timeset(routineN//"_Pia",handle)
    failure=.FALSE.
    ! divide the (ia|jb) integrals by Delta_ij^ab
    DO b=1, my_B_size
      b_global=b+my_B_virtual_start-1
      DO a=1, virtual
        local_ab(a,b)=-local_ab(a,b)/&
                      (Eigenval(homo+a)+Eigenval(homo+b_global)-Eigenval(my_i+iiB-1)-Eigenval(my_j+jjB-1))
      END DO
    END DO
    CALL dgemm('T','N',my_B_size,my_B_size,virtual,1.0_dp,&
                   t_ab(:,:),virtual,local_ab(:,:),virtual,&
                   1.0_dp,mp2_env%ri_grad%P_ab(1:my_B_size,my_B_virtual_start:my_B_virtual_end),my_B_size)
    ! update diagonal part of P_ij
    P_ij_diag=-SUM(local_ab*t_ab)*2.0_dp
    mp2_env%ri_grad%P_ij(my_i+iiB-1,my_i+iiB-1)=mp2_env%ri_grad%P_ij(my_i+iiB-1,my_i+iiB-1)+P_ij_diag
    IF(my_i/=my_j) THEN
      ! CALL dgemm('N','T',virtual,virtual,my_B_size,1.0_dp,&
      !                t_ab(:,:),virtual,local_ab(:,:),virtual,&
      !                1.0_dp,mp2_env%ri_grad%P_ab(:,:),virtual)
      CALL dgemm('N','T',my_B_size,virtual,my_B_size,1.0_dp,&
                     t_ab(my_B_virtual_start:my_B_virtual_end,1:my_B_size),my_B_size,&
                     local_ab(:,:),virtual,&
                     1.0_dp,mp2_env%ri_grad%P_ab(1:my_B_size,1:virtual),my_B_size)
      ! CHECK CORRECTNESS IF BLOCK SIZE greater than 1
      mp2_env%ri_grad%P_ij(my_j+jjB-1,my_j+jjB-1)=mp2_env%ri_grad%P_ij(my_j+jjB-1,my_j+jjB-1)+P_ij_diag
    END IF
    DO proc_shift=1, para_env_sub%num_pe-1
      proc_send=sub_proc_map(para_env_sub%mepos+proc_shift)
      proc_receive=sub_proc_map(para_env_sub%mepos-proc_shift)

      rec_B_size=sizes_B_virtual(proc_receive)
      rec_B_virtual_end=ends_B_virtual(proc_receive)
      rec_B_virtual_start=starts_B_virtual(proc_receive)

      send_B_size=sizes_B_virtual(proc_send)
      send_B_virtual_end=ends_B_virtual(proc_send)
      send_B_virtual_start=starts_B_virtual(proc_send)

      ALLOCATE(external_ab(virtual,rec_B_size),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
      external_ab=0.0_dp

      CALL  mp_sendrecv(local_ab(1:virtual,1:my_B_size),proc_send,&
                        external_ab(1:virtual,1:rec_B_size),proc_receive,&
                        para_env_sub%group)

      CALL dgemm('T','N',my_B_size,rec_B_size,virtual,1.0_dp,&
                     t_ab(:,:),virtual,external_ab(:,:),virtual,&
                     1.0_dp,mp2_env%ri_grad%P_ab(1:my_B_size,rec_B_virtual_start:rec_B_virtual_end),my_B_size)

      DEALLOCATE(external_ab)

      IF(my_i/=my_j) THEN
        ALLOCATE(external_ab(my_B_size,virtual),STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
        external_ab=0.0_dp
        
        ALLOCATE(send_ab(send_B_size,virtual),STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
        send_ab=0.0_dp                

        CALL dgemm('N','T',send_B_size,virtual,my_B_size,1.0_dp,&
                       t_ab(send_B_virtual_start:send_B_virtual_end,1:my_B_size),send_B_size,&
                       local_ab(:,:),virtual,&
                       0.0_dp,send_ab(1:send_B_size,1:virtual),send_B_size)

        CALL  mp_sendrecv(send_ab,proc_send,&
                          external_ab,proc_receive,&
                          para_env_sub%group)

        mp2_env%ri_grad%P_ab(:,:)=mp2_env%ri_grad%P_ab+external_ab

        DEALLOCATE(external_ab)
        DEALLOCATE(send_ab)
      END IF

    END DO
    CALL timestop(handle)

    ! Now, Gamma_P_ia

    CALL timeset(routineN//"_Gamma",handle)
    failure=.FALSE.
    ! CALL dgemm('T','N',my_B_size,my_B_size_beta,dimen_RI,1.0_dp,&
    !                local_i_aL(:,:,iiB),dimen_RI,local_j_aL(:,:,jjB),dimen_RI,&
    !                0.0_dp,local_ab(my_B_virtual_start:my_B_virtual_end,1:my_B_size_beta),my_B_size)
    
    CALL dgemm('N','T',my_B_size,dimen_RI,my_B_size,1.0_dp,&
                   t_ab(my_B_virtual_start:my_B_virtual_end,1:my_B_size),my_B_size,&
                   local_j_aL(1:dimen_RI,1:my_B_size,jjB),dimen_RI,&
                   1.0_dp,Y_i_aP(1:my_B_size,1:dimen_RI,iiB),my_B_size)

    ALLOCATE(external_ab(my_B_size,dimen_RI),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
    external_ab=0.0_dp
    DO proc_shift=1, para_env_sub%num_pe-1
      proc_send=sub_proc_map(para_env_sub%mepos+proc_shift)
      proc_receive=sub_proc_map(para_env_sub%mepos-proc_shift)

      rec_B_size=sizes_B_virtual(proc_receive)
      rec_B_virtual_end=ends_B_virtual(proc_receive)
      rec_B_virtual_start=starts_B_virtual(proc_receive)

      send_B_size=sizes_B_virtual(proc_send)
      send_B_virtual_end=ends_B_virtual(proc_send)
      send_B_virtual_start=starts_B_virtual(proc_send)

      ALLOCATE(send_ab(send_B_size,dimen_RI),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
      send_ab=0.0_dp

      CALL dgemm('N','T',send_B_size,dimen_RI,my_B_size,1.0_dp,&
                     t_ab(send_B_virtual_start:send_B_virtual_end,1:my_B_size),send_B_size,&
                     local_j_aL(1:dimen_RI,1:my_B_size,jjB),dimen_RI,&
                     0.0_dp,send_ab(1:send_B_size,1:dimen_RI),send_B_size)

      CALL  mp_sendrecv(send_ab,proc_send,&
                        external_ab,proc_receive,&
                        para_env_sub%group)

      Y_i_aP(1:my_B_size,1:dimen_RI,iiB)=Y_i_aP(1:my_B_size,1:dimen_RI,iiB)+external_ab

      DEALLOCATE(send_ab)
    END DO 
    DEALLOCATE(external_ab)
 

    IF(my_i/=my_j) THEN
      CALL dgemm('T','T',my_B_size,dimen_RI,my_B_size,1.0_dp,&
                     t_ab(my_B_virtual_start:my_B_virtual_end,1:my_B_size),my_B_size,&
                     local_i_aL(1:dimen_RI,1:my_B_size,iiB),dimen_RI,&
                     1.0_dp,Y_j_aP(1:my_B_size,1:dimen_RI,jjB),my_B_size)

      DO proc_shift=1, para_env_sub%num_pe-1
        proc_send=sub_proc_map(para_env_sub%mepos+proc_shift)
        proc_receive=sub_proc_map(para_env_sub%mepos-proc_shift)

        rec_B_size=sizes_B_virtual(proc_receive)
        rec_B_virtual_end=ends_B_virtual(proc_receive)
        rec_B_virtual_start=starts_B_virtual(proc_receive)

        ALLOCATE(external_ab(dimen_RI,rec_B_size),STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
        external_ab=0.0_dp

        CALL  mp_sendrecv(local_i_aL(1:dimen_RI,1:my_B_size,iiB),proc_send,&
                          external_ab,proc_receive,&
                          para_env_sub%group)

        CALL dgemm('T','T',my_B_size,dimen_RI,rec_B_size,1.0_dp,&
                       t_ab(rec_B_virtual_start:rec_B_virtual_end,1:my_B_size),rec_B_size,&
                       external_ab(1:dimen_RI,1:rec_B_size),dimen_RI,&
                       1.0_dp,Y_j_aP(1:my_B_size,1:dimen_RI,jjB),my_B_size)

        DEALLOCATE(external_ab)
      END DO
    END IF
    
    ! Y_i_aP=0.0_dp
    ! Y_j_aP=0.0_dp

    CALL timestop(handle)

    END SUBROUTINE mp2_update_P_gamma

! *****************************************************************************
!> \brief ...
!> \param mp2_env ...
!> \param error_sub ...
!> \param end_point ...
!> \param iiB ...
!> \param jjb ...
!> \param ij_counter_rec ...
!> \param ij_index ...
!> \param irep ...
!> \param kkk ...
!> \param Lend_pos ...
!> \param lll ...
!> \param Lstart_pos ...
!> \param my_B_size ...
!> \param my_block_size ...
!> \param my_group_L_size ...
!> \param my_i ...
!> \param my_ij_pairs ...
!> \param my_j ...
!> \param ngroup ...
!> \param num_integ_group ...
!> \param proc_receive ...
!> \param proc_send ...
!> \param proc_shift ...
!> \param rec_block_size ...
!> \param rec_i ...
!> \param rec_ij_index ...
!> \param rec_j ...
!> \param send_L_size ...
!> \param start_point ...
!> \param integ_group_pos2color_sub ...
!> \param num_ij_pairs ...
!> \param proc_map ...
!> \param ij_map ...
!> \param ranges_info_array ...
!> \param Y_i_aP ...
!> \param Y_j_aP ...
!> \param para_env_exchange ...
!> \param null_mat_rec ...
!> \param null_mat_send ...
!> \param sizes_array ...
! *****************************************************************************
    SUBROUTINE mp2_redistribute_gamma(mp2_env,error_sub,end_point, iiB, jjb, ij_counter_rec, &
               ij_index, irep, kkk, Lend_pos, lll, Lstart_pos, my_B_size, my_block_size, my_group_L_size, my_i,&
               my_ij_pairs, my_j, ngroup, num_integ_group, proc_receive, proc_send, proc_shift, rec_block_size,&
               rec_i, rec_ij_index, rec_j, send_L_size, start_point, integ_group_pos2color_sub,num_ij_pairs, &
               proc_map, ij_map, ranges_info_array, Y_i_aP, Y_j_aP, para_env_exchange, &
               null_mat_rec, null_mat_send,sizes_array )

    TYPE(mp2_type), POINTER                  :: mp2_env
    TYPE(cp_error_type), INTENT(inout)       :: error_sub
    INTEGER :: end_point, iiB, jjb, ij_counter_rec, ij_index, irep, kkk, &
      Lend_pos, lll, Lstart_pos, my_B_size, my_block_size, my_group_L_size, &
      my_i, my_ij_pairs, my_j, ngroup, num_integ_group, proc_receive, &
      proc_send, proc_shift, rec_block_size, rec_i, rec_ij_index, rec_j, &
      send_L_size, start_point
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: integ_group_pos2color_sub, &
                                                num_ij_pairs, proc_map
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: ij_map
    INTEGER, ALLOCATABLE, DIMENSION(:, :, :) :: ranges_info_array
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: Y_i_aP, Y_j_aP
    TYPE(cp_para_env_type), POINTER          :: para_env_exchange
    REAL(KIND=dp)                            :: null_mat_rec(:,:,:), &
                                                null_mat_send(:,:,:)
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: sizes_array

    CHARACTER(LEN=*), PARAMETER :: routineN = 'mp2_redistribute_gamma', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, handle2, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: BI_C_rec, BI_C_send

    CALL timeset(routineN//"_comm2",handle)
    IF(ij_index<=my_ij_pairs) THEN
      ! somethig to send
      ! start with myself
      CALL timeset(routineN//"_comm2_w",handle2)
      DO irep=0, num_integ_group-1
        Lstart_pos=ranges_info_array(1,irep,para_env_exchange%mepos)
        Lend_pos=ranges_info_array(2,irep,para_env_exchange%mepos)
        start_point=ranges_info_array(3,irep,para_env_exchange%mepos)
        end_point=ranges_info_array(4,irep,para_env_exchange%mepos)
        DO iiB=1, my_block_size
          !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(kkk,lll) &
          !$OMP                           SHARED(start_point,end_point,Lstart_pos,Lend_pos,&
          !$OMP                                  mp2_env,my_i,iiB,my_B_size,Y_i_aP)
          DO kkk=start_point, end_point
            lll=kkk-start_point+Lstart_pos
            mp2_env%ri_grad%Gamma_P_ia(my_i+iiB-1,1:my_B_size,kkk)=&
            mp2_env%ri_grad%Gamma_P_ia(my_i+iiB-1,1:my_B_size,kkk)+&
            Y_i_aP(1:my_B_size,lll,iiB)
          END DO
          !$OMP END PARALLEL DO
        END DO
        DO jjB=1, my_block_size
          !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(kkk,lll) &
          !$OMP                           SHARED(start_point,end_point,Lstart_pos,Lend_pos,&
          !$OMP                                  mp2_env,my_j,jjB,my_B_size,Y_j_aP)
          DO kkk=start_point, end_point
            lll=kkk-start_point+Lstart_pos
            mp2_env%ri_grad%Gamma_P_ia(my_j+jjB-1,1:my_B_size,kkk)=&
            mp2_env%ri_grad%Gamma_P_ia(my_j+jjB-1,1:my_B_size,kkk)+&
            Y_j_aP(1:my_B_size,lll,jjB)
          END DO
          !$OMP END PARALLEL DO
        END DO
        ! DO iiB=1, my_block_size
        !   mp2_env%ri_grad%Gamma_P_ia(my_i+iiB-1,1:my_B_size,start_point:end_point)=&
        !   mp2_env%ri_grad%Gamma_P_ia(my_i+iiB-1,1:my_B_size,start_point:end_point)+&
        !   Y_i_aP(1:my_B_size,Lstart_pos:Lend_pos,iiB)
        ! END DO
        ! DO jjB=1, my_block_size
        !   mp2_env%ri_grad%Gamma_P_ia(my_j+jjB-1,1:my_B_size,start_point:end_point)=&
        !   mp2_env%ri_grad%Gamma_P_ia(my_j+jjB-1,1:my_B_size,start_point:end_point)+&
        !   Y_j_aP(1:my_B_size,Lstart_pos:Lend_pos,jjB)
        ! END DO
      END DO
      CALL timestop(handle2)

      ! Y_i_aP(my_B_size,dimen_RI,block_size)

      DO proc_shift=1, para_env_exchange%num_pe-1
        proc_send=proc_map(para_env_exchange%mepos+proc_shift)
        proc_receive=proc_map(para_env_exchange%mepos-proc_shift)

        send_L_size=sizes_array(proc_send)

        ALLOCATE(BI_C_send(2*my_block_size,my_B_size,send_L_size),STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
        CALL timeset(routineN//"_comm2_w",handle2)
        BI_C_send=0.0_dp
        DO irep=0, num_integ_group-1
          Lstart_pos=ranges_info_array(1,irep,proc_send)
          Lend_pos=ranges_info_array(2,irep,proc_send)
          start_point=ranges_info_array(3,irep,proc_send)
          end_point=ranges_info_array(4,irep,proc_send)
          DO iiB=1, my_block_size
            !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(kkk,lll) &
            !$OMP                           SHARED(start_point,end_point,Lstart_pos,Lend_pos,&
            !$OMP                                  BI_C_send,iiB,my_B_size,Y_i_aP)
            DO kkk=start_point, end_point
              lll=kkk-start_point+Lstart_pos
              BI_C_send(iiB,1:my_B_size,kkk)=Y_i_aP(1:my_B_size,lll,iiB)
            END DO
            !$OMP END PARALLEL DO
          END DO
          DO jjB=1, my_block_size
            !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(kkk,lll) &
            !$OMP                           SHARED(start_point,end_point,Lstart_pos,Lend_pos,my_block_size,&
            !$OMP                                  BI_C_send,jjB,my_B_size,Y_j_aP)
            DO kkk=start_point, end_point
              lll=kkk-start_point+Lstart_pos
              BI_C_send(jjB+my_block_size,1:my_B_size,kkk)=Y_j_aP(1:my_B_size,lll,jjB)
            END DO
            !$OMP END PARALLEL DO
          END DO
          ! DO iiB=1, my_block_size
          !   BI_C_send(iiB,1:my_B_size,start_point:end_point)=&
          !      Y_i_aP(1:my_B_size,Lstart_pos:Lend_pos,iiB)
          ! END DO
          ! DO jjB=1, my_block_size
          !   BI_C_send(jjB+my_block_size,1:my_B_size,start_point:end_point)=&
          !      Y_j_aP(1:my_B_size,Lstart_pos:Lend_pos,jjB)
          ! END DO
        END DO
        CALL timestop(handle2)

        rec_ij_index=num_ij_pairs(proc_receive)

        IF(ij_index<=rec_ij_index) THEN
          ! we know that proc_receive has something to send for us, let's see what
          ij_counter_rec=&
          (ij_index-MIN(1,integ_group_pos2color_sub(proc_receive)))*ngroup+integ_group_pos2color_sub(proc_receive)

          rec_i=ij_map(ij_counter_rec,1)
          rec_j=ij_map(ij_counter_rec,2)
          rec_block_size=ij_map(ij_counter_rec,3)

          ALLOCATE(BI_C_rec(2*rec_block_size,my_B_size,my_group_L_size),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
          BI_C_rec=0.0_dp

          CALL  mp_sendrecv(BI_C_send,proc_send,&
                            BI_C_rec,proc_receive,&
                            para_env_exchange%group)

          CALL timeset(routineN//"_comm2_w",handle2)
          DO irep=0, num_integ_group-1
            Lstart_pos=ranges_info_array(1,irep,para_env_exchange%mepos)
            Lend_pos=ranges_info_array(2,irep,para_env_exchange%mepos)
            start_point=ranges_info_array(3,irep,para_env_exchange%mepos)
            end_point=ranges_info_array(4,irep,para_env_exchange%mepos)
            DO iiB=1, rec_block_size
              !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(kkk) &
              !$OMP                           SHARED(start_point,end_point,&
              !$OMP                                  mp2_env,rec_i,iiB,my_B_size,BI_C_rec)
              DO kkk=start_point, end_point
                mp2_env%ri_grad%Gamma_P_ia(rec_i+iiB-1,1:my_B_size,kkk)=&
                mp2_env%ri_grad%Gamma_P_ia(rec_i+iiB-1,1:my_B_size,kkk)+&
                BI_C_rec(iiB,1:my_B_size,kkk)
              END DO
              !$OMP END PARALLEL DO
            END DO
            DO jjB=1, rec_block_size
              !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(kkk) &
              !$OMP                           SHARED(start_point,end_point,rec_block_size,&
              !$OMP                                  mp2_env,rec_j,jjB,my_B_size,BI_C_rec)
              DO kkk=start_point, end_point
                mp2_env%ri_grad%Gamma_P_ia(rec_j+jjB-1,1:my_B_size,kkk)=&
                mp2_env%ri_grad%Gamma_P_ia(rec_j+jjB-1,1:my_B_size,kkk)+&
                BI_C_rec(jjB+rec_block_size,1:my_B_size,kkk)
              END DO
              !$OMP END PARALLEL DO
            END DO
            ! DO iiB=1, rec_block_size
            !   mp2_env%ri_grad%Gamma_P_ia(rec_i+iiB-1,1:my_B_size,start_point:end_point)=&
            !   mp2_env%ri_grad%Gamma_P_ia(rec_i+iiB-1,1:my_B_size,start_point:end_point)+&
            !   BI_C_rec(iiB,1:my_B_size,start_point:end_point)
            ! END DO
            ! DO jjB=1, rec_block_size
            !   mp2_env%ri_grad%Gamma_P_ia(rec_j+jjB-1,1:my_B_size,start_point:end_point)=&
            !   mp2_env%ri_grad%Gamma_P_ia(rec_j+jjB-1,1:my_B_size,start_point:end_point)+&
            !   BI_C_rec(jjB+rec_block_size,1:my_B_size,start_point:end_point)
            ! END DO
          END DO
          CALL timestop(handle2)

          DEALLOCATE(BI_C_rec)

        ELSE
          ! we have something to send but nothing to receive

          CALL  mp_sendrecv(BI_C_send,proc_send,&
                            null_mat_rec,proc_receive,&
                            para_env_exchange%group)

        END IF

        DEALLOCATE(BI_C_send)
      END DO   

    ELSE
      ! noting to send check if we have to receive
      DO proc_shift=1, para_env_exchange%num_pe-1
        proc_send=proc_map(para_env_exchange%mepos+proc_shift)
        proc_receive=proc_map(para_env_exchange%mepos-proc_shift)
        rec_ij_index=num_ij_pairs(proc_receive)

        IF(ij_index<=rec_ij_index) THEN
          ! we know that proc_receive has something to send for us, let's see what
          ij_counter_rec=&
          (ij_index-MIN(1,integ_group_pos2color_sub(proc_receive)))*ngroup+integ_group_pos2color_sub(proc_receive)

          rec_i=ij_map(ij_counter_rec,1)
          rec_j=ij_map(ij_counter_rec,2)
          rec_block_size=ij_map(ij_counter_rec,3)

          ALLOCATE(BI_C_rec(2*rec_block_size,my_B_size,my_group_L_size),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
          BI_C_rec=0.0_dp

          CALL  mp_sendrecv(null_mat_send,proc_send,&
                            BI_C_rec,proc_receive,&
                            para_env_exchange%group)

          CALL timeset(routineN//"_comm2_w",handle2)
          DO irep=0, num_integ_group-1
            Lstart_pos=ranges_info_array(1,irep,para_env_exchange%mepos)
            Lend_pos=ranges_info_array(2,irep,para_env_exchange%mepos)
            start_point=ranges_info_array(3,irep,para_env_exchange%mepos)
            end_point=ranges_info_array(4,irep,para_env_exchange%mepos)
            DO iiB=1, rec_block_size
              !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(kkk) &
              !$OMP                           SHARED(start_point,end_point,&
              !$OMP                                  mp2_env,rec_i,iiB,my_B_size,BI_C_rec)
              DO kkk=start_point, end_point
                mp2_env%ri_grad%Gamma_P_ia(rec_i+iiB-1,1:my_B_size,kkk)=&
                mp2_env%ri_grad%Gamma_P_ia(rec_i+iiB-1,1:my_B_size,kkk)+&
                BI_C_rec(iiB,1:my_B_size,kkk)
              END DO
              !$OMP END PARALLEL DO
            END DO
            DO jjB=1, rec_block_size
              !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(kkk) &
              !$OMP                           SHARED(start_point,end_point,rec_block_size,&
              !$OMP                                  mp2_env,rec_j,jjB,my_B_size,BI_C_rec)
              DO kkk=start_point, end_point
                mp2_env%ri_grad%Gamma_P_ia(rec_j+jjB-1,1:my_B_size,kkk)=&
                mp2_env%ri_grad%Gamma_P_ia(rec_j+jjB-1,1:my_B_size,kkk)+&
                BI_C_rec(jjB+rec_block_size,1:my_B_size,kkk)
              END DO
              !$OMP END PARALLEL DO
            END DO
            ! DO iiB=1, rec_block_size
            !   mp2_env%ri_grad%Gamma_P_ia(rec_i+iiB-1,1:my_B_size,start_point:end_point)=&
            !   mp2_env%ri_grad%Gamma_P_ia(rec_i+iiB-1,1:my_B_size,start_point:end_point)+&
            !   BI_C_rec(iiB,1:my_B_size,start_point:end_point)
            ! END DO
            ! DO jjB=1, rec_block_size
            !   mp2_env%ri_grad%Gamma_P_ia(rec_j+jjB-1,1:my_B_size,start_point:end_point)=&
            !   mp2_env%ri_grad%Gamma_P_ia(rec_j+jjB-1,1:my_B_size,start_point:end_point)+&
            !   BI_C_rec(jjB+rec_block_size,1:my_B_size,start_point:end_point)
            ! END DO
          END DO
          CALL timestop(handle2)

          DEALLOCATE(BI_C_rec)

        ELSE
          ! nothing to send nothing to receive
          CALL  mp_sendrecv(null_mat_send,proc_send,&
                            null_mat_rec,proc_receive,&
                            para_env_exchange%group)

        END IF
      END DO

    END IF
    CALL timestop(handle)

    END SUBROUTINE mp2_redistribute_gamma

END MODULE mp2_ri_gpw
