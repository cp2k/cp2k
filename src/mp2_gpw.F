!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Rountines to calculate MP2 energy using pw
!> \par History
!>      10.2011 created [Joost VandeVondele and Mauro Del Ben]
! *****************************************************************************
MODULE mp2_gpw
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE cell_types,                      ONLY: cell_type
  USE cp_blacs_env,                    ONLY: BLACS_GRID_SQUARE,&
                                             cp_blacs_env_create,&
                                             cp_blacs_env_release,&
                                             cp_blacs_env_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
  USE cp_dbcsr_interface,              ONLY: &
       array_i1d_obj, array_new, array_nullify, array_release, &
       cp_create_bl_distribution, cp_dbcsr_clear_mempools, cp_dbcsr_copy, &
       cp_dbcsr_create, cp_dbcsr_distribution, cp_dbcsr_distribution_release, &
       cp_dbcsr_filter, cp_dbcsr_get_info, cp_dbcsr_init, cp_dbcsr_init_p, &
       cp_dbcsr_iterator, cp_dbcsr_iterator_blocks_left, &
       cp_dbcsr_iterator_next_block, cp_dbcsr_iterator_start, &
       cp_dbcsr_iterator_stop, cp_dbcsr_multiply, cp_dbcsr_p_type, &
       cp_dbcsr_release, cp_dbcsr_reserve_all_blocks, &
       cp_dbcsr_row_block_sizes, cp_dbcsr_set, cp_dbcsr_type, &
       dbcsr_distribution_mp, dbcsr_distribution_new, dbcsr_distribution_obj, &
       dbcsr_distribution_row_dist, dbcsr_mp_npcols, dbcsr_mp_nprows, &
       dbcsr_type_no_symmetry, dbcsr_type_real_default, dbcsr_type_symmetric
  USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                             cp_dbcsr_dist2d_to_dist
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_info,&
                                             cp_fm_release,&
                                             cp_fm_type
  USE cp_para_env,                     ONLY: cp_para_env_create,&
                                             cp_para_env_release
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE distribution_1d_types,           ONLY: distribution_1d_release,&
                                             distribution_1d_type
  USE distribution_2d_types,           ONLY: distribution_2d_release,&
                                             distribution_2d_type
  USE distribution_methods,            ONLY: distribute_molecules_1d,&
                                             distribute_molecules_2d
  USE f77_blas
  USE input_constants,                 ONLY: use_orb_basis_set
  USE input_section_types,             ONLY: section_vals_val_get
  USE kinds,                           ONLY: dp
  USE machine,                         ONLY: default_output_unit,&
                                             m_flush,&
                                             m_memory
  USE message_passing,                 ONLY: mp_comm_split_direct,&
                                             mp_max,&
                                             mp_min,&
                                             mp_sendrecv,&
                                             mp_sum
  USE molecule_kind_types,             ONLY: molecule_kind_type
  USE molecule_types_new,              ONLY: molecule_type
  USE mp2_cphf,                        ONLY: solve_z_vector_eq
  USE mp2_ri_gpw,                      ONLY: mp2_ri_gpw_compute_en,&
                                             mp2_ri_gpw_compute_in
  USE mp2_ri_grad,                     ONLY: calc_ri_mp2_nonsep
  USE mp2_types,                       ONLY: mp2_type
  USE particle_methods,                ONLY: get_particle_set
  USE particle_types,                  ONLY: particle_type
  USE pw_env_methods,                  ONLY: pw_env_create,&
                                             pw_env_rebuild
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_release,&
                                             pw_env_type
  USE pw_methods,                      ONLY: pw_scale,&
                                             pw_transfer
  USE pw_poisson_methods,              ONLY: pw_poisson_solve
  USE pw_poisson_types,                ONLY: pw_poisson_type
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_p_type,&
                                             pw_release
  USE qs_collocate_density,            ONLY: calculate_wavefunction
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_integrate_potential,          ONLY: integrate_v_rspace
  USE qs_interactions,                 ONLY: init_interaction_radii
  USE qs_kind_types,                   ONLY: qs_kind_type
  USE qs_ks_types,                     ONLY: qs_ks_env_type
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
  USE qs_neighbor_list_types,          ONLY: deallocate_neighbor_list_set,&
                                             neighbor_list_set_p_type
  USE qs_neighbor_lists,               ONLY: atom2d_build,&
                                             atom2d_cleanup,&
                                             build_neighbor_lists,&
                                             local_atoms_type,&
                                             pair_radius_setup
  USE rpa_ri_gpw,                      ONLY: rpa_ri_compute_en
  USE task_list_methods,               ONLY: generate_qs_task_list
  USE task_list_types,                 ONLY: allocate_task_list,&
                                             deallocate_task_list,&
                                             task_list_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: get_limit
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mp2_gpw'

  PUBLIC :: mp2_gpw_main


  CONTAINS

! *****************************************************************************
!> \brief with a big bang to mp2
!> \param qs_env ...
!> \param mp2_env ...
!> \param Emp2 ...
!> \param Emp2_Cou ...
!> \param Emp2_EX ...
!> \param Emp2_S ...
!> \param Emp2_T ...
!> \param mos_mp2 ...
!> \param para_env ...
!> \param unit_nr ...
!> \param calc_forces ...
!> \param calc_ex ...
!> \param error ...
!> \param do_ri_mp2 ...
!> \param do_ri_rpa ...
!> \param do_ri_sos_laplace_mp2 ...
!> \author Mauro Del Ben and Joost VandeVondele
! *****************************************************************************
  SUBROUTINE mp2_gpw_main(qs_env,mp2_env,Emp2,Emp2_Cou,Emp2_EX,Emp2_S,Emp2_T,&
                          mos_mp2,para_env,unit_nr,calc_forces,calc_ex,error,do_ri_mp2,do_ri_rpa,&
                          do_ri_sos_laplace_mp2)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(mp2_type), POINTER                  :: mp2_env
    REAL(KIND=dp)                            :: Emp2, Emp2_Cou, Emp2_EX, &
                                                Emp2_S, Emp2_T
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos_mp2
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER                                  :: unit_nr
    LOGICAL, INTENT(IN)                      :: calc_forces
    LOGICAL                                  :: calc_ex
    TYPE(cp_error_type), INTENT(inout)       :: error
    LOGICAL, OPTIONAL                        :: do_ri_mp2, do_ri_rpa, &
                                                do_ri_sos_laplace_mp2

    CHARACTER(LEN=*), PARAMETER :: routineN = 'mp2_gpw_main', &
      routineP = moduleN//':'//routineN

    INTEGER :: blacs_grid_layout, color_sub, comm_sub, dimen, dimen_RI, &
      handle, homo, homo_beta, i, i_multigrid, local_unit_nr, my_group_L_end, &
      my_group_L_size, my_group_L_start, n_multigrid, natom, nelectron, &
      nelectron_beta, nkind, nmo, nspins, stat
    INTEGER, ALLOCATABLE, DIMENSION(:) :: ends_array, ends_B_virtual, &
      ends_B_virtual_beta, sizes_array, sizes_B_virtual, &
      sizes_B_virtual_beta, starts_array, starts_B_virtual, &
      starts_B_virtual_beta
    INTEGER, DIMENSION(:), POINTER           :: rbs
    LOGICAL :: blacs_repeatable, failure, my_do_ri_mp2, my_do_ri_rpa, &
      my_do_ri_sos_laplace_mp2, skip_load_balance_distributed
    LOGICAL, ALLOCATABLE, DIMENSION(:)       :: orb_present
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: orb_radius
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: pair_radius
    REAL(KIND=dp) :: cutoff_old, Emp2_AB, Emp2_BB, Emp2_Cou_BB, Emp2_d2_AB, &
      Emp2_d_AB, Emp2_EX_BB, eps_gvg_rspace_old, eps_pgf_orb_old, &
      eps_rho_rspace_old, progression_factor, relative_cutoff_old, subcells
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: e_cutoff_old, Eigenval, &
                                                Eigenval_beta
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: BIb_C, BIb_C_beta
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mo_eigenvalues
    TYPE(array_i1d_obj)                      :: row_blk_sizes
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env_sub
    TYPE(cp_dbcsr_p_type)                    :: mat_munu
    TYPE(cp_dbcsr_type), POINTER             :: mo_coeff_o, mo_coeff_o_beta, &
                                                mo_coeff_v, mo_coeff_v_beta
    TYPE(cp_error_type)                      :: error_sub
    TYPE(cp_fm_type), POINTER                :: mo_coeff, mo_coeff_beta
    TYPE(cp_logger_type), POINTER            :: logger, logger_sub
    TYPE(cp_para_env_type), POINTER          :: para_env_sub
    TYPE(dbcsr_distribution_obj), POINTER    :: dbcsr_dist_sub
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(distribution_1d_type), POINTER      :: local_molecules_sub, &
                                                local_particles_sub
    TYPE(distribution_2d_type), POINTER      :: distribution_2d_sub
    TYPE(local_atoms_type), ALLOCATABLE, &
      DIMENSION(:)                           :: atom2d
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb_sub
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env_sub
    TYPE(pw_p_type)                          :: pot_g, rho_g, rho_r
    TYPE(pw_poisson_type), POINTER           :: poisson_env
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(task_list_type), POINTER            :: task_list_sub

    CALL timeset(routineN,handle)
    failure=.FALSE.

    ! check if we want to do ri-mp2
    my_do_ri_mp2=.FALSE.
    IF(PRESENT(do_ri_mp2)) my_do_ri_mp2=do_ri_mp2

    ! check if we want to do ri-rpa
    my_do_ri_rpa=.FALSE.
    IF(PRESENT(do_ri_rpa)) my_do_ri_rpa=do_ri_rpa

    ! check if we want to do ri-sos-laplace-mp2
    my_do_ri_sos_laplace_mp2=.FALSE.
    IF(PRESENT(do_ri_sos_laplace_mp2)) my_do_ri_sos_laplace_mp2=do_ri_sos_laplace_mp2

    ! ... setup needed to be able to qs_integrate in a subgroup.

    ! a para_env
    color_sub=para_env%mepos/mp2_env%mp2_num_proc
    CALL mp_comm_split_direct(para_env%group,comm_sub,color_sub)
    NULLIFY(para_env_sub)
    CALL cp_para_env_create(para_env_sub,comm_sub,error=error)

    ! each of the sub groups might need to generate output
    logger => cp_error_get_logger(error)
    IF (para_env%mepos==para_env%source) THEN
       local_unit_nr=cp_logger_get_default_unit_nr(logger,local=.FALSE.)
    ELSE
       local_unit_nr=default_output_unit
    ENDIF

    IF (unit_nr>0) THEN
      WRITE (UNIT=unit_nr,FMT="(T3,A,T71,F10.1)")&
               "GPW_INFO| Density cutoff [a.u.]:", mp2_env%mp2_gpw%cutoff*0.5_dp
      WRITE (UNIT=unit_nr,FMT="(T3,A,T71,F10.1)")&
            "GPW_INFO| Relative density cutoff [a.u.]:", mp2_env%mp2_gpw%relative_cutoff*0.5_dp
      CALL m_flush(unit_nr)
    ENDIF

    ! a logger
    NULLIFY(logger_sub)
    CALL cp_logger_create(logger_sub,para_env=para_env_sub,&
               default_global_unit_nr=local_unit_nr, close_global_unit_on_dealloc=.FALSE.)
    CALL cp_logger_set(logger_sub,local_filename="MP2_localLog")
    ! set to a custom print level (we could also have a different print level for para_env%source)
    logger_sub%iter_info%print_level=mp2_env%mp2_gpw%print_level

    ! a error
    CALL cp_error_init(error_sub, stop_level=cp_failure_level, logger=logger_sub)

    ! a blacs_env (ignore the globenv stored defaults for now)
    blacs_grid_layout=BLACS_GRID_SQUARE
    blacs_repeatable=.TRUE.
    NULLIFY(blacs_env_sub)
    CALL cp_blacs_env_create(blacs_env_sub,para_env_sub,&
                             blacs_grid_layout,&
                             blacs_repeatable,error=error_sub)

    ! get stuff
    CALL get_qs_env(qs_env,&
                    ks_env=ks_env,&
                    qs_kind_set=qs_kind_set,&
                    cell=cell,&
                    particle_set=particle_set,&
                    atomic_kind_set=atomic_kind_set,&
                    molecule_set=molecule_set,&
                    molecule_kind_set=molecule_kind_set,&
                    dft_control=dft_control,&
                    error=error)

    ! re-init the radii to be able to generate pair lists with MP2-appropriate screening
    eps_pgf_orb_old=dft_control%qs_control%eps_pgf_orb
    eps_rho_rspace_old=dft_control%qs_control%eps_rho_rspace
    eps_gvg_rspace_old=dft_control%qs_control%eps_gvg_rspace
    ! hack hack hack XXXXXXXXXXXXXXX ... to be fixed
    dft_control%qs_control%eps_pgf_orb   =mp2_env%mp2_gpw%eps_grid
    dft_control%qs_control%eps_rho_rspace=mp2_env%mp2_gpw%eps_grid
    dft_control%qs_control%eps_gvg_rspace=mp2_env%mp2_gpw%eps_grid
    CALL init_interaction_radii(dft_control%qs_control,atomic_kind_set, qs_kind_set, error_sub)

    ! get a distribution_1d
    NULLIFY(local_particles_sub,local_molecules_sub)
    CALL distribute_molecules_1d(atomic_kind_set=atomic_kind_set,&
                                 particle_set=particle_set,&
                                 local_particles=local_particles_sub,&
                                 molecule_kind_set=molecule_kind_set,&
                                 molecule_set=molecule_set,&
                                 local_molecules=local_molecules_sub,&
                                 force_env_section=qs_env%input,&
                                 error=error_sub)

    ! get a distribution_2d
    NULLIFY(distribution_2d_sub)
    CALL distribute_molecules_2d(cell=cell,&
                                 atomic_kind_set=atomic_kind_set,&
                                 qs_kind_set=qs_kind_set,&
                                 particle_set=particle_set,&
                                 molecule_kind_set=molecule_kind_set,&
                                 molecule_set=molecule_set,&
                                 distribution_2d=distribution_2d_sub,&
                                 blacs_env=blacs_env_sub,&
                                 force_env_section=qs_env%input, error=error_sub)

    ! Build the sub  orbital-orbital overlap neighbor lists
    NULLIFY(sab_orb_sub)
    CALL section_vals_val_get(qs_env%input,"DFT%SUBCELLS",r_val=subcells,error=error)
    nkind = SIZE(atomic_kind_set)
    ALLOCATE (orb_present(nkind),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (orb_radius(nkind),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    orb_radius(:) = 0.0_dp
    ALLOCATE (pair_radius(nkind,nkind),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (atom2d(nkind),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL atom2d_build(atom2d,orb_radius,orb_present,local_particles_sub,distribution_2d_sub,&
                      atomic_kind_set,qs_kind_set,molecule_set,molecule_only=.FALSE.,dftb=.FALSE.,&
                      particle_set=particle_set,error=error)
    CALL pair_radius_setup(orb_present,orb_present,orb_radius,orb_radius,pair_radius,error)
    CALL build_neighbor_lists(sab_orb_sub,particle_set,atom2d,cell,pair_radius,&
                              mic=.FALSE.,subcells=subcells,molecular=.FALSE.,name="sab_orb_sub",error=error)
    CALL atom2d_cleanup(atom2d,error)
    DEALLOCATE (atom2d,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(orb_present,orb_radius,pair_radius,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ! a dbcsr_dist
    ALLOCATE(dbcsr_dist_sub,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL cp_dbcsr_dist2d_to_dist (distribution_2d_sub, dbcsr_dist_sub, error_sub)

    ! build a dbcsr matrix the hard way
    natom = SIZE(particle_set)
    ALLOCATE (rbs(natom), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL get_particle_set(particle_set, qs_kind_set, nsgf=rbs,error=error)
    CALL array_nullify (row_blk_sizes)
    CALL array_new (row_blk_sizes, rbs, gift=.TRUE.)
    ALLOCATE(mat_munu%matrix, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL cp_dbcsr_init(mat_munu%matrix,error=error_sub)
    CALL cp_dbcsr_create(matrix=mat_munu%matrix,&
               name="(ai|munu)",&
               dist=dbcsr_dist_sub, matrix_type=dbcsr_type_symmetric,&
               row_blk_size=row_blk_sizes, col_blk_size=row_blk_sizes, &
               nblks=0, nze=0, error=error_sub)
    CALL cp_dbcsr_alloc_block_from_nbl(mat_munu%matrix,sab_orb_sub,error_sub)
    CALL array_release (row_blk_sizes)

    ! and the array of mos
    CALL get_qs_env(qs_env=qs_env, dft_control=dft_control, error=error)
    CALL get_mo_set(mo_set=mos_mp2(1)%mo_set, nelectron=nelectron,&
                    eigenvalues=mo_eigenvalues,nmo=nmo,homo=homo,&
                    mo_coeff=mo_coeff,nao=dimen)
    ALLOCATE(Eigenval(dimen),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    Eigenval(:)=mo_eigenvalues(:)

    ! new routine: replicate a full matrix from one para_env to a smaller one
    ! keeping the memory usage as small as possible in this case the 
    ! output the two part of the C matrix (virtual, occupied)
    CALL replicate_mat_to_subgroup(qs_env,mp2_env,para_env,para_env_sub,mo_coeff,dimen,homo,mat_munu,&
                                   mo_coeff_o,mo_coeff_v,blacs_env_sub,unit_nr,error,error_sub)

    ! if open shell case replicate also the coefficient matrix for the beta orbitals 
    nspins=SIZE(mos_mp2)
    IF(nspins==2) THEN
      CALL get_mo_set(mo_set=mos_mp2(2)%mo_set, nelectron=nelectron_beta,&
                      eigenvalues=mo_eigenvalues,homo=homo_beta,&
                      mo_coeff=mo_coeff_beta)
      ALLOCATE(Eigenval_beta(dimen),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      Eigenval_beta(:)=mo_eigenvalues(:)

      CALL replicate_mat_to_subgroup(qs_env,mp2_env,para_env,para_env_sub,mo_coeff_beta,dimen,homo_beta,mat_munu,&
                                     mo_coeff_o_beta,mo_coeff_v_beta,blacs_env_sub,unit_nr,error,error_sub)
    END IF

    ! hack hack hack XXXXXXXXXXXXXXX rebuilds the pw_en with the new cutoffs
    progression_factor=dft_control%qs_control%progression_factor
    n_multigrid=SIZE(dft_control%qs_control%e_cutoff)
    ALLOCATE(e_cutoff_old(n_multigrid), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    e_cutoff_old(:)=dft_control%qs_control%e_cutoff
    cutoff_old=dft_control%qs_control%cutoff

    dft_control%qs_control%cutoff=mp2_env%mp2_gpw%cutoff*0.5_dp
    dft_control%qs_control%e_cutoff(1)=dft_control%qs_control%cutoff
    DO i_multigrid=2, n_multigrid
      dft_control%qs_control%e_cutoff(i_multigrid) = dft_control%qs_control%e_cutoff(i_multigrid-1) &
          / progression_factor
    END DO

    relative_cutoff_old=dft_control%qs_control%relative_cutoff
    dft_control%qs_control%relative_cutoff=mp2_env%mp2_gpw%relative_cutoff*0.5_dp

    ! a pw_env
    NULLIFY(pw_env_sub)
    CALL pw_env_create(pw_env_sub,error_sub)
    CALL pw_env_rebuild(pw_env_sub,qs_env,para_env_sub,error_sub)

    CALL pw_env_get(pw_env_sub, auxbas_pw_pool=auxbas_pw_pool,&
                    poisson_env=poisson_env,error=error_sub)
    ! hack hack hack XXXXXXXXXXXXXXX

    ! now we need a task list, hard code skip_load_balance_distributed
    NULLIFY(task_list_sub)
    skip_load_balance_distributed=dft_control%qs_control%skip_load_balance_distributed
    CALL allocate_task_list(task_list_sub,error_sub)
    CALL generate_qs_task_list(ks_env, task_list_sub, &
                 reorder_rs_grid_ranks=.TRUE., soft_valid=.FALSE., &
                 skip_load_balance_distributed=skip_load_balance_distributed,&
                 pw_env_external=pw_env_sub, sab_orb_external=sab_orb_sub, error=error_sub)

    ! get some of the grids ready
    NULLIFY(rho_r%pw,rho_g%pw,pot_g%pw)
    CALL pw_pool_create_pw(auxbas_pw_pool,rho_r%pw,&
                            use_data=REALDATA3D,&
                            in_space=REALSPACE,error=error_sub)
    CALL pw_pool_create_pw(auxbas_pw_pool,rho_g%pw,&
                            use_data=COMPLEXDATA1D,&
                            in_space=RECIPROCALSPACE,error=error_sub)
    CALL pw_pool_create_pw(auxbas_pw_pool,pot_g%pw,&
                            use_data=COMPLEXDATA1D,&
                            in_space=RECIPROCALSPACE,error=error_sub)

    ! run the FFT once, to set up buffers and to take into account the memory
    rho_r%pw%cr3d = 0.0D0
    CALL pw_transfer(rho_r%pw, rho_g%pw, error=error_sub)

    ! now we're kind of ready to go....
    Emp2_S=0.0_dp
    Emp2_T=0.0_dp
    IF(my_do_ri_mp2.OR.my_do_ri_rpa.OR.my_do_ri_sos_laplace_mp2) THEN
      ! RI-GPW integrals (same stuff for both RPA and MP2)
      IF(nspins==2) THEN
        ! open shell case (RI) here the (ia|K) integrals are computed for both the alpha and beta components
        CALL mp2_ri_gpw_compute_in(BIb_C,ends_array,ends_B_virtual,sizes_array,sizes_B_virtual,starts_array,starts_B_virtual,&
                                   dimen_RI,qs_env,para_env,para_env_sub,color_sub,dft_control,cell,particle_set,&
                                   atomic_kind_set, qs_kind_set,mo_coeff,nmo,homo,rho_r,rho_g,pot_g,&
                                   mat_munu,sab_orb_sub,pw_env_sub,poisson_env,auxbas_pw_pool,task_list_sub,&
                                   mo_coeff_o,mo_coeff_v,mp2_env%mp2_gpw%eps_filter,unit_nr,&
                                   mp2_env%mp2_memory,mp2_env%calc_PQ_cond_num,calc_forces,error,blacs_env_sub,error_sub,&
                                   BIb_C_beta,ends_B_virtual_beta,sizes_B_virtual_beta,starts_B_virtual_beta,&
                                   homo_beta,mo_coeff_o_beta,mo_coeff_v_beta)
      ELSE
        ! closed shell case (RI)
        CALL mp2_ri_gpw_compute_in(BIb_C,ends_array,ends_B_virtual,sizes_array,sizes_B_virtual,starts_array,starts_B_virtual,&
                                   dimen_RI,qs_env,para_env,para_env_sub,color_sub,dft_control,cell,particle_set,&
                                   atomic_kind_set, qs_kind_set,mo_coeff,nmo,homo,rho_r,rho_g,pot_g,&
                                   mat_munu,sab_orb_sub,pw_env_sub,poisson_env,auxbas_pw_pool,task_list_sub,&
                                   mo_coeff_o,mo_coeff_v,mp2_env%mp2_gpw%eps_filter,unit_nr,&
                                   mp2_env%mp2_memory,mp2_env%calc_PQ_cond_num,calc_forces,error,blacs_env_sub,error_sub)
      END IF
    ELSE
      ! Canonical MP2-GPW
      IF(nspins==2) THEN
        ! alpha-alpha and alpha-beta components
        IF (unit_nr>0) WRITE(unit_nr,*)
        IF (unit_nr>0) WRITE(unit_nr,'(T3,A)') 'Alpha (ia|'
        CALL mp2_gpw_compute(Emp2,Emp2_Cou,Emp2_EX,qs_env,para_env,para_env_sub,color_sub,dft_control,cell,particle_set,&
                             atomic_kind_set,qs_kind_set,mo_coeff,Eigenval,nmo,homo,rho_r,rho_g,pot_g,&
                             mat_munu,sab_orb_sub,pw_env_sub,poisson_env,auxbas_pw_pool,task_list_sub,&
                             mo_coeff_o,mo_coeff_v,mp2_env%mp2_gpw%eps_filter,unit_nr,&
                             mp2_env%mp2_memory,calc_ex,blacs_env_sub,error_sub,&
                             homo_beta,mo_coeff_o_beta,mo_coeff_v_beta,Eigenval_beta,Emp2_AB)

        ! beta-beta component
        IF (unit_nr>0) WRITE(unit_nr,*)
        IF (unit_nr>0) WRITE(unit_nr,'(T3,A)') 'Beta (ia|'
        CALL mp2_gpw_compute(Emp2_BB,Emp2_Cou_BB,Emp2_EX_BB,qs_env,para_env,para_env_sub,color_sub,dft_control,cell,particle_set,&
                             atomic_kind_set,qs_kind_set,mo_coeff_beta,Eigenval_beta,nmo,homo_beta,rho_r,rho_g,pot_g,&
                             mat_munu,sab_orb_sub,pw_env_sub,poisson_env,auxbas_pw_pool,task_list_sub,&
                             mo_coeff_o_beta,mo_coeff_v_beta,mp2_env%mp2_gpw%eps_filter,unit_nr,&
                             mp2_env%mp2_memory,calc_ex,blacs_env_sub,error_sub)
        
      ELSE 
        ! closed shell case
        CALL mp2_gpw_compute(Emp2,Emp2_Cou,Emp2_EX,qs_env,para_env,para_env_sub,color_sub,dft_control,cell,particle_set,&
                             atomic_kind_set,qs_kind_set,mo_coeff,Eigenval,nmo,homo,rho_r,rho_g,pot_g,&
                             mat_munu,sab_orb_sub,pw_env_sub,poisson_env,auxbas_pw_pool,task_list_sub,&
                             mo_coeff_o,mo_coeff_v,mp2_env%mp2_gpw%eps_filter,unit_nr,&
                             mp2_env%mp2_memory,calc_ex,blacs_env_sub,error_sub)
      END IF
    END IF

    ! and now free the whole lot
    CALL pw_pool_give_back_pw(auxbas_pw_pool,rho_r%pw,error=error_sub)
    CALL pw_pool_give_back_pw(auxbas_pw_pool,rho_g%pw,error=error_sub)
    CALL pw_pool_give_back_pw(auxbas_pw_pool,pot_g%pw,error=error_sub)

    ! Free possibly large buffers allocated by dbcsr on the GPU,
    ! large hybrid dgemm/pdgemm's coming later will need the space.
    CALL cp_dbcsr_clear_mempools(error_sub)

    ! moved down
    ! CALL deallocate_task_list(task_list_sub,error=error_sub)
    ! CALL pw_env_release(pw_env_sub, error=error_sub)

    IF(calc_forces) THEN
      ! make a copy of mo_coeff_o and mo_coeff_v
      NULLIFY(mp2_env%ri_grad%mo_coeff_o)
      CALL cp_dbcsr_init_p(mp2_env%ri_grad%mo_coeff_o,error=error_sub)
      CALL cp_dbcsr_copy(mp2_env%ri_grad%mo_coeff_o,mo_coeff_o,name="mo_coeff_o",error=error_sub)
      NULLIFY(mp2_env%ri_grad%mo_coeff_v)
      CALL cp_dbcsr_init_p(mp2_env%ri_grad%mo_coeff_v,error=error_sub)
      CALL cp_dbcsr_copy(mp2_env%ri_grad%mo_coeff_v,mo_coeff_v,name="mo_coeff_v",error=error_sub)
      my_group_L_size=sizes_array(color_sub)
      my_group_L_start=starts_array(color_sub)
      my_group_L_end=ends_array(color_sub)
    END IF

    CALL cp_dbcsr_release(mo_coeff_o,error=error_sub)
    DEALLOCATE(mo_coeff_o, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL cp_dbcsr_release(mo_coeff_v,error=error_sub)
    DEALLOCATE(mo_coeff_v, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF(nspins==2) THEN
      CALL cp_dbcsr_release(mo_coeff_o_beta,error=error_sub)
      DEALLOCATE(mo_coeff_o_beta, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      CALL cp_dbcsr_release(mo_coeff_v_beta,error=error_sub)
      DEALLOCATE(mo_coeff_v_beta, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

    IF(.NOT.calc_forces) THEN
      ! release stuff
      CALL cp_dbcsr_release(mat_munu%matrix,error=error_sub)
      DEALLOCATE(mat_munu%matrix, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

      CALL cp_dbcsr_distribution_release(dbcsr_dist_sub)
      DEALLOCATE(dbcsr_dist_sub,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

      DO i=1,SIZE(sab_orb_sub)
         CALL deallocate_neighbor_list_set(sab_orb_sub(i)%neighbor_list_set)
      END DO
      DEALLOCATE(sab_orb_sub,stat=stat)
      CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)

      CALL distribution_2d_release(distribution_2d_sub,error=error_sub)

      CALL distribution_1d_release(local_particles_sub,error=error_sub)
      CALL distribution_1d_release(local_molecules_sub,error=error_sub)
    END IF

    ! decide if to doing RI-RPA or RI-MP2
    IF(my_do_ri_rpa.OR.my_do_ri_sos_laplace_mp2) THEN
      ! RI-RPA
      IF(nspins==2) THEN
        CALL rpa_ri_compute_en(Emp2,mp2_env,BIb_C,qs_env,para_env,para_env_sub,color_sub,&
                               ends_array,ends_B_virtual,sizes_array,sizes_B_virtual,starts_array,starts_B_virtual,&
                               Eigenval,nmo,homo,dimen_RI,unit_nr,error,error_sub,my_do_ri_sos_laplace_mp2,&
                               BIb_C_beta,homo_beta,Eigenval_beta,&
                               ends_B_virtual_beta,sizes_B_virtual_beta, starts_B_virtual_beta)
      ELSE
        CALL rpa_ri_compute_en(Emp2,mp2_env,BIb_C,qs_env,para_env,para_env_sub,color_sub,&
                               ends_array,ends_B_virtual,sizes_array,sizes_B_virtual,starts_array,starts_B_virtual,&
                               Eigenval,nmo,homo,dimen_RI,unit_nr,error,error_sub,my_do_ri_sos_laplace_mp2)
      END IF
    ELSE
      IF(my_do_ri_mp2) THEN
        ! RI-MP2-GPW compute energy
        IF(nspins==2) THEN
          ! alpha-alpha component
          CALL mp2_ri_gpw_compute_en(Emp2,Emp2_Cou,Emp2_EX,BIb_C,qs_env,mp2_env,para_env,para_env_sub,color_sub,&
                                     ends_array,ends_B_virtual,sizes_array,sizes_B_virtual,starts_array,starts_B_virtual,&
                                     Eigenval,nmo,homo,dimen_RI,unit_nr,calc_forces,calc_ex,error,error_sub,&
                                     open_shell_SS=.TRUE.)
 
          ! beta-beta component
          CALL mp2_ri_gpw_compute_en(Emp2_BB,Emp2_Cou_BB,Emp2_EX_BB,BIb_C_beta,qs_env,mp2_env,para_env,para_env_sub,color_sub,&
                                     ends_array,ends_B_virtual_beta,sizes_array,&
                                     sizes_B_virtual_beta,starts_array,starts_B_virtual_beta,&
                                     Eigenval_beta,nmo,homo_beta,dimen_RI,unit_nr,calc_forces,calc_ex,error,error_sub,&
                                     open_shell_SS=.TRUE.)


          ! alpha-beta case
          CALL mp2_ri_gpw_compute_en(Emp2_d_AB,Emp2_AB,Emp2_d2_AB,BIb_C,qs_env,mp2_env,para_env,para_env_sub,color_sub,&
                                     ends_array,ends_B_virtual,sizes_array,sizes_B_virtual,starts_array,starts_B_virtual,&
                                     Eigenval,nmo,homo,dimen_RI,unit_nr,calc_forces,.FALSE.,error,error_sub,&
                                     .FALSE.,BIb_C_beta,homo_beta,Eigenval_beta,&
                                     ends_B_virtual_beta,sizes_B_virtual_beta, starts_B_virtual_beta)

        ELSE
          ! closed shell case
          CALL mp2_ri_gpw_compute_en(Emp2,Emp2_Cou,Emp2_EX,BIb_C,qs_env,mp2_env,para_env,para_env_sub,color_sub,&
                                     ends_array,ends_B_virtual,sizes_array,sizes_B_virtual,starts_array,starts_B_virtual,&
                                     Eigenval,nmo,homo,dimen_RI,unit_nr,calc_forces,calc_ex,error,error_sub)
        END IF
        ! if we need forces time to calculate the MP2 non-separable contribution 
        ! and start coputing the largrangian
        IF(calc_forces) THEN
          ! since we have to compute again integrals reinitialize the stuff we need
          ! get some of the grids ready
          NULLIFY(rho_r%pw,rho_g%pw,pot_g%pw)
          CALL pw_pool_create_pw(auxbas_pw_pool,rho_r%pw,&
                                  use_data=REALDATA3D,&
                                  in_space=REALSPACE,error=error_sub)
          CALL pw_pool_create_pw(auxbas_pw_pool,rho_g%pw,&
                                  use_data=COMPLEXDATA1D,&
                                  in_space=RECIPROCALSPACE,error=error_sub)
          CALL pw_pool_create_pw(auxbas_pw_pool,pot_g%pw,&
                                  use_data=COMPLEXDATA1D,&
                                  in_space=RECIPROCALSPACE,error=error_sub)

          ! the mu_nu matrix (again)
          ! XXXXXXXXXXXXXXXXXXXXXXXXX
          ! ! build a dbcsr matrix the hard way
          ! CALL get_particle_set(particle_set=particle_set,nsgf=rbs,error=error)
          ! CALL array_nullify (row_blk_sizes)
          ! CALL array_new (row_blk_sizes, rbs, gift=.TRUE.)
          ! ALLOCATE(mat_munu%matrix, STAT=stat)
          ! CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ! CALL cp_dbcsr_init(mat_munu%matrix,error=error_sub)
          ! CALL cp_dbcsr_create(matrix=mat_munu%matrix,&
          !            name="(ai|munu)",&
          !            dist=dbcsr_dist_sub, &
          !            matrix_type=dbcsr_type_symmetric,&
          !            ! matrix_type=dbcsr_type_no_symmetry,&
          !            row_blk_size=row_blk_sizes, col_blk_size=row_blk_sizes, &
          !            nblks=0, nze=0, error=error_sub)
          ! CALL cp_dbcsr_alloc_block_from_nbl(mat_munu%matrix,sab_orb_sub,error_sub)
          ! CALL array_release (row_blk_sizes)

          CALL  calc_ri_mp2_nonsep(qs_env,mp2_env,para_env,para_env_sub,color_sub,dft_control,cell,particle_set,&
                                   atomic_kind_set,qs_kind_set,mo_coeff,nmo,homo,dimen_RI,Eigenval,&
                                   my_group_L_start,my_group_L_end,my_group_L_size,rho_r,rho_g,pot_g,&
                                   mat_munu,sab_orb_sub,pw_env_sub,poisson_env,auxbas_pw_pool,task_list_sub,&
                                   unit_nr,blacs_env_sub,error,error_sub)

          ! release
          CALL pw_pool_give_back_pw(auxbas_pw_pool,rho_r%pw,error=error_sub)
          CALL pw_pool_give_back_pw(auxbas_pw_pool,rho_g%pw,error=error_sub)
          CALL pw_pool_give_back_pw(auxbas_pw_pool,pot_g%pw,error=error_sub)

          CALL cp_dbcsr_release(mp2_env%ri_grad%mo_coeff_o,error=error_sub)
          DEALLOCATE(mp2_env%ri_grad%mo_coeff_o, STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

          CALL cp_dbcsr_release(mp2_env%ri_grad%mo_coeff_v,error=error_sub)
          DEALLOCATE(mp2_env%ri_grad%mo_coeff_v, STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

          CALL cp_dbcsr_release(mat_munu%matrix,error=error_sub)
          DEALLOCATE(mat_munu%matrix, STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

          CALL cp_dbcsr_distribution_release(dbcsr_dist_sub)
          DEALLOCATE(dbcsr_dist_sub,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

          DO i=1,SIZE(sab_orb_sub)
             CALL deallocate_neighbor_list_set(sab_orb_sub(i)%neighbor_list_set)
          END DO
          DEALLOCATE(sab_orb_sub,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)

          CALL distribution_2d_release(distribution_2d_sub,error=error_sub)

          CALL distribution_1d_release(local_particles_sub,error=error_sub)
          CALL distribution_1d_release(local_molecules_sub,error=error_sub)

        END IF
      END IF

      IF(nspins==2) THEN
        ! make order on the MP2 energy contributions
        Emp2_Cou=Emp2_Cou*0.25_dp
        Emp2_EX=Emp2_EX*0.5_dp

        Emp2_Cou_BB=Emp2_Cou_BB*0.25_dp
        Emp2_EX_BB=Emp2_EX_BB*0.5_dp

        Emp2_S=Emp2_AB
        Emp2_T=Emp2_Cou+Emp2_Cou_BB+Emp2_EX+Emp2_EX_BB

        Emp2_Cou=Emp2_Cou+Emp2_Cou_BB+Emp2_AB
        Emp2_EX=Emp2_EX+Emp2_EX_BB
        Emp2=Emp2_EX+Emp2_Cou
      END IF

    END IF

    !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXx
    ! moved from above
    CALL deallocate_task_list(task_list_sub,error=error_sub)

    CALL pw_env_release(pw_env_sub, error=error_sub)

    ! CALL cp_dbcsr_distribution_release(dbcsr_dist_sub)
    ! DEALLOCATE(dbcsr_dist_sub,STAT=stat)
    ! CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ! DO i=1,SIZE(sab_orb_sub)
    !    CALL deallocate_neighbor_list_set(sab_orb_sub(i)%neighbor_list_set)
    ! END DO
    ! DEALLOCATE(sab_orb_sub,stat=stat)
    ! CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)

    ! CALL distribution_2d_release(distribution_2d_sub,error=error_sub)

    ! CALL distribution_1d_release(local_particles_sub,error=error_sub)
    ! CALL distribution_1d_release(local_molecules_sub,error=error_sub)
    !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXx

    ! re-init the radii to be able to generate pair lists with MP2-appropriate screening
    dft_control%qs_control%eps_pgf_orb=eps_pgf_orb_old
    dft_control%qs_control%eps_rho_rspace=eps_rho_rspace_old
    dft_control%qs_control%eps_gvg_rspace=eps_gvg_rspace_old
    CALL init_interaction_radii(dft_control%qs_control,atomic_kind_set, qs_kind_set, error_sub)

    ! restore the initial value of the cutoff
    dft_control%qs_control%e_cutoff=e_cutoff_old
    dft_control%qs_control%cutoff=cutoff_old
    dft_control%qs_control%relative_cutoff=relative_cutoff_old

    CALL cp_blacs_env_release(blacs_env_sub, error=error_sub)

    CALL cp_error_dealloc_ref(error_sub,error=error)

    CALL cp_logger_release(logger_sub)

    CALL cp_para_env_release(para_env_sub,error=error)

    ! finally solve the z-vector equation if forces are required
    IF(calc_forces) THEN
      CALL solve_z_vector_eq(qs_env,mp2_env,para_env,dft_control,cell,particle_set,&
                             atomic_kind_set,mo_coeff,nmo,homo,Eigenval,unit_nr,error)
    END IF   

    DEALLOCATE(Eigenval,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF(nspins==2) THEN
      DEALLOCATE(Eigenval_beta,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

    CALL timestop(handle)

  END SUBROUTINE mp2_gpw_main

! *****************************************************************************
!> \brief ...
!> \param Emp2 ...
!> \param Emp2_Cou ...
!> \param Emp2_EX ...
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
!> \param Eigenval ...
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
!> \param eps_filter ...
!> \param unit_nr ...
!> \param mp2_memory ...
!> \param calc_ex ...
!> \param blacs_env_sub ...
!> \param error_sub ...
!> \param homo_beta ...
!> \param mo_coeff_o_beta ...
!> \param mo_coeff_v_beta ...
!> \param Eigenval_beta ...
!> \param Emp2_AB ...
! *****************************************************************************
  SUBROUTINE mp2_gpw_compute(Emp2,Emp2_Cou,Emp2_EX,qs_env,para_env,para_env_sub,color_sub,dft_control,&
                 cell,particle_set,atomic_kind_set,qs_kind_set,mo_coeff,Eigenval,nmo,homo,&
                 rho_r,rho_g,pot_g,mat_munu,sab_orb_sub,pw_env_sub,&
                 poisson_env,auxbas_pw_pool,task_list_sub,mo_coeff_o,mo_coeff_v,eps_filter,unit_nr,&
                 mp2_memory,calc_ex,blacs_env_sub,error_sub,homo_beta,mo_coeff_o_beta,mo_coeff_v_beta,Eigenval_beta,Emp2_AB)

    REAL(KIND=dp)                            :: Emp2, Emp2_Cou, Emp2_EX
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
    REAL(KIND=dp), DIMENSION(:)              :: Eigenval
    INTEGER                                  :: nmo, homo
    TYPE(pw_p_type)                          :: rho_r, rho_g, pot_g
    TYPE(cp_dbcsr_p_type)                    :: mat_munu
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb_sub
    TYPE(pw_env_type), POINTER               :: pw_env_sub
    TYPE(pw_poisson_type), POINTER           :: poisson_env
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(task_list_type), POINTER            :: task_list_sub
    TYPE(cp_dbcsr_type), POINTER             :: mo_coeff_o, mo_coeff_v
    REAL(KIND=dp)                            :: eps_filter
    INTEGER                                  :: unit_nr
    REAL(KIND=dp)                            :: mp2_memory
    LOGICAL                                  :: calc_ex
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env_sub
    TYPE(cp_error_type), INTENT(inout)       :: error_sub
    INTEGER, OPTIONAL                        :: homo_beta
    TYPE(cp_dbcsr_type), OPTIONAL, POINTER   :: mo_coeff_o_beta, &
                                                mo_coeff_v_beta
    REAL(KIND=dp), DIMENSION(:), OPTIONAL    :: Eigenval_beta
    REAL(KIND=dp), OPTIONAL                  :: Emp2_AB

    CHARACTER(LEN=*), PARAMETER :: routineN = 'mp2_gpw_compute', &
      routineP = moduleN//':'//routineN

    INTEGER :: a, a_group_counter, b, b_global, b_group_counter, blk, col, &
      col_offset, col_size, color_counter, comm_exchange, EX_end, &
      EX_end_send, EX_start, EX_start_send, group_counter, handle, handle2, &
      handle3, i, i_counter, i_group_counter, index_proc_shift, j, &
      max_b_size, max_batch_size_A, max_batch_size_I, max_row_col_local, &
      mepos_in_EX_group, my_A_batch_size, my_A_virtual_end, &
      my_A_virtual_start, my_B_size, my_B_virtual_end, my_B_virtual_start, &
      my_I_batch_size, my_I_occupied_end, my_I_occupied_start, my_q_position, &
      ncol_local, nfullcols_total, nfullrows_total, ngroup, nrow_local, one, p
    INTEGER :: p_best, proc_receive, proc_send, q, q_best, row, row_offset, &
      row_size, size_EX, size_EX_send, size_of_exchange_group, stat, &
      sub_sub_color, virtual, virtual_beta, wfn_calc, wfn_calc_best
    INTEGER, ALLOCATABLE, DIMENSION(:) :: proc_map, sub_proc_map, &
      vector_B_sizes, vector_batch_A_size_group, vector_batch_I_size_group
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: color_array, &
                                                exchange_group_sizes, &
                                                local_col_row_info
    INTEGER, DIMENSION(:), POINTER           :: col_indices, row_indices
    LOGICAL                                  :: do_alpha_beta, failure
    REAL(KIND=dp)                            :: mem_min, mem_real, mem_try, &
                                                pair_energy, wfn_size
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: my_Cocc, my_Cvirt
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: BIb_C, BIb_Ex, BIb_send
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: data_block
    TYPE(array_i1d_obj)                      :: col_blk_size, col_dist, &
                                                row_blk_size, row_dist
    TYPE(cp_dbcsr_iterator)                  :: iter
    TYPE(cp_dbcsr_type)                      :: matrix_ia_jb, &
                                                matrix_ia_jb_beta, &
                                                matrix_ia_jnu, &
                                                matrix_ia_jnu_beta
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct
    TYPE(cp_fm_type), POINTER                :: fm_BIb_jb
    TYPE(cp_para_env_type), POINTER          :: para_env_exchange
    TYPE(dbcsr_distribution_obj)             :: dist
    TYPE(pw_p_type)                          :: psi_a
    TYPE(pw_p_type), ALLOCATABLE, &
      DIMENSION(:)                           :: psi_i

    CALL timeset(routineN,handle)
    failure=.FALSE.

    do_alpha_beta=.FALSE.
    IF(PRESENT(homo_beta).AND.&
       PRESENT(mo_coeff_o_beta).AND.&
       PRESENT(mo_coeff_v_beta).AND.&
       PRESENT(Eigenval_beta).AND.&
       PRESENT(Emp2_AB)) do_alpha_beta=.TRUE.

    ! initialize and create the matrix (ia|jnu)
    CALL cp_dbcsr_init(matrix_ia_jnu,error=error_sub)
    CALL cp_dbcsr_create(matrix_ia_jnu,template=mo_coeff_o,error=error_sub)

    ! Allocate Sparse matrices: (ia|jb)
    CALL cp_create_bl_distribution (col_dist, col_blk_size, nmo-homo, &
          dbcsr_mp_npcols(dbcsr_distribution_mp(cp_dbcsr_distribution(mo_coeff_v))))
    CALL cp_create_bl_distribution (row_dist, row_blk_size, homo, &
          dbcsr_mp_nprows(dbcsr_distribution_mp(cp_dbcsr_distribution(mo_coeff_o))))
    CALL dbcsr_distribution_new (dist, dbcsr_distribution_mp(cp_dbcsr_distribution(mo_coeff_o)),&
          row_dist,col_dist)
    CALL cp_dbcsr_init(matrix_ia_jb,error=error_sub)
    CALL cp_dbcsr_create(matrix_ia_jb,"matrix_ia_jb",dist,dbcsr_type_no_symmetry,&
         row_blk_size,col_blk_size,0,0,dbcsr_type_real_default,error=error_sub)
    CALL cp_dbcsr_distribution_release (dist)
    CALL array_release (col_blk_size)
    CALL array_release (col_dist)
    CALL array_release (row_blk_size)
    CALL array_release (row_dist)

    ! set all to zero in such a way that the memory is actually allocated
    CALL cp_dbcsr_set(matrix_ia_jnu,0.0_dp,error=error_sub)
    CALL cp_dbcsr_set(matrix_ia_jb,0.0_dp,error=error_sub)
    CALL cp_dbcsr_set(mat_munu%matrix,0.0_dp,error=error_sub)

    IF(calc_ex) THEN
      ! create the analogous of matrix_ia_jb in fm type
      NULLIFY(fm_BIb_jb)
      NULLIFY(fm_struct)
      CALL cp_dbcsr_get_info(matrix_ia_jb,nfullrows_total=nfullrows_total,nfullcols_total=nfullcols_total)
      CALL cp_fm_struct_create(fm_struct,context=blacs_env_sub,nrow_global=nfullrows_total,&
                               ncol_global=nfullcols_total,para_env=para_env_sub,error=error_sub)
      CALL cp_fm_create(fm_BIb_jb,fm_struct,name="fm_BIb_jb",error=error_sub)

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
    END IF

    IF(do_alpha_beta) THEN
      ! initialize and create the matrix (ia|jnu)
      CALL cp_dbcsr_init(matrix_ia_jnu_beta,error=error_sub)
      CALL cp_dbcsr_create(matrix_ia_jnu_beta,template=mo_coeff_o_beta,error=error_sub)

      ! Allocate Sparse matrices: (ia|jb)
      CALL cp_create_bl_distribution (col_dist, col_blk_size, nmo-homo_beta, &
            dbcsr_mp_npcols(dbcsr_distribution_mp(cp_dbcsr_distribution(mo_coeff_v_beta))))
      CALL cp_create_bl_distribution (row_dist, row_blk_size, homo_beta, &
            dbcsr_mp_nprows(dbcsr_distribution_mp(cp_dbcsr_distribution(mo_coeff_o_beta))))
      CALL dbcsr_distribution_new (dist, dbcsr_distribution_mp(cp_dbcsr_distribution(mo_coeff_o_beta)),&
            row_dist,col_dist)
      CALL cp_dbcsr_init(matrix_ia_jb_beta,error=error_sub)
      CALL cp_dbcsr_create(matrix_ia_jb_beta,"matrix_ia_jb_beta",dist,dbcsr_type_no_symmetry,&
           row_blk_size,col_blk_size,0,0,dbcsr_type_real_default,error=error_sub)
      CALL cp_dbcsr_distribution_release (dist)
      CALL array_release (col_blk_size)
      CALL array_release (col_dist)
      CALL array_release (row_blk_size)
      CALL array_release (row_dist)
      virtual_beta=nmo-homo_beta

      CALL cp_dbcsr_set(matrix_ia_jnu_beta,0.0_dp,error=error_sub)
      CALL cp_dbcsr_set(matrix_ia_jb_beta,0.0_dp,error=error_sub)
    END IF

    mem_real=m_memory()
    mem_real=(mem_real+1024*1024-1)/(1024*1024)
    ! mp_min .... a hack.. it should be mp_max, but as it turns out, on some processes the previously freed memory (hfx)
    ! has not been given back to the OS yet. 
    CALL mp_min(mem_real,para_env%group)

    virtual=nmo-homo

    wfn_size=REAL(SIZE(rho_r%pw%cr3d),KIND=dp)
    CALL mp_max(wfn_size,para_env%group)

    ngroup=para_env%num_pe/para_env_sub%num_pe

    ! calculate the minimal memory required per MPI task (p=occupied division,q=virtual division)
    p_best=ngroup
    q_best=1
    mem_min=HUGE(0)
    DO p=1,ngroup
       q=ngroup/p
       IF (p*q.NE.ngroup) CYCLE

       CALL estimate_memory_usage(wfn_size, p, q, para_env_sub%num_pe, nmo, virtual, homo, calc_ex, mem_try)

       IF(mem_try<=mem_min) THEN
         mem_min=mem_try
         p_best=p
         q_best=q
       END IF
    END DO
    IF (unit_nr>0) WRITE(unit_nr,'(T3,A,T69,F9.2,A3)') 'Minimum required memory per MPI process for MP2:',&
                                                       mem_min, ' MB'

    mem_real=mp2_memory-mem_real
    mem_real=MAX(mem_real,mem_min)
    IF (unit_nr>0) WRITE(unit_nr,'(T3,A,T69,F9.2,A3)') 'Available memory per MPI process for MP2:',&
                                                    mem_real, ' MB'

    wfn_calc_best=HUGE(wfn_calc_best)
    DO p=1,ngroup
       q=ngroup/p
       IF (p*q.NE.ngroup) CYCLE

       CALL estimate_memory_usage(wfn_size, p, q, para_env_sub%num_pe, nmo, virtual, homo, calc_ex, mem_try)

       IF(mem_try>mem_real) CYCLE
       wfn_calc=((homo+p-1)/p)+((virtual+q-1)/q)
       IF (wfn_calc<wfn_calc_best) THEN
          wfn_calc_best=wfn_calc
          p_best=p
          q_best=q
       ENDIF
    ENDDO

    max_batch_size_I=(homo+p_best-1)/p_best
    max_batch_size_A=(virtual+q_best-1)/q_best

    IF (unit_nr>0) THEN
      WRITE (UNIT=unit_nr,FMT="(T3,A,T77,i4)")&
               "MP2_GPW| max. batch size for the occupied states:", max_batch_size_I
      WRITE (UNIT=unit_nr,FMT="(T3,A,T77,i4)")&
            "MP2_GPW| max. batch size for the virtual states:",max_batch_size_A
    ENDIF

    ALLOCATE(vector_batch_I_size_group(0:p_best-1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
    ALLOCATE(vector_batch_A_size_group(0:q_best-1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)

    vector_batch_I_size_group=max_batch_size_I
    IF(SUM(vector_batch_I_size_group)/=homo) THEN
      one=1
      IF(SUM(vector_batch_I_size_group)>homo) one=-1
        i=-1
        DO
          i=i+1
          vector_batch_I_size_group(i)=vector_batch_I_size_group(i)+one
          IF(SUM(vector_batch_I_size_group)==homo) EXIT
          IF(i==p_best-1) i=-1
        END DO
    END IF

    vector_batch_A_size_group=max_batch_size_A
    IF(SUM(vector_batch_A_size_group)/=virtual) THEN
      one=1
      IF(SUM(vector_batch_A_size_group)>virtual) one=-1
        i=-1
        DO
          i=i+1
          vector_batch_A_size_group(i)=vector_batch_A_size_group(i)+one
          IF(SUM(vector_batch_A_size_group)==virtual) EXIT
          IF(i==q_best-1) i=-1
        END DO
    END IF

    !XXXXXXXXXXXXX inverse group distribution
    group_counter=0
    a_group_counter=0
    my_A_virtual_start=1
    DO j=0, q_best-1
      my_I_occupied_start=1
      i_group_counter=0
      DO i=0, p_best-1
        group_counter=group_counter+1
        IF(color_sub==group_counter-1) EXIT
        my_I_occupied_start=my_I_occupied_start+vector_batch_I_size_group(i)
        i_group_counter=i_group_counter+1
      END DO
      my_q_position=j
      IF(color_sub==group_counter-1) EXIT
      my_A_virtual_start=my_A_virtual_start+vector_batch_A_size_group(j)
      a_group_counter=a_group_counter+1
    END DO
    !XXXXXXXXXXXXX inverse group distribution

    my_I_occupied_end=my_I_occupied_start+vector_batch_I_size_group(i_group_counter)-1
    my_I_batch_size=vector_batch_I_size_group(i_group_counter)
    my_A_virtual_end=my_A_virtual_start+vector_batch_A_size_group(a_group_counter)-1
    my_A_batch_size=vector_batch_A_size_group(a_group_counter)

    DEALLOCATE(vector_batch_I_size_group)
    DEALLOCATE(vector_batch_A_size_group)

    ! replicate on a local array on proc 0 the occupied and virtual wavevectior
    ! needed for the calculation of the WF's by calculate_wavefunction
    ! (external vector)
    CALL grep_occ_virt_wavefunc(qs_env,para_env_sub,nmo,homo,virtual,&
                                my_I_occupied_start,my_I_occupied_end,my_I_batch_size,&
                                my_A_virtual_start,my_A_virtual_end,my_A_batch_size,&
                                mo_coeff_o,mo_coeff_v,my_Cocc,my_Cvirt,unit_nr,error_sub)


    ! divide the b states in the sub_group in such a way to create
    ! b_start and b_end for each proc inside the sub_group
    max_b_size=(virtual+para_env_sub%num_pe-1)/para_env_sub%num_pe
    ALLOCATE(vector_B_sizes(0:para_env_sub%num_pe-1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
    vector_B_sizes=max_b_size
    IF(SUM(vector_B_sizes)/=virtual) THEN
      one=1
      IF(SUM(vector_B_sizes)>virtual) one=-1
      i=-1
      DO
        i=i+1
        vector_B_sizes(i)=vector_B_sizes(i)+one
        IF(SUM(vector_B_sizes)==virtual) EXIT
        IF(i==para_env_sub%num_pe-1) i=-1
      END DO
    END IF
    ! now give to each proc its b_start and b_end
    b_group_counter=0
    my_B_virtual_start=1
    DO j=0, para_env_sub%num_pe-1
      b_group_counter=b_group_counter+1
      IF(b_group_counter-1==para_env_sub%mepos) EXIT
      my_B_virtual_start=my_B_virtual_start+vector_B_sizes(j)
    END DO
    my_B_virtual_end=my_B_virtual_start+vector_B_sizes(para_env_sub%mepos)-1
    my_B_size=vector_B_sizes(para_env_sub%mepos)

    DEALLOCATE(vector_B_sizes)

    ! create an array containing a different "color" for each pair of
    ! A_start and B_start, communication will take place only among
    ! those proc that have the same A_start and B_start
    ALLOCATE(color_array(0:para_env_sub%num_pe-1,0:q_best-1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
    color_array=0
    color_counter=0
    DO j=0, q_best-1
      DO i=0, para_env_sub%num_pe-1
        color_counter=color_counter+1
        color_array(i,j)=color_counter
      END DO
    END DO
    sub_sub_color=color_array(para_env_sub%mepos,my_q_position)

    DEALLOCATE(color_array)

    ! now create a group that contains all the proc that have the same 2 virtual starting points
    ! in this way it is possible to sum the common integrals needed for the full MP2 energy
    ! in mp_comm_split_direct the color is give by my_a_virtual_start and my_b_virtual_start
    CALL mp_comm_split_direct(para_env%group,comm_exchange,sub_sub_color)
    NULLIFY(para_env_exchange)
    CALL cp_para_env_create(para_env_exchange,comm_exchange,error=error_sub)

    ! crate the proc maps
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

    ! create an array containing the information for communication
    ALLOCATE(exchange_group_sizes(0:para_env_exchange%num_pe-1,3))
    exchange_group_sizes=0
    exchange_group_sizes(para_env_exchange%mepos,1)=my_I_occupied_start
    exchange_group_sizes(para_env_exchange%mepos,2)=my_I_occupied_end
    exchange_group_sizes(para_env_exchange%mepos,3)=my_I_batch_size
    CALL mp_sum(exchange_group_sizes,para_env_exchange%group)
    mepos_in_EX_group=para_env_exchange%mepos
    size_of_exchange_group=para_env_exchange%num_pe

    NULLIFY(psi_a%pw)
    CALL pw_pool_create_pw(auxbas_pw_pool,psi_a%pw,&
                            use_data=REALDATA3D,&
                            in_space=REALSPACE,error=error_sub)


    ALLOCATE(psi_i(my_I_occupied_start:my_I_occupied_end),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
    DO i=my_I_occupied_start, my_I_occupied_end
      NULLIFY(psi_i(i)%pw)
      CALL pw_pool_create_pw(auxbas_pw_pool,psi_i(i)%pw,&
                            use_data=REALDATA3D,&
                            in_space=REALSPACE,error=error_sub)
      CALL calculate_wavefunction(mo_coeff,i,psi_i(i),rho_g, atomic_kind_set,&
                        qs_kind_set,cell,dft_control,particle_set, &
                        pw_env_sub,external_vector=my_Cocc(:,i-my_I_occupied_start+1),error=error_sub)
    END DO

    Emp2=0.0_dp
    Emp2_Cou=0.0_dp
    Emp2_EX=0.0_dp
    IF(do_alpha_beta) Emp2_AB=0.0_dp
    IF(calc_ex) THEN
      ALLOCATE(BIb_C(my_B_size,homo,my_I_batch_size),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
    END IF

    CALL timeset(routineN//"_loop",handle2)
    DO a=homo+my_A_virtual_start, homo+my_A_virtual_end

      IF(calc_ex) BIb_C=0.0_dp

      ! psi_a
      CALL calculate_wavefunction(mo_coeff,a,psi_a,rho_g, atomic_kind_set,&
                      qs_kind_set,cell,dft_control,particle_set, &
                      pw_env_sub,external_vector=my_Cvirt(:,a-(homo+my_A_virtual_start)+1),error=error_sub)
      i_counter=0
      DO i=my_I_occupied_start, my_I_occupied_end
        i_counter=i_counter+1

        ! potential
        CALL timeset(routineN//"_pot",handle3)
        rho_r%pw%cr3d = psi_i(i)%pw%cr3d * psi_a%pw%cr3d
        CALL pw_transfer(rho_r%pw, rho_g%pw, error=error_sub)
        CALL pw_poisson_solve(poisson_env,rho_g%pw, pair_energy, pot_g%pw,error=error_sub)
        CALL pw_transfer(pot_g%pw, rho_r%pw, error=error_sub)
        CALL pw_scale(rho_r%pw,rho_r%pw%pw_grid%dvol, error=error_sub)
        CALL timestop(handle3)


        ! and finally (ia|munu)
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

        CALL timeset(routineN//"_E_Cou",handle3)
        CALL cp_dbcsr_iterator_start(iter, matrix_ia_jb)
        DO WHILE (cp_dbcsr_iterator_blocks_left (iter))
           CALL cp_dbcsr_iterator_next_block(iter, row, col, data_block,blk,&
                row_size=row_size, col_size=col_size, &
                row_offset=row_offset, col_offset=col_offset)
           DO b=1,col_size
           DO j=1,row_size
              ! Compute the coulomb MP2 energy
              Emp2_Cou=Emp2_Cou-2.0_dp*data_block(j,b)**2/&
                       (Eigenval(a)+Eigenval(homo+col_offset+b-1)-Eigenval(i)-Eigenval(row_offset+j-1))
           ENDDO
           ENDDO
        ENDDO
        CALL cp_dbcsr_iterator_stop(iter)
        IF(do_alpha_beta) THEN
          ! Compute the coulomb only= SO = MP2 alpha-beta  MP2 energy component
          CALL cp_dbcsr_iterator_start(iter, matrix_ia_jb_beta)
          DO WHILE (cp_dbcsr_iterator_blocks_left (iter))
             CALL cp_dbcsr_iterator_next_block(iter, row, col, data_block,blk,&
                  row_size=row_size, col_size=col_size, &
                  row_offset=row_offset, col_offset=col_offset)
             DO b=1,col_size
             DO j=1,row_size
                ! Compute the coulomb MP2 energy alpha beta case
                Emp2_AB=Emp2_AB-data_block(j,b)**2/&
                         (Eigenval(a)+Eigenval_beta(homo_beta+col_offset+b-1)-Eigenval(i)-Eigenval_beta(row_offset+j-1))
             ENDDO
             ENDDO
          ENDDO
          CALL cp_dbcsr_iterator_stop(iter)         
        END IF        
        CALL timestop(handle3)

        ! now collect my local data from all the other members of the group
        ! b_start, b_end
        IF(calc_ex) THEN
          CALL timeset(routineN//"_E_Ex_1",handle3)
          CALL copy_dbcsr_to_fm(matrix_ia_jb, fm_BIb_jb, error=error_sub)
          CALL grep_my_integrals(para_env_sub,fm_BIb_jb,BIb_C(1:my_B_size,1:homo,i_counter),max_row_col_local,&
                                 homo,virtual,sub_proc_map,local_col_row_info,&
                                 my_B_virtual_end,my_B_virtual_start,my_B_size,&
                                 error_sub)
          CALL timestop(handle3)
        END IF

      END DO

      IF(calc_ex) THEN
        CALL timeset(routineN//"_E_Ex_2",handle3)
        ! calculate the contribution to MP2 energy for my local data
        DO i=1,my_I_batch_size
          DO j=my_I_occupied_start,my_I_occupied_end
            DO b=1,my_B_size
              b_global=b-1+my_B_virtual_start
              Emp2_EX=Emp2_EX+BIb_C(b,j,i)*BIb_C(b,i+my_I_occupied_start-1,j-my_I_occupied_start+1)&
                    /(Eigenval(a)+Eigenval(homo+b_global)-Eigenval(i+my_I_occupied_start-1)-Eigenval(j))
            END DO
          END DO
        END DO

        ! start communicating and collecting exchange contributions from
        ! other processes in my exchange group
        DO index_proc_shift=1, size_of_exchange_group-1
          proc_send=proc_map(mepos_in_EX_group+index_proc_shift)
          proc_receive=proc_map(mepos_in_EX_group-index_proc_shift)

          EX_start=exchange_group_sizes(proc_receive,1)
          EX_end=exchange_group_sizes(proc_receive,2)
          size_EX=exchange_group_sizes(proc_receive,3)

          ALLOCATE(BIb_EX(my_B_size,my_I_batch_size,size_EX),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
          BIb_EX=0.0_dp

          EX_start_send=exchange_group_sizes(proc_send,1)
          EX_end_send=exchange_group_sizes(proc_send,2)
          size_EX_send=exchange_group_sizes(proc_send,3)

          ALLOCATE(BIb_send(my_B_size,size_EX_send,my_I_batch_size),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
          BIb_send(1:my_B_size,1:size_EX_send,1:my_I_batch_size)=BIb_C(1:my_B_size,EX_start_send:EX_end_send,1:my_I_batch_size)

          ! send and receive the exchange array
          CALL mp_sendrecv(BIb_send,proc_send,BIb_EX,proc_receive,para_env_exchange%group)

          DO i=1,my_I_batch_size
            DO j=1,size_EX
              DO b=1, my_B_size
                b_global=b-1+my_B_virtual_start
                Emp2_EX=Emp2_EX+BIb_C(b,j+EX_start-1,i)*BIb_EX(b,i,j)&
                      /(Eigenval(a)+Eigenval(homo+b_global)-Eigenval(i+my_I_occupied_start-1)-Eigenval(j+EX_start-1))
              END DO
            END DO
          END DO

          DEALLOCATE(BIb_EX)
          DEALLOCATE(BIb_send)

        END DO
        CALL timestop(handle3)
      END IF

    ENDDO
    CALL timestop(handle2)

    CALL mp_sum(Emp2_Cou,para_env%group)
    CALL mp_sum(Emp2_EX,para_env%group)
    Emp2=Emp2_Cou+Emp2_EX
    IF(do_alpha_beta) CALL mp_sum(Emp2_AB,para_env%group)

    DEALLOCATE(my_Cocc)
    DEALLOCATE(my_Cvirt)

    IF(calc_ex) THEN
      CALL cp_fm_release(fm_BIb_jb, error=error_sub)
      DEALLOCATE(local_col_row_info)
      DEALLOCATE(BIb_C)
    END IF
    DEALLOCATE(proc_map)
    DEALLOCATE(sub_proc_map)
    DEALLOCATE(exchange_group_sizes)

    CALL cp_para_env_release(para_env_exchange,error=error_sub)

    CALL cp_dbcsr_release(matrix_ia_jnu,error=error_sub)
    CALL cp_dbcsr_release(matrix_ia_jb,error=error_sub)
    IF(do_alpha_beta) THEN
      CALL cp_dbcsr_release(matrix_ia_jnu_beta,error=error_sub)
      CALL cp_dbcsr_release(matrix_ia_jb_beta,error=error_sub)
    END IF

    DO i=my_I_occupied_start, my_I_occupied_end
      CALL pw_release(psi_i(i)%pw,error=error_sub)
    END DO
    DEALLOCATE(psi_i)

    CALL pw_pool_give_back_pw(auxbas_pw_pool,psi_a%pw,error=error_sub)

    CALL timestop(handle)

  END SUBROUTINE mp2_gpw_compute

! *****************************************************************************
!> \brief ...
!> \param wfn_size ...
!> \param p ...
!> \param q ...
!> \param num_w ...
!> \param nmo ...
!> \param virtual ...
!> \param homo ...
!> \param calc_ex ...
!> \param mem_try ...
! *****************************************************************************
  SUBROUTINE estimate_memory_usage(wfn_size,p,q,num_w,nmo, virtual, homo,calc_ex,mem_try)
    REAL(KIND=dp)                            :: wfn_size
    INTEGER                                  :: p, q, num_w, nmo, virtual, &
                                                homo
    LOGICAL                                  :: calc_ex
    REAL(KIND=dp)                            :: mem_try

    mem_try=0.0_dp
    ! integrals
    mem_try=mem_try+virtual*REAL(homo,KIND=dp)**2/(p*num_w)
    ! array for the coefficient matrix and wave vectors
    mem_try=mem_try+REAL(homo,KIND=dp)*nmo/p+&
                    REAL(virtual,KIND=dp)*nmo/q+&
                    2.0_dp*MAX(REAL(homo,KIND=dp)*nmo/p,REAL(virtual,KIND=dp)*nmo/q)
    ! temporary array for MO integrals and MO integrals to be exchanged 
    IF(calc_ex) THEN
      mem_try=mem_try+2.0_dp*MAX(virtual*REAL(homo,KIND=dp)*MIN(1,num_w-1)/num_w,&
                                 virtual*REAL(homo,KIND=dp)**2/(p*p*num_w))
    ELSE
      mem_try=mem_try+2.0_dp*virtual*REAL(homo,KIND=dp)
    END IF
    ! wfn
    mem_try=mem_try+((homo+p-1)/p)*wfn_size
    ! Mb
    mem_try=mem_try*8.0D+00/1024.0D+00**2

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
!> \param qs_env ...
!> \param mp2_env ...
!> \param para_env ...
!> \param para_env_sub ...
!> \param mo_coeff ...
!> \param dimen ...
!> \param homo ...
!> \param mat_munu ...
!> \param mo_coeff_o ...
!> \param mo_coeff_v ...
!> \param blacs_env_sub ...
!> \param unit_nr ...
!> \param error ...
!> \param error_sub ...
! *****************************************************************************
  SUBROUTINE replicate_mat_to_subgroup(qs_env,mp2_env,para_env,para_env_sub,mo_coeff,dimen,homo,mat_munu,&
                                       mo_coeff_o,mo_coeff_v,blacs_env_sub,unit_nr,error,error_sub)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(mp2_type), POINTER                  :: mp2_env
    TYPE(cp_para_env_type), POINTER          :: para_env, para_env_sub
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    INTEGER                                  :: dimen, homo
    TYPE(cp_dbcsr_p_type)                    :: mat_munu
    TYPE(cp_dbcsr_type), POINTER             :: mo_coeff_o, mo_coeff_v
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env_sub
    INTEGER                                  :: unit_nr
    TYPE(cp_error_type), INTENT(inout)       :: error, error_sub

    CHARACTER(LEN=*), PARAMETER :: routineN = 'replicate_mat_to_subgroup', &
      routineP = moduleN//':'//routineN

    INTEGER :: blk, col, col_offset, col_size, handle, i, i_global, iiB, &
      iproc, itmp(2), j, j_global, jjB, max_row_col_local, my_mu_end, &
      my_mu_size, my_mu_start, ncol_local, ncol_rec, nrow_local, nrow_rec, &
      proc_receive, proc_receive_static, proc_send, proc_send_static, &
      proc_shift, rec_mu_end, rec_mu_size, rec_mu_start, row, row_offset, &
      row_size, stat, virtual
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: ends_array, proc_map, &
                                                sizes_array, starts_array
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: local_col_row_info, &
                                                rec_col_row_info
    INTEGER, DIMENSION(:), POINTER           :: col_indices, col_indices_rec, &
                                                row_indices, row_indices_rec
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: C, Cocc, Cvirt, rec_C
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: data_block, local_C, &
                                                local_C_internal
    TYPE(array_i1d_obj)                      :: col_blk_size, col_dist
    TYPE(cp_dbcsr_iterator)                  :: iter
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_coeff
    TYPE(dbcsr_distribution_obj)             :: dist

    CALL timeset(routineN,handle)
    failure=.FALSE.

    ALLOCATE(sizes_array(0:para_env_sub%num_pe-1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(starts_array(0:para_env_sub%num_pe-1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    starts_array=0
    ALLOCATE(ends_array(0:para_env_sub%num_pe-1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ends_array=0

    DO iproc=0,para_env_sub%num_pe-1
       itmp=get_limit(dimen,para_env_sub%num_pe,iproc) 
       starts_array(iproc)=itmp(1)
       ends_array(iproc)=itmp(2)
       sizes_array(iproc)=itmp(2)-itmp(1)+1
    ENDDO
 
    my_mu_size=sizes_array(para_env_sub%mepos)
    my_mu_start=starts_array(para_env_sub%mepos)
    my_mu_end=ends_array(para_env_sub%mepos)

    ! local storage for the C matrix
    ALLOCATE(C(my_mu_size,dimen),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    C=0.0_dp

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

    CALL cp_fm_get_info(matrix=mo_coeff,&
                        matrix_struct=fm_struct_coeff,&
                        nrow_local=nrow_local,&
                        ncol_local=ncol_local,&
                        row_indices=row_indices,&
                        col_indices=col_indices,&
                        local_data=local_C_internal,&
                        error=error)

    ALLOCATE(local_C(nrow_local,ncol_local),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    local_C=local_C_internal(1:nrow_local,1:ncol_local)
    NULLIFY(local_C_internal) 

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

    ! accumulate data on C buffer starting from myself
    DO iiB=1, nrow_local
      i_global=row_indices(iiB)
      IF(i_global>=my_mu_start.AND.i_global<=my_mu_end) THEN
        DO jjB=1, ncol_local
          j_global=col_indices(jjB)
          C(i_global-my_mu_start+1,j_global)=local_C(iiB,jjB)
        END DO
      END IF
    END DO

    ! start ring communication for collecting the data from the other
    proc_send_static=proc_map(para_env%mepos+1)
    proc_receive_static=proc_map(para_env%mepos-1)
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

      ALLOCATE(rec_C(nrow_rec,ncol_rec),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      rec_C=0.0_dp

      ! then send and receive the real data
      CALL  mp_sendrecv(local_C,proc_send_static,rec_C,proc_receive_static,para_env%group)

      ! accumulate the received data on C buffer 
      DO iiB=1, nrow_rec
        i_global=row_indices_rec(iiB)
        IF(i_global>=my_mu_start.AND.i_global<=my_mu_end) THEN
          DO jjB=1, ncol_rec
            j_global=col_indices_rec(jjB)
            C(i_global-my_mu_start+1,j_global)=rec_C(iiB,jjB)
          END DO
        END IF
      END DO

      local_col_row_info(:,:)=rec_col_row_info
      DEALLOCATE(local_C)
      ALLOCATE(local_C(nrow_rec,ncol_rec),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      local_C=rec_C

      DEALLOCATE(col_indices_rec)
      DEALLOCATE(row_indices_rec)
      DEALLOCATE(rec_C)
    END DO

    DEALLOCATE(local_col_row_info)
    DEALLOCATE(rec_col_row_info)
    DEALLOCATE(proc_map)

    ! proc_map, for the sub_group
    ALLOCATE(proc_map(-para_env_sub%num_pe:2*para_env_sub%num_pe-1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DO iiB=0,para_env_sub%num_pe-1
      proc_map(iiB)=iiB
      proc_map(-iiB-1)=para_env_sub%num_pe-iiB-1
      proc_map(para_env_sub%num_pe+iiB)=iiB
    END DO

    ! split the C matrix into occupied and virtual
    ALLOCATE(Cocc(my_mu_size,homo),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    Cocc(1:my_mu_size,1:homo)=C(1:my_mu_size,1:homo)
 
    virtual=dimen-homo
    ALLOCATE(Cvirt(my_mu_size,virtual),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    Cvirt(1:my_mu_size,1:virtual)=C(1:my_mu_size,homo+1:dimen)

    DEALLOCATE(C)

    ! the occupied
    CALL cp_create_bl_distribution (col_dist, col_blk_size, homo, &
         dbcsr_mp_npcols(dbcsr_distribution_mp(cp_dbcsr_distribution(mat_munu%matrix))))
    CALL dbcsr_distribution_new (dist, dbcsr_distribution_mp (cp_dbcsr_distribution(mat_munu%matrix)),&
         dbcsr_distribution_row_dist(cp_dbcsr_distribution(mat_munu%matrix)), col_dist)
    NULLIFY(mo_coeff_o)
    CALL cp_dbcsr_init_p(mo_coeff_o,error=error_sub)
    CALL cp_dbcsr_create(mo_coeff_o, "mo_coeff_o", dist, dbcsr_type_no_symmetry,&
         cp_dbcsr_row_block_sizes(mat_munu%matrix), col_blk_size, 0, 0, dbcsr_type_real_default,&
         error=error_sub)
    CALL cp_dbcsr_distribution_release (dist)
    CALL array_release (col_blk_size)
    CALL array_release (col_dist)
    CALL cp_dbcsr_reserve_all_blocks(mo_coeff_o,error)

    ! accumulate data on mo_coeff_o starting from myself
    CALL cp_dbcsr_iterator_start(iter, mo_coeff_o)
    DO WHILE (cp_dbcsr_iterator_blocks_left (iter))
       CALL cp_dbcsr_iterator_next_block(iter, row, col, data_block,blk,&
            row_size=row_size, col_size=col_size, &
            row_offset=row_offset, col_offset=col_offset)
       DO i=1,row_size
         i_global=row_offset+i-1
         IF(i_global>=my_mu_start.AND.i_global<=my_mu_end) THEN
           DO j=1,col_size
             j_global=col_offset+j-1
             data_block(i,j)=Cocc(i_global-my_mu_start+1,col_offset+j-1)
           ENDDO
         END IF
       ENDDO
    ENDDO
    CALL cp_dbcsr_iterator_stop(iter)

    ! start ring communication in the subgroup for collecting the data from the other
    ! proc (occupied)
    DO proc_shift=1, para_env_sub%num_pe-1
      proc_send=proc_map(para_env_sub%mepos+proc_shift)
      proc_receive=proc_map(para_env_sub%mepos-proc_shift)

      rec_mu_start=starts_array(proc_receive)
      rec_mu_end=ends_array(proc_receive)
      rec_mu_size=sizes_array(proc_receive)

      ALLOCATE(rec_C(rec_mu_size,homo),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      rec_C=0.0_dp

      ! then send and receive the real data
      CALL  mp_sendrecv(Cocc,proc_send,rec_C,proc_receive,para_env_sub%group) 

      ! accumulate data on mo_coeff_o the data received from proc_rec
      CALL cp_dbcsr_iterator_start(iter, mo_coeff_o)
      DO WHILE (cp_dbcsr_iterator_blocks_left (iter))
         CALL cp_dbcsr_iterator_next_block(iter, row, col, data_block,blk,&
              row_size=row_size, col_size=col_size, &
              row_offset=row_offset, col_offset=col_offset)
         DO i=1,row_size
           i_global=row_offset+i-1
           IF(i_global>=rec_mu_start.AND.i_global<=rec_mu_end) THEN
             DO j=1,col_size
               j_global=col_offset+j-1
               data_block(i,j)=rec_C(i_global-rec_mu_start+1,col_offset+j-1)
             ENDDO
           END IF
         ENDDO
      ENDDO
      CALL cp_dbcsr_iterator_stop(iter)

      DEALLOCATE(rec_C)

    END DO
    CALL cp_dbcsr_filter(mo_coeff_o,mp2_env%mp2_gpw%eps_filter,error=error_sub)

    DEALLOCATE(Cocc)

    ! the virtual
    CALL cp_create_bl_distribution (col_dist, col_blk_size, virtual, &
         dbcsr_mp_npcols(dbcsr_distribution_mp(cp_dbcsr_distribution(mat_munu%matrix))))
    CALL dbcsr_distribution_new (dist, dbcsr_distribution_mp (cp_dbcsr_distribution(mat_munu%matrix)),&
         dbcsr_distribution_row_dist(cp_dbcsr_distribution(mat_munu%matrix)), col_dist)
    NULLIFY(mo_coeff_v)
    CALL cp_dbcsr_init_p(mo_coeff_v,error=error_sub)
    CALL cp_dbcsr_create(mo_coeff_v, "mo_coeff_v", dist, dbcsr_type_no_symmetry,&
         cp_dbcsr_row_block_sizes(mat_munu%matrix), col_blk_size, 0, 0, dbcsr_type_real_default,&
         error=error_sub)
    CALL cp_dbcsr_distribution_release (dist)
    CALL array_release (col_blk_size)
    CALL array_release (col_dist)
    CALL cp_dbcsr_reserve_all_blocks(mo_coeff_v,error)

    ! accumulate data on mo_coeff_v starting from myself
    CALL cp_dbcsr_iterator_start(iter, mo_coeff_v)
    DO WHILE (cp_dbcsr_iterator_blocks_left (iter))
       CALL cp_dbcsr_iterator_next_block(iter, row, col, data_block,blk,&
            row_size=row_size, col_size=col_size, &
            row_offset=row_offset, col_offset=col_offset)
       DO i=1,row_size
         i_global=row_offset+i-1
         IF(i_global>=my_mu_start.AND.i_global<=my_mu_end) THEN
           DO j=1,col_size
             j_global=col_offset+j-1
             data_block(i,j)=Cvirt(i_global-my_mu_start+1,col_offset+j-1)
           ENDDO
         END IF
       ENDDO
    ENDDO
    CALL cp_dbcsr_iterator_stop(iter)

    ! start ring communication in the subgroup for collecting the data from the other
    ! proc (virtual)
    DO proc_shift=1, para_env_sub%num_pe-1
      proc_send=proc_map(para_env_sub%mepos+proc_shift)
      proc_receive=proc_map(para_env_sub%mepos-proc_shift)

      rec_mu_start=starts_array(proc_receive)
      rec_mu_end=ends_array(proc_receive)
      rec_mu_size=sizes_array(proc_receive)

      ALLOCATE(rec_C(rec_mu_size,virtual),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      rec_C=0.0_dp

      ! then send and receive the real data
      CALL  mp_sendrecv(Cvirt,proc_send,rec_C,proc_receive,para_env_sub%group)

      ! accumulate data on mo_coeff_v the data received from proc_rec
      CALL cp_dbcsr_iterator_start(iter, mo_coeff_v)
      DO WHILE (cp_dbcsr_iterator_blocks_left (iter))
         CALL cp_dbcsr_iterator_next_block(iter, row, col, data_block,blk,&
              row_size=row_size, col_size=col_size, &
              row_offset=row_offset, col_offset=col_offset)
         DO i=1,row_size
           i_global=row_offset+i-1
           IF(i_global>=rec_mu_start.AND.i_global<=rec_mu_end) THEN
             DO j=1,col_size
               j_global=col_offset+j-1
               data_block(i,j)=rec_C(i_global-rec_mu_start+1,col_offset+j-1)
             ENDDO
           END IF
         ENDDO
      ENDDO
      CALL cp_dbcsr_iterator_stop(iter)

      DEALLOCATE(rec_C)

    END DO
    CALL cp_dbcsr_filter(mo_coeff_v,mp2_env%mp2_gpw%eps_filter,error=error_sub)

    DEALLOCATE(Cvirt)

    DEALLOCATE(proc_map)
    DEALLOCATE(sizes_array)
    DEALLOCATE(starts_array)
    DEALLOCATE(ends_array)
    DEALLOCATE(local_C)

    CALL timestop(handle)

  END SUBROUTINE

! *****************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param para_env_sub ...
!> \param dimen ...
!> \param homo ...
!> \param virtual ...
!> \param my_I_occupied_start ...
!> \param my_I_occupied_end ...
!> \param my_I_batch_size ...
!> \param my_A_virtual_start ...
!> \param my_A_virtual_end ...
!> \param my_A_batch_size ...
!> \param mo_coeff_o ...
!> \param mo_coeff_v ...
!> \param my_Cocc ...
!> \param my_Cvirt ...
!> \param unit_nr ...
!> \param error_sub ...
! *****************************************************************************
  SUBROUTINE grep_occ_virt_wavefunc(qs_env,para_env_sub,dimen,homo,virtual,&
                                    my_I_occupied_start,my_I_occupied_end,my_I_batch_size,&
                                    my_A_virtual_start,my_A_virtual_end,my_A_batch_size,&
                                    mo_coeff_o,mo_coeff_v,my_Cocc,my_Cvirt,unit_nr,error_sub)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_para_env_type), POINTER          :: para_env_sub
    INTEGER :: dimen, homo, virtual, my_I_occupied_start, my_I_occupied_end, &
      my_I_batch_size, my_A_virtual_start, my_A_virtual_end, my_A_batch_size
    TYPE(cp_dbcsr_type), POINTER             :: mo_coeff_o, mo_coeff_v
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: my_Cocc, my_Cvirt
    INTEGER                                  :: unit_nr
    TYPE(cp_error_type), INTENT(inout)       :: error_sub

    CHARACTER(LEN=*), PARAMETER :: routineN = 'grep_occ_virt_wavefunc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, col, col_offset, &
                                                col_size, handle, i, &
                                                i_global, j, j_global, row, &
                                                row_offset, row_size, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: data_block
    TYPE(cp_dbcsr_iterator)                  :: iter

    CALL timeset(routineN,handle)
    failure=.FALSE.

    ALLOCATE(my_Cocc(dimen,my_I_batch_size),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
    my_Cocc=0.0_dp

    ALLOCATE(my_Cvirt(dimen,my_A_batch_size),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
    my_Cvirt=0.0_dp

    ! accumulate data from mo_coeff_o into Cocc
    CALL cp_dbcsr_iterator_start(iter, mo_coeff_o)
    DO WHILE (cp_dbcsr_iterator_blocks_left (iter))
       CALL cp_dbcsr_iterator_next_block(iter, row, col, data_block,blk,&
            row_size=row_size, col_size=col_size, &
            row_offset=row_offset, col_offset=col_offset)
      DO j=1,col_size
        j_global=col_offset+j-1
        IF(j_global>=my_I_occupied_start.AND.j_global<=my_I_occupied_end) THEN
          DO i=1,row_size
            i_global=row_offset+i-1
            my_Cocc(i_global,j_global-my_I_occupied_start+1)=data_block(i,j)
          END DO
        END IF
      END DO
    ENDDO
    CALL cp_dbcsr_iterator_stop(iter)   

    CALL mp_sum(my_Cocc,para_env_sub%group)

    ! accumulate data from mo_coeff_o into Cocc
    CALL cp_dbcsr_iterator_start(iter, mo_coeff_v)
    DO WHILE (cp_dbcsr_iterator_blocks_left (iter))
       CALL cp_dbcsr_iterator_next_block(iter, row, col, data_block,blk,&
            row_size=row_size, col_size=col_size, &
            row_offset=row_offset, col_offset=col_offset)
      DO j=1,col_size
        j_global=col_offset+j-1
        IF(j_global>=my_A_virtual_start.AND.j_global<=my_A_virtual_end) THEN
          DO i=1,row_size
            i_global=row_offset+i-1
            my_Cvirt(i_global,j_global-my_A_virtual_start+1)=data_block(i,j)
          END DO
        END IF
      END DO
    ENDDO
    CALL cp_dbcsr_iterator_stop(iter)

    CALL mp_sum(my_Cvirt,para_env_sub%group) 

    CALL timestop(handle)

  END SUBROUTINE grep_occ_virt_wavefunc

END MODULE mp2_gpw
