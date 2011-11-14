!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2011  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Rountines to calculate MP2 energy using pw 
!> \par History
!>      10.2011 created [Joost VandeVondele and Mauro Del Ben]
! *****************************************************************************
MODULE mp2_gpw
  USE array_types,                     ONLY: array_i1d_obj,&
                                             array_new,&
                                             array_nullify,&
                                             array_release
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE cell_types,                      ONLY: cell_type
  USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                             cp_blacs_env_release
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_interface,              ONLY: &
       cp_create_bl_distribution, cp_dbcsr_create, cp_dbcsr_distribution, &
       cp_dbcsr_distribution_release, cp_dbcsr_filter, cp_dbcsr_init, &
       cp_dbcsr_init_p, cp_dbcsr_iterator_blocks_left, &
       cp_dbcsr_iterator_next_block, cp_dbcsr_iterator_start, &
       cp_dbcsr_iterator_stop, cp_dbcsr_multiply, cp_dbcsr_release, &
       cp_dbcsr_row_block_sizes, cp_dbcsr_set
  USE cp_dbcsr_operations,             ONLY: copy_fm_to_dbcsr,&
                                             cp_dbcsr_alloc_block_from_nbl,&
                                             cp_dbcsr_dist2d_to_dist
  USE cp_dbcsr_types,                  ONLY: cp_dbcsr_iterator,&
                                             cp_dbcsr_p_type,&
                                             cp_dbcsr_type
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_submatrix,&
                                             cp_fm_release,&
                                             cp_fm_set_all,&
                                             cp_fm_type
  USE cp_para_env,                     ONLY: cp_para_env_create,&
                                             cp_para_env_release
  USE cp_para_types,                   ONLY: cp_blacs_env_type,&
                                             cp_para_env_type
  USE cp_subsys_types,                 ONLY: cp_subsys_type
  USE dbcsr_methods,                   ONLY: dbcsr_distribution_mp,&
                                             dbcsr_distribution_new,&
                                             dbcsr_distribution_row_dist,&
                                             dbcsr_mp_npcols,&
                                             dbcsr_mp_nprows
  USE dbcsr_types,                     ONLY: dbcsr_distribution_obj,&
                                             dbcsr_type_no_symmetry,&
                                             dbcsr_type_real_default,&
                                             dbcsr_type_symmetric
  USE distribution_1d_types,           ONLY: distribution_1d_release,&
                                             distribution_1d_type
  USE distribution_2d_types,           ONLY: distribution_2d_release,&
                                             distribution_2d_type
  USE distribution_methods,            ONLY: distribute_molecules_1d,&
                                             distribute_molecules_2d
  USE f77_blas
  USE input_constants,                 ONLY: BLACS_GRID_SQUARE,&
                                             use_orb_basis_set
  USE input_section_types,             ONLY: section_vals_val_get
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_comm_split_direct,&
                                             mp_sendrecv,&
                                             mp_sum
  USE molecule_kind_types,             ONLY: molecule_kind_type
  USE molecule_types_new,              ONLY: molecule_type
  USE mp2_types,                       ONLY: mp2_type
  USE particle_types,                  ONLY: get_particle_set,&
                                             particle_type
  USE pw_env_methods,                  ONLY: pw_env_create
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
                                             pw_p_type
  USE qs_collocate_density,            ONLY: calculate_wavefunction
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_integrate_potential,          ONLY: integrate_v_rspace
  USE qs_interactions,                 ONLY: init_interaction_radii
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type,&
                                             mo_set_type
  USE qs_neighbor_list_types,          ONLY: deallocate_neighbor_list_set,&
                                             neighbor_list_set_p_type
  USE qs_neighbor_lists,               ONLY: atom2d_build,&
                                             atom2d_cleanup,&
                                             build_neighbor_lists,&
                                             local_atoms_type,&
                                             pair_radius_setup
  USE task_list_methods,               ONLY: generate_qs_task_list
  USE task_list_types,                 ONLY: allocate_task_list,&
                                             deallocate_task_list,&
                                             task_list_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mp2_gpw'

  PUBLIC :: mp2_gpw_main


  CONTAINS

! *****************************************************************************
!> \brief with a big bang to mp2 
!> \author Mauro Del Ben and Joost VandeVondele
! *****************************************************************************
  SUBROUTINE mp2_gpw_main(qs_env,mp2_env,para_env,unit_nr,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(mp2_type), POINTER                  :: mp2_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER                                  :: unit_nr
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'mp2_gpw_main', &
      routineP = moduleN//':'//routineN

    INTEGER :: blacs_grid_layout, blk, col, col_offset, col_size, color_sub, &
      comm_sub, dimen, handle, homo, i, j, local_unit_nr, nao, natom, &
      nelectron, nkind, nmo, row, row_offset, row_size, stat
    INTEGER, DIMENSION(:), POINTER           :: rbs
    LOGICAL                                  :: blacs_repeatable, failure, &
                                                skip_load_balance_distributed
    LOGICAL, ALLOCATABLE, DIMENSION(:)       :: orb_present
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: orb_radius
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: pair_radius
    REAL(KIND=dp)                            :: eps_gvg_rspace_old, &
                                                eps_pgf_orb_old, &
                                                eps_rho_rspace_old, subcells
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: Auto
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: C
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mo_eigenvalues
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: data_block
    TYPE(array_i1d_obj)                      :: col_blk_size, col_dist, &
                                                row_blk_sizes
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env_sub
    TYPE(cp_dbcsr_iterator)                  :: iter
    TYPE(cp_dbcsr_p_type)                    :: mat_munu
    TYPE(cp_dbcsr_type), POINTER             :: mo_coeff_o, mo_coeff_v
    TYPE(cp_error_type)                      :: error_sub
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: fm_mo_coeff_o, fm_mo_coeff_v, &
                                                mo_coeff
    TYPE(cp_logger_type), POINTER            :: logger, logger_sub
    TYPE(cp_para_env_type), POINTER          :: para_env_sub
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(dbcsr_distribution_obj)             :: dist
    TYPE(dbcsr_distribution_obj), POINTER    :: dbcsr_dist_sub
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(distribution_1d_type), POINTER      :: local_molecules_sub, &
                                                local_particles_sub
    TYPE(distribution_2d_type), POINTER      :: distribution_2d_sub
    TYPE(local_atoms_type), ALLOCATABLE, &
      DIMENSION(:)                           :: atom2d
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(mo_set_type), POINTER               :: mo_set
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
    TYPE(task_list_type), POINTER            :: task_list_sub

    CALL timeset(routineN,handle)
    failure=.FALSE.

    ! ... setup needed to be able to qs_integrate in a subgroup.
    logger => cp_error_get_logger(error)
    local_unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)

    ! a para_env
    color_sub=para_env%mepos/mp2_env%mp2_num_proc
    CALL mp_comm_split_direct(para_env%group,comm_sub,color_sub)
    NULLIFY(para_env_sub)
    CALL cp_para_env_create(para_env_sub,comm_sub,error=error)

    ! a logger
    NULLIFY(logger_sub)
    CALL cp_logger_create(logger_sub,para_env=para_env_sub,&
               default_global_unit_nr=local_unit_nr, close_global_unit_on_dealloc=.FALSE.)
    CALL cp_logger_set(logger_sub,local_filename="MP2_localLog")

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
    subsys=>qs_env%subsys
    particle_set => subsys%particles%els
    atomic_kind_set => subsys%atomic_kinds%els
    molecule_set => subsys%molecules_new%els
    molecule_kind_set => subsys%molecule_kinds_new%els

    ! re-init the radii to be able to generate pair lists with MP2-appropriate screening
    eps_pgf_orb_old=qs_env%dft_control%qs_control%eps_pgf_orb
    eps_rho_rspace_old=qs_env%dft_control%qs_control%eps_rho_rspace
    eps_gvg_rspace_old=qs_env%dft_control%qs_control%eps_gvg_rspace
    ! hack hack hack XXXXXXXXXXXXXXX
    qs_env%dft_control%qs_control%eps_pgf_orb   =mp2_env%mp2_gpw%eps_grid
    qs_env%dft_control%qs_control%eps_rho_rspace=mp2_env%mp2_gpw%eps_grid
    qs_env%dft_control%qs_control%eps_gvg_rspace=mp2_env%mp2_gpw%eps_grid
    CALL init_interaction_radii(qs_env%dft_control%qs_control,atomic_kind_set, error_sub)

    ! get a distribution_1d
    NULLIFY(local_particles_sub,local_molecules_sub)
    CALL distribute_molecules_1d(particle_kind_set=atomic_kind_set,&
                                 particle_set=particle_set,&
                                 local_particles=local_particles_sub,&
                                 molecule_kind_set=molecule_kind_set,&
                                 molecule_set=molecule_set,&
                                 local_molecules=local_molecules_sub,&
                                 force_env_section=qs_env%input,&
                                 error=error_sub)

    ! get a distribution_2d
    cell=>qs_env%cell
    NULLIFY(distribution_2d_sub)
    CALL distribute_molecules_2d(cell=cell,&
                                 particle_kind_set=atomic_kind_set,&
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
                      atomic_kind_set,molecule_set,molecule_only=.FALSE.,dftb=.FALSE.,&
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
    CALL get_particle_set(particle_set=particle_set,nsgf=rbs,error=error)
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

    ! a pw_env
    NULLIFY(pw_env_sub)
    CALL pw_env_create(pw_env_sub,qs_env,para_env_sub,error_sub)
    CALL pw_env_get(pw_env_sub, auxbas_pw_pool=auxbas_pw_pool,&
                    poisson_env=poisson_env,error=error_sub)

    ! now we need a task list, hard code skip_load_balance_distributed
    NULLIFY(task_list_sub)
    skip_load_balance_distributed=qs_env%dft_control%qs_control%skip_load_balance_distributed
    CALL allocate_task_list(task_list_sub,error_sub)
    CALL generate_qs_task_list(qs_env, task_list_sub, &
                 reorder_rs_grid_ranks=.TRUE., soft_valid=.FALSE., &
                 skip_load_balance_distributed=skip_load_balance_distributed,&
                 pw_env_external=pw_env_sub, sab_orb_external=sab_orb_sub, error=error_sub)

    ! get some of the grids ready
    NULLIFY(rho_r%pw,rho_g%pw,pot_g%pw)
    CALL pw_pool_create_pw(auxbas_pw_pool,rho_r%pw,&
                            use_data=REALDATA3D,&
                            in_space=REALSPACE,error=error_sub)
    !CALL pw_pool_create_pw(auxbas_pw_pool,psi_a%pw,&
    !                        use_data=REALDATA3D,&
    !                        in_space=REALSPACE,error=error_sub)
    !CALL pw_pool_create_pw(auxbas_pw_pool,psi_i%pw,&
    !                        use_data=REALDATA3D,&
    !                        in_space=REALSPACE,error=error_sub)
    CALL pw_pool_create_pw(auxbas_pw_pool,rho_g%pw,&
                            use_data=COMPLEXDATA1D,&
                            in_space=RECIPROCALSPACE,error=error_sub)
    CALL pw_pool_create_pw(auxbas_pw_pool,pot_g%pw,&
                            use_data=COMPLEXDATA1D,&
                            in_space=RECIPROCALSPACE,error=error_sub)

    ! and the array of mos
    CALL get_qs_env(qs_env=qs_env, mos=mos, dft_control=dft_control, error=error)
    CALL get_mo_set( mo_set=mos(1)%mo_set, nelectron=nelectron,&
                      eigenvalues=mo_eigenvalues,nmo=nmo,homo=homo,&
                      mo_coeff=mo_coeff,nao=dimen)
    ALLOCATE(C(dimen,dimen),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(Auto(dimen),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL cp_fm_get_submatrix(mo_coeff,C, 1, 1, dimen, dimen, .FALSE., error)
    Auto(:)=mo_eigenvalues(:)

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

    ! a very ugly trick to create a dbcsr matrix with all blocks ...
    CALL cp_fm_struct_create(fm_struct_tmp,para_env=para_env_sub,context=blacs_env_sub, &
                             nrow_global=dimen,ncol_global=homo,error=error_sub)
    CALL cp_fm_create(fm_mo_coeff_o, fm_struct_tmp,name="occ",error=error_sub)
    CALL cp_fm_struct_release(fm_struct_tmp,error=error_sub)
    CALL cp_fm_set_all(fm_mo_coeff_o,0.0_dp,error=error_sub)
    CALL copy_fm_to_dbcsr(fm_mo_coeff_o,mo_coeff_o, error=error_sub)
    CALL cp_fm_release(fm_mo_coeff_o, error=error_sub)

    CALL cp_dbcsr_iterator_start(iter, mo_coeff_o)
    DO WHILE (cp_dbcsr_iterator_blocks_left (iter))
       CALL cp_dbcsr_iterator_next_block(iter, row, col, data_block,blk,&
            row_size=row_size, col_size=col_size, &
            row_offset=row_offset, col_offset=col_offset)
       DO j=1,col_size
       DO i=1,row_size
          data_block(i,j)=C(row_offset+i-1,col_offset+j-1)
       ENDDO
       ENDDO
    ENDDO
    CALL cp_dbcsr_iterator_stop(iter)
    CALL cp_dbcsr_filter(mo_coeff_o,mp2_env%mp2_gpw%eps_filter,error=error_sub)


    ! the virtual
    CALL cp_create_bl_distribution (col_dist, col_blk_size, nmo-homo, &
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

    ! a very ugly trick to create a dbcsr matrix with all blocks ...
    CALL cp_fm_struct_create(fm_struct_tmp,para_env=para_env_sub,context=blacs_env_sub, &
                             nrow_global=dimen,ncol_global=nmo-homo,error=error_sub)
    CALL cp_fm_create(fm_mo_coeff_v, fm_struct_tmp,name="virt",error=error_sub)
    CALL cp_fm_struct_release(fm_struct_tmp,error=error_sub)
    CALL cp_fm_set_all(fm_mo_coeff_v,0.0_dp,error=error_sub)
    CALL copy_fm_to_dbcsr(fm_mo_coeff_v,mo_coeff_v, error=error_sub)
    CALL cp_fm_release(fm_mo_coeff_v, error=error_sub)

    CALL cp_dbcsr_iterator_start(iter, mo_coeff_v)
    DO WHILE (cp_dbcsr_iterator_blocks_left (iter))
       CALL cp_dbcsr_iterator_next_block(iter, row, col, data_block,blk,&
            row_size=row_size, col_size=col_size, &
            row_offset=row_offset, col_offset=col_offset)
       DO j=1,col_size
       DO i=1,row_size
          data_block(i,j)=C(row_offset+i-1,homo+col_offset+j-1)
       ENDDO
       ENDDO
    ENDDO
    CALL cp_dbcsr_iterator_stop(iter)
    CALL cp_dbcsr_filter(mo_coeff_v,mp2_env%mp2_gpw%eps_filter,error=error_sub)

    ! now we're kind of ready to go....
    CALL mp2_gpw_compute(qs_env,para_env,para_env_sub,color_sub,dft_control,cell,particle_set,&
                          atomic_kind_set,mo_coeff,C,auto,nmo,homo,rho_r,rho_g,pot_g,&
                          mat_munu,sab_orb_sub,pw_env_sub,poisson_env,auxbas_pw_pool,task_list_sub,&
                          mo_coeff_o,mo_coeff_v,mp2_env%mp2_gpw%eps_filter,unit_nr,error_sub)

    ! and now free the whole lot 
    DEALLOCATE(C,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(Auto,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL pw_pool_give_back_pw(auxbas_pw_pool,rho_r%pw,error=error_sub)
    CALL pw_pool_give_back_pw(auxbas_pw_pool,rho_g%pw,error=error_sub)
    CALL pw_pool_give_back_pw(auxbas_pw_pool,pot_g%pw,error=error_sub)
    !CALL pw_pool_give_back_pw(auxbas_pw_pool,psi_i%pw,error=error_sub)
    !CALL pw_pool_give_back_pw(auxbas_pw_pool,psi_a%pw,error=error_sub)

    CALL deallocate_task_list(task_list_sub,error=error_sub)

    CALL pw_env_release(pw_env_sub, error=error_sub)

    CALL cp_dbcsr_release(mo_coeff_o,error=error_sub)
    DEALLOCATE(mo_coeff_o, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL cp_dbcsr_release(mo_coeff_v,error=error_sub)
    DEALLOCATE(mo_coeff_v, STAT=stat)
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

    ! re-init the radii to be able to generate pair lists with MP2-appropriate screening
    qs_env%dft_control%qs_control%eps_pgf_orb=eps_pgf_orb_old
    qs_env%dft_control%qs_control%eps_rho_rspace=eps_rho_rspace_old
    qs_env%dft_control%qs_control%eps_gvg_rspace=eps_gvg_rspace_old
    CALL init_interaction_radii(qs_env%dft_control%qs_control,atomic_kind_set, error_sub)


    CALL cp_blacs_env_release(blacs_env_sub, error=error_sub)

    CALL cp_error_dealloc_ref(error_sub,error=error)

    CALL cp_logger_release(logger_sub)

    CALL cp_para_env_release(para_env_sub,error=error)
 
    CALL timestop(handle)

  END SUBROUTINE mp2_gpw_main

  SUBROUTINE mp2_gpw_compute(qs_env,para_env,para_env_sub,color_sub,dft_control,&
                 cell,particle_set,atomic_kind_set,mo_coeff,C,auto,nmo,homo,&
                 rho_r,rho_g,pot_g,mat_munu,sab_orb_sub,pw_env_sub,&
                 poisson_env,auxbas_pw_pool,task_list_sub,mo_coeff_o,mo_coeff_v,eps_filter,unit_nr,&
                 error_sub)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_para_env_type), POINTER          :: para_env, para_env_sub
    INTEGER                                  :: color_sub
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(cell_type), POINTER                 :: cell
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    REAL(KIND=dp), DIMENSION(:, :)           :: C
    REAL(KIND=dp), DIMENSION(:)              :: Auto
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
    TYPE(cp_error_type), INTENT(inout)       :: error_sub

    CHARACTER(LEN=*), PARAMETER :: routineN = 'mp2_gpw_compute', &
      routineP = moduleN//':'//routineN

    INTEGER :: a, a_group_counter, b, blk, col, col_offset, col_size, &
      comm_exchange, EX_end, EX_end_send, EX_start, EX_start_send, &
      group_counter, handle, i, i_counter, i_group_counter, index_proc_shift, &
      j, max_batch_size_A, max_batch_size_I, mepos_in_EX_group, &
      my_A_batch_size, my_A_virtual_end, my_A_virtual_start, my_I_batch_size, &
      my_I_occupied_end, my_I_occupied_start, ngroup, one, p, p_best, &
      proc_receive, proc_send, q, q_best, row, row_offset, row_size, size_EX, &
      size_EX_send, size_of_exchange_group, stat, virtual, wfn_calc, &
      wfn_calc_best
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: proc_map, &
                                                vector_batch_A_size_group, &
                                                vector_batch_I_size_group
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: exchange_group_sizes
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: E_mp2, pair_energy
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: BIb_C, BIb_Ex, BIb_send
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: data_block
    TYPE(array_i1d_obj)                      :: col_blk_size, col_dist, &
                                                row_blk_size, row_dist
    TYPE(cp_dbcsr_iterator)                  :: iter
    TYPE(cp_dbcsr_type)                      :: matrix_ia_jb, matrix_ia_jnu
    TYPE(cp_para_env_type), POINTER          :: para_env_exchange
    TYPE(dbcsr_distribution_obj)             :: dist
    TYPE(pw_p_type)                          :: psi_a
    TYPE(pw_p_type), ALLOCATABLE, &
      DIMENSION(:)                           :: psi_i

    CALL timeset(routineN,handle)
    failure=.FALSE.

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

    virtual=nmo-homo

    ngroup=para_env%num_pe/para_env_sub%num_pe
    wfn_calc_best=HUGE(wfn_calc_best)
    DO p=1,ngroup
       DO q=1,ngroup
          IF (p*q.NE.ngroup) CYCLE
          wfn_calc=((homo+p-1)/p)+((virtual+q-1)/q)
          IF (wfn_calc<wfn_calc_best) THEN
             wfn_calc_best=wfn_calc
             p_best=p
             q_best=q
          ENDIF
       ENDDO
    ENDDO

    max_batch_size_I=(homo+p_best-1)/p_best
    max_batch_size_A=(virtual+q_best-1)/q_best
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

    ! give to each group its (ia| block
    group_counter=0
    i_group_counter=0
    my_I_occupied_start=1
    DO i=0, p_best-1
      my_A_virtual_start=1
      a_group_counter=0
      DO j=0, q_best-1
        group_counter=group_counter+1
        IF(color_sub==group_counter-1) EXIT
        my_A_virtual_start=my_A_virtual_start+vector_batch_A_size_group(j)
        a_group_counter=a_group_counter+1
      END DO
      IF(color_sub==group_counter-1) EXIT
      my_I_occupied_start=my_I_occupied_start+vector_batch_I_size_group(i)
      i_group_counter=i_group_counter+1
    END DO
    my_I_occupied_end=my_I_occupied_start+vector_batch_I_size_group(i_group_counter)-1
    my_I_batch_size=vector_batch_I_size_group(i_group_counter)
    my_A_virtual_end=my_A_virtual_start+vector_batch_A_size_group(a_group_counter)-1
    my_A_batch_size=vector_batch_A_size_group(a_group_counter)

    DEALLOCATE(vector_batch_I_size_group)
    DEALLOCATE(vector_batch_A_size_group)

    ! now create a group that contains all the proc that have the same virtual starting point
    ! in this way it is possible to sum the common integrals needed for the full MP2 energy
    ! in mp_comm_split_direct the color is give by my_a_virtual_start
    ! only proc with rank 0 in the sub_group will be responsible for communication in the exchange group
    ! then setup all is needed for communication
    IF (para_env_sub%mepos==0) THEN
      CALL mp_comm_split_direct(para_env%group,comm_exchange,my_A_virtual_start)
      NULLIFY(para_env_exchange)
      CALL cp_para_env_create(para_env_exchange,comm_exchange,error=error_sub)

      ! crate the proc map
      ALLOCATE(proc_map(-para_env_exchange%num_pe:2*para_env_exchange%num_pe-1))
      DO i=0,para_env_exchange%num_pe-1
        proc_map(i)=i
        proc_map(-i-1)=para_env_exchange%num_pe-i-1
        proc_map(para_env_exchange%num_pe+i)=i
      END DO

      ! create an array containing the information for communication
      ALLOCATE(exchange_group_sizes(0:para_env_exchange%num_pe-1,3))
      exchange_group_sizes=0
      exchange_group_sizes(para_env_exchange%mepos,1)=my_I_occupied_start
      exchange_group_sizes(para_env_exchange%mepos,2)=my_I_occupied_end
      exchange_group_sizes(para_env_exchange%mepos,3)=my_I_batch_size
      CALL mp_sum(exchange_group_sizes,para_env_exchange%group)
      ! broadcast the array to all processes of the sub_group
      mepos_in_EX_group=para_env_exchange%mepos
      size_of_exchange_group=para_env_exchange%num_pe
      CALL mp_bcast(mepos_in_EX_group,0,para_env_sub%group)
      CALL mp_bcast(size_of_exchange_group,0,para_env_sub%group)
      CALL mp_bcast(exchange_group_sizes,0,para_env_sub%group)
      CALL mp_bcast(proc_map,0,para_env_sub%group)

    ELSE
      ! nmo+para_env%mepos should be replaced with MPI_UNDEFINED (add in message_passing.F)
      CALL mp_comm_split_direct(para_env%group,comm_exchange,nmo+para_env%mepos)

      ! get the info for allocating array received from communication
      CALL mp_bcast(mepos_in_EX_group,0,para_env_sub%group)
      CALL mp_bcast(size_of_exchange_group,0,para_env_sub%group)
      ALLOCATE(exchange_group_sizes(0:size_of_exchange_group-1,3)) 
      exchange_group_sizes=0
      CALL mp_bcast(exchange_group_sizes,0,para_env_sub%group)
      ALLOCATE(proc_map(-size_of_exchange_group:2*size_of_exchange_group-1))
      CALL mp_bcast(proc_map,0,para_env_sub%group)
      
    END IF

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
      CALL calculate_wavefunction(mo_coeff,i,psi_i(i),rho_g, &
                        atomic_kind_set,cell,dft_control,particle_set, &
                        pw_env_sub,external_vector=C(:,i),error=error_sub)
    END DO

    E_mp2=0.0_dp
    ALLOCATE(BIb_C(virtual,homo,my_I_batch_size))
    DO a=homo+my_A_virtual_start, homo+my_A_virtual_end
    ! psi_a
      BIb_C=0.0_dp
      CALL calculate_wavefunction(mo_coeff,a,psi_a,rho_g, &
                      atomic_kind_set,cell,dft_control,particle_set, &
                      pw_env_sub,external_vector=C(:,a),error=error_sub)
      i_counter=0
      DO i=my_I_occupied_start, my_I_occupied_end
        i_counter=i_counter+1
        !ia_counter=ia_counter+1
        !IF(MOD(ia_counter,para_env%num_pe/para_env_sub%num_pe).NE.para_env%mepos/para_env_sub%num_pe) CYCLE
        !! psi_i
        !CALL calculate_wavefunction(mo_coeff,i,psi_i,rho_g, &
        !                atomic_kind_set,cell,dft_control,particle_set, &
        !                pw_env_sub,external_vector=C(:,i),error=error_sub)

       ! potential
       rho_r%pw%cr3d = psi_i(i)%pw%cr3d * psi_a%pw%cr3d
       CALL pw_transfer(rho_r%pw, rho_g%pw, error=error_sub)
       CALL pw_poisson_solve(poisson_env,rho_g%pw, pair_energy, pot_g%pw,error=error_sub)
       CALL pw_transfer(pot_g%pw, rho_r%pw, error=error_sub)
       CALL pw_scale(rho_r%pw,rho_r%pw%pw_grid%dvol, error=error_sub)

       ! and finally (ia|munu)
       CALL cp_dbcsr_set(mat_munu%matrix,0.0_dp,error=error_sub)
       CALL integrate_v_rspace(rho_r, h=mat_munu,qs_env=qs_env,calculate_forces=.FALSE.,compute_tau=.FALSE.,gapw=.FALSE.,&
          basis_set_id=use_orb_basis_set, &
          sab_orb_external=sab_orb_sub, pw_env_external=pw_env_sub, task_list_external=task_list_sub, error=error_sub)
       
       ! multiply and goooooooo ...
       CALL cp_dbcsr_multiply("N", "N", 1.0_dp, mat_munu%matrix, mo_coeff_o, &
                               0.0_dp, matrix_ia_jnu, filter_eps=eps_filter, error=error_sub)
       CALL cp_dbcsr_multiply("T", "N", 1.0_dp, matrix_ia_jnu, mo_coeff_v, &
                               0.0_dp, matrix_ia_jb, filter_eps=eps_filter, error=error_sub)

       ! Coulomb only MP2 .....
       CALL cp_dbcsr_iterator_start(iter, matrix_ia_jb)
       DO WHILE (cp_dbcsr_iterator_blocks_left (iter))
          CALL cp_dbcsr_iterator_next_block(iter, row, col, data_block,blk,&
               row_size=row_size, col_size=col_size, &
               row_offset=row_offset, col_offset=col_offset)
          DO b=1,col_size
          DO j=1,row_size
             ! This was for coulom only
             ! E_mp2=E_mp2-data_block(j,b)**2/(Auto(a)+Auto(homo+col_offset+b-1)-Auto(i)-Auto(row_offset+j-1))
             !
             ! now fill the matrix
             BIb_C(col_offset+b-1,row_offset+j-1,i_counter)=data_block(j,b)
          ENDDO
          ENDDO
       ENDDO
       CALL cp_dbcsr_iterator_stop(iter)

      END DO

      ! now mp_sum the matrix BIb in such a way that all 
      ! proc have the data for contracting the FULL mp2 energy
      CALL mp_sum(BIb_C,para_env_sub%group)
      
      ! calculate the contribution to MP2 energy for my local data
      i_counter=0
      DO i=1,my_I_batch_size
        i_counter=i_counter+1
        IF(MOD(i_counter,para_env_sub%num_pe)/=para_env_sub%mepos) CYCLE
        DO j=my_I_occupied_start,my_I_occupied_end
          DO b=1, virtual
            E_mp2=E_mp2-BIb_C(b,j,i) * &
                  (2.0_dp*BIb_C(b,j,i)-&
                  BIb_C(b,i+my_I_occupied_start-1,j-my_I_occupied_start+1))&
                  /(Auto(a)+Auto(homo+b)-Auto(i+my_I_occupied_start-1)-Auto(j))
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

        ALLOCATE(BIb_EX(virtual,my_I_batch_size,size_EX))
        BIb_EX=0.0_dp

        IF(para_env_sub%mepos==0) THEN
          EX_start_send=exchange_group_sizes(proc_send,1)
          EX_end_send=exchange_group_sizes(proc_send,2)
          size_EX_send=exchange_group_sizes(proc_send,3)

          ALLOCATE(BIb_send(virtual,size_EX_send,my_I_batch_size))
          BIb_send(1:virtual,1:size_EX_send,1:my_I_batch_size)=BIb_C(1:virtual,EX_start_send:EX_end_send,1:my_I_batch_size)
      
          ! send and receive the exchange array and broadcast to the other members
          ! of the sub_group
          CALL mp_sendrecv(BIb_send,proc_send,BIb_EX,proc_receive,para_env_exchange%group)
          CALL mp_bcast(BIb_EX,0,para_env_sub%group)
        ELSE
          CALL mp_bcast(BIb_EX,0,para_env_sub%group)
        END IF

        i_counter=0
        DO i=1,my_I_batch_size
          i_counter=i_counter+1
          IF(MOD(i_counter,para_env_sub%num_pe)/=para_env_sub%mepos) CYCLE
          DO j=1,size_EX
            DO b=1, virtual
              E_mp2=E_mp2-BIb_C(b,j+EX_start-1,i) * &
                    (2.0_dp*BIb_C(b,j+EX_start-1,i)-&
                    BIb_EX(b,i,j))&
                    /(Auto(a)+Auto(homo+b)-Auto(i+my_I_occupied_start-1)-Auto(j+EX_start-1))
            END DO
          END DO
        END DO
       

        DEALLOCATE(BIb_EX)
        IF(para_env_sub%mepos==0) DEALLOCATE(BIb_send)

      END DO

    ENDDO

    CALL mp_sum(E_mp2,para_env%group)
    IF(para_env%mepos==0) THEN
      IF (unit_nr>0) WRITE(unit_nr,*)
      IF (unit_nr>0) WRITE(unit_nr,*) 'MP2 Energy', E_mp2
    END IF
    ! let's all dump integrals
    ! CALL cp_dbcsr_print(matrix_ia_jb, nodata=.FALSE., matlab_format=.TRUE., &
    !     variable_name="(ia|jb)", unit_nr=1000+para_env%mepos, error=error_sub)

    DEALLOCATE(BIb_C)
    DEALLOCATE(proc_map)
    DEALLOCATE(exchange_group_sizes)

    IF (para_env_sub%mepos==0) THEN
      CALL cp_para_env_release(para_env_exchange,error=error_sub)
    ENDIF


    CALL cp_dbcsr_release(matrix_ia_jnu,error=error_sub)
    CALL cp_dbcsr_release(matrix_ia_jb,error=error_sub)

    DO i=my_I_occupied_start, my_I_occupied_end  
      CALL pw_pool_give_back_pw(auxbas_pw_pool,psi_i(i)%pw,error=error_sub)
    END DO 
    DEALLOCATE(psi_i)
    CALL pw_pool_give_back_pw(auxbas_pw_pool,psi_a%pw,error=error_sub)

    CALL timestop(handle)

  END SUBROUTINE mp2_gpw_compute

END MODULE mp2_gpw

