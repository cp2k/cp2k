!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Interface 
!> \par History
!>       2011.05 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
MODULE almo_scf_special
  USE almo_scf_methods,                ONLY: almo_scf_get_t_blk,&
                                             almo_scf_ortho_blk,&
                                             almo_scf_p_get_t_blk,&
                                             get_sigma_or_s
  USE almo_scf_types,                  ONLY: almo_mat_dim_aobasis,&
                                             almo_mat_dim_occ,&
                                             almo_mat_dim_virt,&
                                             almo_mat_dim_virt_disc,&
                                             almo_mat_dim_virt_full,&
                                             almo_objectM1_type
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cell_methods,                    ONLY: cell_type,&
                                             pbc
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
  USE cp_dbcsr_interface,              ONLY: &
       array_data, array_i1d_obj, array_new, array_release, &
       cp_dbcsr_allocate_matrix_set, cp_dbcsr_complete_redistribute, &
       cp_dbcsr_copy, cp_dbcsr_copy_into_existing, cp_dbcsr_create, &
       cp_dbcsr_desymmetrize, cp_dbcsr_distribution, &
       cp_dbcsr_distribution_release, cp_dbcsr_filter, cp_dbcsr_finalize, &
       cp_dbcsr_get_block_p, cp_dbcsr_get_info, cp_dbcsr_get_num_blocks, &
       cp_dbcsr_get_stored_coordinates, cp_dbcsr_init, cp_dbcsr_multiply, &
       cp_dbcsr_nblkcols_total, cp_dbcsr_nblkrows_total, cp_dbcsr_p_type, &
       cp_dbcsr_release, cp_dbcsr_reserve_block2d, cp_dbcsr_row_block_sizes, &
       cp_dbcsr_set, cp_dbcsr_type, cp_dbcsr_work_create, &
       dbcsr_distribution_col_dist, dbcsr_distribution_mp, &
       dbcsr_distribution_new, dbcsr_distribution_obj, &
       dbcsr_distribution_row_dist, dbcsr_mp_group, dbcsr_mp_mynode, &
       dbcsr_mp_numnodes, dbcsr_type_no_symmetry
  USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_units,                        ONLY: cp_unit_to_cp2k
  USE input_constants,                 ONLY: almo_constraint_ao_overlap,&
                                             almo_constraint_block_diagonal,&
                                             almo_constraint_distance,&
                                             almo_domain_layout_molecular,&
                                             almo_mat_distr_atomic,&
                                             almo_mat_distr_molecular,&
                                             do_bondparm_covalent,&
                                             do_bondparm_vdw
  USE input_section_types,             ONLY: section_vals_type
  USE iterate_matrix,                  ONLY: matrix_sqrt_Newton_Schulz
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_allgather
  USE molecule_types_new,              ONLY: get_molecule_set_info,&
                                             molecule_type
  USE particle_types,                  ONLY: particle_type
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_initial_guess,                ONLY: calculate_atomic_block_dm,&
                                             calculate_mopac_dm
  USE qs_ks_methods,                   ONLY: qs_ks_create,&
                                             qs_ks_did_change,&
                                             qs_ks_update_qs_env
  USE qs_ks_types,                     ONLY: qs_ks_env_type,&
                                             qs_ks_release
  USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                             neighbor_list_iterate,&
                                             neighbor_list_iterator_create,&
                                             neighbor_list_iterator_p_type,&
                                             neighbor_list_iterator_release,&
                                             neighbor_list_set_p_type
  USE qs_rho_methods,                  ONLY: qs_rho_update_rho
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'almo_scf_special'

  PUBLIC :: almo_level2_spec2_1, almo_levelX_spec6_0, almo_levelX_spec5_0,&
            almo_levelX_spec3_0,&
            almo_level2_spec1_0, almo_levelX_spec2_0

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param matrix_new ...
!> \param matrix_qs ...
!> \param main_objectM1 ...
!> \param name_new ...
!> \param size_keys ...
!> \param symmetry_new ...
!> \param spin_key ...
!> \param init_domains ...
!> \param error ...
!> \par History
!>       2011.05 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE almo_levelX_spec6_0(matrix_new,matrix_qs,main_objectM1,&
                                name_new,size_keys,symmetry_new,&
                                spin_key,init_domains,error)

    TYPE(cp_dbcsr_type)                      :: matrix_new, matrix_qs
    TYPE(almo_objectM1_type), INTENT(IN)     :: main_objectM1
    CHARACTER(len=*), INTENT(IN)             :: name_new
    INTEGER, DIMENSION(2), INTENT(IN)        :: size_keys
    CHARACTER, INTENT(IN)                    :: symmetry_new
    INTEGER, INTENT(IN)                      :: spin_key
    LOGICAL, INTENT(IN)                      :: init_domains
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'almo_levelX_spec6_0', &
      routineP = moduleN//':'//routineN

    INTEGER :: dimen, handle, hold, iatom, iblock_col, iblock_row, imol, &
      mynode, natoms, nblkrows_tot, nlength, nmols, row
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: block_sizes_new, &
                                                distr_new_array
    INTEGER, DIMENSION(:), POINTER           :: blk_distr, blk_sizes
    LOGICAL                                  :: active, failure, tr
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: p_new_block
    TYPE(array_i1d_obj)                      :: col_distr_new, col_sizes_new, &
                                                row_distr_new, row_sizes_new
    TYPE(dbcsr_distribution_obj)             :: dist_new, dist_qs

    CALL timeset(routineN,handle)
    dist_qs = cp_dbcsr_distribution(matrix_qs)
    
    natoms=main_objectM1%natoms
    nmols=main_objectM1%nmolecules
    
    DO dimen=1,2
    
     IF (dimen==1) THEN 
      blk_distr => array_data(dbcsr_distribution_row_dist(dist_qs))
     ELSE
      blk_distr => array_data(dbcsr_distribution_col_dist(dist_qs))
     ENDIF

     IF (size_keys(dimen)==almo_mat_dim_aobasis) THEN 
    
      blk_sizes => array_data(cp_dbcsr_row_block_sizes(matrix_qs))
    
      IF (main_objectM1%mat_distr_aos==almo_mat_distr_atomic) THEN 
       ALLOCATE(block_sizes_new(natoms),distr_new_array(natoms))
       block_sizes_new(:)=blk_sizes(:)
       distr_new_array(:)=blk_distr(:)
      ELSE IF (main_objectM1%mat_distr_aos==almo_mat_distr_molecular) THEN
       ALLOCATE(block_sizes_new(nmols),distr_new_array(nmols))
       block_sizes_new(:)=0
       DO iatom=1,natoms
        block_sizes_new(main_objectM1%domain_index_of_atom(iatom)) = &
          block_sizes_new(main_objectM1%domain_index_of_atom(iatom)) + &
          blk_sizes(iatom)
       ENDDO
       DO imol=1,nmols
        distr_new_array(imol) = &
                  blk_distr(main_objectM1%first_atom_of_domain(imol))
       ENDDO
      ELSE
       CPErrorMessage(cp_failure_level,routineP,"Illegal distribution",error)
       CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
      ENDIF
    
     ELSE 
     
      IF (size_keys(dimen)==almo_mat_dim_occ .OR.&
          size_keys(dimen)==almo_mat_dim_virt .OR. &
          size_keys(dimen)==almo_mat_dim_virt_disc .OR. &
          size_keys(dimen)==almo_mat_dim_virt_full) THEN 
    
       IF (main_objectM1%mat_distr_mos==almo_mat_distr_atomic) THEN
        nlength=natoms
        ALLOCATE(block_sizes_new(nlength))
        IF (size_keys(dimen)==almo_mat_dim_occ) THEN
        ELSE IF (size_keys(dimen)==almo_mat_dim_virt) THEN
        ENDIF
       ELSE IF (main_objectM1%mat_distr_mos==almo_mat_distr_molecular) THEN
        nlength=nmols
        ALLOCATE(block_sizes_new(nlength))
        IF (size_keys(dimen)==almo_mat_dim_occ) THEN
         block_sizes_new(:)=main_objectM1%nocc_of_domain(:,spin_key)
        ELSE IF (size_keys(dimen)==almo_mat_dim_virt_disc) THEN
         block_sizes_new(:)=main_objectM1%nvirt_disc_of_domain(:,spin_key)
        ELSE IF (size_keys(dimen)==almo_mat_dim_virt_full) THEN
         block_sizes_new(:)=main_objectM1%nvirt_full_of_domain(:,spin_key)
        ELSE IF (size_keys(dimen)==almo_mat_dim_virt) THEN
         block_sizes_new(:)=main_objectM1%nvirt_of_domain(:,spin_key)
        ENDIF
       ELSE
        CPErrorMessage(cp_failure_level,routineP,"Illegal distribution",error)
        CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
       ENDIF
      
      ELSE
 
       CPErrorMessage(cp_failure_level,routineP,"Illegal dimension",error)
       CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)

      ENDIF 
      
      ALLOCATE(distr_new_array(nlength))
      IF (main_objectM1%mat_distr_mos==almo_mat_distr_atomic) THEN 
       distr_new_array(:)=blk_distr(:)
      ELSE IF (main_objectM1%mat_distr_mos==almo_mat_distr_molecular) THEN
       DO imol=1,nmols
        distr_new_array(imol) = &
                  blk_distr(main_objectM1%first_atom_of_domain(imol))
       ENDDO
      ENDIF
    
     ENDIF 
     
     IF (dimen==1) THEN 
      CALL array_new(row_sizes_new, block_sizes_new, lb=1)
      CALL array_new(row_distr_new, distr_new_array, lb=1)
     ELSE 
      CALL array_new(col_sizes_new, block_sizes_new, lb=1)
      CALL array_new(col_distr_new, distr_new_array, lb=1)
     ENDIF
     DEALLOCATE(block_sizes_new,distr_new_array)
    
    ENDDO 

    CALL dbcsr_distribution_new(dist_new,&
     dbcsr_distribution_mp(dist_qs),row_distr_new,col_distr_new)
    CALL array_release(col_distr_new)
    CALL array_release(row_distr_new)
    
    CALL cp_dbcsr_init(matrix_new,error=error)
    CALL cp_dbcsr_create (matrix_new, name_new,&
     dist_new, symmetry_new,&
     row_sizes_new, col_sizes_new, error=error)
    CALL array_release(col_sizes_new)
    CALL array_release(row_sizes_new)
    CALL cp_dbcsr_distribution_release(dist_new)

    IF (init_domains) THEN
  
       mynode = dbcsr_mp_mynode(dbcsr_distribution_mp(cp_dbcsr_distribution(matrix_new)))
       CALL cp_dbcsr_work_create(matrix_new, work_mutable=.TRUE., error=error)

       nblkrows_tot = cp_dbcsr_nblkrows_total(matrix_new)
       DO row = 1, nblkrows_tot
          tr = .FALSE.
          iblock_row = row
          iblock_col = row
          CALL cp_dbcsr_get_stored_coordinates(matrix_new, iblock_row, iblock_col, tr, hold)
   
          IF(hold.EQ.mynode) THEN
          
             active=.TRUE.

             IF (active) THEN
                NULLIFY (p_new_block)
                CALL cp_dbcsr_reserve_block2d(matrix_new, iblock_row, iblock_col, p_new_block)
                CPPostcondition(ASSOCIATED(p_new_block),cp_failure_level,routineP,error,failure)
                p_new_block(:,:) = 1.0_dp
             ENDIF

          ENDIF 
       ENDDO

    ENDIF 

    CALL cp_dbcsr_finalize(matrix_new, error=error)

    CALL timestop (handle)

  END SUBROUTINE almo_levelX_spec6_0
  
! *****************************************************************************
!> \brief ...
!> \param matrix_qs ...
!> \param matrix_almo ...
!> \param main_objectM1 ...
!> \param keep_sparsity ...
!> \param error ...
!> \par History
!>       2011.06 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE almo_levelX_spec5_0(matrix_qs,matrix_almo,main_objectM1,&
                               keep_sparsity,error)
    TYPE(cp_dbcsr_type)                      :: matrix_qs, matrix_almo
    TYPE(almo_objectM1_type)                 :: main_objectM1
    LOGICAL, INTENT(IN)                      :: keep_sparsity
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'almo_levelX_spec5_0', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure
    TYPE(cp_dbcsr_type)                      :: matrix_qs_nosym

    CALL timeset(routineN,handle)

    SELECT CASE(main_objectM1%mat_distr_aos)
    CASE(almo_mat_distr_atomic)
       CALL cp_dbcsr_copy(matrix_almo,matrix_qs,&
            keep_sparsity=keep_sparsity,error=error)
    CASE(almo_mat_distr_molecular)
       CALL cp_dbcsr_init (matrix_qs_nosym, error=error)
       CALL cp_dbcsr_create (matrix_qs_nosym, template=matrix_qs,&
            matrix_type=dbcsr_type_no_symmetry, error=error)
       CALL cp_dbcsr_desymmetrize (matrix_qs, matrix_qs_nosym,&
            error=error)

       CALL cp_dbcsr_set(matrix_almo, 0.0_dp, error=error)
       CALL cp_dbcsr_complete_redistribute(matrix_qs_nosym, matrix_almo,&
               keep_sparsity=keep_sparsity,error=error);
       CALL cp_dbcsr_release (matrix_qs_nosym, error=error)

    CASE DEFAULT 
        CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
    END SELECT

    CALL timestop(handle)

  END SUBROUTINE almo_levelX_spec5_0

! *****************************************************************************
!> \brief ...
!> \param matrix_almo ...
!> \param matrix_qs ...
!> \param main_objectM1 ...
!> \param error ...
!> \par History
!>       2011.06 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE almo_levelX_spec4_1(matrix_almo,matrix_qs,main_objectM1,error)
    TYPE(cp_dbcsr_type)                      :: matrix_almo, matrix_qs
    TYPE(almo_objectM1_type), INTENT(IN)     :: main_objectM1
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'almo_levelX_spec4_1', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure

    CALL timeset(routineN,handle)
    SELECT CASE(main_objectM1%mat_distr_aos)
    CASE(almo_mat_distr_atomic)
          CALL cp_dbcsr_copy_into_existing (matrix_qs, matrix_almo, error=error)
    CASE(almo_mat_distr_molecular)
       CALL cp_dbcsr_set(matrix_qs,0.0_dp,error=error)
          CALL cp_dbcsr_complete_redistribute(matrix_almo, matrix_qs, keep_sparsity=.TRUE., error=error) 
    CASE DEFAULT 
       CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
    END SELECT

    CALL timestop(handle)

  END SUBROUTINE almo_levelX_spec4_1

! *****************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param main_objectM1 ...
!> \param error ...
!> \par History
!>       2011.05 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE almo_level2_spec2_1(qs_env,main_objectM1,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(almo_objectM1_type)                 :: main_objectM1
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'almo_level2_spec2_1', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, nspin, unit_nr
    INTEGER, DIMENSION(2)                    :: nelectron_spin
    LOGICAL                                  :: has_unit_metric
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_s
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(section_vals_type), POINTER         :: input

    CALL timeset(routineN,handle)

    logger => cp_error_get_logger(error)
    IF (logger%para_env%mepos==logger%para_env%source) THEN
       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
    ELSE
       unit_nr=-1
    ENDIF
    
    CALL get_qs_env(qs_env,&
                    dft_control=dft_control,&
                    matrix_s=matrix_s,&
                    matrix_ks=matrix_ks,&
                    ks_env=ks_env,&
                    input=input,&
                    energy=energy,&
                    atomic_kind_set=atomic_kind_set,&
                    particle_set=particle_set,&
                    has_unit_metric=has_unit_metric,&
                    para_env=para_env,&
                    nelectron_spin=nelectron_spin,&
                    error=error)

    nspin=dft_control%nspins

    IF (.NOT.ASSOCIATED(matrix_ks)) THEN
       CALL cp_dbcsr_allocate_matrix_set(matrix_ks,nspin,error)
       DO ispin=1,nspin
          ALLOCATE(matrix_ks(ispin)%matrix)
          CALL cp_dbcsr_init(matrix_ks(ispin)%matrix,error=error)
          CALL cp_dbcsr_create(matrix_ks(ispin)%matrix,&
               template=matrix_s(1)%matrix,error=error)
          CALL cp_dbcsr_alloc_block_from_nbl(matrix_ks(ispin)%matrix,&
               qs_env%sab_orb,error=error)
          CALL cp_dbcsr_set(matrix_ks(ispin)%matrix,0.0_dp,error=error)
       ENDDO
       CALL set_qs_env(qs_env=qs_env,matrix_ks=matrix_ks,error=error)
    ENDIF
    IF (.NOT.ASSOCIATED(ks_env)) THEN
       CALL qs_ks_create(ks_env,qs_env=qs_env,error=error)
       CALL set_qs_env(qs_env, ks_env=ks_env,error=error)
       CALL qs_ks_release(ks_env,error=error)
    END IF

    DO ispin=1,nspin
       CALL cp_dbcsr_set(qs_env%rho%rho_ao(ispin)%matrix,0.0_dp,error=error)
       IF (dft_control%qs_control%dftb .OR. dft_control%qs_control%semi_empirical) THEN
          CALL calculate_mopac_dm(qs_env%rho%rho_ao(ispin)%matrix,&
               matrix_s(1)%matrix, has_unit_metric, &
               dft_control,particle_set, atomic_kind_set,&
               ispin, nspin, nelectron_spin(ispin),&
               para_env, error)
       ELSE
          CALL calculate_atomic_block_dm(qs_env%rho%rho_ao(ispin)%matrix,&
               matrix_s(1)%matrix, particle_set, atomic_kind_set, &
               ispin, nspin, nelectron_spin(ispin), unit_nr, error)
       ENDIF

       CALL almo_levelX_spec5_0(qs_env%rho%rho_ao(ispin)%matrix,&
                              main_objectM1%dpp_b(ispin),main_objectM1,&
                              .FALSE.,error)
       CALL cp_dbcsr_filter(main_objectM1%dpp_b(ispin),&
               main_objectM1%eps_filter,error=error)

    ENDDO

    CALL almo_scf_get_t_blk(main_objectM1,error)
    CALL almo_scf_ortho_blk(main_objectM1,error)
    CALL almo_scf_p_get_t_blk(main_objectM1,&
            use_sigma_inv_guess=.FALSE.,error=error)

    DO ispin=1,nspin
       CALL almo_levelX_spec4_1(main_objectM1%dpp(ispin),&
                              qs_env%rho%rho_ao(ispin)%matrix,&
                              main_objectM1,error)
    ENDDO

    CALL qs_rho_update_rho(qs_env%rho,qs_env=qs_env, error=error)
    CALL qs_ks_did_change(qs_env%ks_env,rho_changed=.TRUE., error=error)
    CALL qs_ks_update_qs_env(qs_env%ks_env,qs_env=qs_env,&
         calculate_forces=.FALSE., just_energy=.FALSE.,error=error)

    IF (unit_nr>0) THEN
       WRITE(unit_nr,'(T2,A,F20.9)') "Energy of the initial guess:",energy%total
       WRITE(unit_nr,'()') 
    ENDIF

    CALL timestop(handle)

  END SUBROUTINE almo_level2_spec2_1

! *****************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param main_objectM1 ...
!> \param energy_new ...
!> \param error ...
!> \par History
!>       2011.05 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE almo_levelX_spec3_0(qs_env,main_objectM1,energy_new,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(almo_objectM1_type)                 :: main_objectM1
    REAL(KIND=dp)                            :: energy_new
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'almo_levelX_spec3_0', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, nspin
    TYPE(qs_energy_type), POINTER            :: energy

    CALL timeset(routineN,handle)

    nspin=main_objectM1%nspins
    CALL get_qs_env(qs_env, energy=energy, error=error)

    DO ispin=1,nspin
       CALL almo_levelX_spec4_1(main_objectM1%dpp(ispin),&
                              qs_env%rho%rho_ao(ispin)%matrix,&
                              main_objectM1,&
                              error)
    END DO

    CALL qs_rho_update_rho(qs_env%rho,qs_env=qs_env, error=error)
    CALL qs_ks_did_change(qs_env%ks_env,rho_changed=.TRUE., error=error)
    CALL qs_ks_update_qs_env(qs_env%ks_env,qs_env=qs_env,&
         calculate_forces=.FALSE.,just_energy=.FALSE.,&
         print_active=.TRUE., error=error)
    energy_new=energy%total

    CALL timestop(handle)

  END SUBROUTINE almo_levelX_spec3_0

! *****************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param energy ...
!> \param error ...
!> \par History
!>       2013.03 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE almo_levelX_spec2_0(qs_env,energy,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    REAL(KIND=dp)                            :: energy
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'almo_levelX_spec2_0', &
      routineP = moduleN//':'//routineN

    TYPE(qs_energy_type), POINTER            :: qs_energy

    CALL get_qs_env(qs_env, energy=qs_energy, error=error)
    qs_energy%total=energy

  END SUBROUTINE almo_levelX_spec2_0

! *****************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param main_objectM1 ...
!> \param error ...
!> \par History
!>       2011.11 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE almo_levelX_spec1_0(qs_env,main_objectM1,error)
    
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(almo_objectM1_type)                 :: main_objectM1
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'almo_levelX_spec1_0', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, unit_nr
    LOGICAL                                  :: failure
    TYPE(cp_dbcsr_type)                      :: matrix_sigma_sqrt, &
                                                matrix_sigma_sqrt_inv, &
                                                matrix_tmp1
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)
    failure = .FALSE.

    logger => cp_error_get_logger(error)
    IF (logger%para_env%mepos==logger%para_env%source) THEN
       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
    ELSE
       unit_nr=-1
    ENDIF

    DO ispin=1,main_objectM1%nspins

       CALL cp_dbcsr_init(matrix_sigma_sqrt,error=error)
       CALL cp_dbcsr_init(matrix_sigma_sqrt_inv,error=error)
       CALL cp_dbcsr_create(matrix_sigma_sqrt,&
               template=main_objectM1%rem(ispin),&
               matrix_type=dbcsr_type_no_symmetry,error=error)
       CALL cp_dbcsr_create(matrix_sigma_sqrt_inv,&
               template=main_objectM1%rem(ispin),&
               matrix_type=dbcsr_type_no_symmetry,error=error)
       CALL get_sigma_or_s(bra=main_objectM1%enter(ispin),&
               ket=main_objectM1%enter(ispin),&
               overlap=main_objectM1%rem(ispin),&
               metric=main_objectM1%matrix_s(1),&
               retain_overlap_sparsity=.FALSE.,&
               eps_filter=main_objectM1%eps_filter,&
               error=error)
       CALL matrix_sqrt_Newton_Schulz(matrix_sigma_sqrt,&
                                      matrix_sigma_sqrt_inv,&
                                      main_objectM1%rem(ispin),&
                                      main_objectM1%eps_filter, &
                                      3, 1.0E-4_dp, 40, error=error)
       CALL cp_dbcsr_release(matrix_sigma_sqrt,error=error)
       CALL cp_dbcsr_init(matrix_tmp1,error=error)
       CALL cp_dbcsr_create(matrix_tmp1,&
               template=main_objectM1%enter(ispin),&
               matrix_type=dbcsr_type_no_symmetry,error=error)
       CALL cp_dbcsr_multiply("N","N",1.0_dp,main_objectM1%enter(ispin),&
               matrix_sigma_sqrt_inv,&
               0.0_dp,matrix_tmp1,filter_eps=main_objectM1%eps_filter,&
               error=error)
       CALL copy_dbcsr_to_fm(matrix_tmp1,&
               qs_env%mos(ispin)%mo_set%mo_coeff,error=error)
       CALL almo_levelX_spec4_1(main_objectM1%dpp(ispin),&
                              qs_env%rho%rho_ao(ispin)%matrix,&
                              main_objectM1,error)

       CALL cp_dbcsr_release(matrix_sigma_sqrt_inv,error=error)
       CALL cp_dbcsr_release(matrix_tmp1,error=error)
    ENDDO

    CALL timestop(handle)

  END SUBROUTINE almo_levelX_spec1_0

! *****************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param main_objectM1 ...
!> \param error ...
!> \par History
!>       2011.11 created [Rustam Z. Khaliullin]
!> \author Rustam Z. Khaliullin
! *****************************************************************************
  SUBROUTINE almo_level2_spec1_0(qs_env,main_objectM1,error)
    
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(almo_objectM1_type), INTENT(INOUT)  :: main_objectM1
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'almo_level2_spec1_0', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: sym
    INTEGER :: col, contact_atom_1, contact_atom_2, domain_col, &
      domain_map_local_entries, domain_row, global_entries, &
      global_list_length, grid1, GroupID, handle, hold, iatom, iatom2, &
      iblock_col, iblock_row, idomain, idomain2, ientry, igrid, ineig, &
      ineighbor, iNode, inode2, ipair, ispin, jatom, jatom2, jdomain2, &
      local_list_length, max_domain_neighbors, max_neig, mynode, &
      nblkcols_tot, nblkrows_tot, nblks, ndomains, neig_temp, nnode2, nNodes, &
      row, unit_nr
    INTEGER, ALLOCATABLE, DIMENSION(:) :: current_number_neighbors, &
      domain_entries_cpu, domain_map_global, domain_map_local, &
      first_atom_of_molecule, global_list, last_atom_of_molecule, &
      list_length_cpu, list_offset_cpu, local_list, offset_for_cpu
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: domain_grid, &
                                                domain_neighbor_list, &
                                                domain_neighbor_list_excessive
    LOGICAL                                  :: already_listed, block_active, &
                                                failure, found, &
                                                max_neig_fails, tr
    REAL(KIND=dp) :: contact1_radius, contact2_radius, distance, &
      distance_squared, overlap, r0, r1, s0, s1, trial_distance_squared
    REAL(KIND=dp), DIMENSION(3)              :: rab
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: p_new_block
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    TYPE(cp_dbcsr_type)                      :: matrix_s_sym
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: nl_iterator, nl_iterator2
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_almo
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

    CALL timeset(routineN,handle)
       
    logger => cp_error_get_logger(error)
    IF (logger%para_env%mepos==logger%para_env%source) THEN
       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
    ELSE
       unit_nr=-1
    ENDIF

    ndomains=main_objectM1%ndomains

    CALL get_qs_env(qs_env=qs_env,&
      particle_set=particle_set,&
      molecule_set=molecule_set,&
      cell=cell,&
      matrix_s=matrix_s,&
      sab_almo=sab_almo,&
      error=error)

    IF (main_objectM1%domain_layout_mos==almo_domain_layout_molecular .OR. &
       main_objectM1%domain_layout_aos==almo_domain_layout_molecular) THEN
       ALLOCATE(first_atom_of_molecule(main_objectM1%nmolecules))
       ALLOCATE(last_atom_of_molecule(main_objectM1%nmolecules))
       CALL get_molecule_set_info(molecule_set,&
               mol_to_first_atom=first_atom_of_molecule,&
               mol_to_last_atom=last_atom_of_molecule,&
               error=error)
    ENDIF

    CALL cp_dbcsr_init(matrix_s_sym,error=error)
    CALL cp_dbcsr_create(matrix_s_sym,&
            template=main_objectM1%matrix_s(1),&
            matrix_type=dbcsr_type_no_symmetry,error=error)
    CALL cp_dbcsr_get_info(main_objectM1%matrix_s(1),&
            matrix_type=sym)
    IF (sym.eq.dbcsr_type_no_symmetry) THEN
       CALL cp_dbcsr_copy(matrix_s_sym,main_objectM1%matrix_s(1),error=error)
    ELSE
       CALL cp_dbcsr_desymmetrize(main_objectM1%matrix_s(1),&
               matrix_s_sym,error=error)
    ENDIF

    ALLOCATE(main_objectM1%quench_t(main_objectM1%nspins))
    ALLOCATE(main_objectM1%domain_map(main_objectM1%nspins))

    ispin=1

       CALL almo_levelX_spec6_0(matrix_new=main_objectM1%quench_t(ispin),&
            matrix_qs=matrix_s(1)%matrix,&
            main_objectM1=main_objectM1,&
            name_new="T_QUENCHER",&
            size_keys=(/almo_mat_dim_aobasis,almo_mat_dim_occ/),&
            symmetry_new=dbcsr_type_no_symmetry,&
            spin_key=ispin,&
            init_domains=.FALSE.,&
            error=error)
    
       CALL cp_dbcsr_work_create(main_objectM1%quench_t(ispin),&
               work_mutable=.TRUE., error=error)

       nblkrows_tot = cp_dbcsr_nblkrows_total(main_objectM1%quench_t(ispin))
       nblkcols_tot = cp_dbcsr_nblkcols_total(main_objectM1%quench_t(ispin))

       mynode = dbcsr_mp_mynode(dbcsr_distribution_mp(&
          cp_dbcsr_distribution(main_objectM1%quench_t(ispin))))
       nNodes = dbcsr_mp_numnodes(dbcsr_distribution_mp(&
          cp_dbcsr_distribution(main_objectM1%quench_t(ispin))))
       GroupID = dbcsr_mp_group(dbcsr_distribution_mp(&
          cp_dbcsr_distribution(main_objectM1%quench_t(ispin))))

       local_list_length=0
       CALL neighbor_list_iterator_create(nl_iterator,sab_almo)
       DO WHILE (neighbor_list_iterate(nl_iterator)==0)
          CALL get_iterator_info(nl_iterator,&
                  iatom=iatom2,jatom=jatom2,inode=inode2,nnode=nnode2)
          IF (inode2==1) THEN
             local_list_length=local_list_length+nnode2
          END IF
       END DO
       CALL neighbor_list_iterator_release(nl_iterator)

       ALLOCATE(local_list(2*local_list_length))
       local_list(:)=0
       local_list_length=0
       CALL neighbor_list_iterator_create(nl_iterator2,sab_almo)
       DO WHILE (neighbor_list_iterate(nl_iterator2)==0)
          CALL get_iterator_info(nl_iterator2,&
                  iatom=iatom2,jatom=jatom2)
          local_list(2*local_list_length+1)=iatom2
          local_list(2*local_list_length+2)=jatom2
          local_list_length=local_list_length+1
       ENDDO 
       CALL neighbor_list_iterator_release(nl_iterator2)

       ALLOCATE(list_length_cpu(nNodes),list_offset_cpu(nNodes))
       CALL mp_allgather(2*local_list_length,list_length_cpu,GroupID)
       
       list_offset_cpu(1)=0
       DO iNode=2,nNodes
          list_offset_cpu(iNode)=list_offset_cpu(iNode-1) + &
             list_length_cpu(iNode-1)
       ENDDO
       global_list_length=list_offset_cpu(nNodes)+list_length_cpu(nNodes)
       global_list_length=global_list_length/2
 
       ALLOCATE(global_list(2*global_list_length))
       CALL mp_allgather(local_list,global_list,&
               list_length_cpu,list_offset_cpu,GroupID)
       DEALLOCATE(list_length_cpu,list_offset_cpu)
       DEALLOCATE(local_list)
       
       ALLOCATE(current_number_neighbors(main_objectM1%ndomains))
       current_number_neighbors(:)=0
       DO ipair=1, global_list_length
          iatom2=global_list(2*(ipair-1)+1)
          jatom2=global_list(2*(ipair-1)+2)
          idomain2=main_objectM1%domain_index_of_atom(iatom2)
          jdomain2=main_objectM1%domain_index_of_atom(jatom2)
          current_number_neighbors(idomain2)=current_number_neighbors(idomain2)+1
          IF (idomain2.ne.jdomain2) THEN
             current_number_neighbors(jdomain2)=current_number_neighbors(jdomain2)+1
          ENDIF
       ENDDO 
 
       max_domain_neighbors=MAXVAL(current_number_neighbors)
       ALLOCATE(domain_neighbor_list_excessive(ndomains,max_domain_neighbors))
       current_number_neighbors(:)=1
       DO ipair=1, ndomains
          domain_neighbor_list_excessive(ipair,1)=ipair
       ENDDO
       DO ipair=1, global_list_length
          iatom2=global_list(2*(ipair-1)+1)
          jatom2=global_list(2*(ipair-1)+2)
          idomain2=main_objectM1%domain_index_of_atom(iatom2)
          jdomain2=main_objectM1%domain_index_of_atom(jatom2)
          already_listed=.FALSE.
          DO ineighbor=1,current_number_neighbors(idomain2)
             IF (domain_neighbor_list_excessive(idomain2,ineighbor).eq.jdomain2) THEN
                already_listed=.TRUE.
                EXIT
             ENDIF
          ENDDO
          IF (.NOT.already_listed) THEN
             current_number_neighbors(idomain2)=current_number_neighbors(idomain2)+1
             domain_neighbor_list_excessive(idomain2,current_number_neighbors(idomain2))=jdomain2
             IF (idomain2.ne.jdomain2) THEN
                current_number_neighbors(jdomain2)=current_number_neighbors(jdomain2)+1
                domain_neighbor_list_excessive(jdomain2,current_number_neighbors(jdomain2))=idomain2
             ENDIF
          ENDIF
       ENDDO
       DEALLOCATE(global_list)
 
       max_domain_neighbors=MAXVAL(current_number_neighbors)
       ALLOCATE(domain_neighbor_list(ndomains,max_domain_neighbors))
       domain_neighbor_list(:,:)=0
       domain_neighbor_list(:,:)=domain_neighbor_list_excessive(:,1:max_domain_neighbors)
       DEALLOCATE(domain_neighbor_list_excessive)

       ALLOCATE(main_objectM1%domain_map(ispin)%index1(ndomains))
       ALLOCATE(main_objectM1%domain_map(ispin)%pairs(max_domain_neighbors*ndomains,2))
       main_objectM1%domain_map(ispin)%pairs(:,:)=0
       main_objectM1%domain_map(ispin)%index1(:)=0
       domain_map_local_entries=0
       
       
       DO row = 1, nblkrows_tot
          DO col = 1, current_number_neighbors(row)
             tr = .FALSE.
             iblock_row = row
             iblock_col = domain_neighbor_list(row,col)
             CALL cp_dbcsr_get_stored_coordinates(main_objectM1%quench_t(ispin),&
                     iblock_row, iblock_col, tr, hold)

             IF(hold.EQ.mynode) THEN

                domain_row=main_objectM1%domain_index_of_ao_block(iblock_row)
                domain_col=main_objectM1%domain_index_of_mo_block(iblock_col)
                
                SELECT CASE (main_objectM1%constraint_type)
                CASE (almo_constraint_block_diagonal)
               
                   block_active=.FALSE.
                   IF (main_objectM1%domain_layout_mos==almo_domain_layout_molecular) THEN
   
                      IF (main_objectM1%domain_layout_aos==almo_domain_layout_molecular) THEN
   
                         IF (domain_row==domain_col) THEN
                            block_active=.TRUE.
                         ENDIF
   
                      ELSE 
   
                         CPErrorMessage(cp_failure_level,routineP,"Illegal: atomic domains and molecular groups",error)
                         CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
   
                      ENDIF
                   
                   ELSE 
                      
                      IF (main_objectM1%domain_layout_aos==almo_domain_layout_molecular) THEN
   
                         CPErrorMessage(cp_failure_level,routineP,"Illegal: molecular domains and atomic groups",error)
                         CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
   
                      ELSE 
                         
                         IF (domain_row==domain_col) THEN
                            block_active=.TRUE.
                         ENDIF
                      
                      ENDIF
                   
                   ENDIF 
                   
                   IF ( block_active ) THEN
                      
                      NULLIFY (p_new_block)
                      CALL cp_dbcsr_reserve_block2d(main_objectM1%quench_t(ispin),&
                              iblock_row, iblock_col, p_new_block)
                      CPPostcondition(ASSOCIATED(p_new_block),cp_failure_level,routineP,error,failure)
                      p_new_block(:,:) = 1.0_dp
                      
                      IF (domain_map_local_entries.ge.max_domain_neighbors*main_objectM1%ndomains) THEN
                         CPErrorMessage(cp_failure_level,routineP,"weird... max_domain_neighbors is exceeded",error)
                         CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
                      ENDIF
                      main_objectM1%domain_map(ispin)%pairs(domain_map_local_entries+1,1)=iblock_row
                      main_objectM1%domain_map(ispin)%pairs(domain_map_local_entries+1,2)=iblock_col
                      domain_map_local_entries=domain_map_local_entries+1
                   
                   ENDIF

                CASE (almo_constraint_ao_overlap)

                   IF (main_objectM1%domain_layout_mos==almo_domain_layout_molecular) THEN
   
                      IF (main_objectM1%domain_layout_aos==almo_domain_layout_molecular) THEN
   
                         CALL cp_dbcsr_get_block_p(matrix_s_sym,&
                                 iblock_row, iblock_col, p_new_block, found)
                         IF (found) THEN
                            overlap=MAXVAL(ABS(p_new_block))
                         ELSE
                            overlap=0.0_dp
                         ENDIF
   
                      ELSE 
   
                         CPErrorMessage(cp_failure_level,routineP,"atomic domains and molecular groups - NYI",error)
                         CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
   
                      ENDIF
                   
                   ELSE 
                      
                      IF (main_objectM1%domain_layout_aos==almo_domain_layout_molecular) THEN
   
                         CPErrorMessage(cp_failure_level,routineP,"molecular domains and atomic groups - NYI",error)
                         CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
   
                      ELSE 
                         
                         CALL cp_dbcsr_get_block_p(matrix_s_sym,&
                                 iblock_row, iblock_col, p_new_block, found)
                         IF (found) THEN
                            overlap=MAXVAL(ABS(p_new_block))
                         ELSE
                            overlap=0.0_dp
                         ENDIF
                      
                      ENDIF
                   
                   ENDIF 
                   
                   s0=-LOG10(ABS(main_objectM1%quencher_s0))
                   s1=-LOG10(ABS(main_objectM1%quencher_s1))
                   IF (overlap.eq.0.0_dp) THEN
                      overlap=-LOG10(ABS(main_objectM1%eps_filter))+100.0_dp
                   ELSE
                      overlap=-LOG10(overlap)
                   ENDIF
                   IF ( s0.lt.0.0_dp ) THEN
                      CALL stop_program(routineN,moduleN,__LINE__,&
                              "S0 is less than zero")
                   ENDIF
                   IF ( s1.le.0.0_dp ) THEN
                      CALL stop_program(routineN,moduleN,__LINE__,&
                              "S1 is less than or equal to zero")
                   ENDIF
                   IF ( s0.ge.s1 ) THEN
                      CALL stop_program(routineN,moduleN,__LINE__,&
                              "S0 is greater than or equal to S1")
                   ENDIF

                   IF ( overlap.lt.s1 ) THEN
                      NULLIFY (p_new_block)
                      CALL cp_dbcsr_reserve_block2d(main_objectM1%quench_t(ispin),&
                              iblock_row, iblock_col, p_new_block)
                      CPPostcondition(ASSOCIATED(p_new_block),cp_failure_level,routineP,error,failure)
                      IF ( overlap.le.s0 ) THEN
                         p_new_block(:,:) = 1.0_dp
                      ELSE
                         p_new_block(:,:) = 1.0_dp/(1.0_dp+EXP(-(s0-s1)/(s0-overlap)-(s0-s1)/(overlap-s1)))
                      ENDIF
                      
                      IF (ABS(p_new_block(1,1)).gt.ABS(main_objectM1%eps_filter)) THEN
                         IF (domain_map_local_entries.ge.max_domain_neighbors*main_objectM1%ndomains) THEN
                            CPErrorMessage(cp_failure_level,routineP,"weird... max_domain_neighbors is exceeded",error)
                            CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
                         ENDIF
                         main_objectM1%domain_map(ispin)%pairs(domain_map_local_entries+1,1)=iblock_row
                         main_objectM1%domain_map(ispin)%pairs(domain_map_local_entries+1,2)=iblock_col
                         domain_map_local_entries=domain_map_local_entries+1
                      ENDIF
                   
                   ENDIF

                CASE (almo_constraint_distance)

                   IF (main_objectM1%domain_layout_mos==almo_domain_layout_molecular) THEN
   
                      IF (main_objectM1%domain_layout_aos==almo_domain_layout_molecular) THEN
   
                         IF (domain_row==domain_col) THEN
                            distance=0.0_dp
                            contact_atom_1=first_atom_of_molecule(domain_row)
                            contact_atom_2=first_atom_of_molecule(domain_col)
                         ELSE
                            distance_squared=1.0E+100_dp
                            contact_atom_1=-1
                            contact_atom_2=-1
                            DO iatom=first_atom_of_molecule(domain_row),last_atom_of_molecule(domain_row)
                               DO jatom=first_atom_of_molecule(domain_col),last_atom_of_molecule(domain_col)
                                  rab(:) = pbc(particle_set(iatom)%r(:),particle_set(jatom)%r(:),cell)
                                  trial_distance_squared = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
                                  IF (trial_distance_squared.lt.distance_squared) THEN
                                     distance_squared=trial_distance_squared
                                     contact_atom_1=iatom
                                     contact_atom_2=jatom
                                  ENDIF
                               ENDDO 
                            ENDDO 
                            CPPrecondition(contact_atom_1.gt.0,cp_failure_level,routineP,error,failure)
                            distance=SQRT(distance_squared)
                         ENDIF
   
                      ELSE 
   
                         CPErrorMessage(cp_failure_level,routineP,"atomic domains and molecular groups - NYI",error)
                         CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
   
                      ENDIF
                   
                   ELSE 
                      
                      IF (main_objectM1%domain_layout_aos==almo_domain_layout_molecular) THEN
   
                         CPErrorMessage(cp_failure_level,routineP,"molecular domains and atomic groups - NYI",error)
                         CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
   
                      ELSE 
                         
                         rab(:) = pbc(particle_set(domain_row)%r(:),particle_set(domain_col)%r(:),cell)
                         distance = SQRT(rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3))
                         contact_atom_1=domain_row
                         contact_atom_2=domain_col
                      
                      ENDIF
                   
                   ENDIF 
                   
                   IF(main_objectM1%quencher_radius_type==do_bondparm_covalent) THEN
                      CALL get_atomic_kind(atomic_kind=particle_set(contact_atom_1)%atomic_kind,&
                              rcov=contact1_radius)
                      CALL get_atomic_kind(atomic_kind=particle_set(contact_atom_2)%atomic_kind,&
                              rcov=contact2_radius)
                   ELSE IF(main_objectM1%quencher_radius_type==do_bondparm_vdw) THEN
                      CALL get_atomic_kind(atomic_kind=particle_set(contact_atom_1)%atomic_kind,&
                              rvdw=contact1_radius)
                      CALL get_atomic_kind(atomic_kind=particle_set(contact_atom_2)%atomic_kind,&
                              rvdw=contact2_radius)
                   ELSE
                      CALL stop_program(routineN,moduleN,__LINE__,&
                                        "Illegal quencher_radius_type")
                   END IF
                   contact1_radius = cp_unit_to_cp2k(contact1_radius,"angstrom",error=error)
                   contact2_radius = cp_unit_to_cp2k(contact2_radius,"angstrom",error=error)

                   r0=main_objectM1%quencher_r0_factor*(contact1_radius+contact2_radius)
                   r1=main_objectM1%quencher_r1_factor*(contact1_radius+contact2_radius)
                   
                   IF ( r0.lt.0.0_dp ) THEN
                      CALL stop_program(routineN,moduleN,__LINE__,&
                              "R0 is less than zero")
                   ENDIF
                   IF ( r1.le.0.0_dp ) THEN
                      CALL stop_program(routineN,moduleN,__LINE__,&
                              "R1 is less than or equal to zero")
                   ENDIF
                   IF ( r0.gt.r1 ) THEN
                      CALL stop_program(routineN,moduleN,__LINE__,&
                              "R0 is greater than or equal to R1")
                   ENDIF
                   
                   IF ( distance.lt.r1 ) THEN
                      NULLIFY (p_new_block)
                      CALL cp_dbcsr_reserve_block2d(main_objectM1%quench_t(ispin),&
                              iblock_row, iblock_col, p_new_block)
                      CPPostcondition(ASSOCIATED(p_new_block),cp_failure_level,routineP,error,failure)
                      IF ( distance.le.r0 ) THEN
                         p_new_block(:,:) = 1.0_dp
                      ELSE
CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
                         p_new_block(:,:) = 1.0_dp/(1.0_dp+EXP((r1-r0)/(r0-distance)+(r1-r0)/(r1-distance)))
                      ENDIF
                      
                      IF (ABS(p_new_block(1,1)).gt.ABS(main_objectM1%eps_filter)) THEN
                         IF (domain_map_local_entries.ge.max_domain_neighbors*main_objectM1%ndomains) THEN
                            CPErrorMessage(cp_failure_level,routineP,"weird... max_domain_neighbors is exceeded",error)
                            CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
                         ENDIF
                         main_objectM1%domain_map(ispin)%pairs(domain_map_local_entries+1,1)=iblock_row
                         main_objectM1%domain_map(ispin)%pairs(domain_map_local_entries+1,2)=iblock_col
                         domain_map_local_entries=domain_map_local_entries+1
                      ENDIF

                   ENDIF

                CASE DEFAULT
                   CALL stop_program(routineN,moduleN,__LINE__,&
                         "Illegal constraint type")
                END SELECT

             ENDIF 

          ENDDO
       ENDDO 
       
       DEALLOCATE(domain_neighbor_list)
       DEALLOCATE(current_number_neighbors)

       CALL cp_dbcsr_finalize(main_objectM1%quench_t(ispin),error=error)
       CALL cp_dbcsr_filter(main_objectM1%quench_t(ispin),&
               main_objectM1%eps_filter,error=error)

       nblks=cp_dbcsr_get_num_blocks(main_objectM1%quench_t(ispin))
       IF (nblks.ne.domain_map_local_entries) THEN
          CPErrorMessage(cp_failure_level,routineP,"number of blocks is wrong",error)
          CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
       ENDIF

       ALLOCATE(domain_entries_cpu(nNodes),offset_for_cpu(nNodes))
       CALL mp_allgather(2*domain_map_local_entries,domain_entries_cpu,GroupID)
   
       offset_for_cpu(1)=0
       DO iNode=2,nNodes
          offset_for_cpu(iNode)=offset_for_cpu(iNode-1) + &
             domain_entries_cpu(iNode-1)
       ENDDO
       global_entries=offset_for_cpu(nNodes)+domain_entries_cpu(nNodes)
       
       ALLOCATE(domain_map_global(global_entries))
       ALLOCATE(domain_map_local(2*domain_map_local_entries))
       DO ientry=1,domain_map_local_entries
          domain_map_local(2*(ientry-1)+1)=main_objectM1%domain_map(ispin)%pairs(ientry,1)
          domain_map_local(2*ientry)=main_objectM1%domain_map(ispin)%pairs(ientry,2)
       ENDDO
       CALL mp_allgather(domain_map_local,domain_map_global,&
               domain_entries_cpu,offset_for_cpu,GroupID)
       DEALLOCATE(domain_entries_cpu,offset_for_cpu)
       DEALLOCATE(domain_map_local)
       
       DEALLOCATE(main_objectM1%domain_map(ispin)%index1)
       DEALLOCATE(main_objectM1%domain_map(ispin)%pairs)
       ALLOCATE(main_objectM1%domain_map(ispin)%index1(ndomains))
       ALLOCATE(main_objectM1%domain_map(ispin)%pairs(global_entries/2,2))
       main_objectM1%domain_map(ispin)%pairs(:,:)=0
       main_objectM1%domain_map(ispin)%index1(:)=0
   
       max_neig=max_domain_neighbors
       max_neig_fails=.TRUE.
       max_neig_loop: DO WHILE (max_neig_fails)
          ALLOCATE(domain_grid(main_objectM1%ndomains,0:max_neig))
          domain_grid(:,:)=0
          domain_grid(:,0)=1
          global_entries=global_entries/2
          DO ientry=1,global_entries
             grid1=domain_map_global(2*ientry)
             ineig=domain_map_global(2*(ientry-1)+1)
             DO igrid=1,domain_grid(grid1,0)
                IF (ineig.lt.domain_grid(grid1,igrid)) THEN
                   neig_temp=ineig
                   ineig=domain_grid(grid1,igrid)
                   domain_grid(grid1,igrid)=neig_temp
                ELSE
                   IF (domain_grid(grid1,igrid).eq.0) THEN
                      domain_grid(grid1,igrid)=ineig
                      domain_grid(grid1,0)=domain_grid(grid1,0)+1
                      IF (domain_grid(grid1,0).gt.max_neig) THEN
                         DEALLOCATE(domain_grid)
                         max_neig=max_neig*2
                         CYCLE max_neig_loop
                      ENDIF
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
          max_neig_fails=.FALSE.
       ENDDO max_neig_loop
       DEALLOCATE(domain_map_global)
       
       ientry=1
       DO idomain=1,main_objectM1%ndomains
          DO ineig=1,domain_grid(idomain,0)-1
             main_objectM1%domain_map(ispin)%pairs(ientry,1)=domain_grid(idomain,ineig)
             main_objectM1%domain_map(ispin)%pairs(ientry,2)=idomain
             ientry=ientry+1
          ENDDO
          main_objectM1%domain_map(ispin)%index1(idomain)=ientry
       ENDDO
       DEALLOCATE(domain_grid)

    IF (main_objectM1%nspins.eq.2) THEN
       CALL cp_dbcsr_copy(main_objectM1%quench_t(2),&
               main_objectM1%quench_t(1),error=error)
       main_objectM1%domain_map(2)%pairs(:,:)=&
          main_objectM1%domain_map(1)%pairs(:,:)
       main_objectM1%domain_map(2)%index1(:)=&
          main_objectM1%domain_map(1)%index1(:)
    ENDIF

    CALL cp_dbcsr_release(matrix_s_sym,error=error)

    IF (main_objectM1%domain_layout_mos==almo_domain_layout_molecular .OR. &
       main_objectM1%domain_layout_aos==almo_domain_layout_molecular) THEN
       DEALLOCATE(first_atom_of_molecule)
       DEALLOCATE(last_atom_of_molecule)
    ENDIF

    CALL timestop(handle)
  
  END SUBROUTINE almo_level2_spec1_0

END MODULE almo_scf_special

