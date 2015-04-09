!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Interface between ALMO SCF and QS
!> \par History
!>       2011.05 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
MODULE almo_scf_qs
  USE almo_scf_methods,                ONLY: almo_scf_p_blk_to_t_blk,&
                                             almo_scf_t_blk_to_p,&
                                             almo_scf_t_blk_to_t_blk_orthonormal
  USE almo_scf_types,                  ONLY: almo_mat_dim_aobasis,&
                                             almo_mat_dim_occ,&
                                             almo_mat_dim_virt,&
                                             almo_mat_dim_virt_disc,&
                                             almo_mat_dim_virt_full,&
                                             almo_scf_env_type
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cell_types,                      ONLY: cell_type,&
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
       cp_dbcsr_get_stored_coordinates, cp_dbcsr_init, &
       cp_dbcsr_nblkcols_total, cp_dbcsr_nblkrows_total, cp_dbcsr_p_type, &
       cp_dbcsr_release, cp_dbcsr_reserve_block2d, cp_dbcsr_row_block_sizes, &
       cp_dbcsr_set, cp_dbcsr_type, cp_dbcsr_work_create, &
       dbcsr_distribution_col_dist, dbcsr_distribution_mp, &
       dbcsr_distribution_new, dbcsr_distribution_obj, &
       dbcsr_distribution_row_dist, dbcsr_mp_group, dbcsr_mp_mynode, &
       dbcsr_mp_numnodes, dbcsr_type_no_symmetry
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_units,                        ONLY: cp_unit_to_cp2k
  USE input_constants,                 ONLY: &
       almo_constraint_ao_overlap, almo_constraint_block_diagonal, &
       almo_constraint_distance, almo_domain_layout_molecular, &
       almo_mat_distr_atomic, almo_mat_distr_molecular, atomic_guess, &
       do_bondparm_covalent, do_bondparm_vdw, molecular_guess
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_allgather
  USE molecule_types_new,              ONLY: get_molecule_set_info,&
                                             molecule_type
  USE mscfg_types,                     ONLY: get_matrix_from_submatrices,&
                                             molecular_scf_guess_env_type
  USE particle_types,                  ONLY: particle_type
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_initial_guess,                ONLY: calculate_atomic_block_dm,&
                                             calculate_mopac_dm
  USE qs_kind_types,                   ONLY: qs_kind_type
  USE qs_ks_methods,                   ONLY: qs_ks_update_qs_env
  USE qs_ks_types,                     ONLY: qs_ks_did_change,&
                                             qs_ks_env_type,&
                                             set_ks_env
  USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                             neighbor_list_iterate,&
                                             neighbor_list_iterator_create,&
                                             neighbor_list_iterator_p_type,&
                                             neighbor_list_iterator_release,&
                                             neighbor_list_set_p_type
  USE qs_rho_methods,                  ONLY: qs_rho_update_rho
  USE qs_rho_types,                    ONLY: qs_rho_get,&
                                             qs_rho_type
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'almo_scf_qs'

  PUBLIC :: almo_scf_init_qs, matrix_almo_create, matrix_qs_to_almo,&
            almo_scf_dm_to_ks,&
            almo_scf_construct_quencher, almo_scf_update_ks_energy

CONTAINS

! *****************************************************************************
!> \brief create the ALMO matrix templates
!> \param matrix_new ...
!> \param matrix_qs ...
!> \param almo_scf_env ...
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
  SUBROUTINE matrix_almo_create(matrix_new,matrix_qs,almo_scf_env,&
                                name_new,size_keys,symmetry_new,&
                                spin_key,init_domains,error)

    TYPE(cp_dbcsr_type)                      :: matrix_new, matrix_qs
    TYPE(almo_scf_env_type), INTENT(IN)      :: almo_scf_env
    CHARACTER(len=*), INTENT(IN)             :: name_new
    INTEGER, DIMENSION(2), INTENT(IN)        :: size_keys
    CHARACTER, INTENT(IN)                    :: symmetry_new
    INTEGER, INTENT(IN)                      :: spin_key
    LOGICAL, INTENT(IN)                      :: init_domains
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'matrix_almo_create', &
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

! dimension size: AO, MO, etc
!                 almo_mat_dim_aobasis - no. of AOs,
!                 almo_mat_dim_occ     - no. of occupied MOs
!                 almo_mat_dim_domains - no. of domains
! symmetry type: dbcsr_type_no_symmetry, dbcsr_type_symmetric,
!  dbcsr_type_antisymmetric, dbcsr_type_hermitian, dbcsr_type_antihermitian
!  (see dbcsr_lib/dbcsr_types.F for other values)
! spin_key: either 1 or 2 (0 is allowed for matrics in the AO basis)
!    TYPE(cp_dbcsr_iterator)                  :: iter
!    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: allones
!-----------------------------------------------------------------------

    CALL timeset(routineN,handle)

    ! RZK-warning The structure of the matrices can be optimized:
    ! 1. Diagonal matrices must be distributed evenly over the processes.
    !    This can be achived by distributing cpus: 012012-rows and 001122-cols
    !    block_diagonal_flag is introduced but not used
    ! 2. Multiplication of diagonally dominant matrices will be faster
    !    if the diagonal blocks are local to the same processes.
    ! 3. Systems of molecules of drastically different sizes might need
    !    better distribution.

    ! obtain distribution from the qs matrix - it might be useful
    ! to get the structure of the AO dimensions
    dist_qs = cp_dbcsr_distribution(matrix_qs)
    
    natoms=almo_scf_env%natoms
    nmols=almo_scf_env%nmolecules
    
    DO dimen=1,2 ! 1 - row, 2 - column dimension
    
     ! distribution pattern is the same for all matrix types (ao, occ, virt)
     IF (dimen==1) THEN !rows
      blk_distr => array_data(dbcsr_distribution_row_dist(dist_qs))
     ELSE !columns
      blk_distr => array_data(dbcsr_distribution_col_dist(dist_qs))
     ENDIF

     IF (size_keys(dimen)==almo_mat_dim_aobasis) THEN ! this dimension is AO
    
      ! structure of an AO dimension can be copied from matrix_qs
      blk_sizes => array_data(cp_dbcsr_row_block_sizes(matrix_qs))
    
      ! atomic clustering of AOs 
      IF (almo_scf_env%mat_distr_aos==almo_mat_distr_atomic) THEN 
       ALLOCATE(block_sizes_new(natoms),distr_new_array(natoms))
       block_sizes_new(:)=blk_sizes(:)
       distr_new_array(:)=blk_distr(:)
      ! molecular clustering of AOs
      ELSE IF (almo_scf_env%mat_distr_aos==almo_mat_distr_molecular) THEN
       ALLOCATE(block_sizes_new(nmols),distr_new_array(nmols))
       block_sizes_new(:)=0
       DO iatom=1,natoms
        block_sizes_new(almo_scf_env%domain_index_of_atom(iatom)) = &
          block_sizes_new(almo_scf_env%domain_index_of_atom(iatom)) + &
          blk_sizes(iatom)
       ENDDO
       DO imol=1,nmols
        distr_new_array(imol) = &
                  blk_distr(almo_scf_env%first_atom_of_domain(imol))
       ENDDO
      ELSE
       CPErrorMessage(cp_failure_level,routineP,"Illegal distribution",error)
       CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
      ENDIF
    
     ELSE ! this dimension is not AO
     
      IF (size_keys(dimen)==almo_mat_dim_occ .OR.&
          size_keys(dimen)==almo_mat_dim_virt .OR. &
          size_keys(dimen)==almo_mat_dim_virt_disc .OR. &
          size_keys(dimen)==almo_mat_dim_virt_full) THEN ! this dim is MO
    
       ! atomic clustering of MOs
       IF (almo_scf_env%mat_distr_mos==almo_mat_distr_atomic) THEN
        nlength=natoms
        ALLOCATE(block_sizes_new(nlength))
        IF (size_keys(dimen)==almo_mat_dim_occ) THEN
         ! currently distributing atomic distr of mos is not allowed
         ! RZK-warning define nocc_of_atom and nvirt_atom to implement it
         !block_sizes_new(:)=almo_scf_env%nocc_of_atom(:,spin_key)
        ELSE IF (size_keys(dimen)==almo_mat_dim_virt) THEN
         !block_sizes_new(:)=almo_scf_env%nvirt_of_atom(:,spin_key)
        ENDIF
       ! molecular clustering of MOs
       ELSE IF (almo_scf_env%mat_distr_mos==almo_mat_distr_molecular) THEN
        nlength=nmols
        ALLOCATE(block_sizes_new(nlength))
        IF (size_keys(dimen)==almo_mat_dim_occ) THEN
         block_sizes_new(:)=almo_scf_env%nocc_of_domain(:,spin_key)
        ELSE IF (size_keys(dimen)==almo_mat_dim_virt_disc) THEN
         block_sizes_new(:)=almo_scf_env%nvirt_disc_of_domain(:,spin_key)
        ELSE IF (size_keys(dimen)==almo_mat_dim_virt_full) THEN
         block_sizes_new(:)=almo_scf_env%nvirt_full_of_domain(:,spin_key)
        ELSE IF (size_keys(dimen)==almo_mat_dim_virt) THEN
         block_sizes_new(:)=almo_scf_env%nvirt_of_domain(:,spin_key)
        ENDIF
       ELSE
        CPErrorMessage(cp_failure_level,routineP,"Illegal distribution",error)
        CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
       ENDIF
      
      ELSE
 
       CPErrorMessage(cp_failure_level,routineP,"Illegal dimension",error)
       CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)

      ENDIF ! end choosing dim size (occ, virt)
      
      ! distribution for MOs is copied from AOs
      ALLOCATE(distr_new_array(nlength))
      ! atomic clustering  
      IF (almo_scf_env%mat_distr_mos==almo_mat_distr_atomic) THEN 
       distr_new_array(:)=blk_distr(:)
      ! molecular clustering
      ELSE IF (almo_scf_env%mat_distr_mos==almo_mat_distr_molecular) THEN
       DO imol=1,nmols
        distr_new_array(imol) = &
                  blk_distr(almo_scf_env%first_atom_of_domain(imol))
       ENDDO
      ENDIF
      ! old distribution - removed because it might fail for direct diag
      !allocate(distr_new_array(nlength))
      !if (dimen==1) then !rows
      ! nprocesses = dbcsr_mp_nprows(dbcsr_distribution_mp(dist_qs))
      ! distr_new_array=(/ (mod(i-1,nprocesses), i=1,nlength) /)
      !else !columns
      ! nprocesses = dbcsr_mp_npcols(dbcsr_distribution_mp(dist_qs))
      ! distr_new_array=(/ (mod(i-1,nprocesses), i=1,nlength) /)
      !endif
    
     ENDIF ! end choosing dimension size (AOs vs .NOT.AOs)
     
     ! create final arrays
     IF (dimen==1) THEN !rows
      CALL array_new(row_sizes_new, block_sizes_new, lb=1)
      CALL array_new(row_distr_new, distr_new_array, lb=1)
     ELSE !columns
      CALL array_new(col_sizes_new, block_sizes_new, lb=1)
      CALL array_new(col_distr_new, distr_new_array, lb=1)
     ENDIF
     DEALLOCATE(block_sizes_new,distr_new_array)
    
    ENDDO ! both rows and columns are done

    ! Create the distribution
    CALL dbcsr_distribution_new(dist_new,&
     dbcsr_distribution_mp(dist_qs),row_distr_new,col_distr_new)
    CALL array_release(col_distr_new)
    CALL array_release(row_distr_new)
    
    ! Create the matrix
    CALL cp_dbcsr_init(matrix_new,error=error)
    CALL cp_dbcsr_create (matrix_new, name_new,&
     dist_new, symmetry_new,&
     row_sizes_new, col_sizes_new, error=error)
    CALL array_release(col_sizes_new)
    CALL array_release(row_sizes_new)
    CALL cp_dbcsr_distribution_release(dist_new)

    ! fill out reqired blocks with 1.0_dp to tell the dbcsr library
    ! which blocks to keep
    IF (init_domains) THEN
  
       mynode = dbcsr_mp_mynode(dbcsr_distribution_mp(cp_dbcsr_distribution(matrix_new)))
       CALL cp_dbcsr_work_create(matrix_new, work_mutable=.TRUE., error=error)

       ! startQQQ - this part of the code scales quadratically
       ! therefore it is replaced with a less general but linear scaling algorithm below
       ! the quadratic algorithm is kept to be re-written later
       !QQQnblkrows_tot = cp_dbcsr_nblkrows_total(matrix_new)
       !QQQnblkcols_tot = cp_dbcsr_nblkcols_total(matrix_new)
       !QQQDO row = 1, nblkrows_tot
       !QQQ   DO col = 1, nblkcols_tot
       !QQQ      tr = .FALSE.
       !QQQ      iblock_row = row
       !QQQ      iblock_col = col
       !QQQ      CALL cp_dbcsr_get_stored_coordinates(matrix_new, iblock_row, iblock_col, tr, hold)
   
       !QQQ      IF(hold.EQ.mynode) THEN
       !QQQ      
       !QQQ         ! RZK-warning replace with a function which says if this 
       !QQQ         ! distribution block is active or not
       !QQQ         ! Translate indeces of distribution blocks to domain blocks
       !QQQ         if (size_keys(1)==almo_mat_dim_aobasis) then
       !QQQ           domain_row=almo_scf_env%domain_index_of_ao_block(iblock_row)
       !QQQ         else if (size_keys(2)==almo_mat_dim_occ .OR. &
       !QQQ                  size_keys(2)==almo_mat_dim_virt .OR. &
       !QQQ                  size_keys(2)==almo_mat_dim_virt_disc .OR. &
       !QQQ                  size_keys(2)==almo_mat_dim_virt_full) then 
       !QQQ           domain_row=almo_scf_env%domain_index_of_mo_block(iblock_row)
       !QQQ         else
       !QQQ           CPErrorMessage(cp_failure_level,routineP,"Illegal dimension",error)
       !QQQ           CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
       !QQQ         endif
   
       !QQQ         if (size_keys(2)==almo_mat_dim_aobasis) then
       !QQQ           domain_col=almo_scf_env%domain_index_of_ao_block(iblock_col)
       !QQQ         else if (size_keys(2)==almo_mat_dim_occ .OR. &
       !QQQ                  size_keys(2)==almo_mat_dim_virt .OR. &
       !QQQ                  size_keys(2)==almo_mat_dim_virt_disc .OR. &
       !QQQ                  size_keys(2)==almo_mat_dim_virt_full) then 
       !QQQ           domain_col=almo_scf_env%domain_index_of_mo_block(iblock_col)
       !QQQ         else
       !QQQ           CPErrorMessage(cp_failure_level,routineP,"Illegal dimension",error)
       !QQQ           CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
       !QQQ         endif

       !QQQ         ! Finds if we need this block
       !QQQ         ! only the block-diagonal constraint is implemented here
       !QQQ         active=.false.
       !QQQ         if (domain_row==domain_col) active=.true.

       !QQQ         IF (active) THEN
       !QQQ            NULLIFY (p_new_block)
       !QQQ            CALL cp_dbcsr_reserve_block2d(matrix_new, iblock_row, iblock_col, p_new_block)
       !QQQ            CPPostcondition(ASSOCIATED(p_new_block),cp_failure_level,routineP,error,failure)
       !QQQ            p_new_block(:,:) = 1.0_dp
       !QQQ         ENDIF

       !QQQ      ENDIF ! mynode
       !QQQ   ENDDO
       !QQQENDDO
       ! endQQQ - end of the quadratic part
       ! start linear-scaling replacement:
       ! works only for molecular blocks AND molecular distributions
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

          ENDIF ! mynode
       ENDDO
       ! end lnear-scaling replacement

    ENDIF ! init_domains

    CALL cp_dbcsr_finalize(matrix_new, error=error)

    CALL timestop (handle)

  END SUBROUTINE matrix_almo_create
  
! *****************************************************************************
!> \brief convert between two types of matrices: QS style to ALMO style
!> \param matrix_qs ...
!> \param matrix_almo ...
!> \param almo_scf_env ...
!> \param keep_sparsity ...
!> \param error ...
!> \par History
!>       2011.06 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE matrix_qs_to_almo(matrix_qs,matrix_almo,almo_scf_env,&
                               keep_sparsity,error)
    TYPE(cp_dbcsr_type)                      :: matrix_qs, matrix_almo
    TYPE(almo_scf_env_type)                  :: almo_scf_env
    LOGICAL, INTENT(IN)                      :: keep_sparsity
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'matrix_qs_to_almo', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure
    TYPE(cp_dbcsr_type)                      :: matrix_qs_nosym

    CALL timeset(routineN,handle)
    !RZK-warning if it's not a N(AO)xN(AO) matrix then stop

    SELECT CASE(almo_scf_env%mat_distr_aos)
    CASE(almo_mat_distr_atomic)
       ! automatic data_type conversion
       CALL cp_dbcsr_copy(matrix_almo,matrix_qs,&
            keep_sparsity=keep_sparsity,error=error)
    CASE(almo_mat_distr_molecular)
       ! desymmetrize the qs matrix
       CALL cp_dbcsr_init (matrix_qs_nosym, error=error)
       CALL cp_dbcsr_create (matrix_qs_nosym, template=matrix_qs,&
            matrix_type=dbcsr_type_no_symmetry, error=error)
       CALL cp_dbcsr_desymmetrize (matrix_qs, matrix_qs_nosym,&
            error=error)

       ! perform the magic complete_redistribute
       ! before calling complete_redistribute set all blocks to zero
       ! otherwise the non-zero elements of the redistributed matrix,
       ! which are in zero-blocks of the original matrix, will remain
       ! in the final redistributed matrix. this is a bug in
       ! complete_redistribute. RZK-warning it should be later corrected by calling
       ! dbcsr_set to 0.0 from within complete_redistribute
       CALL cp_dbcsr_set(matrix_almo, 0.0_dp, error=error)
       CALL cp_dbcsr_complete_redistribute(matrix_qs_nosym, matrix_almo,&
               keep_sparsity=keep_sparsity,error=error);
       CALL cp_dbcsr_release (matrix_qs_nosym, error=error)

    CASE DEFAULT 
        CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
    END SELECT

    CALL timestop(handle)

  END SUBROUTINE matrix_qs_to_almo

! *****************************************************************************
!> \brief convert between two types of matrices: ALMO style to QS style
!> \param matrix_almo ...
!> \param matrix_qs ...
!> \param almo_scf_env ...
!> \param error ...
!> \par History
!>       2011.06 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE matrix_almo_to_qs(matrix_almo,matrix_qs,almo_scf_env,error)
    TYPE(cp_dbcsr_type)                      :: matrix_almo, matrix_qs
    TYPE(almo_scf_env_type), INTENT(IN)      :: almo_scf_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'matrix_almo_to_qs', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure

    CALL timeset(routineN,handle)
    ! RZK-warning if it's not a N(AO)xN(AO) matrix then stop

!    IF (ls_mstruct%single_precision) THEN
!       CALL cp_dbcsr_init (matrix_tmp, error=error)
!       CALL cp_dbcsr_create (matrix_tmp, template=matrix_ls,&
!            data_type=dbcsr_type_real_8, error=error)
!       CALL cp_dbcsr_copy (matrix_tmp, matrix_ls, error=error)
!    ENDIF

    SELECT CASE(almo_scf_env%mat_distr_aos)
    CASE(almo_mat_distr_atomic)
!       IF (ls_mstruct%single_precision) THEN
!          CALL cp_dbcsr_copy_into_existing (matrix_qs, matrix_tmp, error=error)
!       ELSE
          CALL cp_dbcsr_copy_into_existing (matrix_qs, matrix_almo, error=error)
!       ENDIF
    CASE(almo_mat_distr_molecular)
       CALL cp_dbcsr_set(matrix_qs,0.0_dp,error=error)
!       IF (ls_mstruct%single_precision) THEN
!          CALL cp_dbcsr_complete_redistribute(matrix_tmp, matrix_qs, keep_sparsity=.TRUE., error=error) 
!       ELSE
          CALL cp_dbcsr_complete_redistribute(matrix_almo, matrix_qs, keep_sparsity=.TRUE., error=error) 
!       ENDIF
    CASE DEFAULT 
       CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
    END SELECT

!    IF (ls_mstruct%single_precision) THEN
!       CALL cp_dbcsr_release(matrix_tmp,error=error)
!    ENDIF

    CALL timestop(handle)

  END SUBROUTINE matrix_almo_to_qs

! *****************************************************************************
!> \brief Initialization of QS and ALMOs
!>        Some parts can be factored-out since they are common 
!>        for the other SCF methods
!> \param qs_env ...
!> \param almo_scf_env ...
!> \param error ...
!> \par History
!>       2011.05 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE almo_scf_init_qs(qs_env,almo_scf_env,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(almo_scf_env_type)                  :: almo_scf_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'almo_scf_init_qs', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, nspin, unit_nr
    INTEGER, DIMENSION(2)                    :: nelectron_spin
    LOGICAL                                  :: failure, has_unit_metric
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_s, rho_ao
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(molecular_scf_guess_env_type), &
      POINTER                                :: mscfg_env
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(section_vals_type), POINTER         :: input

    CALL timeset(routineN,handle)
    NULLIFY(rho, sab_orb, rho_ao)
    ! get a useful output_unit
    logger => cp_error_get_logger(error)
    IF (logger%para_env%mepos==logger%para_env%source) THEN
       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
    ELSE
       unit_nr=-1
    ENDIF

    ! get basic quantities from the qs_env
    CALL get_qs_env(qs_env,&
                    dft_control=dft_control,&
                    matrix_s=matrix_s,&
                    matrix_ks=matrix_ks,&
                    ks_env=ks_env,&
                    input=input,&
                    energy=energy,&
                    atomic_kind_set=atomic_kind_set,&
                    qs_kind_set=qs_kind_set,&
                    particle_set=particle_set,&
                    has_unit_metric=has_unit_metric,&
                    para_env=para_env,&
                    nelectron_spin=nelectron_spin,&
                    mscfg_env=mscfg_env,&
                    rho=rho,&
                    sab_orb=sab_orb,&
                    error=error)

    CALL qs_rho_get(rho, rho_ao=rho_ao, error=error)
    CPPostcondition(ASSOCIATED(mscfg_env),cp_failure_level,routineP,error,failure)

    nspin=dft_control%nspins

    ! create matrix_ks if necessary
    IF (.NOT.ASSOCIATED(matrix_ks)) THEN
       CALL cp_dbcsr_allocate_matrix_set(matrix_ks,nspin,error)
       DO ispin=1,nspin
          ALLOCATE(matrix_ks(ispin)%matrix)
          CALL cp_dbcsr_init(matrix_ks(ispin)%matrix,error=error)
          CALL cp_dbcsr_create(matrix_ks(ispin)%matrix,&
               template=matrix_s(1)%matrix,error=error)
          CALL cp_dbcsr_alloc_block_from_nbl(matrix_ks(ispin)%matrix,sab_orb,error=error)
          CALL cp_dbcsr_set(matrix_ks(ispin)%matrix,0.0_dp,error=error)
       ENDDO
       CALL set_ks_env(ks_env,matrix_ks=matrix_ks,error=error)
    ENDIF

    ! create an initial guess
    SELECT CASE (almo_scf_env%almo_scf_guess)
    CASE(molecular_guess)
       
       DO ispin=1,nspin
         
          ! the calculations on "isolated" molecules has already been done
          ! all we need to do is convert the MOs of molecules into
          ! the ALMO matrix taking into account different distributions 
          CALL get_matrix_from_submatrices(mscfg_env,&
                  almo_scf_env%matrix_t_blk(ispin), ispin, error)
          CALL cp_dbcsr_filter(almo_scf_env%matrix_t_blk(ispin),&
                  almo_scf_env%eps_filter,error=error)

       ENDDO
       
    CASE (atomic_guess)

       DO ispin=1,nspin
          CALL cp_dbcsr_set(rho_ao(ispin)%matrix,0.0_dp,error=error)
          IF (dft_control%qs_control%dftb .OR. dft_control%qs_control%semi_empirical) THEN
             CALL calculate_mopac_dm(rho_ao(ispin)%matrix,&
                  matrix_s(1)%matrix, has_unit_metric, &
                  dft_control, particle_set, atomic_kind_set, qs_kind_set,&
                  nspin, nelectron_spin(ispin),&
                  para_env, error)
          ELSE
             CALL calculate_atomic_block_dm(rho_ao(ispin)%matrix,&
                  matrix_s(1)%matrix, particle_set, atomic_kind_set, qs_kind_set, &
                  ispin, nspin, nelectron_spin(ispin), unit_nr, error)
          ENDIF

          ! copy the atomic-block dm into matrix_p_blk
          CALL matrix_qs_to_almo(rho_ao(ispin)%matrix,&
                                 almo_scf_env%matrix_p_blk(ispin),almo_scf_env,&
                                 .FALSE.,error)
          CALL cp_dbcsr_filter(almo_scf_env%matrix_p_blk(ispin),&
                  almo_scf_env%eps_filter,error=error)

       ENDDO
       
       ! obtain orbitals from density matrix (currently ALMO SCF needs orbitals)
       ! RZK-warning currently the orbitals are generated
       ! by diagonalizing p blocks and taking only the
       ! eigenvectors with highest occupation
       ! essentially producing the ionic orbitals
       ! while it works for simple systems this guess
       ! will lead to convergence problems with complex molecules
       CALL almo_scf_p_blk_to_t_blk(almo_scf_env,error)
    
    END SELECT

    CALL almo_scf_t_blk_to_t_blk_orthonormal(almo_scf_env,error)
    CALL almo_scf_t_blk_to_p(almo_scf_env,&
            use_sigma_inv_guess=.FALSE.,error=error)

    DO ispin=1,nspin
       CALL matrix_almo_to_qs(almo_scf_env%matrix_p(ispin),&
                              rho_ao(ispin)%matrix,&
                              almo_scf_env,error)
    ENDDO

    CALL qs_rho_update_rho(rho,qs_env=qs_env, error=error)
    CALL qs_ks_did_change(qs_env%ks_env,rho_changed=.TRUE., error=error)
    CALL qs_ks_update_qs_env(qs_env, calculate_forces=.FALSE.,&
                             just_energy=.FALSE.,error=error)

    IF (unit_nr>0) THEN
       IF (almo_scf_env%almo_scf_guess.eq.molecular_guess) THEN
          WRITE(unit_nr,'(T2,A38,F40.10)') "Single-molecule energy:", &
            SUM(mscfg_env%energy_of_frag)
       ENDIF
       WRITE(unit_nr,'(T2,A38,F40.10)') "Energy of the initial guess:",energy%total
       WRITE(unit_nr,'()') 
    ENDIF

    CALL timestop(handle)

  END SUBROUTINE almo_scf_init_qs

! *****************************************************************************
!> \brief use the density matrix in almo_scf_env
!>        to compute the new energy and KS matrix
!> \param qs_env ...
!> \param almo_scf_env ...
!> \param energy_new ...
!> \param error ...
!> \par History
!>       2011.05 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE almo_scf_dm_to_ks(qs_env,almo_scf_env,energy_new,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(almo_scf_env_type)                  :: almo_scf_env
    REAL(KIND=dp)                            :: energy_new
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'almo_scf_dm_to_ks', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, nspin
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: rho_ao
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_rho_type), POINTER               :: rho

    NULLIFY(rho, rho_ao)
    CALL timeset(routineN,handle)

    nspin=almo_scf_env%nspins
    CALL get_qs_env(qs_env, energy=energy, rho=rho, error=error)
    CALL qs_rho_get(rho, rho_ao=rho_ao, error=error)

    ! set the new density matrix 
    DO ispin=1,nspin
       CALL matrix_almo_to_qs(almo_scf_env%matrix_p(ispin),&
                              rho_ao(ispin)%matrix,&
                              almo_scf_env,&
                              error)
    END DO

    ! compute the corresponding KS matrix and new energy
    CALL qs_rho_update_rho(rho,qs_env=qs_env, error=error)
    CALL qs_ks_did_change(qs_env%ks_env,rho_changed=.TRUE., error=error)
    CALL qs_ks_update_qs_env(qs_env, calculate_forces=.FALSE., just_energy=.FALSE.,&
                             print_active=.TRUE., error=error)
    energy_new=energy%total

    CALL timestop(handle)

  END SUBROUTINE almo_scf_dm_to_ks

! *****************************************************************************
!> \brief update qs_env total energy
!> \param qs_env ...
!> \param energy ...
!> \param error ...
!> \par History
!>       2013.03 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE almo_scf_update_ks_energy(qs_env,energy,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    REAL(KIND=dp)                            :: energy
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'almo_scf_update_ks_energy', &
      routineP = moduleN//':'//routineN

    TYPE(qs_energy_type), POINTER            :: qs_energy

    CALL get_qs_env(qs_env, energy=qs_energy, error=error)
    qs_energy%total=energy

  END SUBROUTINE almo_scf_update_ks_energy

! *****************************************************************************
!> \brief Creates the matrix that imposes absolute locality on MOs 
!> \param qs_env ...
!> \param almo_scf_env ...
!> \param error ...
!> \par History
!>       2011.11 created [Rustam Z. Khaliullin]
!> \author Rustam Z. Khaliullin
! *****************************************************************************
  SUBROUTINE almo_scf_construct_quencher(qs_env,almo_scf_env,error)
    
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(almo_scf_env_type), INTENT(INOUT)   :: almo_scf_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'almo_scf_construct_quencher', &
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
                                                delayed_increment, failure, &
                                                found, max_neig_fails, tr
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
       
    ! get a useful output_unit
    logger => cp_error_get_logger(error)
    IF (logger%para_env%mepos==logger%para_env%source) THEN
       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
    ELSE
       unit_nr=-1
    ENDIF

    ndomains=almo_scf_env%ndomains

    CALL get_qs_env(qs_env=qs_env,&
      particle_set=particle_set,&
      molecule_set=molecule_set,&
      cell=cell,&
      matrix_s=matrix_s,&
      sab_almo=sab_almo,&
      error=error)

    ! if we are dealing with molecules get info about them
    IF (almo_scf_env%domain_layout_mos==almo_domain_layout_molecular .OR. &
       almo_scf_env%domain_layout_aos==almo_domain_layout_molecular) THEN
       ALLOCATE(first_atom_of_molecule(almo_scf_env%nmolecules))
       ALLOCATE(last_atom_of_molecule(almo_scf_env%nmolecules))
       CALL get_molecule_set_info(molecule_set,&
               mol_to_first_atom=first_atom_of_molecule,&
               mol_to_last_atom=last_atom_of_molecule,&
               error=error)
    ENDIF

    ! create a symmetrized copy of the ao overlap
    CALL cp_dbcsr_init(matrix_s_sym,error=error)
    CALL cp_dbcsr_create(matrix_s_sym,&
            template=almo_scf_env%matrix_s(1),&
            matrix_type=dbcsr_type_no_symmetry,error=error)
    CALL cp_dbcsr_get_info(almo_scf_env%matrix_s(1),&
            matrix_type=sym)
    IF (sym.eq.dbcsr_type_no_symmetry) THEN
       CALL cp_dbcsr_copy(matrix_s_sym,almo_scf_env%matrix_s(1),error=error)
    ELSE
       CALL cp_dbcsr_desymmetrize(almo_scf_env%matrix_s(1),&
               matrix_s_sym,error=error)
    ENDIF

    ALLOCATE(almo_scf_env%quench_t(almo_scf_env%nspins))
    ALLOCATE(almo_scf_env%domain_map(almo_scf_env%nspins))

    !DO ispin=1,almo_scf_env%nspins
    ispin=1

       ! create the sparsity template for the occupied orbitals
       CALL matrix_almo_create(matrix_new=almo_scf_env%quench_t(ispin),&
            matrix_qs=matrix_s(1)%matrix,&
            almo_scf_env=almo_scf_env,&
            name_new="T_QUENCHER",&
            size_keys=(/almo_mat_dim_aobasis,almo_mat_dim_occ/),&
            symmetry_new=dbcsr_type_no_symmetry,&
            spin_key=ispin,&
            init_domains=.FALSE.,&
            error=error)
    
       ! initialize distance quencher
       CALL cp_dbcsr_work_create(almo_scf_env%quench_t(ispin),&
               work_mutable=.TRUE., error=error)

       nblkrows_tot = cp_dbcsr_nblkrows_total(almo_scf_env%quench_t(ispin))
       nblkcols_tot = cp_dbcsr_nblkcols_total(almo_scf_env%quench_t(ispin))

       mynode = dbcsr_mp_mynode(dbcsr_distribution_mp(&
          cp_dbcsr_distribution(almo_scf_env%quench_t(ispin))))
       nNodes = dbcsr_mp_numnodes(dbcsr_distribution_mp(&
          cp_dbcsr_distribution(almo_scf_env%quench_t(ispin))))
       GroupID = dbcsr_mp_group(dbcsr_distribution_mp(&
          cp_dbcsr_distribution(almo_scf_env%quench_t(ispin))))

       ! create global atom neighbor list from the local lists
       ! first, calculate number of local pairs
       local_list_length=0
       CALL neighbor_list_iterator_create(nl_iterator,sab_almo)
       DO WHILE (neighbor_list_iterate(nl_iterator)==0)
          ! nnode - total number of neighbors for iatom
          ! inode - current neighbor count
          CALL get_iterator_info(nl_iterator,&
                  iatom=iatom2,jatom=jatom2,inode=inode2,nnode=nnode2)
          !WRITE(*,*) "GET INFO: ",iatom2, jatom2, inode2, nnode2
          IF (inode2==1) THEN
             local_list_length=local_list_length+nnode2
          END IF
       END DO
       CALL neighbor_list_iterator_release(nl_iterator)

       ! second, extract the local list to an array
       ALLOCATE(local_list(2*local_list_length))
       local_list(:)=0
       local_list_length=0
       CALL neighbor_list_iterator_create(nl_iterator2,sab_almo)
       DO WHILE (neighbor_list_iterate(nl_iterator2)==0)
          CALL get_iterator_info(nl_iterator2,&
                  iatom=iatom2,jatom=jatom2)
!WRITE(*,*) "GET INFO: ",iatom2, jatom2
          local_list(2*local_list_length+1)=iatom2
          local_list(2*local_list_length+2)=jatom2
          local_list_length=local_list_length+1
       ENDDO ! end loop over pairs of atoms
       CALL neighbor_list_iterator_release(nl_iterator2)

       ! third, communicate local length to the other nodes
       ALLOCATE(list_length_cpu(nNodes),list_offset_cpu(nNodes))
       CALL mp_allgather(2*local_list_length,list_length_cpu,GroupID)
       
       ! fourth, create a global list
       list_offset_cpu(1)=0
       DO iNode=2,nNodes
          list_offset_cpu(iNode)=list_offset_cpu(iNode-1) + &
             list_length_cpu(iNode-1)
       ENDDO
       global_list_length=list_offset_cpu(nNodes)+list_length_cpu(nNodes)
 
       ! fifth, communicate all list data
       ALLOCATE(global_list(global_list_length))
!WRITE(*,*) "LENGTH: ", list_length_cpu
!WRITE(*,*) "OFFSET: ", list_offset_cpu
!WRITE(*,*) "LOCAL:  ", local_list
       CALL mp_allgather(local_list,global_list,&
               list_length_cpu,list_offset_cpu,GroupID)
       DEALLOCATE(list_length_cpu,list_offset_cpu)
       DEALLOCATE(local_list)
       
       ! calculate maximum number of atoms surrounding the domain
       ALLOCATE(current_number_neighbors(almo_scf_env%ndomains))
       current_number_neighbors(:)=0
       global_list_length=global_list_length/2
       DO ipair=1, global_list_length
          iatom2=global_list(2*(ipair-1)+1)
          jatom2=global_list(2*(ipair-1)+2)
          idomain2=almo_scf_env%domain_index_of_atom(iatom2)
          jdomain2=almo_scf_env%domain_index_of_atom(jatom2)
          ! add to the list
          current_number_neighbors(idomain2)=current_number_neighbors(idomain2)+1
          ! add j,i with i,j
          IF (idomain2.ne.jdomain2) THEN
             current_number_neighbors(jdomain2)=current_number_neighbors(jdomain2)+1
          ENDIF
       ENDDO
       max_domain_neighbors=MAXVAL(current_number_neighbors)

       ! use the global atom neighbor list to create a global domain neighbor list
       ALLOCATE(domain_neighbor_list_excessive(ndomains,max_domain_neighbors))
       current_number_neighbors(:)=1
       DO ipair=1, ndomains
          domain_neighbor_list_excessive(ipair,1)=ipair
       ENDDO
       DO ipair=1, global_list_length
          iatom2=global_list(2*(ipair-1)+1)
          jatom2=global_list(2*(ipair-1)+2)
          idomain2=almo_scf_env%domain_index_of_atom(iatom2)
          jdomain2=almo_scf_env%domain_index_of_atom(jatom2)
          already_listed=.FALSE.
          DO ineighbor=1,current_number_neighbors(idomain2)
             IF (domain_neighbor_list_excessive(idomain2,ineighbor).eq.jdomain2) THEN
                already_listed=.TRUE.
                EXIT
             ENDIF
          ENDDO
          IF (.NOT.already_listed) THEN
             ! add to the list
             current_number_neighbors(idomain2)=current_number_neighbors(idomain2)+1
             domain_neighbor_list_excessive(idomain2,current_number_neighbors(idomain2))=jdomain2
             ! add j,i with i,j
             IF (idomain2.ne.jdomain2) THEN
                current_number_neighbors(jdomain2)=current_number_neighbors(jdomain2)+1
                domain_neighbor_list_excessive(jdomain2,current_number_neighbors(jdomain2))=idomain2
             ENDIF
          ENDIF
       ENDDO ! end loop over pairs of atoms
       DEALLOCATE(global_list)

       max_domain_neighbors=MAXVAL(current_number_neighbors)
       ALLOCATE(domain_neighbor_list(ndomains,max_domain_neighbors))
       domain_neighbor_list(:,:)=0
       domain_neighbor_list(:,:)=domain_neighbor_list_excessive(:,1:max_domain_neighbors)
       DEALLOCATE(domain_neighbor_list_excessive)

       ALLOCATE(almo_scf_env%domain_map(ispin)%index1(ndomains))
       ALLOCATE(almo_scf_env%domain_map(ispin)%pairs(max_domain_neighbors*ndomains,2))
       almo_scf_env%domain_map(ispin)%pairs(:,:)=0
       almo_scf_env%domain_map(ispin)%index1(:)=0
       domain_map_local_entries=0
   
       
       ! RZK-warning intermediate [0,1] quencher values are ill-defined
       ! for molecules (not continuous and conceptually inadequate)
       
       ! O(N) loop over domain pairs
       DO row = 1, nblkrows_tot
          DO col = 1, current_number_neighbors(row)
          !DO col = 1, nblkcols_tot
             tr = .FALSE.
             iblock_row = row
             iblock_col = domain_neighbor_list(row,col)
             !iblock_col = col
!IF (unit_nr>0) WRITE(*,*) iblock_row, iblock_col
             CALL cp_dbcsr_get_stored_coordinates(almo_scf_env%quench_t(ispin),&
                     iblock_row, iblock_col, tr, hold)

             IF(hold.EQ.mynode) THEN

                ! Translate indices of distribution blocks to indices of domain blocks
                ! Rows are AOs
                domain_row=almo_scf_env%domain_index_of_ao_block(iblock_row)
                ! Columns are electrons (i.e. MOs)
                domain_col=almo_scf_env%domain_index_of_mo_block(iblock_col)
                
                SELECT CASE (almo_scf_env%constraint_type)
                CASE (almo_constraint_block_diagonal)
               
                   block_active=.FALSE.
                   ! type of electron groups
                   IF (almo_scf_env%domain_layout_mos==almo_domain_layout_molecular) THEN
   
                      ! type of ao domains
                      IF (almo_scf_env%domain_layout_aos==almo_domain_layout_molecular) THEN
   
                         ! ao domains are molecular / electron groups are molecular
                         IF (domain_row==domain_col) THEN
                            block_active=.TRUE.
                         ENDIF
   
                      ELSE ! ao domains are atomic
   
                         ! ao domains are atomic / electron groups are molecular
                         CPErrorMessage(cp_failure_level,routineP,"Illegal: atomic domains and molecular groups",error)
                         CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
   
                      ENDIF
                   
                   ELSE ! electron groups are atomic
                      
                      ! type of ao domains
                      IF (almo_scf_env%domain_layout_aos==almo_domain_layout_molecular) THEN
   
                         ! ao domains are molecular / electron groups are atomic
                         CPErrorMessage(cp_failure_level,routineP,"Illegal: molecular domains and atomic groups",error)
                         CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
   
                      ELSE 
                         
                         ! ao domains are atomic / electron groups are atomic
                         IF (domain_row==domain_col) THEN
                            block_active=.TRUE.
                         ENDIF
                      
                      ENDIF
                   
                   ENDIF ! end type of electron groups
                   
                   IF ( block_active ) THEN
                      
                      NULLIFY (p_new_block)
                      CALL cp_dbcsr_reserve_block2d(almo_scf_env%quench_t(ispin),&
                              iblock_row, iblock_col, p_new_block)
                      CPPostcondition(ASSOCIATED(p_new_block),cp_failure_level,routineP,error,failure)
                      p_new_block(:,:) = 1.0_dp
                      
                      IF (domain_map_local_entries.ge.max_domain_neighbors*almo_scf_env%ndomains) THEN
                         CPErrorMessage(cp_failure_level,routineP,"weird... max_domain_neighbors is exceeded",error)
                         CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
                      ENDIF
                      almo_scf_env%domain_map(ispin)%pairs(domain_map_local_entries+1,1)=iblock_row
                      almo_scf_env%domain_map(ispin)%pairs(domain_map_local_entries+1,2)=iblock_col
                      domain_map_local_entries=domain_map_local_entries+1
                   
                   ENDIF

                CASE (almo_constraint_ao_overlap)

                   ! type of electron groups
                   IF (almo_scf_env%domain_layout_mos==almo_domain_layout_molecular) THEN
   
                      ! type of ao domains
                      IF (almo_scf_env%domain_layout_aos==almo_domain_layout_molecular) THEN
   
                         ! ao domains are molecular / electron groups are molecular
                         
                         ! compute the maximum overlap between the atoms of the two molecules 
                         CALL cp_dbcsr_get_block_p(matrix_s_sym,&
                                 iblock_row, iblock_col, p_new_block, found)
                         IF (found) THEN
                         !   CPErrorMessage(cp_failure_level,routineP,"S block not found",error)
                         !   CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
                            overlap=MAXVAL(ABS(p_new_block))
                         ELSE
                            overlap=0.0_dp
                         ENDIF
   
                      ELSE ! ao domains are atomic
   
                         ! ao domains are atomic / electron groups are molecular
                         ! overlap_between_atom_and_molecule(atom=domain_row,molecule=domain_col)
                         CPErrorMessage(cp_failure_level,routineP,"atomic domains and molecular groups - NYI",error)
                         CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
   
                      ENDIF
                   
                   ELSE ! electron groups are atomic
                      
                      ! type of ao domains
                      IF (almo_scf_env%domain_layout_aos==almo_domain_layout_molecular) THEN
   
                         ! ao domains are molecular / electron groups are atomic
                         ! overlap_between_atom_and_molecule(atom=domain_col,molecule=domain_row)
                         CPErrorMessage(cp_failure_level,routineP,"molecular domains and atomic groups - NYI",error)
                         CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
   
                      ELSE 
                         
                         ! ao domains are atomic / electron groups are atomic
                         ! compute max overlap between atoms: domain_row and domain_col
                         CALL cp_dbcsr_get_block_p(matrix_s_sym,&
                                 iblock_row, iblock_col, p_new_block, found)
                         IF (found) THEN
                         !   CPErrorMessage(cp_failure_level,routineP,"S block not found",error)
                         !   CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
                            overlap=MAXVAL(ABS(p_new_block))
                         ELSE
                            overlap=0.0_dp
                         ENDIF
                      
                      ENDIF
                   
                   ENDIF ! end type of electron groups
                   
                   s0=-LOG10(ABS(almo_scf_env%quencher_s0))
                   s1=-LOG10(ABS(almo_scf_env%quencher_s1))
                   IF (overlap.eq.0.0_dp) THEN
                      overlap=-LOG10(ABS(almo_scf_env%eps_filter))+100.0_dp
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

                   ! Fill in non-zero blocks if AOs are close to the electron center
                   IF ( overlap.lt.s1 ) THEN
                      NULLIFY (p_new_block)
                      CALL cp_dbcsr_reserve_block2d(almo_scf_env%quench_t(ispin),&
                              iblock_row, iblock_col, p_new_block)
                      CPPostcondition(ASSOCIATED(p_new_block),cp_failure_level,routineP,error,failure)
                      IF ( overlap.le.s0 ) THEN
                         p_new_block(:,:) = 1.0_dp
                         !WRITE(*,'(A15,2I7,3F8.3,E11.3)') "INTRA-BLOCKS: ",&
                         !   iblock_col, iblock_row, s0, s1, overlap, p_new_block(1,1)
                      ELSE
                         p_new_block(:,:) = 1.0_dp/(1.0_dp+EXP(-(s0-s1)/(s0-overlap)-(s0-s1)/(overlap-s1)))
                         !WRITE(*,'(A15,2I7,3F8.3,E11.3)') "INTER-BLOCKS: ",&
                         !   iblock_col, iblock_row, s0, s1, overlap, p_new_block(1,1)
                      ENDIF
                      
                      IF (ABS(p_new_block(1,1)).gt.ABS(almo_scf_env%eps_filter)) THEN
                         IF (domain_map_local_entries.ge.max_domain_neighbors*almo_scf_env%ndomains) THEN
                            CPErrorMessage(cp_failure_level,routineP,"weird... max_domain_neighbors is exceeded",error)
                            CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
                         ENDIF
                         almo_scf_env%domain_map(ispin)%pairs(domain_map_local_entries+1,1)=iblock_row
                         almo_scf_env%domain_map(ispin)%pairs(domain_map_local_entries+1,2)=iblock_col
                         domain_map_local_entries=domain_map_local_entries+1
                      ENDIF
                   
                   ENDIF

                CASE (almo_constraint_distance)

                   ! type of electron groups
                   IF (almo_scf_env%domain_layout_mos==almo_domain_layout_molecular) THEN
   
                      ! type of ao domains
                      IF (almo_scf_env%domain_layout_aos==almo_domain_layout_molecular) THEN
   
                         ! ao domains are molecular / electron groups are molecular
                         
                         ! compute distance between molecules: domain_row and domain_col
                         ! distance between molecules is defined as the smallest 
                         ! distance among all atom pairs
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
                               ENDDO ! jatom
                            ENDDO ! iatom
                            CPPrecondition(contact_atom_1.gt.0,cp_failure_level,routineP,error,failure)
                            distance=SQRT(distance_squared)
                         ENDIF
   
                      ELSE ! ao domains are atomic
   
                         ! ao domains are atomic / electron groups are molecular
                         !distance_between_atom_and_molecule(atom=domain_row,molecule=domain_col)
                         CPErrorMessage(cp_failure_level,routineP,"atomic domains and molecular groups - NYI",error)
                         CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
   
                      ENDIF
                   
                   ELSE ! electron groups are atomic
                      
                      ! type of ao domains
                      IF (almo_scf_env%domain_layout_aos==almo_domain_layout_molecular) THEN
   
                         ! ao domains are molecular / electron groups are atomic
                         !distance_between_atom_and_molecule(atom=domain_col,molecule=domain_row)
                         CPErrorMessage(cp_failure_level,routineP,"molecular domains and atomic groups - NYI",error)
                         CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
   
                      ELSE 
                         
                         ! ao domains are atomic / electron groups are atomic
                         ! compute distance between atoms: domain_row and domain_col
                         rab(:) = pbc(particle_set(domain_row)%r(:),particle_set(domain_col)%r(:),cell)
                         distance = SQRT(rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3))
                         contact_atom_1=domain_row
                         contact_atom_2=domain_col
                      
                      ENDIF
                   
                   ENDIF ! end type of electron groups
                   
                   ! get atomic radii to compute distance cutoff threshold
                   IF(almo_scf_env%quencher_radius_type==do_bondparm_covalent) THEN
                      CALL get_atomic_kind(atomic_kind=particle_set(contact_atom_1)%atomic_kind,&
                              rcov=contact1_radius)
                      CALL get_atomic_kind(atomic_kind=particle_set(contact_atom_2)%atomic_kind,&
                              rcov=contact2_radius)
                   ELSE IF(almo_scf_env%quencher_radius_type==do_bondparm_vdw) THEN
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

!RZK-warning the procedure is faulty for molecules: the closest contacts should be found using
! the element specific radii

                   ! compute inner and outer cutoff radii
                   r0=almo_scf_env%quencher_r0_factor*(contact1_radius+contact2_radius)
                      !+almo_scf_env%quencher_r0_shift
                   r1=almo_scf_env%quencher_r1_factor*(contact1_radius+contact2_radius)
                      !+almo_scf_env%quencher_r1_shift
                   
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
                   
                   ! Fill in non-zero blocks if AOs are close to the electron center
                   IF ( distance.lt.r1 ) THEN
                      NULLIFY (p_new_block)
                      CALL cp_dbcsr_reserve_block2d(almo_scf_env%quench_t(ispin),&
                              iblock_row, iblock_col, p_new_block)
                      CPPostcondition(ASSOCIATED(p_new_block),cp_failure_level,routineP,error,failure)
                      IF ( distance.le.r0 ) THEN
                         p_new_block(:,:) = 1.0_dp
                         !WRITE(*,'(A15,2I7,5F8.3,E11.3)') "INTRA-BLOCKS: ",&
                         !   iblock_col, iblock_row, contact1_radius,&
                         !   contact2_radius, r0, r1, distance, p_new_block(1,1)
                      ELSE
! remove the intermediate values from the quencher temporarily
CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
                         p_new_block(:,:) = 1.0_dp/(1.0_dp+EXP((r1-r0)/(r0-distance)+(r1-r0)/(r1-distance)))
                         !WRITE(*,'(A15,2I7,5F8.3,E11.3)') "INTER-BLOCKS: ",&
                         !   iblock_col, iblock_row, contact1_radius,&
                         !   contact2_radius, r0, r1, distance, p_new_block(1,1)
                      ENDIF
                      
                      IF (ABS(p_new_block(1,1)).gt.ABS(almo_scf_env%eps_filter)) THEN
                         IF (domain_map_local_entries.ge.max_domain_neighbors*almo_scf_env%ndomains) THEN
                            CPErrorMessage(cp_failure_level,routineP,"weird... max_domain_neighbors is exceeded",error)
                            CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
                         ENDIF
                         almo_scf_env%domain_map(ispin)%pairs(domain_map_local_entries+1,1)=iblock_row
                         almo_scf_env%domain_map(ispin)%pairs(domain_map_local_entries+1,2)=iblock_col
                         domain_map_local_entries=domain_map_local_entries+1
                      ENDIF

                   ENDIF

                CASE DEFAULT
                   CALL stop_program(routineN,moduleN,__LINE__,&
                         "Illegal constraint type")
                END SELECT

             ENDIF ! mynode

          ENDDO
       ENDDO ! end O(N) loop over pairs
       
       DEALLOCATE(domain_neighbor_list)
       DEALLOCATE(current_number_neighbors)

       CALL cp_dbcsr_finalize(almo_scf_env%quench_t(ispin),error=error)
       !CALL cp_dbcsr_scale(almo_scf_env%quench_t(ispin),&
       !        almo_scf_env%envelope_amplitude,error=error)
       CALL cp_dbcsr_filter(almo_scf_env%quench_t(ispin),&
               almo_scf_env%eps_filter,error=error)

       ! check that both domain_map and quench_t have the same number of entries 
       nblks=cp_dbcsr_get_num_blocks(almo_scf_env%quench_t(ispin))
       IF (nblks.ne.domain_map_local_entries) THEN
          CPErrorMessage(cp_failure_level,routineP,"number of blocks is wrong",error)
          CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
       ENDIF

       ! communicate local parts of the domain map
       !nNodes = dbcsr_mp_numnodes(dbcsr_distribution_mp(&
       !   cp_dbcsr_distribution(almo_scf_env%quench_t(ispin))))
       !GroupID = dbcsr_mp_group(dbcsr_distribution_mp(&
       !   cp_dbcsr_distribution(almo_scf_env%quench_t(ispin))))
   
       ! first, communicate map sizes on the other nodes
       ALLOCATE(domain_entries_cpu(nNodes),offset_for_cpu(nNodes))
       CALL mp_allgather(2*domain_map_local_entries,domain_entries_cpu,GroupID)

       ! second, create 
       offset_for_cpu(1)=0
       DO iNode=2,nNodes
          offset_for_cpu(iNode)=offset_for_cpu(iNode-1) + &
             domain_entries_cpu(iNode-1)
       ENDDO
       global_entries=offset_for_cpu(nNodes)+domain_entries_cpu(nNodes)
       
       ! communicate all entries
       ALLOCATE(domain_map_global(global_entries))
       ALLOCATE(domain_map_local(2*domain_map_local_entries))
       DO ientry=1,domain_map_local_entries
          domain_map_local(2*(ientry-1)+1)=almo_scf_env%domain_map(ispin)%pairs(ientry,1)
          domain_map_local(2*ientry)=almo_scf_env%domain_map(ispin)%pairs(ientry,2)
       ENDDO
       CALL mp_allgather(domain_map_local,domain_map_global,&
               domain_entries_cpu,offset_for_cpu,GroupID)
       DEALLOCATE(domain_entries_cpu,offset_for_cpu)
       DEALLOCATE(domain_map_local)

       DEALLOCATE(almo_scf_env%domain_map(ispin)%index1)
       DEALLOCATE(almo_scf_env%domain_map(ispin)%pairs)
       ALLOCATE(almo_scf_env%domain_map(ispin)%index1(ndomains))
       ALLOCATE(almo_scf_env%domain_map(ispin)%pairs(global_entries/2,2))
       almo_scf_env%domain_map(ispin)%pairs(:,:)=0
       almo_scf_env%domain_map(ispin)%index1(:)=0
       
       ! unpack the received data into a local variable
       ! since we do not know the maximum global number of neighbors
       ! try one. if fails increase the maximum number and try again
       ! until it succeeds
       max_neig=max_domain_neighbors
       max_neig_fails=.TRUE.
       max_neig_loop: DO WHILE (max_neig_fails)
          ALLOCATE(domain_grid(almo_scf_env%ndomains,0:max_neig))
          domain_grid(:,:)=0
          ! init the number of collected neighbors
          domain_grid(:,0)=1
          ! loop over the records
          global_entries=global_entries/2
          DO ientry=1,global_entries
             ! get the center
             grid1=domain_map_global(2*ientry)
             ! get the neighbor
             ineig=domain_map_global(2*(ientry-1)+1)
             ! check boundaries
             IF (domain_grid(grid1,0).gt.max_neig) THEN
                ! this neighbor will overstep the boundaries
                ! stop the trial and increase the max number of neigbors
                DEALLOCATE(domain_grid)
                max_neig=max_neig*2
                CYCLE max_neig_loop
             ENDIF
             ! for the current center loop over the collected neighbors
             ! to insert the current record in a numerical order
             delayed_increment=.FALSE.
             DO igrid=1,domain_grid(grid1,0)
                ! compare the current neighbor with that already in the 'book'
                IF (ineig.lt.domain_grid(grid1,igrid)) THEN
                   ! if this one is smaller then insert it here and pick up the one 
                   ! from the book to continue inserting
                   neig_temp=ineig
                   ineig=domain_grid(grid1,igrid)
                   domain_grid(grid1,igrid)=neig_temp
                ELSE
                   IF (domain_grid(grid1,igrid).eq.0) THEN
                      ! got the empty slot now - insert the record
                      domain_grid(grid1,igrid)=ineig
                      ! increase the record counter but do it outside the loop
                      delayed_increment=.TRUE.
                   ENDIF
                ENDIF
             ENDDO
             IF (delayed_increment) THEN
                domain_grid(grid1,0)=domain_grid(grid1,0)+1
             ELSE
                ! should not be here - all records must be inserted
                CPErrorMessage(cp_failure_level,routineP,"all records must be inserted",error)
                CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
             ENDIF
          ENDDO
          max_neig_fails=.FALSE.
       ENDDO max_neig_loop
       DEALLOCATE(domain_map_global)

       ientry=1
       DO idomain=1,almo_scf_env%ndomains
          DO ineig=1,domain_grid(idomain,0)-1
             almo_scf_env%domain_map(ispin)%pairs(ientry,1)=domain_grid(idomain,ineig)
             almo_scf_env%domain_map(ispin)%pairs(ientry,2)=idomain
             ientry=ientry+1
          ENDDO
          almo_scf_env%domain_map(ispin)%index1(idomain)=ientry
       ENDDO
       DEALLOCATE(domain_grid)

    !ENDDO ! ispin
    IF (almo_scf_env%nspins.eq.2) THEN
       CALL cp_dbcsr_copy(almo_scf_env%quench_t(2),&
               almo_scf_env%quench_t(1),error=error)
       almo_scf_env%domain_map(2)%pairs(:,:)=&
          almo_scf_env%domain_map(1)%pairs(:,:)
       almo_scf_env%domain_map(2)%index1(:)=&
          almo_scf_env%domain_map(1)%index1(:)
    ENDIF

    CALL cp_dbcsr_release(matrix_s_sym,error=error)

    IF (almo_scf_env%domain_layout_mos==almo_domain_layout_molecular .OR. &
       almo_scf_env%domain_layout_aos==almo_domain_layout_molecular) THEN
       DEALLOCATE(first_atom_of_molecule)
       DEALLOCATE(last_atom_of_molecule)
    ENDIF

    !mynode = dbcsr_mp_mynode(dbcsr_distribution_mp(&
    !   cp_dbcsr_distribution(almo_scf_env%quench_t(ispin))))
    !nblkrows_tot = cp_dbcsr_nblkrows_total(almo_scf_env%quench_t(ispin))
    !nblkcols_tot = cp_dbcsr_nblkcols_total(almo_scf_env%quench_t(ispin))
    !DO row = 1, nblkrows_tot
    !   DO col = 1, nblkcols_tot
    !      tr = .FALSE.
    !      iblock_row = row
    !      iblock_col = col
    !      CALL cp_dbcsr_get_stored_coordinates(almo_scf_env%quench_t(ispin),&
    !              iblock_row, iblock_col, tr, hold)
    !      CALL cp_dbcsr_get_block_p(almo_scf_env%quench_t(ispin),&
    !              row, col, p_new_block, found)
    !      write(*,*) "RST_NOTE:", mynode, row, col, hold, found  
    !   ENDDO
    !ENDDO

    CALL timestop(handle)
  
  END SUBROUTINE almo_scf_construct_quencher

END MODULE almo_scf_qs

