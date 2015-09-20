!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief  Routines that link DBCSR and CP2K concepts together
!> \author Ole Schuett
!> \par History
!>         01.2014 created
! *****************************************************************************
MODULE cp_dbcsr_cp2k_link
  USE ao_util,                         ONLY: exp_radius
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: gto_basis_set_p_type,&
                                             gto_basis_set_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_interface,              ONLY: &
       accdrv_avoid_after_busy, accdrv_binning_binsize, accdrv_binning_nbins, &
       accdrv_do_inhomogenous, accdrv_min_flop_process, accdrv_min_flop_sort, &
       accdrv_posterior_buffers, accdrv_posterior_streams, &
       accdrv_priority_buffers, accdrv_priority_streams, cp_dbcsr_finalize, &
       cp_dbcsr_get_block_p, cp_dbcsr_get_matrix_type, cp_dbcsr_has_symmetry, &
       cp_dbcsr_reserve_blocks, cp_dbcsr_set, cp_dbcsr_type, dbcsr_error_set, &
       dbcsr_error_stop, dbcsr_error_type, &
       dbcsr_get_conf_avg_elements_images, dbcsr_get_conf_comm_thread_load, &
       dbcsr_get_conf_mm_driver, dbcsr_get_conf_mm_stacksize, &
       dbcsr_get_conf_mpi_mem, dbcsr_get_conf_nstacks, &
       dbcsr_get_conf_randmat_seed, dbcsr_get_conf_use_comm_thread, &
       dbcsr_get_conf_use_mpi_filtering, dbcsr_get_conf_use_mpi_rma, &
       dbcsr_set_conf_avg_elements_images, dbcsr_set_conf_comm_thread_load, &
       dbcsr_set_conf_external_timing, dbcsr_set_conf_max_ele_block, &
       dbcsr_set_conf_mm_driver, dbcsr_set_conf_mm_stacksize, &
       dbcsr_set_conf_nstacks, dbcsr_set_conf_randmat_seed, &
       dbcsr_set_conf_use_comm_thread, dbcsr_set_conf_use_mpi_filtering, &
       dbcsr_set_conf_use_mpi_rma, dbcsr_type_no_symmetry, has_acc, has_mpi, &
       mm_driver_acc, mm_driver_blas, mm_driver_matmul, mm_driver_smm, &
       mm_name_acc, mm_name_blas, mm_name_matmul, mm_name_smm, multrec_limit
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp,&
                                             real_4
  USE orbital_pointers,                ONLY: nso
  USE qs_integral_utils,               ONLY: basis_set_list_setup
  USE qs_kind_types,                   ONLY: qs_kind_type
  USE qs_ks_types,                     ONLY: get_ks_env,&
                                             qs_ks_env_type
  USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                             get_neighbor_list_set_p,&
                                             neighbor_list_iterate,&
                                             neighbor_list_iterator_create,&
                                             neighbor_list_iterator_p_type,&
                                             neighbor_list_iterator_release,&
                                             neighbor_list_set_p_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_dbcsr_cp2k_link'


  PUBLIC :: cp_dbcsr_config
  PUBLIC :: cp_dbcsr_print_config
  PUBLIC :: cp_dbcsr_alloc_block_from_nbl
  PUBLIC :: cp_dbcsr_to_csr_screening

  PRIVATE

CONTAINS


! *****************************************************************************
!> \brief Configures options for DBCSR
!> \param root_section ...
! *****************************************************************************
  SUBROUTINE cp_dbcsr_config(root_section)
    TYPE(section_vals_type), POINTER         :: root_section

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_config', &
      routineP = moduleN//':'//routineN

    INTEGER :: avg_elements_images, comm_thread_load, error_handle, &
      max_ele_block, mm_driver, mm_ss, randmat_seed
    INTEGER, DIMENSION(3)                    :: nstacks
    LOGICAL                                  :: use_comm_thread, &
                                                use_mpi_filtering, use_mpi_rma
    TYPE(dbcsr_error_type)                   :: dbcsr_error
    TYPE(section_vals_type), POINTER         :: dbcsr_section

! need to call dbcsr_set_conf_external_timing() before dbcsr_error_set()
! otherwise this would lead to mismatched calls to timeset/timestop

    CALL dbcsr_set_conf_external_timing(timeset, timestop, error=dbcsr_error)

    CALL dbcsr_error_set (routineN, error_handle, dbcsr_error)
    dbcsr_section => section_vals_get_subs_vals(root_section,&
      "GLOBAL%DBCSR")
    CALL section_vals_val_get(dbcsr_section,&
         "mm_stack_size", i_val=mm_ss)
    CALL section_vals_val_get(dbcsr_section,&
         "mm_driver", i_val=mm_driver)
    CALL section_vals_val_get(dbcsr_section,&
         "MAX_ELEMENTS_PER_BLOCK", i_val=max_ele_block)
    CALL section_vals_val_get(dbcsr_section,&
         "avg_elements_images", i_val=avg_elements_images)
    CALL section_vals_val_get(dbcsr_section,&
         "n_size_mnk_stacks", i_val=nstacks(1))
    nstacks(2:3) = nstacks(1)
    CALL section_vals_val_get(dbcsr_section,&
         "use_mpi_filtering", l_val=use_mpi_filtering)
    CALL section_vals_val_get(dbcsr_section,&
         "use_mpi_rma", l_val=use_mpi_rma)
    CALL section_vals_val_get(dbcsr_section,&
         "use_comm_thread", l_val=use_comm_thread)
    CALL section_vals_val_get(dbcsr_section,&
         "comm_thread_load", i_val=comm_thread_load)
    CALL section_vals_val_get(dbcsr_section,&
         "multrec_limit", i_val=multrec_limit)
    CALL section_vals_val_get(dbcsr_section,&
         "randmat_seed", i_val=randmat_seed)

    CALL section_vals_val_get(dbcsr_section,&
         "ACC%priority_streams", i_val=accdrv_priority_streams)
    CALL section_vals_val_get(dbcsr_section,&
         "ACC%priority_buffers", i_val=accdrv_priority_buffers)
    CALL section_vals_val_get(dbcsr_section,&
         "ACC%posterior_streams", i_val=accdrv_posterior_streams)
    CALL section_vals_val_get(dbcsr_section,&
         "ACC%posterior_buffers", i_val=accdrv_posterior_buffers)

    CALL section_vals_val_get(dbcsr_section,&
         "ACC%min_flop_process", i_val=accdrv_min_flop_process)
    CALL section_vals_val_get(dbcsr_section,&
         "ACC%min_flop_sort", i_val=accdrv_min_flop_sort)
    CALL section_vals_val_get(dbcsr_section,&
         "ACC%process_inhomogenous", l_val=accdrv_do_inhomogenous)
    CALL section_vals_val_get(dbcsr_section,&
         "ACC%avoid_after_busy", l_val=accdrv_avoid_after_busy)
    CALL section_vals_val_get(dbcsr_section,&
         "ACC%binning_nbins", i_val=accdrv_binning_nbins)
    CALL section_vals_val_get(dbcsr_section,&
         "ACC%binning_binsize", i_val=accdrv_binning_binsize)

    CALL dbcsr_set_conf_mm_driver (mm_driver, error=dbcsr_error)
    CALL dbcsr_set_conf_max_ele_block (max_ele_block, error=dbcsr_error)
    IF(mm_ss > 0) &
    CALL dbcsr_set_conf_mm_stacksize (mm_ss, error=dbcsr_error)
    IF(avg_elements_images>0) &
         CALL dbcsr_set_conf_avg_elements_images(avg_elements_images, error=dbcsr_error)
    IF(nstacks(1) > 0) &
    CALL dbcsr_set_conf_nstacks(nstacks, error=dbcsr_error)
    CALL dbcsr_set_conf_use_mpi_filtering(use_mpi_filtering, error=dbcsr_error)
    CALL dbcsr_set_conf_use_mpi_rma(use_mpi_rma, error=dbcsr_error)
    CALL dbcsr_set_conf_use_comm_thread(use_comm_thread, error=dbcsr_error)
    CALL dbcsr_set_conf_comm_thread_load(comm_thread_load, error=dbcsr_error)
    CALL dbcsr_set_conf_randmat_seed(randmat_seed, error=dbcsr_error)

    CALL dbcsr_error_stop (error_handle, dbcsr_error)
  END SUBROUTINE cp_dbcsr_config


! *****************************************************************************
!> \brief Prints configuration for DBCSR
!> \param unit_nr ...
! *****************************************************************************
  SUBROUTINE cp_dbcsr_print_config(unit_nr)
    INTEGER, INTENT(IN), OPTIONAL            :: unit_nr

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_print_config', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length)     :: mm_name
    INTEGER                                  :: unit_num
    INTEGER, DIMENSION(3)                    :: n_mnk_stacks
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    logger => cp_get_default_logger()
    unit_num = cp_logger_get_default_io_unit(logger)
    IF (PRESENT (unit_nr)) unit_num=unit_nr

    IF (unit_num>0 .AND. cp_logger_would_log(logger, cp_note_level)) THEN
       SELECT CASE (dbcsr_get_conf_mm_driver())
       CASE(mm_driver_blas)   ; mm_name = mm_name_blas
       CASE(mm_driver_matmul) ; mm_name = mm_name_matmul
       CASE(mm_driver_smm)    ; mm_name = mm_name_smm
       CASE(mm_driver_acc)   ; mm_name = mm_name_acc
       END SELECT
       WRITE(UNIT=unit_num, FMT='(1X,A,T41,A40)')&
        "DBCSR| Multiplication driver", ADJUSTR(mm_name(1:40))

       WRITE(UNIT=unit_num, FMT='(1X,A,T70,I11)')&
        "DBCSR| Multrec recursion limit", multrec_limit
       WRITE(UNIT=unit_num, FMT='(1X,A,T70,I11)')&
        "DBCSR| Multiplication stack size", dbcsr_get_conf_mm_stacksize()

       IF (dbcsr_get_conf_avg_elements_images()>0) THEN
          WRITE(UNIT=unit_num, FMT='(1X,A,T70,I11)')&
               "DBCSR| Average elements for images", dbcsr_get_conf_avg_elements_images()
       ELSE
          WRITE(UNIT=unit_num, FMT='(1X,A,T72,A)')&
               "DBCSR| Maximum elements for images","UNLIMITED"
       ENDIF
       WRITE(UNIT=unit_num, FMT='(1X,A,T70,I11)')&
        "DBCSR| Randmat seed", dbcsr_get_conf_randmat_seed()
       
       CALL dbcsr_get_conf_nstacks (n_mnk_stacks, error=dbcsr_error)
       IF (ALL(n_mnk_stacks .EQ. n_mnk_stacks(1))) THEN
           WRITE(UNIT=unit_num, FMT='(1X,A,T70,I11)')&
            "DBCSR| Multiplication size stacks", n_mnk_stacks(1)
       ELSE
           WRITE(UNIT=unit_num, FMT='(1X,A,T70,I11)')&
            "DBCSR| Multiplication size m stacks", n_mnk_stacks(1)
           WRITE(UNIT=unit_num, FMT='(1X,A,T70,I11)')&
            "DBCSR| Multiplication size n stacks", n_mnk_stacks(2)
           WRITE(UNIT=unit_num, FMT='(1X,A,T70,I11)')&
            "DBCSR| Multiplication size l stacks", n_mnk_stacks(3)
       ENDIF

       IF (has_mpi) THEN
           WRITE(UNIT=unit_num, FMT='(1X,A,T80,L1)')&
            "DBCSR| Use MPI memory allocation", dbcsr_get_conf_mpi_mem ()
           WRITE(UNIT=unit_num, FMT='(1X,A,T80,L1)')&
            "DBCSR| Use RMA algorithm", dbcsr_get_conf_use_mpi_rma ()
           WRITE(UNIT=unit_num, FMT='(1X,A,T80,L1)')&
             "DBCSR| Use Communication thread",  dbcsr_get_conf_use_comm_thread()
           WRITE(UNIT=unit_num, FMT='(1X,A,T70,I11)')&
             "DBCSR| Communication thread load", dbcsr_get_conf_comm_thread_load()
           IF (dbcsr_get_conf_use_mpi_rma ()) THEN
              WRITE(UNIT=unit_num, FMT='(1X,A,T80,L1)')&
                   "DBCSR| Use MPI filtering", dbcsr_get_conf_use_mpi_filtering ()
           ENDIF
       ENDIF

       IF (has_acc) THEN
           WRITE(UNIT=unit_num, FMT='(1X,A,T70,I11)')&
            "DBCSR| ACC: Number of priority stack-buffers", accdrv_priority_buffers
           WRITE(UNIT=unit_num, FMT='(1X,A,T70,I11)')&
            "DBCSR| ACC: Number of posterior stack-buffers", accdrv_posterior_buffers
           WRITE(UNIT=unit_num, FMT='(1X,A,T70,I11)')&
            "DBCSR| ACC: Number of priority streams", accdrv_priority_streams
           WRITE(UNIT=unit_num, FMT='(1X,A,T70,I11)')&
            "DBCSR| ACC: Number of posterior streams", accdrv_posterior_streams
           WRITE(UNIT=unit_num, FMT='(1X,A,T80,L1)')&
            "DBCSR| ACC: Avoid driver after busy ", accdrv_avoid_after_busy
           WRITE(UNIT=unit_num, FMT='(1X,A,T80,L1)')&
            "DBCSR| ACC: Process inhomogenous stacks", accdrv_do_inhomogenous
           WRITE(UNIT=unit_num, FMT='(1X,A,T70,I11)')&
            "DBCSR| ACC: Min. flop for processing", accdrv_min_flop_process
           WRITE(UNIT=unit_num, FMT='(1X,A,T70,I11)')&
            "DBCSR| ACC: Min. flop for sorting", accdrv_min_flop_sort
           WRITE(UNIT=unit_num, FMT='(1X,A,T70,I11)')&
            "DBCSR| ACC: Number of binning bins", accdrv_binning_nbins
           WRITE(UNIT=unit_num, FMT='(1X,A,T70,I11)')&
            "DBCSR| ACC: Size of binning bins", accdrv_binning_binsize
       END IF
    ENDIF

  END SUBROUTINE cp_dbcsr_print_config



! *****************************************************************************
!> \brief allocate the blocks of a dbcsr based on the neighbor list
!> \param matrix        the matrix
!> \param sab_orb       the corresponding neighbor list
!> \par History
!>      11.2009 created vw
!>      01.2014 moved here from cp_dbcsr_operations (Ole Schuett)
!> \author vw
!> \note
! *****************************************************************************

  SUBROUTINE cp_dbcsr_alloc_block_from_nbl(matrix,sab_orb)

    TYPE(cp_dbcsr_type)                      :: matrix
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'cp_dbcsr_alloc_block_from_nbl', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=1)                         :: symmetry
    INTEGER                                  :: blk_cnt, handle, iatom, icol, &
                                                inode, irow, jatom, last_jatom
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: cols, rows, tmp
    LOGICAL                                  :: failure, new_atom_b
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: nl_iterator

    CALL timeset(routineN,handle)

    failure = .FALSE.

    symmetry = cp_dbcsr_get_matrix_type(matrix)

    CPPrecondition(ASSOCIATED(sab_orb),cp_failure_level,routineP,failure)

    CALL cp_dbcsr_finalize (matrix)
    ALLOCATE (rows(1), cols(1))
    blk_cnt = 0

    CALL neighbor_list_iterator_create(nl_iterator,sab_orb)
    DO WHILE (neighbor_list_iterate(nl_iterator)==0)
       CALL get_iterator_info(nl_iterator,iatom=iatom,jatom=jatom,inode=inode)
       IF(inode==1) last_jatom = 0
       IF (jatom /= last_jatom) THEN
          new_atom_b = .TRUE.
          last_jatom = jatom
       ELSE
          new_atom_b = .FALSE.
          CYCLE
       END IF
       IF (blk_cnt+1 .GT. SIZE(rows)) THEN
          ALLOCATE (tmp (blk_cnt))
          tmp(:) = rows(:)
          DEALLOCATE (rows)
          ALLOCATE (rows((blk_cnt+1)*2))
          rows(1:blk_cnt) = tmp(1:blk_cnt)
          tmp(:) = cols(:)
          DEALLOCATE (cols)
          ALLOCATE (cols((blk_cnt+1)*2))
          cols(1:blk_cnt) = tmp(1:blk_cnt)
          DEALLOCATE (tmp)
       ENDIF
       blk_cnt = blk_cnt+1
       IF(symmetry==dbcsr_type_no_symmetry) THEN
          rows(blk_cnt) = iatom
          cols(blk_cnt) = jatom
       ELSE
          IF(iatom<=jatom) THEN
             irow = iatom
             icol = jatom
          ELSE
             irow = jatom
             icol = iatom
          END IF
          rows(blk_cnt) = irow
          cols(blk_cnt) = icol
       END IF

    END DO
    CALL neighbor_list_iterator_release(nl_iterator)

    !
    CALL cp_dbcsr_reserve_blocks (matrix, rows(1:blk_cnt), cols(1:blk_cnt))
    DEALLOCATE (rows)
    DEALLOCATE (cols)
    CALL cp_dbcsr_finalize( matrix)

    CALL timestop(handle)

  END SUBROUTINE cp_dbcsr_alloc_block_from_nbl

! *****************************************************************************
!> \brief Apply distance screening to refine sparsity pattern of matrices in CSR
!>        format (using eps_pgf_orb). Currently this is used for the external 
!>        library PEXSI.
!> \param ks_env ...
!> \param[in, out] csr_sparsity DBCSR matrix defining CSR sparsity pattern. 
!>                              This matrix must be initialized and allocated
!>                              with exactly the same DBCSR sparsity pattern as 
!>                              the DBCSR matrix that is used to create the CSR
!>                              matrix. It must have symmetric DBCSR format and 
!>                              must not be filtered.
!> \par History
!>      02.2015 created [Patrick Seewald]
!> \author Patrick Seewald
! *****************************************************************************
  SUBROUTINE cp_dbcsr_to_csr_screening(ks_env, csr_sparsity)
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: csr_sparsity

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_to_csr_screening', &
      routineP = moduleN//':'//routineN

    INTEGER :: atom_a, atom_b, handle, iatom, icol, ikind, ipgf, irow, iset, &
      isgf, ishell, jatom, jkind, jpgf, jset, jsgf, jshell, natom, nkind, &
      nset_a, nset_b
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind
    INTEGER, DIMENSION(:), POINTER           :: npgf_a, npgf_b, nshell_a, &
                                                nshell_b
    INTEGER, DIMENSION(:, :), POINTER        :: l_a, l_b
    LOGICAL                                  :: do_symmetric, failure, found
    REAL(KIND=dp)                            :: dab, eps_pgf_orb, r_a, r_b
    REAL(KIND=dp), DIMENSION(3)              :: rab
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: rpgfa, rpgfb, zet_a, zet_b
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: gcc_a, gcc_b
    REAL(KIND=real_4), DIMENSION(:, :), &
      POINTER                                :: screen_blk
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(gto_basis_set_p_type), &
      DIMENSION(:), POINTER                  :: basis_set_list_a, &
                                                basis_set_list_b
    TYPE(gto_basis_set_type), POINTER        :: basis_set_a, basis_set_b
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: nl_iterator
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: neighbour_list
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set

    NULLIFY(screen_blk, atomic_kind_set, basis_set_list_a, &
            basis_set_list_b, basis_set_a, basis_set_b, nl_iterator, &
            qs_kind_set, dft_control)

    CALL timeset(routineN,handle)

    CPPrecondition(cp_dbcsr_has_symmetry(csr_sparsity),cp_failure_level,routineP,failure)

    CALL get_ks_env(ks_env,&
                    sab_orb=neighbour_list,&
                    atomic_kind_set=atomic_kind_set,&
                    natom=natom,&
                    qs_kind_set=qs_kind_set,&
                    dft_control=dft_control)
 
    eps_pgf_orb = dft_control%qs_control%eps_pgf_orb

    nkind = SIZE(qs_kind_set)
    ALLOCATE (atom_of_kind(natom))
    CALL get_atomic_kind_set(atomic_kind_set,atom_of_kind=atom_of_kind)
    CPPrecondition(SIZE(neighbour_list) > 0,cp_failure_level,routineP,failure)
    CALL get_neighbor_list_set_p(neighbor_list_sets=neighbour_list,symmetric=do_symmetric)
    CPPrecondition(do_symmetric,cp_failure_level,routineP,failure)
    ALLOCATE (basis_set_list_a(nkind),basis_set_list_b(nkind))
    CALL basis_set_list_setup(basis_set_list_a,"ORB",qs_kind_set)
    CALL basis_set_list_setup(basis_set_list_b,"ORB",qs_kind_set)

    ! csr_sparsity can obtain values 0 (if zero element) or 1 (if non-zero element)
    CALL cp_dbcsr_set (csr_sparsity, 0.0)

    CALL neighbor_list_iterator_create(nl_iterator, neighbour_list)

    ! Iterate over interacting pairs of atoms corresponding to non-zero
    ! DBCSR blocks
    DO WHILE (neighbor_list_iterate(nl_iterator)==0)
      CALL get_iterator_info(nl_iterator, &
                             ikind=ikind, jkind=jkind,&
                             iatom=iatom, jatom=jatom, &
                             r=rab)

      basis_set_a => basis_set_list_a(ikind)%gto_basis_set
      IF (.NOT.ASSOCIATED(basis_set_a)) CYCLE
      basis_set_b => basis_set_list_b(jkind)%gto_basis_set
      IF (.NOT.ASSOCIATED(basis_set_b)) CYCLE

      atom_a = atom_of_kind(iatom)
      atom_b = atom_of_kind(jatom)

      nset_a = basis_set_a%nset
      nset_b = basis_set_b%nset
      npgf_a => basis_set_a%npgf
      npgf_b => basis_set_b%npgf
      nshell_a => basis_set_a%nshell
      nshell_b => basis_set_b%nshell

      l_a => basis_set_a%l
      l_b => basis_set_b%l
      gcc_a => basis_set_a%gcc
      gcc_b => basis_set_b%gcc
      zet_a => basis_set_a%zet
      zet_b => basis_set_b%zet

      rpgfa =>  basis_set_a%pgf_radius
      rpgfb =>  basis_set_b%pgf_radius

      IF (iatom <= jatom) THEN
        irow = iatom
        icol = jatom
      ELSE
        irow = jatom
        icol = iatom
      END IF

      CALL cp_dbcsr_get_block_p(matrix=csr_sparsity,row=irow,col=icol, &
                                block=screen_blk,found=found)

      CPPostcondition(found,cp_failure_level,routineP,failure)

      ! Distance between atoms a and b
      dab = SQRT(rab(1)**2+rab(2)**2+rab(3)**2)
  
      ! iterate over pairs of primitive GTOs i,j, get their radii r_i, r_j according 
      ! to eps_pgf_orb. Define all matrix elements as non-zero to which a 
      ! contribution from two Gaussians i,j exists with r_i + r_j >= dab.
  
      isgf=0
      DO iset=1,nset_a
        DO ishell=1,nshell_a(iset)
          jsgf=0
          DO jset=1,nset_b
            DO jshell=1,nshell_b(jset)
              gto_loop: DO ipgf=1,npgf_a(iset)
                DO jpgf=1,npgf_b(jset)
                  IF (rpgfa(ipgf,iset)+rpgfb(jpgf,jset) .GE. dab) THEN
                    ! more selective screening with radius calculated for each primitive GTO
                    r_a = exp_radius(l_a(ishell,iset), &
                                     zet_a(ipgf,iset), &
                                     eps_pgf_orb, &
                                     gcc_a(ipgf,ishell,iset))
                    r_b = exp_radius(l_b(jshell,jset), &
                                     zet_b(jpgf,jset), &
                                     eps_pgf_orb, &
                                     gcc_b(jpgf,jshell,jset))
                    IF(r_a+r_b .GE. dab) THEN
                      IF(irow .EQ. iatom) THEN 
                        screen_blk(isgf+1:isgf+nso(l_a(ishell,iset)), &
                                   jsgf+1:jsgf+nso(l_b(jshell,jset))) = 1.0_dp
                      ELSE
                        screen_blk(jsgf+1:jsgf+nso(l_b(jshell,jset)), &
                                   isgf+1:isgf+nso(l_a(ishell,iset))) = 1.0_dp
                      ENDIF
                      EXIT gto_loop
                    ENDIF
                  ENDIF
                ENDDO
              ENDDO gto_loop
              jsgf=jsgf+nso(l_b(jshell,jset))
            ENDDO
          ENDDO
          isgf=isgf+nso(l_a(ishell,iset))
        ENDDO
      ENDDO
    ENDDO

    CALL neighbor_list_iterator_release(nl_iterator)
    DEALLOCATE (basis_set_list_a,basis_set_list_b)

    CALL timestop(handle)
  END SUBROUTINE cp_dbcsr_to_csr_screening

END MODULE cp_dbcsr_cp2k_link
