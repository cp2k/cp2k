!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief  Routines that link DBCSR and CP2K concepts together
!> \author Ole Schuett
!> \par History
!>         01.2014 created
! *****************************************************************************
MODULE cp_dbcsr_cp2k_link
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_finalize,&
                                             cp_dbcsr_get_matrix_type,&
                                             cp_dbcsr_reserve_blocks
  USE cp_dbcsr_types,                  ONLY: cp_dbcsr_type
  USE dbcsr_config,                    ONLY: &
       cudadrv_avoid_after_busy, cudadrv_binning_binsize, &
       cudadrv_binning_nbins, cudadrv_do_inhomogenous, &
       cudadrv_min_flop_process, cudadrv_min_flop_sort, &
       cudadrv_posterior_buffers, cudadrv_posterior_streams, &
       cudadrv_priority_buffers, cudadrv_priority_streams, &
       dbcsr_get_conf_combtypes, dbcsr_get_conf_comm_thread_load, &
       dbcsr_get_conf_mm_driver, dbcsr_get_conf_mm_stacksize, &
       dbcsr_get_conf_mpi_mem, dbcsr_get_conf_nstacks, &
       dbcsr_get_conf_subcomm, dbcsr_get_conf_use_comm_thread, &
       dbcsr_set_conf_combtypes, dbcsr_set_conf_comm_thread_load, &
       dbcsr_set_conf_external_timing, dbcsr_set_conf_max_ele_block, &
       dbcsr_set_conf_mm_driver, dbcsr_set_conf_mm_stacksize, &
       dbcsr_set_conf_mpi_mem, dbcsr_set_conf_nstacks, &
       dbcsr_set_conf_subcomm, dbcsr_set_conf_use_comm_thread, has_cuda, &
       has_mpi, mm_driver_blas, mm_driver_cuda, mm_driver_matmul, &
       mm_driver_plasma, mm_driver_smm, mm_name_blas, mm_name_cuda, &
       mm_name_matmul, mm_name_plasma, mm_name_smm, multrec_limit
  USE dbcsr_dist_operations
  USE dbcsr_error_handling,            ONLY: dbcsr_error_set,&
                                             dbcsr_error_stop,&
                                             dbcsr_error_type
  USE dbcsr_io
  USE dbcsr_types,                     ONLY: dbcsr_type_no_symmetry
  USE f77_blas
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length
  USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
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

  PRIVATE

CONTAINS


! *****************************************************************************
!> \brief Configures options for DBCSR
! *****************************************************************************
  SUBROUTINE cp_dbcsr_config(root_section, error)
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_config', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: comm_thread_load, &
                                                error_handle, max_ele_block, &
                                                mm_driver, mm_ss
    INTEGER, DIMENSION(3)                    :: nstacks
    LOGICAL                                  :: use_combtypes, &
                                                use_comm_thread, use_mpi_mem, &
                                                use_subcomms
    TYPE(dbcsr_error_type)                   :: dbcsr_error
    TYPE(section_vals_type), POINTER         :: dbcsr_section

! need to call dbcsr_set_conf_external_timing() before dbcsr_error_set()
! otherwise this would lead to mismatched calls to timeset/timestop

    CALL dbcsr_set_conf_external_timing(timeset, timestop, error=dbcsr_error)

    CALL dbcsr_error_set (routineN, error_handle, dbcsr_error)
    dbcsr_section => section_vals_get_subs_vals(root_section,&
      "GLOBAL%DBCSR",error=error)
    CALL section_vals_val_get(dbcsr_section,&
         "subcommunicators", l_val=use_subcomms, error=error)
    CALL section_vals_val_get(dbcsr_section,&
         "combined_types", l_val=use_combtypes, error=error)
    CALL section_vals_val_get(dbcsr_section,&
         "use_mpi_allocation", l_val=use_mpi_mem, error=error)
    CALL section_vals_val_get(dbcsr_section,&
         "mm_stack_size", i_val=mm_ss, error=error)
    CALL section_vals_val_get(dbcsr_section,&
         "mm_driver", i_val=mm_driver, error=error)
    CALL section_vals_val_get(dbcsr_section,&
         "MAX_ELEMENTS_PER_BLOCK", i_val=max_ele_block, error=error)
    CALL section_vals_val_get(dbcsr_section,&
         "n_size_mnk_stacks", i_val=nstacks(1), error=error)
    nstacks(2:3) = nstacks(1)
    CALL section_vals_val_get(dbcsr_section,&
         "use_comm_thread", l_val=use_comm_thread, error=error)
    CALL section_vals_val_get(dbcsr_section,&
         "comm_thread_load", i_val=comm_thread_load, error=error)
    CALL section_vals_val_get(dbcsr_section,&
         "multrec_limit", i_val=multrec_limit, error=error)

    CALL section_vals_val_get(dbcsr_section,&
         "CUDA%priority_streams", i_val=cudadrv_priority_streams, error=error)
    CALL section_vals_val_get(dbcsr_section,&
         "CUDA%priority_buffers", i_val=cudadrv_priority_buffers, error=error)
    CALL section_vals_val_get(dbcsr_section,&
         "CUDA%posterior_streams", i_val=cudadrv_posterior_streams, error=error)
    CALL section_vals_val_get(dbcsr_section,&
         "CUDA%posterior_buffers", i_val=cudadrv_posterior_buffers, error=error)

    CALL section_vals_val_get(dbcsr_section,&
         "CUDA%min_flop_process", i_val=cudadrv_min_flop_process, error=error)
    CALL section_vals_val_get(dbcsr_section,&
         "CUDA%min_flop_sort", i_val=cudadrv_min_flop_sort, error=error)
    CALL section_vals_val_get(dbcsr_section,&
         "CUDA%process_inhomogenous", l_val=cudadrv_do_inhomogenous, error=error)
    CALL section_vals_val_get(dbcsr_section,&
         "CUDA%avoid_after_busy", l_val=cudadrv_avoid_after_busy, error=error)
    CALL section_vals_val_get(dbcsr_section,&
         "CUDA%binning_nbins", i_val=cudadrv_binning_nbins, error=error)
    CALL section_vals_val_get(dbcsr_section,&
         "CUDA%binning_binsize", i_val=cudadrv_binning_binsize, error=error)

    CALL dbcsr_set_conf_mm_driver (mm_driver, error=dbcsr_error)
    CALL dbcsr_set_conf_subcomm (use_subcomms, error=dbcsr_error)
    CALL dbcsr_set_conf_combtypes (use_combtypes, error=dbcsr_error)
    CALL dbcsr_set_conf_mpi_mem (use_mpi_mem, error=dbcsr_error)
    CALL dbcsr_set_conf_max_ele_block (max_ele_block, error=dbcsr_error)
    IF(mm_ss > 0) &
    CALL dbcsr_set_conf_mm_stacksize (mm_ss, error=dbcsr_error)
    IF(nstacks(1) > 0) &
    CALL dbcsr_set_conf_nstacks(nstacks, error=dbcsr_error)

    CALL dbcsr_set_conf_use_comm_thread(use_comm_thread, error=dbcsr_error)
    CALL dbcsr_set_conf_comm_thread_load(comm_thread_load, error=dbcsr_error)

    CALL dbcsr_error_stop (error_handle, dbcsr_error)
  END SUBROUTINE cp_dbcsr_config


! *****************************************************************************
!> \brief Prints configuration for DBCSR
! *****************************************************************************
  SUBROUTINE cp_dbcsr_print_config(unit_nr, error)
    INTEGER, INTENT(IN), OPTIONAL            :: unit_nr
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_print_config', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length)     :: mm_name
    INTEGER                                  :: unit_num
    INTEGER, DIMENSION(3)                    :: n_mnk_stacks
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    logger => cp_error_get_logger(error)
    unit_num = cp_logger_get_default_io_unit(logger)
    IF (PRESENT (unit_nr)) unit_num=unit_nr

    IF (unit_num>0 .AND. cp_logger_would_log(logger, cp_note_level)) THEN
       SELECT CASE (dbcsr_get_conf_mm_driver())
       CASE(mm_driver_blas)   ; mm_name = mm_name_blas
       CASE(mm_driver_matmul) ; mm_name = mm_name_matmul
       CASE(mm_driver_smm)    ; mm_name = mm_name_smm
       CASE(mm_driver_plasma) ; mm_name = mm_name_plasma
       CASE(mm_driver_cuda)   ; mm_name = mm_name_cuda
       END SELECT
       WRITE(UNIT=unit_num, FMT='(1X,A,T41,A40)')&
        "DBCSR| Multiplication driver", ADJUSTR(mm_name(1:40))

       WRITE(UNIT=unit_num, FMT='(1X,A,T70,I11)')&
        "DBCSR| Multrec recursion limit", multrec_limit
       WRITE(UNIT=unit_num, FMT='(1X,A,T70,I11)')&
        "DBCSR| Multiplication stack size", dbcsr_get_conf_mm_stacksize()

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
            "DBCSR| Use subcommunicators",  dbcsr_get_conf_subcomm ()
           WRITE(UNIT=unit_num, FMT='(1X,A,T80,L1)')&
            "DBCSR| Use MPI combined types" ,  dbcsr_get_conf_combtypes ()
           WRITE(UNIT=unit_num, FMT='(1X,A,T80,L1)')&
            "DBCSR| Use MPI memory allocation", dbcsr_get_conf_mpi_mem ()
           WRITE(UNIT=unit_num, FMT='(1X,A,T80,L1)')&
            "DBCSR| Use Communication thread",  dbcsr_get_conf_use_comm_thread()
           WRITE(UNIT=unit_num, FMT='(1X,A,T70,I11)')&
            "DBCSR| Communication thread load", dbcsr_get_conf_comm_thread_load()
       ENDIF

       IF (has_cuda) THEN
           WRITE(UNIT=unit_num, FMT='(1X,A,T70,I11)')&
            "DBCSR| CUDA: Number of priority stack-buffers", cudadrv_priority_buffers
           WRITE(UNIT=unit_num, FMT='(1X,A,T70,I11)')&
            "DBCSR| CUDA: Number of posterior stack-buffers", cudadrv_posterior_buffers
           WRITE(UNIT=unit_num, FMT='(1X,A,T70,I11)')&
            "DBCSR| CUDA: Number of priority streams", cudadrv_priority_streams
           WRITE(UNIT=unit_num, FMT='(1X,A,T70,I11)')&
            "DBCSR| CUDA: Number of posterior streams", cudadrv_posterior_streams
           WRITE(UNIT=unit_num, FMT='(1X,A,T80,L1)')&
            "DBCSR| CUDA: Avoid driver after busy ", cudadrv_avoid_after_busy
           WRITE(UNIT=unit_num, FMT='(1X,A,T80,L1)')&
            "DBCSR| CUDA: Process inhomogenous stacks", cudadrv_do_inhomogenous
           WRITE(UNIT=unit_num, FMT='(1X,A,T70,I11)')&
            "DBCSR| CUDA: Min. flop for processing", cudadrv_min_flop_process
           WRITE(UNIT=unit_num, FMT='(1X,A,T70,I11)')&
            "DBCSR| CUDA: Min. flop for sorting", cudadrv_min_flop_sort
           WRITE(UNIT=unit_num, FMT='(1X,A,T70,I11)')&
            "DBCSR| CUDA: Number of binning bins", cudadrv_binning_nbins
           WRITE(UNIT=unit_num, FMT='(1X,A,T70,I11)')&
            "DBCSR| CUDA: Size of binning bins", cudadrv_binning_binsize
       END IF
    ENDIF

  END SUBROUTINE cp_dbcsr_print_config



! *****************************************************************************
!> \brief allocate the blocks of a dbcsr based on the neighbor list
!> \note
!> \param matrix        the matrix
!> \param sab_orb       the corresponding neighbor list
!> \param error
!> \par History
!>      11.2009 created vw
!>      01.2014 moved here from cp_dbcsr_operations (Ole Schuett)
!> \author vw
! *****************************************************************************

  SUBROUTINE cp_dbcsr_alloc_block_from_nbl(matrix,sab_orb,error)

    TYPE(cp_dbcsr_type)                      :: matrix
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb
    TYPE(cp_error_type), INTENT(inout)       :: error

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

    CPPrecondition(ASSOCIATED(sab_orb),cp_failure_level,routineP,error,failure)

    CALL cp_dbcsr_finalize (matrix, error=error)
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
    CALL cp_dbcsr_reserve_blocks (matrix, rows(1:blk_cnt), cols(1:blk_cnt),&
         error=error)
    DEALLOCATE (rows)
    DEALLOCATE (cols)
    CALL cp_dbcsr_finalize( matrix, error=error )

    CALL timestop(handle)

  END SUBROUTINE cp_dbcsr_alloc_block_from_nbl

END MODULE cp_dbcsr_cp2k_link
