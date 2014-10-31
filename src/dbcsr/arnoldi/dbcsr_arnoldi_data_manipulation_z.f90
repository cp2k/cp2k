  SUBROUTINE  setup_arnoldi_data_z (arnoldi_data, matrix, max_iter)
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data
    TYPE(dbcsr_obj_type_p), DIMENSION(:)     :: matrix
    INTEGER                                  :: max_iter

    CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_arnoldi_data_z', &
      routineP = moduleN//':'//routineN

    INTEGER                                           :: nrow_local
    TYPE(arnoldi_data_z), POINTER           :: ar_data

    ALLOCATE(ar_data)
    CALL dbcsr_get_info(matrix=matrix(1)%matrix, nfullrows_local=nrow_local)
    ALLOCATE(ar_data%f_vec(nrow_local))
    ALLOCATE(ar_data%Hessenberg(max_iter+1, max_iter))
    ALLOCATE(ar_data%local_history(nrow_local, max_iter))

    ALLOCATE(ar_data%evals(max_iter))
    ALLOCATE(ar_data%revec(max_iter, max_iter))

    CALL set_data_z(arnoldi_data,ar_data)

  END SUBROUTINE setup_arnoldi_data_z

  SUBROUTINE deallocate_arnoldi_data_z (arnoldi_data)
    TYPE(dbcsr_arnoldi_data)                     :: arnoldi_data

    CHARACTER(LEN=*), PARAMETER :: routineN = 'deallocate_arnoldi_data_z', &
      routineP = moduleN//':'//routineN

    TYPE(arnoldi_data_z), POINTER            :: ar_data

    ar_data=>get_data_z(arnoldi_data)
    IF(ASSOCIATED(ar_data%f_vec))DEALLOCATE(ar_data%f_vec)
    IF(ASSOCIATED(ar_data%Hessenberg))DEALLOCATE(ar_data%Hessenberg)
    IF(ASSOCIATED(ar_data%local_history))DEALLOCATE(ar_data%local_history)
    IF(ASSOCIATED(ar_data%evals))DEALLOCATE(ar_data%evals)
    IF(ASSOCIATED(ar_data%revec))DEALLOCATE(ar_data%revec)
    DEALLOCATE(ar_data)

  END SUBROUTINE deallocate_arnoldi_data_z
