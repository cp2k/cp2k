!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2013  Urban Borstnik and the CP2K developers group   !
!-----------------------------------------------------------------------------!

  PURE FUNCTION dev_mem_get_alloc_z (dev_mem) RESULT (data_size)
    TYPE(dbcsr_cuda_mem_type_c8), INTENT(IN) :: dev_mem
    INTEGER                                  :: data_size

    data_size = dev_mem%len
  END FUNCTION dev_mem_get_alloc_z

  PURE FUNCTION dev_mem_get_size_z (dev_mem) RESULT (data_size)
    TYPE(dbcsr_cuda_mem_type_c8), INTENT(IN) :: dev_mem
    INTEGER                                  :: data_size

    data_size = dev_mem%req_len
  END FUNCTION dev_mem_get_size_z
