!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief type to store parallelization informations (at the moment assumes 1d
!>      position and uses mpi).
!>      type def in this separate module to avoid circular ref.
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
MODULE cp_para_types
  

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_para_types'

  PUBLIC :: cp_para_env_type, cp_para_env_p_type, cp_para_cart_type

! *****************************************************************************
!> \brief stores all the informations relevant to an mpi environement
!> \param mepos rank of the actual processor
!> \param group id of the actual mpi communicator (I know that it is confusing,
!>        but this is not the mpi group)
!> \param num_pe number of processors in the communicator
!> \param source rank of a special processor (for example the one for i-o,
!>        or the master
!> \param owns_group if it owns the group (and thus should free it when
!>        this object is deallocated)
!> \param ref_count the reference count, when it is zero this object gets
!>        deallocated
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  TYPE cp_para_env_type
     LOGICAL :: owns_group, ionode
     INTEGER :: mepos, source, num_pe, ref_count
     INTEGER :: group
  END TYPE cp_para_env_type

! *****************************************************************************
!> \brief represent a pointer to a para env (to build arrays)
!> \param para_env the pointer to the para_env
!> \par History
!>      07.2003 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  TYPE cp_para_env_p_type
     TYPE(cp_para_env_type), POINTER :: para_env
  END TYPE cp_para_env_p_type

! *****************************************************************************
!> \brief represent a multidimensional parallel environment
!> \param mepos the position of the actual processor
!> \param group id of the actual group (context, communicator)
!> \param num_pe number of processors in the group in each dimension
!> \param source id of a special processor (for example the one for i-o,
!>        or the master
!> \param owns_group if it owns the group (and thus should free it when
!>        this object is deallocated)
!> \param ref_count the reference count, when it is zero this object gets
!>        deallocated
!> \note
!>      not yet implemented for mpi
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  TYPE cp_para_cart_type
     LOGICAL :: owns_group
     INTEGER :: ndims, rank, ntask
     INTEGER, DIMENSION(:), POINTER :: mepos, source, num_pe
     LOGICAL, DIMENSION(:), POINTER :: periodic
     INTEGER :: group, ref_count
  END TYPE cp_para_cart_type

END MODULE cp_para_types
