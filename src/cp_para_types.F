!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/cp_para_types [1.0] *
!!
!!   NAME
!!     cp_para_types
!!
!!   FUNCTION
!!     type to store parallelization informations (at the moment assumes 1d
!!     position and uses mpi).
!!     type def in this separate module to avoid circular ref.
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE cp_para_types
  USE kinds,                           ONLY: dp
  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='cp_para_types'

  PUBLIC :: cp_para_env_type, cp_blacs_env_type, cp_para_env_p_type
!***
!****************************************************************************

!!****s* cp_para_types/cp_para_env_type [1.0] *
!!
!!   NAME
!!     cp_para_env_type
!!
!!   FUNCTION
!!     stores all the informations relevant to an mpi environement
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - mepos: rank of the actual processor
!!     - group: id of the actual mpi communicator (I know that it is confusing, 
!!       but this is not the mpi group)
!!     - num_pe: number of processors in the communicator
!!     - source: rank of a special processor (for example the one for i-o,
!!       or the master
!!     - owns_group: if it owns the group (and thus should free it when
!!       this object is deallocated)
!!     - ref_count: the reference count, when it is zero this object gets
!!       deallocated
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!   SOURCE
  !***************************************************************************
  TYPE cp_para_env_type
     LOGICAL :: owns_group
     INTEGER :: mepos, group, source, num_pe, ref_count
  END TYPE cp_para_env_type
!!***
!****************************************************************************

!!****s* cp_para_types/cp_para_env_p_type [1.0] *
!!
!!   NAME
!!     cp_para_env_p_type
!!
!!   FUNCTION
!!     represent a pointer to a para env (to build arrays)
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - para_env: the pointer to the para_env
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     07.2003 created [fawzi]
!!
!!   SOURCE
  !***************************************************************************
  TYPE cp_para_env_p_type
     TYPE(cp_para_env_type), POINTER :: para_env
  END TYPE cp_para_env_p_type
!!***
  !****************************************************************************

!!****s* cp_para_types/cp_para_cart_type [1.0] *
!!
!!   NAME
!!     cp_para_cart_type
!!
!!   FUNCTION
!!     represent a multidimensional parallel environment
!!
!!   NOTES
!!     not yet implemented for mpi
!!
!!   ATTRIBUTES
!!     - mepos: the position of the actual processor
!!     - group: id of the actual group (context, communicator)
!!     - num_pe: number of processors in the group in each dimension
!!     - source: id of a special processor (for example the one for i-o,
!!       or the master
!!     - owns_group: if it owns the group (and thus should free it when
!!       this object is deallocated)
!!     - ref_count: the reference count, when it is zero this object gets
!!       deallocated
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!   SOURCE
  !***************************************************************************
  TYPE cp_para_cart_type
     LOGICAL :: owns_group
     INTEGER, DIMENSION(:), POINTER :: mepos, source, num_pe
     INTEGER :: group, kind, ref_count
  END TYPE cp_para_cart_type
!!***
  !****************************************************************************

!!****s* cp_para_types/cp_blacs_env_type [1.0] *
!!
!!   NAME
!!     cp_blacs_env_type
!!
!!   FUNCTION
!!     represent a blacs multidimensional parallel environment
!!     (for the mpi corrispective see cp_paratypes/cp_para_cart_type)
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - mepos: the position of the actual processor
!!     - group: id of the actual group (context, communicator)
!!     - num_pe: number of processors in the group in each dimension
!!     - source: id of a special processor (for example the one for i-o,
!!       or the master
!!     - owns_group: if it owns the group (and thus should free it when
!!       this object is deallocated)
!!     - ref_count: the reference count, when it is zero this object gets
!!       deallocated
!!     - my_pid: process id of the actual processor (global)
!!     - n_pid: number of process ids (this is always global, not referred
!!       to the actual number of processors in the blacs grid)
!!     - the para env associated (and compatible) with this blacs_env
!!     - blacs2mpi: maps mepos(1)-mepos(2) of blacs to its mpi rank
!!     - mpi2blacs(i,rank): maps the mpi rank to the mepos(i)
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!   SOURCE
  !***************************************************************************
  TYPE cp_blacs_env_type
     LOGICAL :: owns_group
     INTEGER, DIMENSION(2) :: mepos, source, num_pe
     INTEGER :: group, my_pid, n_pid, ref_count
     TYPE(cp_para_env_type), POINTER :: para_env
     INTEGER, DIMENSION(:,:), POINTER :: blacs2mpi
     INTEGER, DIMENSION(:,:), POINTER :: mpi2blacs
  END TYPE cp_blacs_env_type
!!***
  !****************************************************************************

END MODULE cp_para_types
