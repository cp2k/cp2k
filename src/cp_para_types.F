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
module cp_para_types
  use kinds, only: wp=>dp
  implicit none
  private

  logical, private, parameter :: debug_this_module=.true.
  character(len=*), private, parameter :: moduleN='cp_para_types'

  public :: cp_para_env_type, cp_blacs_env_type
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
!!     - mepos: id of the actual processor
!!     - group: id of the actual mpi group
!!     - num_pe: number of processors in the group
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
  type cp_para_env_type
     logical :: owns_group
     integer :: mepos, group, source, num_pe, ref_count
  end type cp_para_env_type
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
  type cp_para_cart_type
     logical :: owns_group
     integer, dimension(:), pointer :: mepos, source, num_pe
     integer :: group, kind, ref_count
  end type cp_para_cart_type
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
  !!
  !!   AUTHOR
  !!     Fawzi Mohamed
  !!
  !!   MODIFICATION HISTORY
  !!     08.2002 created [fawzi]
  !!
  !!   SOURCE
  !***************************************************************************
  type cp_blacs_env_type
     logical :: owns_group
     integer, dimension(2) :: mepos, source, num_pe
     integer :: group, my_pid, n_pid, ref_count
  end type cp_blacs_env_type
  !!***
  !****************************************************************************

end module cp_para_types
