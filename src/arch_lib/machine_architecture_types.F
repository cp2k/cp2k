!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2012 Christiane Ribeiro and the CP2K developers group!
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Machine Hardware support for cp2k
!           Data structures for machine representation 
!> \author  Christiane Pousa Ribeiro
!> \date    2011-08-08
!> \version 1.0
!>
!> <b>Modification history:</b>
!> - Created 2011-08-09
! *****************************************************************************
MODULE machine_architecture_types
#if !defined (__HAS_NO_ISO_C_BINDING)
  USE ISO_C_BINDING
#endif  
  USE ma_dummy_c_bindings

 IMPLICIT NONE
  
 PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'machine_architecture_types' 

 PUBLIC :: arch_topology, thread_inf, machine_output

!Interface with the Machine Abstraction  
#if !defined (__HAS_NO_ISO_C_BINDING)
 TYPE, BIND(C) :: arch_topology 
#else
 TYPE :: arch_topology 
#endif
   INTEGER(C_INT)         :: nnodes, nsockets
   INTEGER(C_INT)         :: ncores, npus
   INTEGER(C_INT)         :: ngpus, ncaches
   INTEGER(C_INT)         :: nshared_caches, nsiblings
   INTEGER(C_INT)          :: nmemcontroller, nnetcards                 
  END TYPE arch_topology

!Interface with the console output
#if !defined (__HAS_NO_ISO_C_BINDING)
 TYPE, BIND(C) :: machine_output
   CHARACTER(C_CHAR)           :: console_output(8192)
#else
 TYPE :: machine_output
   CHARACTER                   :: console_output(8192)
#endif
   INTEGER(C_INT)              :: len
 END TYPE machine_output  

!Information of the threads
 TYPE :: thread_inf
   INTEGER                :: core, node, id_omp, id_real 
 END TYPE thread_inf 


 ! For machine architecture internal use 
 PUBLIC                   :: threads, mpi
 INTEGER, PARAMETER       :: threads = 0
 INTEGER, PARAMETER       :: mpi     = 1 

 ! For manual setting of all affinity
 PUBLIC :: manual
 INTEGER, PARAMETER       :: manual = -1

 ! MPI/Thread scheduling policies
 PUBLIC  :: def, nosched, linear, scatter, group
 
 INTEGER, PARAMETER       :: nosched = 0
 INTEGER, PARAMETER       :: linear  = 1
 INTEGER, PARAMETER       :: scatter = 2
 INTEGER, PARAMETER       :: group = 3 
 INTEGER, PARAMETER       :: def = 4 


!Information of the mempory
 TYPE :: mem_inf
   INTEGER                :: mempol, node
 END TYPE mem_inf

 ! For machine architecture internal use
 ! Memory binding policies
 PUBLIC  :: none_pol, os, local, interleave

 INTEGER, PARAMETER       :: os         = 0
 INTEGER, PARAMETER       :: local      = 1
 INTEGER, PARAMETER       :: interleave = 2
 INTEGER, PARAMETER       :: none_pol   = 4


#if defined (__HWLOC) 
  LOGICAL,  PUBLIC, PARAMETER        :: has_ma = .TRUE.
#elif defined (__LIBNUMA)
  LOGICAL,  PUBLIC, PARAMETER        :: has_ma = .TRUE.
#else
  LOGICAL,  PUBLIC, PARAMETER        :: has_ma = .FALSE.
#endif
#if defined (__parallel)
  LOGICAL, PARAMETER, PUBLIC                   :: has_MPI = .TRUE.
#else
  LOGICAL, PARAMETER, PUBLIC                   :: has_MPI = .FALSE.
#endif


 PUBLIC :: ma_mp_type, ma_process
 
! *****************************************************************************
!> \brief A processor (process) grid distribution
!> \var myproc                my processor (process) number
!> \var numproc               number of processors (processes)
!> \var mp_group              message-passing group ID
!> \var refcount              reference counter
! *****************************************************************************
  TYPE ma_mp_type
     INTEGER                            :: myproc, numproc
     INTEGER                            :: mp_group
     INTEGER                            :: refcount
     INTEGER                            :: myid
     LOGICAL                            :: all_proc
  END TYPE ma_mp_type

! *****************************************************************************
!> \brief The process information 
!> \var mp_info            my processor information of the parallel environment
!> \var threads_info       my threads information
!> \var nr_threads         my number of threads
!> \var core               core where the process run
!> \var node               NUMA node where the process run  
! *****************************************************************************
  TYPE ma_process
     TYPE(ma_mp_type)                            :: mp_info
     TYPE(thread_inf), DIMENSION(:), ALLOCATABLE :: threads_info
     INTEGER                                     :: nr_threads
     INTEGER                                     :: nr_neighbors 
     INTEGER                                     :: core, node
     INTEGER                                     :: gpu
  END TYPE ma_process

END MODULE machine_architecture_types
