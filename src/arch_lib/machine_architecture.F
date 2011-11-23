!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2011 Christiane Ribeiro and the CP2K developers group!
!-----------------------------------------------------------------------------!
!> \brief   Machine Hardware support for cp2k
!           Interface to the HWLOC/LIBNUMA library
!> \author  Christiane Pousa Ribeiro
!> \date    2011-08-08
!> \version 1.0
!>
!> <b>Modification history:</b>
!> - Created 2011-08-08
!> - Modified 2011-08-09 - Adding interfaces to C functions
MODULE machine_architecture

#if defined (__HAS_NO_ISO_C_BINDING)
#if defined (__HWLOC) || defined(__LIBNUMA)
  !Compile-time error. These two flags are incompatible
  STOP "The machine abstraction needs ISO_C_BINDING. Compile without __HAS_NO_ISO_C_BINDING"  
#endif
#endif

#if !defined (__HAS_NO_ISO_C_BINDING)
  USE ISO_C_BINDING
#endif 

  USE machine_architecture_types
  USE machine_architecture_utils
  USE ma_dummy_c_bindings
  USE ma_kinds 
  USE ma_message_passing,           ONLY: mp_allgather,mp_sum, mp_proc_name
  USE ma_errors,                    ONLY: ma_error_stop, ma_error_allocation, &
                                          ma_error_finalize
  USE machine, ONLY: m_getpid
 
  IMPLICIT NONE

  PRIVATE 

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'machine_architecture'

  PUBLIC :: ma_init_machine, ma_finalize_machine

  PUBLIC :: ma_get_id

  PUBLIC :: ma_get_nnodes, ma_get_nsockets, ma_get_ncores, ma_get_ncaches 
  PUBLIC :: ma_get_nshared_caches, ma_get_nsiblings
  
  PUBLIC :: ma_get_proc_core
  PUBLIC :: ma_get_mynode, ma_get_mycore
  PUBLIC :: ma_get_proc_node, ma_set_proc_node
  
  PUBLIC :: ma_show_topology, ma_show_machine_branch
  PUBLIC :: ma_show_machine_full

#if defined (__HWLOC)
  PUBLIC :: ma_hw_set_mempol, ma_hw_get_mempol
  PUBLIC :: ma_get_thread_id, ma_set_proc_core
  PUBLIC :: ma_set_core, ma_get_core
  PUBLIC :: ma_set_proc_cores, ma_set_procnode
  PUBLIC :: ma_get_core_node, ma_set_thread_cores
  PUBLIC :: ma_get_cuDev, ma_get_gpu_node
  PUBLIC :: ma_get_netDev
#endif

#if defined (__LIBNUMA)
  PUBLIC :: ma_linux_set_mempol, ma_linux_get_mempol
  PUBLIC :: ma_linux_get_core, ma_linux_set_core
  PUBLIC :: ma_linux_set_proccore, ma_linux_get_proccore
  PUBLIC :: ma_linux_get_thread_id
  PUBLIC :: ma_linux_ncores, ma_linux_get_mynode
  PUBLIC :: ma_set_proc_cores, ma_linux_set_procnode
  PUBLIC :: ma_get_cuDev, ma_get_gpu_node
  PUBLIC :: ma_linux_get_core_node, ma_linux_nnodes
#endif

  ! These are for Machine architecture internal use.
  !
  PUBLIC :: topology 
  PUBLIC :: console_out
  PUBLIC :: ma_hwloc, ma_libnuma
  PUBLIC :: ma_interface, ma_int_hwloc, ma_int_libnuma, ma_int_none
  
  ! For the machine representation
  TYPE(arch_topology)                            :: topology 
  TYPE(machine_output)                           :: console_out

#if defined (__LIBNUMA)
  LOGICAL, PARAMETER :: ma_libnuma = .TRUE.
#else
  LOGICAL, PARAMETER :: ma_libnuma = .FALSE.
#endif
#if defined (__HWLOC)
  LOGICAL, PARAMETER :: ma_hwloc = .TRUE.
#else
  LOGICAL, PARAMETER :: ma_hwloc = .FALSE.
#endif

  CHARACTER, PARAMETER :: ma_int_hwloc = 'H'
  CHARACTER, PARAMETER :: ma_int_libnuma = 'N'
  CHARACTER, PARAMETER :: ma_int_none = '0'
  CHARACTER, SAVE :: ma_interface = ma_int_none

!
! Interfaces to deal with HWLOC
!
#if defined (__HWLOC) 
  INTERFACE 
     FUNCTION ma_init_interface(topo) RESULT (istat) BIND(C, name="hw_topology_init") 
       USE ISO_C_BINDING
       USE machine_architecture_types,   ONLY: arch_topology 
    TYPE(arch_topology)                      :: topo
    INTEGER(KIND=C_INT)                      :: istat

     END FUNCTION ma_init_interface
  END INTERFACE
 
 INTERFACE 
     FUNCTION ma_final_interface(dtopo) RESULT (istat) BIND(C, name="hw_topology_destroy") 
       USE ISO_C_BINDING
       USE machine_architecture_types,   ONLY: arch_topology 
    TYPE(arch_topology), BIND(C)             :: dtopo
    INTEGER(KIND=C_INT)                      :: istat

     END FUNCTION ma_final_interface
  END INTERFACE

 INTERFACE 
     SUBROUTINE ma_show_phys_topo(output_phys) BIND(C, name="hw_phys_pu_topology") 
       USE ISO_C_BINDING
       USE machine_architecture_types,  ONLY: machine_output
    TYPE(machine_output), BIND(C)            :: output_phys

     END SUBROUTINE ma_show_phys_topo
  END INTERFACE


  INTERFACE 
     SUBROUTINE ma_show_branch(output_branch) BIND(C, name="hw_high_level_show") 
       USE ISO_C_BINDING
       USE machine_architecture_types,  ONLY: machine_output
    TYPE(machine_output), BIND(C)            :: output_branch

     END SUBROUTINE ma_show_branch
  END INTERFACE

  INTERFACE 
     SUBROUTINE ma_show_topo(output_topo) BIND(C, name="hw_machine_topology") 
       USE ISO_C_BINDING
       USE machine_architecture_types,  ONLY: machine_output
    TYPE(machine_output), BIND(C)            :: output_topo

     END SUBROUTINE ma_show_topo
  END INTERFACE

  INTERFACE 
    FUNCTION ma_get_cuDev(coreId, myRank, nMPIs) RESULT (cuDev) BIND(C, name="hw_my_gpu") 
    USE ISO_C_BINDING
    INTEGER(KIND=C_INT), VALUE               :: coreId, myRank, nMPIs
    INTEGER(KIND=C_INT)                      :: cuDev

    END FUNCTION ma_get_cuDev
  END INTERFACE 

 INTERFACE 
  FUNCTION ma_get_netDev(coreId) RESULT (netDev) BIND(C, name="hw_my_netCard") 
    USE ISO_C_BINDING
    INTEGER(KIND=C_INT), VALUE               :: coreId
    INTEGER(KIND=C_INT)                      :: netDev

    END FUNCTION ma_get_netDev
  END INTERFACE 

  INTERFACE 
     FUNCTION ma_get_node() RESULT (node) BIND(C, name="hw_my_node") 
       USE ISO_C_BINDING
    INTEGER(KIND=C_INT)                      :: node

     END FUNCTION ma_get_node
  END INTERFACE

  INTERFACE
     FUNCTION ma_get_core_node(core) RESULT (node) BIND(C,name="hw_get_myNode")
       USE ISO_C_BINDING
    INTEGER(KIND=C_INT), VALUE               :: core
    INTEGER(KIND=C_INT)                      :: node

     END FUNCTION ma_get_core_node
  END INTERFACE                     

  INTERFACE 
     FUNCTION ma_get_procnode() RESULT (node) BIND(C, name="hw_get_proc_node") 
       USE ISO_C_BINDING
    INTEGER(KIND=C_INT)                      :: node

     END FUNCTION ma_get_procnode
  END INTERFACE

  INTERFACE
    SUBROUTINE ma_set_procnode(node) BIND(C, name="hw_set_proc_node")
       USE ISO_C_BINDING
    INTEGER(KIND=C_INT), VALUE               :: node

    END SUBROUTINE ma_set_procnode
  END INTERFACE

 INTERFACE
    SUBROUTINE ma_set_proc_cores(dist, core) BIND(C, name="hw_set_proc_cores")
       USE ISO_C_BINDING
    INTEGER(KIND=C_INT), VALUE               :: dist, core

    END SUBROUTINE ma_set_proc_cores
  END INTERFACE

  INTERFACE
    SUBROUTINE ma_set_thread_cores(dist, core) BIND(C, name="hw_set_thread_cores")
       USE ISO_C_BINDING
    INTEGER(KIND=C_INT), VALUE               :: dist, core

    END SUBROUTINE ma_set_thread_cores
  END INTERFACE

  INTERFACE 
      SUBROUTINE ma_set_core(cpu) BIND(C, name="hw_set_my_core") 
       USE ISO_C_BINDING
    INTEGER(KIND=C_INT), VALUE               :: cpu

     END SUBROUTINE ma_set_core
  END INTERFACE

  INTERFACE 
      FUNCTION ma_get_core() RESULT (core) BIND(C, name="hw_my_core") 
       USE ISO_C_BINDING
    INTEGER(KIND=C_INT)                      :: core

     END FUNCTION ma_get_core
  END INTERFACE

  INTERFACE 
      SUBROUTINE ma_set_proc_core(cpu) BIND(C, name="hw_set_proc_core") 
       USE ISO_C_BINDING
    INTEGER(KIND=C_INT), VALUE               :: cpu

     END SUBROUTINE ma_set_proc_core
  END INTERFACE

  INTERFACE 
      FUNCTION ma_get_proccore(flag) RESULT (core) BIND(C, name="hw_get_proc_core") 
       USE ISO_C_BINDING
    INTEGER(KIND=C_INT), VALUE               :: flag
    INTEGER(KIND=C_INT)                      :: core

     END FUNCTION ma_get_proccore
  END INTERFACE

  INTERFACE 
      FUNCTION ma_get_thread_id() RESULT (id) BIND(C, name="hw_get_myid") 
       USE ISO_C_BINDING
    INTEGER(KIND=C_INT)                      :: id

     END FUNCTION ma_get_thread_id
  END INTERFACE

  INTERFACE
     FUNCTION ma_hw_set_mempol(mempol, node) RESULT (error) BIND(C, name="hw_set_mempol")
       USE ISO_C_BINDING
    INTEGER(KIND=C_INT), VALUE               :: mempol, node
    INTEGER(KIND=C_INT)                      :: error

     END FUNCTION ma_hw_set_mempol
  END INTERFACE 

  INTERFACE
      SUBROUTINE ma_hw_get_mempol(node, mempol) BIND(C, name="hw_get_mempol")
       USE ISO_C_BINDING
    INTEGER(KIND=C_INT)                      :: node, mempol

     END SUBROUTINE ma_hw_get_mempol
  END INTERFACE

  INTERFACE
   FUNCTION ma_get_gpu_node (gpu) RESULT (node)  BIND(C, name="hw_get_gpu_node")
    USE ISO_C_BINDING 
    INTEGER(KIND=C_INT), VALUE               :: gpu
    INTEGER(KIND=C_INT)                      :: node

   END FUNCTION ma_get_gpu_node
  END INTERFACE

#endif

!
! Interfaces to deal with LIBNUMA
!
#if defined (__LIBNUMA)
  INTERFACE
     FUNCTION ma_init_interface_linux(topo) RESULT (istat) BIND(C, name="linux_topology_init")
       USE ISO_C_BINDING
       USE machine_architecture_types,   ONLY: arch_topology
    TYPE(arch_topology)                      :: topo
    INTEGER(KIND=C_INT)                      :: istat

     END FUNCTION ma_init_interface_linux
  END INTERFACE

  INTERFACE
    FUNCTION ma_get_cuDev(coreId, myRank, numMPI) RESULT (cuDev) BIND(C,name="linux_my_gpu")
    USE ISO_C_BINDING
    INTEGER(KIND=C_INT), VALUE               :: coreId, myRank, numMPI
    INTEGER(KIND=C_INT)                      :: cuDev

    END FUNCTION ma_get_cuDev
  END INTERFACE 

  INTERFACE
      SUBROUTINE ma_linux_get_mempol(node, mempol) BIND(C, name="linux_get_mempol")
       USE ISO_C_BINDING
    INTEGER(KIND=C_INT)                      :: node, mempol

     END SUBROUTINE ma_linux_get_mempol
  END INTERFACE

  INTERFACE
     FUNCTION ma_linux_set_mempol(mempol,node) RESULT (error) BIND(C, name="linux_set_mempol")
       USE ISO_C_BINDING
    INTEGER(KIND=C_INT), VALUE               :: mempol, node
    INTEGER(KIND=C_INT)                      :: error

     END FUNCTION ma_linux_set_mempol
  END INTERFACE
  
  INTERFACE
      FUNCTION ma_linux_get_core() RESULT (core) BIND(C, name="linux_my_core")
       USE ISO_C_BINDING
    INTEGER(KIND=C_INT)                      :: core

     END FUNCTION ma_linux_get_core
  END INTERFACE

  INTERFACE
      SUBROUTINE ma_linux_set_core(cpu) BIND(C, name="linux_set_my_core")
       USE ISO_C_BINDING
    INTEGER(KIND=C_INT), VALUE               :: cpu

     END SUBROUTINE ma_linux_set_core
  END INTERFACE
  
  INTERFACE
      SUBROUTINE ma_linux_set_proccore(cpu) BIND(C, name="linux_set_proc_core")
       USE ISO_C_BINDING
    INTEGER(KIND=C_INT), VALUE               :: cpu

     END SUBROUTINE ma_linux_set_proccore
  END INTERFACE

  INTERFACE
      FUNCTION ma_linux_get_proccore() RESULT (core) BIND(C, name="linux_proc_core")
       USE ISO_C_BINDING
    INTEGER(KIND=C_INT)                      :: core

     END FUNCTION ma_linux_get_proccore
  END INTERFACE  

  INTERFACE
      FUNCTION ma_linux_get_thread_id() RESULT (id) BIND(C, name="linux_get_myid")
       USE ISO_C_BINDING
    INTEGER(KIND=C_INT)                      :: id

     END FUNCTION ma_linux_get_thread_id
  END INTERFACE

  INTERFACE
      FUNCTION ma_linux_nnodes() RESULT (nnodes) BIND(C, name="linux_get_nnodes")
       USE ISO_C_BINDING
    INTEGER(KIND=C_INT)                      :: nnodes

     END FUNCTION ma_linux_nnodes
  END INTERFACE

  INTERFACE
      FUNCTION ma_linux_ncores() RESULT (ncores) BIND(C, name="linux_get_ncores")
       USE ISO_C_BINDING
    INTEGER(KIND=C_INT)                      :: ncores

     END FUNCTION ma_linux_ncores
  END INTERFACE

  INTERFACE
      FUNCTION ma_linux_get_mynode() RESULT (node) BIND(C, name="linux_get_nodeid")
       USE ISO_C_BINDING
    INTEGER(KIND=C_INT)                      :: node

     END FUNCTION ma_linux_get_mynode
  END INTERFACE


 INTERFACE
      FUNCTION ma_linux_get_core_node(core) RESULT (node) BIND(C, name="linux_get_nodeid_cpu")
       USE ISO_C_BINDING
    INTEGER(KIND=C_INT), VALUE               :: core
    INTEGER(KIND=C_INT)                      :: node

     END FUNCTION ma_linux_get_core_node
  END INTERFACE

 INTERFACE


    SUBROUTINE ma_set_proc_cores(dist, core) BIND(C, name="linux_set_proc_cores")
       USE ISO_C_BINDING
    INTEGER(KIND=C_INT), VALUE               :: dist, core

    END SUBROUTINE ma_set_proc_cores
  END INTERFACE
 
  INTERFACE
    SUBROUTINE ma_linux_set_procnode(node) BIND(C, name="linux_set_proc_node")
       USE ISO_C_BINDING
    INTEGER(KIND=C_INT), VALUE               :: node

    END SUBROUTINE ma_linux_set_procnode
  END INTERFACE

  INTERFACE
   FUNCTION ma_get_gpu_node (gpu) RESULT (node)  BIND(C, name="linux_get_gpu_node")
    USE ISO_C_BINDING 
    INTEGER(KIND=C_INT), VALUE               :: gpu
    INTEGER(KIND=C_INT)                      :: node

   END FUNCTION ma_get_gpu_node
  END INTERFACE

#endif


CONTAINS

! *****************************************************************************
!> \brief Initialize the machine abstraction library
!> \var ma_env            The parallel environment
!> \var num_threads       The number of threads of each process
! *****************************************************************************
  INTEGER FUNCTION ma_init_machine(core,node) RESULT(istat)
  INTEGER                          :: core,node

   istat = -1
#if defined (__HWLOC) || defined (__LIBNUMA)
   topology%nnodes = -1
#if defined (__HWLOC)  && !defined (__LIBNUMA) 
   istat = ma_init_interface(topology)
   core = ma_get_mycore()
   node = ma_get_mynode()
#elif defined (__LIBNUMA) && !defined (__HWLOC)
   istat = ma_init_interface_linux(topology) 
   core = ma_linux_get_core()
   node = ma_linux_get_mynode()
#endif 
#endif
  END FUNCTION ma_init_machine 

! *****************************************************************************
!> \brief Finalize the machine abstraction library
! *****************************************************************************
   SUBROUTINE ma_finalize_machine()
   INTEGER  :: istat
#if defined (__HWLOC) && !defined (__LIBNUMA)
   istat = ma_final_interface(topology)
   IF( istat .NE. 0) THEN
     CALL ma_error_stop(ma_error_finalize)
   ENDIF
#endif 
  END SUBROUTINE ma_finalize_machine

! *****************************************************************************
!> \brief Get the process ID
! *****************************************************************************
INTEGER FUNCTION ma_get_id() RESULT (pid)
  INTEGER                                  :: getpid

  !pid = getpid()
  CALL m_getpid(pid)
END FUNCTION ma_get_id

! *****************************************************************************
!> \brief Prints in text format the machine abstraction - only one branch
!> \var unit_num the output console
! *****************************************************************************
  SUBROUTINE ma_show_machine_branch(unit_num)
  INTEGER :: unit_num
#if defined (__HWLOC)
  
  CALL ma_show_branch(console_out)
  WRITE(unit_num,*) console_out%console_output(1:console_out%len)
#endif
  END SUBROUTINE ma_show_machine_branch


! *****************************************************************************
!> \brief Prints in text format the machine abstraction - the full machine node
!> \var unit_num the output console
! *****************************************************************************
   SUBROUTINE ma_show_machine_full(unit_num)
   INTEGER :: unit_num
#if defined (__HWLOC) 
   CALL ma_show_phys_topo(console_out)
   WRITE(unit_num,*)console_out%console_output(1:console_out%len)
   CALL ma_show_topo(console_out)
   WRITE(unit_num,*)console_out%console_output(1:console_out%len)
#endif 
  END SUBROUTINE ma_show_machine_full

! *****************************************************************************
! Functions to extract the machine information
! *****************************************************************************
  INTEGER FUNCTION ma_get_mycore() RESULT (core)
   core = -1  
#if defined (__HWLOC)  && !defined (__LIBNUMA)
   core = ma_get_core()
#endif 
#if defined (__LIBNUMA) && !defined (__HWLOC)
  core = ma_linux_get_core()
#endif
  END FUNCTION ma_get_mycore 

  FUNCTION ma_get_mynode() RESULT (node)
    INTEGER                                  :: node

    node  = -1 
#if defined (__HWLOC) && !defined (__LIBNUMA)
   node = ma_get_node()
#endif 
#if defined (__LIBNUMA) && !defined (__HWLOC)
    node = ma_linux_get_mynode()  
#endif
  END FUNCTION ma_get_mynode

  INTEGER FUNCTION ma_get_proc_core(flag) RESULT (core)
   INTEGER  :: flag
   core = -1  
#if defined (__HWLOC)  && !defined (__LIBNUMA)
   core = ma_get_proccore(flag)
#endif 
#if defined (__LIBNUMA) && !defined (__HWLOC)
   core = ma_linux_get_proccore() 
#endif
  END FUNCTION ma_get_proc_core 

  FUNCTION ma_get_proc_node() RESULT (node)
    INTEGER                                  :: node

    node  = -1
#if defined (__HWLOC) && !defined (__LIBNUMA)
   node = ma_get_procnode()
#endif 
#if defined (__LIBNUMA) && !defined (__HWLOC)
   node = ma_linux_get_mynode()
#endif 
  END FUNCTION ma_get_proc_node

SUBROUTINE ma_set_proc_node(node) 
    INTEGER                                  :: node
#if defined (__HWLOC) && !defined (__LIBNUMA)
   CALL ma_set_procnode(node)
#endif 
#if defined (__LIBNUMA) && !defined (__HWLOC)
   CALL ma_linux_set_procnode(node)
#endif 
  END SUBROUTINE ma_set_proc_node

 INTEGER FUNCTION ma_get_nnodes()
     ma_get_nnodes = topology%nnodes 
  END FUNCTION ma_get_nnodes 

 INTEGER FUNCTION ma_get_nsockets()
     ma_get_nsockets = topology%nsockets 
  END FUNCTION ma_get_nsockets 

  INTEGER FUNCTION ma_get_ncores()
     ma_get_ncores = topology%ncores 
  END FUNCTION ma_get_ncores

 INTEGER FUNCTION ma_get_ncaches()
     ma_get_ncaches = topology%ncaches 
  END FUNCTION ma_get_ncaches

 INTEGER FUNCTION ma_get_nshared_caches()
     ma_get_nshared_caches = topology%nshared_caches 
  END FUNCTION ma_get_nshared_caches

 INTEGER FUNCTION ma_get_nsiblings()
     ma_get_nsiblings = topology%nsiblings 
  END FUNCTION ma_get_nsiblings

! *****************************************************************************
!> \brief Prints in text format the machine abstraction - the full machine node
!> \var unit_num the output console
! *****************************************************************************
  SUBROUTINE ma_show_topology (unit_num)
    INTEGER                                  :: unit_num

    IF (topology%nnodes .NE. -1) THEN       
       
     WRITE(unit_num, "(T2,A)")     "MACHINE| Machine Architecture"
    
     IF (topology%nnodes .EQ. 0) THEN       
       WRITE(unit_num, "(T2,A)")  "MACHINE| No NUMA nodes"
     ELSE
      WRITE(unit_num, "(T2,A,27X,I3)")    "MACHINE| Number of NUMA nodes                    ",&
                                topology%nnodes
     ENDIF
     WRITE(unit_num, "(T2,A,27X,I3)")      "MACHINE| Number of sockets                       ",&
                                topology%nsockets
     WRITE(unit_num, "(T2,A,27X,I3)")      "MACHINE| Number of memory controllers per socket ",&
                                topology%nmemcontroller; 
     WRITE(unit_num, "(T2,A,27X,I3)")      "MACHINE| Number of cores                         ",&
                                topology%ncores
     WRITE(unit_num, "(T2,A,27X,I3)")      "MACHINE| Number of cores per socket              ",&
                                topology%ncores/topology%nsockets
     WRITE(unit_num, "(T2,A,27X,I3)")      "MACHINE| Number of siblings per core             ",&
                                topology%nsiblings
#if defined (__CUDAPW) || defined (__DBCSR_CUDA)
     WRITE(unit_num, "(T2,A,27X,I3)")      "MACHINE| Number of graphics processing units     ",&
                                topology%ngpus
#endif
     WRITE(unit_num, "(T2,A,27X,I3)")      "MACHINE| Number of cache levels                  ",&
                                topology%ncaches
     WRITE(unit_num, "(T2,A,27X,I3)")      "MACHINE| Number of shared caches                 ",&
                                topology%nshared_caches
     WRITE(unit_num, "()")     
    ENDIF
  END SUBROUTINE ma_show_topology

END MODULE machine_architecture 
