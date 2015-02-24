!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Machine Hardware support for cp2k
!           Interface to the HWLOC/LIBNUMA library
!> \author  Christiane Pousa Ribeiro
!> \date    2011-08-08
!> \version 1.0
!>
!> <b>Modification history:</b>
!> - Created 2011-08-08
!> - Modified 2011-08-09 - Adding interfaces to C functions
! *****************************************************************************
MODULE machine_architecture

  USE ISO_C_BINDING
  USE ma_errors,                       ONLY: ma_error_allocation,&
                                             ma_error_finalize,&
                                             ma_error_stop
  USE ma_kinds,                        ONLY: default_string_length
  USE machine,                         ONLY: m_getpid
  USE machine_architecture_types,      ONLY: arch_topology,&
                                             group,&
                                             ma_mp_type,&
                                             machine_output
  USE machine_architecture_utils,      ONLY: ascii_to_string,&
                                             string_to_ascii
  USE message_passing,                 ONLY: mp_environ,&
                                             mp_proc_name,&
                                             mp_sum

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'machine_architecture'

  PUBLIC :: ma_init_machine, ma_finalize_machine

  PUBLIC :: ma_get_id, ma_get_thread_id

  PUBLIC :: ma_get_nnodes, ma_get_ncores
  PUBLIC :: ma_get_nmachines

  PUBLIC :: ma_get_proc_core, ma_my_first_core
  PUBLIC :: ma_get_mynode, ma_get_mycore
  PUBLIC :: ma_get_proc_node, ma_set_proc_node

  PUBLIC :: ma_show_topology, ma_show_machine_branch
  PUBLIC :: ma_show_machine_full

  PUBLIC :: ma_get_netDev, ma_get_nnetDev
  PUBLIC :: ma_get_node_netDev, ma_get_core_node
  PUBLIC :: ma_hw_set_mempol

  PUBLIC :: ma_set_proc_cores, ma_set_thread_cores
  PUBLIC :: ma_set_proc_core, ma_set_core
  PUBLIC :: ma_hw_get_mempol
  PUBLIC :: ma_set_thread_allnodes, ma_set_thread_node
  PUBLIC :: ma_set_first_core

  ! These are for Machine architecture internal use.
  !
  PUBLIC :: topology
  PUBLIC :: ma_hwloc, ma_libnuma
  PUBLIC :: ma_interface, ma_int_hwloc, ma_int_libnuma, ma_int_none

  ! For the machine representation
  TYPE(arch_topology)                            :: topology
#if defined (__HWLOC)
  TYPE(machine_output)                           :: console_out
#endif

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
    TYPE(arch_topology)                      :: dtopo
    INTEGER(KIND=C_INT)                      :: istat

     END FUNCTION ma_final_interface
  END INTERFACE

 INTERFACE
     SUBROUTINE ma_show_phys_topo(output_phys) BIND(C, name="hw_phys_pu_topology")
       USE ISO_C_BINDING
       USE machine_architecture_types,  ONLY: machine_output
    TYPE(machine_output)                     :: output_phys

     END SUBROUTINE ma_show_phys_topo
  END INTERFACE


  INTERFACE
     SUBROUTINE ma_show_branch(output_branch) BIND(C, name="hw_high_level_show")
       USE ISO_C_BINDING
       USE machine_architecture_types,  ONLY: machine_output
    TYPE(machine_output)                     :: output_branch

     END SUBROUTINE ma_show_branch
  END INTERFACE

  INTERFACE
     SUBROUTINE ma_show_topo(output_topo) BIND(C, name="hw_machine_topology")
       USE ISO_C_BINDING
       USE machine_architecture_types,  ONLY: machine_output
    TYPE(machine_output)                     :: output_topo

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
  FUNCTION ma_get_netDev(coreId, myRank) RESULT (netDev) BIND(C, name="hw_my_netCard")
    USE ISO_C_BINDING
    INTEGER(KIND=C_INT), VALUE               :: coreId, myRank
    INTEGER(KIND=C_INT)                      :: netDev

    END FUNCTION ma_get_netDev
  END INTERFACE

 INTERFACE
  FUNCTION ma_get_nnetDev(node) RESULT (nnetDev) BIND(C, name="hw_get_nnetcards")
    USE ISO_C_BINDING
    INTEGER(KIND=C_INT), VALUE               :: node
    INTEGER(KIND=C_INT)                      :: nnetDev

    END FUNCTION ma_get_nnetDev
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
    SUBROUTINE ma_set_thread_node(node) BIND(C, name="hw_set_thread_node")
       USE ISO_C_BINDING
    INTEGER(KIND=C_INT), VALUE               :: node

    END SUBROUTINE ma_set_thread_node
  END INTERFACE

  INTERFACE
    SUBROUTINE ma_set_thread_allnodes() BIND(C, name="hw_set_thread_allnodes")
       USE ISO_C_BINDING
    END SUBROUTINE ma_set_thread_allnodes
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
      SUBROUTINE ma_set_first_core(node,proc) BIND(C, name="hw_set_first_core_node")
       USE ISO_C_BINDING
    INTEGER(KIND=C_INT), VALUE               :: node, proc

     END SUBROUTINE ma_set_first_core
  END INTERFACE

  INTERFACE
      SUBROUTINE ma_set_core(cpu) BIND(C, name="hw_set_my_core")
       USE ISO_C_BINDING
    INTEGER(KIND=C_INT), VALUE               :: cpu

     END SUBROUTINE ma_set_core
  END INTERFACE

  INTERFACE
      FUNCTION ma_get_node_netDev() RESULT (node) BIND(C, name="hw_get_node_net")
       USE ISO_C_BINDING
    INTEGER(KIND=C_INT)                      :: node

     END FUNCTION ma_get_node_netDev
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
      FUNCTION ma_get_first_core() RESULT (core) BIND(C, name="hw_get_firstcore")
       USE ISO_C_BINDING
    INTEGER(KIND=C_INT)                      :: core

     END FUNCTION ma_get_first_core
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
!> \brief Empty routines needed to compile without HWLOC
!> \param core ...
!> \retval node ...
! *****************************************************************************
#if !defined (__HWLOC)
 FUNCTION ma_get_core_node(core) RESULT (node)
    INTEGER                                  :: core, node

    node = -1
    STOP "ma_get_core_node: HWLOC not compiled in"
 END FUNCTION ma_get_core_node

! *****************************************************************************
!> \brief ...
!> \param coreId ...
!> \param myRank ...
!> \retval netDev ...
! *****************************************************************************
 FUNCTION ma_get_netDev(coreId, myRank) RESULT (netDev)
    INTEGER                                  :: coreId, myRank, netDev

    netDev = -1
    STOP "ma_get_netDev: HWLOC not compiled in"
 END FUNCTION ma_get_netDev

! *****************************************************************************
!> \brief ...
!> \param node ...
!> \retval nnetDev ...
! *****************************************************************************
 FUNCTION ma_get_nnetDev(node) RESULT (nnetDev)
    INTEGER                                  :: node, nnetDev

    nnetDev = -1
    STOP "ma_get_nnetDev: HWLOC not compiled in"
 END FUNCTION ma_get_nnetDev

! *****************************************************************************
!> \brief ...
!> \retval node ...
! *****************************************************************************
 FUNCTION ma_get_node_netDev() RESULT (node)
    INTEGER                                  :: node

    node = -1
    STOP "ma_get_node_netDev: HWLOC not compiled in"
 END FUNCTION ma_get_node_netDev

! *****************************************************************************
!> \brief ...
!> \retval id ...
! *****************************************************************************
 FUNCTION ma_get_thread_id() RESULT (id)
    INTEGER                                  :: id

    id = -1
    STOP "ma_get_thread_id: HWLOC not compiled in"
 END FUNCTION ma_get_thread_id

! *****************************************************************************
!> \brief ...
!> \param mempol ...
!> \param node ...
!> \retval error ...
! *****************************************************************************
 FUNCTION ma_hw_set_mempol(mempol, node) RESULT (error)
    INTEGER                                  :: mempol, node, error

    error = -1
    STOP "ma_hw_set_mempol: HWLOC not compiled in"
 END FUNCTION ma_hw_set_mempol

! *****************************************************************************
!> \brief ...
!> \param mempol ...
!> \param node ...
! *****************************************************************************
 SUBROUTINE ma_hw_get_mempol(mempol, node)
    INTEGER                                  :: mempol, node

    STOP "ma_hw_get_mempol: HWLOC not compiled in"
 END SUBROUTINE ma_hw_get_mempol

! *****************************************************************************
!> \brief ...
!> \param cpu ...
! *****************************************************************************
 SUBROUTINE ma_set_core(cpu)
    INTEGER                                  :: cpu

    STOP "ma_set_core: HWLOC not compiled in"
 END SUBROUTINE ma_set_core

! *****************************************************************************
!> \brief ...
!> \param node ...
!> \param proc ...
! *****************************************************************************
 SUBROUTINE ma_set_first_core(node, proc)
    INTEGER                                  :: node, proc

    STOP "ma_set_first_core: HWLOC not compiled in"
 END SUBROUTINE ma_set_first_core

! *****************************************************************************
!> \brief ...
!> \param dist ...
!> \param core ...
! *****************************************************************************
 SUBROUTINE ma_set_proc_cores(dist, core)
    INTEGER                                  :: dist, core

    STOP "ma_set_proc_cores: HWLOC not compiled in"
 END SUBROUTINE ma_set_proc_cores

! *****************************************************************************
!> \brief ...
!> \param cpu ...
! *****************************************************************************
 SUBROUTINE ma_set_proc_core(cpu)
    INTEGER                                  :: cpu

    STOP "ma_set_proc_core: HWLOC not compiled in"
 END SUBROUTINE ma_set_proc_core

! *****************************************************************************
!> \brief ...
! *****************************************************************************
 SUBROUTINE ma_set_thread_allnodes()

    STOP "ma_set_thread_allnodes: HWLOC not compiled in"
 END SUBROUTINE ma_set_thread_allnodes

! *****************************************************************************
!> \brief ...
!> \param dist ...
!> \param core ...
! *****************************************************************************
 SUBROUTINE ma_set_thread_cores(dist, core)
    INTEGER                                  :: dist, core

    STOP "ma_set_thread_cores: HWLOC not compiled in"
 END SUBROUTINE ma_set_thread_cores

! *****************************************************************************
!> \brief ...
!> \param node ...
! *****************************************************************************
 SUBROUTINE ma_set_thread_node(node)
    INTEGER                                  :: node

    STOP "ma_set_thread_node: HWLOC not compiled in"
 END SUBROUTINE ma_set_thread_node

#endif


! *****************************************************************************
!> \brief Initialize the machine abstraction library
!> \param core ...
!> \param node ...
!> \retval istat ...
! *****************************************************************************
  FUNCTION ma_init_machine(core,node) RESULT(istat)
    INTEGER                                  :: core, node, istat

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
#if defined (__HWLOC) && !defined (__LIBNUMA)
   INTEGER  :: istat
   istat = ma_final_interface(topology)
   IF( istat .NE. 0) THEN
     CALL ma_error_stop(ma_error_finalize)
   ENDIF
#endif
  END SUBROUTINE ma_finalize_machine

! *****************************************************************************
!> \brief Count the number of machines used in an execution
!> \param group ...
!> \retval nmachines ...
! *****************************************************************************
  FUNCTION ma_get_nmachines (group) RESULT (nmachines)
    INTEGER                                  :: group
    INTEGER                                  :: nmachines
#if defined (__parallel)
    TYPE(ma_mp_type)                         :: ma_env

    CHARACTER(LEN=default_string_length)     :: host_name, string
    INTEGER                                  :: istat, jpe, nr_occu, nr_occu_e
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: all_host
    INTEGER                                  :: numnode, mynode

      CALL mp_environ(numnode, mynode, group)
      ma_env%numproc = numnode
      ma_env%mp_group = group
      CALL m_getpid(ma_env%myid)
      ma_env%myproc = mynode
      ma_env%all_proc = .TRUE.

      ALLOCATE (all_host(default_string_length,ma_env%numproc),STAT=istat)
      IF (istat /= 0) CALL ma_error_stop(ma_error_allocation)

      all_host(:,:) = 0
      CALL mp_proc_name(host_name)

      CALL string_to_ascii(host_name,all_host(:,ma_env%myproc+1))
      CALL mp_sum(all_host,ma_env%mp_group)

      nr_occu = 0
      nr_occu_e = 0
      DO jpe=1,ma_env%numproc
         CALL ascii_to_string(all_host(:,jpe),string)
         IF  (string .NE. host_name) THEN
            nr_occu = nr_occu +1
         ELSE
            nr_occu_e = nr_occu_e + 1
         ENDIF
      END DO
  
  nmachines = (nr_occu / nr_occu_e) + 1
#else
  nmachines = 0 
#endif
   
  END FUNCTION ma_get_nmachines

! *****************************************************************************
!> \brief Get the first core where a MPI process is running
!> \retval core ...
! *****************************************************************************
  FUNCTION ma_my_first_core() RESULT (core)
    INTEGER                                  :: core

    core = 0
#if defined (__HWLOC)
  core = ma_get_first_core()
#endif
  END FUNCTION ma_my_first_core


! *****************************************************************************
!> \brief Get the process ID
!> \retval pid ...
! *****************************************************************************
FUNCTION ma_get_id() RESULT (pid)
    INTEGER                                  :: pid

  CALL m_getpid(pid)
END FUNCTION ma_get_id

! *****************************************************************************
!> \brief Prints in text format the machine abstraction - only one branch
!> \param unit_num the output console 
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
!> \param unit_num the output console 
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
!> \brief ...
!> \retval core ...
! *****************************************************************************
  FUNCTION ma_get_mycore() RESULT (core)
    INTEGER                                  :: core

   core = -1
#if defined (__HWLOC)  && !defined (__LIBNUMA)
   core = ma_get_core()
#endif
#if defined (__LIBNUMA) && !defined (__HWLOC)
  core = ma_linux_get_core()
#endif
  END FUNCTION ma_get_mycore

! *****************************************************************************
!> \brief ...
!> \retval node ...
! *****************************************************************************
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

! *****************************************************************************
!> \brief ...
!> \param flag ...
!> \retval core ...
! *****************************************************************************
  FUNCTION ma_get_proc_core(flag) RESULT (core)
    INTEGER                                  :: flag, core

   core = -1
#if defined (__HWLOC)  && !defined (__LIBNUMA)
   core = ma_get_proccore(flag)
#endif
#if defined (__LIBNUMA) && !defined (__HWLOC)
   core = ma_linux_get_proccore()
#endif
  END FUNCTION ma_get_proc_core

! *****************************************************************************
!> \brief ...
!> \retval node ...
! *****************************************************************************
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

! *****************************************************************************
!> \brief ...
!> \param node ...
! *****************************************************************************
SUBROUTINE ma_set_proc_node(node)
    INTEGER                                  :: node
#if defined (__HWLOC) && !defined (__LIBNUMA)
   CALL ma_set_procnode(node)
#endif
#if defined (__LIBNUMA) && !defined (__HWLOC)
   CALL ma_linux_set_procnode(node)
#endif
  END SUBROUTINE ma_set_proc_node

! *****************************************************************************
!> \brief ...
!> \retval nnodes ...
! *****************************************************************************
  FUNCTION ma_get_nnodes() RESULT (nnodes)
    INTEGER                                  :: nnodes

     nnodes = topology%nnodes
  END FUNCTION ma_get_nnodes

! *****************************************************************************
!> \brief ...
!> \retval nsockets ...
! *****************************************************************************
  FUNCTION ma_get_nsockets() RESULT (nsockets)
    INTEGER                                  :: nsockets

     nsockets = topology%nsockets
  END FUNCTION ma_get_nsockets

! *****************************************************************************
!> \brief ...
!> \retval ncores ...
! *****************************************************************************
  FUNCTION ma_get_ncores() RESULT (ncores)
    INTEGER                                  :: ncores

     ncores = topology%ncores
  END FUNCTION ma_get_ncores

! *****************************************************************************
!> \brief ...
!> \retval ncaches ...
! *****************************************************************************
  FUNCTION ma_get_ncaches() RESULT (ncaches)
    INTEGER                                  :: ncaches

     ncaches = topology%ncaches
  END FUNCTION ma_get_ncaches

! *****************************************************************************
!> \brief ...
!> \retval nshared_caches ...
! *****************************************************************************
  FUNCTION ma_get_nshared_caches() RESULT (nshared_caches)
    INTEGER                                  :: nshared_caches

     nshared_caches = topology%nshared_caches
  END FUNCTION ma_get_nshared_caches

! *****************************************************************************
!> \brief ...
!> \retval nsiblings ...
! *****************************************************************************
  FUNCTION ma_get_nsiblings() RESULT (nsiblings)
    INTEGER                                  :: nsiblings

     nsiblings = topology%nsiblings
  END FUNCTION ma_get_nsiblings

! *****************************************************************************
!> \brief ...
!> \retval netcards ...
! *****************************************************************************
  FUNCTION ma_get_nnet_cards() RESULT (netcards)
    INTEGER                                  :: netcards

    netcards = topology%nnetcards
  END FUNCTION ma_get_nnet_cards

! *****************************************************************************
!> \brief Prints in text format the machine abstraction - the full machine node
!> \param unit_num the output console 
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
#if defined (__PW_CUDA) || defined (__DBCSR_ACC)
     WRITE(unit_num, "(T2,A,27X,I3)")      "MACHINE| Number of graphics processing units     ",&
                                topology%ngpus
#endif
     WRITE(unit_num, "(T2,A,27X,I3)")      "MACHINE| Number of cache levels                  ",&
                                topology%ncaches
     WRITE(unit_num, "(T2,A,27X,I3)")      "MACHINE| Number of shared caches                 ",&
                                topology%nshared_caches
     IF (topology%nnetcards .GT. 0) THEN
        WRITE(unit_num, "(T2,A,27X,I3)")   "MACHINE| Number of network cards                 ",&
                                topology%nnetcards
     ENDIF
     WRITE(unit_num, "()")
    ENDIF
  END SUBROUTINE ma_show_topology

END MODULE machine_architecture
