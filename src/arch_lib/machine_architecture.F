!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2011 Christiane Ribeiro and the CP2K developers group!
!-----------------------------------------------------------------------------!
!> \brief   Machine Hardware support for cp2k
!           Interface to the HWLOC library
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

#if defined (__HWLOC) || defined(__LIBNUMA)
  USE ISO_C_BINDING
#endif 

  USE machine_architecture_types
  USE machine_architecture_utils
  USE ma_config,                    ONLY: print_thread
  USE ma_dummy_c_bindings
  USE ma_message_passing,           ONLY: mp_allgather
  
  PRIVATE 

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'machine_architecture'

  PUBLIC :: ma_init, ma_finalize
  PUBLIC :: ma_error

  PUBLIC :: ma_get_nnodes, ma_get_nsockets, ma_get_ncores, ma_get_ncaches 
  PUBLIC :: ma_get_nshared_caches, ma_get_nsiblings
  
  PUBLIC :: ma_get_mynode, ma_get_mycore
 
  PUBLIC :: ma_sched_threads, ma_no_place, ma_linear_place

  PUBLIC :: ma_verify_place
  PUBLIC :: ma_running_on, ma_thread_running_on, ma_current_thread_run 

  PUBLIC :: ma_show_topology, ma_show_machine_branch
  PUBLIC :: ma_show_machine_full 

  PUBLIC :: ma_thread_run, ma_get_thread_run, ma_show_threads

#if defined (__HWLOC)
  PUBLIC :: ma_set_mempol, ma_get_mempol, ma_string_mempol
#endif

#if defined (__LIBNUMA)
  PUBLIC :: ma_linux_set_mempol, ma_linux_get_mempol
  PUBLIC :: ma_linux_get_core, ma_linux_set_core
  PUBLIC :: ma_linux_set_proccore, ma_linux_get_proccore
  PUBLIC :: ma_linux_get_my_id, ma_linux_nnodes
  PUBLIC :: ma_linux_ncores, ma_linux_get_mynode
#endif
   

  ! These are for Machine architecture internal use.
  !
  PUBLIC :: thread_mapping
  PUBLIC :: topology 
  PUBLIC :: console_output
  PUBLIC :: my_info
  PUBLIC :: ma_hwloc, ma_libnuma
  PUBLIC :: ma_interface, ma_int_hwloc, ma_int_libnuma, ma_int_none
  
  ! For the machine representation
  TYPE(arch_topology)                            :: topology 
  TYPE(machine_output)                           :: console_output

  TYPE(ma_process)                               :: my_info
  ! The global view of the thread mapping on the machine
  TYPE(thread_inf), DIMENSION(:,:), ALLOCATABLE  :: thread_mapping 

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


#if defined (__HWLOC) 

!
! interfaces to deal with HWLOC
!
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
     FUNCTION ma_get_node() RESULT (node) BIND(C, name="hw_my_node") 
       USE ISO_C_BINDING
    INTEGER(KIND=C_INT)                      :: node

     END FUNCTION ma_get_node
  END INTERFACE

  INTERFACE 
     FUNCTION ma_get_procnode() RESULT (node) BIND(C, name="hw_proc_node") 
       USE ISO_C_BINDING
    INTEGER(KIND=C_INT)                      :: node

     END FUNCTION ma_get_procnode
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
      SUBROUTINE ma_set_proccore(cpu) BIND(C, name="hw_set_proc_core") 
       USE ISO_C_BINDING
    INTEGER(KIND=C_INT), VALUE               :: cpu

     END SUBROUTINE ma_set_proccore
  END INTERFACE

  INTERFACE 
      FUNCTION ma_get_proccore() RESULT (core) BIND(C, name="hw_proc_core") 
       USE ISO_C_BINDING
    INTEGER(KIND=C_INT)                      :: core

     END FUNCTION ma_get_proccore
  END INTERFACE

  INTERFACE 
      FUNCTION ma_get_my_id() RESULT (id) BIND(C, name="hw_get_myid") 
       USE ISO_C_BINDING
    INTEGER(KIND=C_INT)                      :: id

     END FUNCTION ma_get_my_id
  END INTERFACE

  INTERFACE
     SUBROUTINE ma_set_mempol(mempol) BIND(C, name="hw_set_mempol")
       USE ISO_C_BINDING
    INTEGER(KIND=C_INT), VALUE               :: mempol

     END SUBROUTINE ma_set_mempol
  END INTERFACE 

  INTERFACE
      SUBROUTINE ma_get_mempol(node, mempol) BIND(C, name="hw_get_mempol")
       USE ISO_C_BINDING
    INTEGER(KIND=C_INT)                      :: node, mempol

     END SUBROUTINE ma_get_mempol
  END INTERFACE

#endif

#if defined (__LIBNUMA)

  INTERFACE
      SUBROUTINE ma_linux_get_mempol(node, mempol) BIND(C, name="linux_get_mempol")
       USE ISO_C_BINDING
    INTEGER(KIND=C_INT)                      :: node, mempol

     END SUBROUTINE ma_linux_get_mempol
  END INTERFACE

  INTERFACE
     SUBROUTINE ma_linux_set_mempol(mempol) BIND(C, name="linux_set_mempol")
       USE ISO_C_BINDING
    INTEGER(KIND=C_INT), VALUE               :: mempol

     END SUBROUTINE ma_linux_set_mempol
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
      FUNCTION ma_linux_get_my_id() RESULT (id) BIND(C, name="linux_get_myid")
       USE ISO_C_BINDING
    INTEGER(KIND=C_INT)                      :: id

     END FUNCTION ma_linux_get_my_id
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

#endif

CONTAINS

! *****************************************************************************
!> \brief To convert a memory policy descriptor to a string
!> \var mempol          Memory policy descriptor
!> \var mempolicy       String to hold the memory policy name
! *****************************************************************************
SUBROUTINE ma_string_mempol(mempol, mempolicy)
    INTEGER, INTENT(IN)                      :: mempol
    CHARACTER(LEN=10)                        :: mempolicy

   SELECT CASE (mempol)
   CASE (os)
     mempolicy = 'DEFAULT'
   CASE (local)
     mempolicy = 'LOCAL'
   CASE (interleave)
     mempolicy = 'INTERLEAVE'
   CASE DEFAULT 
     mempolicy = 'NONE'
   END SELECT
END SUBROUTINE ma_string_mempol

! *****************************************************************************
!> \brief To show the library erros
!> \var message         The error message  
!> \var unit_num        The output console
! *****************************************************************************
  SUBROUTINE ma_error (message, unit_num)
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: message
    INTEGER                                  :: unit_num

    IF (PRESENT(message)) THEN
       WRITE(unit_num,"(T2,A)")        "MACHINE| "//TRIM(message)
       WRITE(unit_num,"(T2,A)")        "MACHINE| Abnormal program termination, stopped due to HWLOC"
    ENDIF  
    STOP 1
  END SUBROUTINE ma_error

! *****************************************************************************
!> \brief Initialize the machine abstraction library
!> \var ma_env            The parallel environment
!> \var num_threads       The number of threads of each process
! *****************************************************************************
  INTEGER FUNCTION ma_init(ma_env, num_threads) RESULT(istat)
   INTEGER                          :: num_threads
   TYPE(ma_mp_type)                 :: ma_env

   topology%nnodes = -1
   istat = -1
#if defined (__HWLOC) || defined (__LIBNUMA)
   ALLOCATE(thread_mapping(ma_env%numproc,num_threads))
   CALL ma_create_process(my_info,ma_env,num_threads)
#endif 
#if defined (__HWLOC) && !defined (__LIBNUMA)
   istat = ma_init_interface(topology)
#endif 
#if defined (__LIBNUMA) && !defined (__HWLOC)
   topology%nnodes = ma_linux_nnodes()
   topology%ncores = ma_linux_ncores() 
   topology%npus = ma_linux_ncores()   
   topology%nsockets = 0
   topology%ncaches = 0
   topology%nshared_caches = 0
   topology%nsiblings = 0
   topology%nmemcontroller = 0
#endif 
  END FUNCTION ma_init 

! *****************************************************************************
!> \brief Finalize the machine abstraction library
! *****************************************************************************
   SUBROUTINE ma_finalize()
   INTEGER                :: error
#if defined (__HWLOC) && !defined (__LIBNUMA)
   error = ma_final_interface(topology)
#endif 
#if defined (__LIBNUMA) || defined (__HWLOC)
   DEALLOCATE(thread_mapping)
#endif
  END SUBROUTINE ma_finalize

! *****************************************************************************
!> \brief Prints in text format the machine abstraction - only one branch
!> \var unit_num the output console
! *****************************************************************************
  SUBROUTINE ma_show_machine_branch(unit_num)
  INTEGER :: unit_num
#if defined (__HWLOC)
  CALL ma_show_branch(console_out)
  WRITE(unit_num,*)console_out%console_output(1:console_out%len)
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
!> \brief Get the threads mapping on the machine cores
!> \var id thread OMP id
! *****************************************************************************
  SUBROUTINE ma_no_place(id)
  INTEGER :: id,core
#if defined (__HWLOC)  && !defined (__LIBNUMA)
  my_info%threads_info(id+1)%core = ma_get_mycore()  
  my_info%threads_info(id+1)%node = ma_get_mynode() 
  my_info%threads_info(id+1)%id_omp = id
  my_info%threads_info(id+1)%id_real = ma_get_my_id()  
#endif
#if defined (__LIBNUMA) && !defined (__HWLOC)
  core = ma_linux_get_core()
  my_info%threads_info(id+1)%core = core
  my_info%threads_info(id+1)%node = ma_linux_get_mynode()
  my_info%threads_info(id+1)%id_omp = id
  my_info%threads_info(id+1)%id_real = ma_linux_get_my_id()
#endif
  END SUBROUTINE ma_no_place

! *****************************************************************************
!> \brief Set the threads mapping on the machine cores - round-robin
!> \var id thread OMP id
! *****************************************************************************
  SUBROUTINE ma_linear_place(id)
  INTEGER :: id, core
#if defined (__HWLOC) && !defined (__LIBNUMA)
  core =  MOD(my_info%mp_info%myproc*my_info%nr_threads+id,ma_get_ncores())
  CALL ma_set_core(core)
  my_info%threads_info(id+1)%core = core  
  my_info%threads_info(id+1)%node = ma_get_mynode() 
  my_info%threads_info(id+1)%id_omp = id
  my_info%threads_info(id+1)%id_real = ma_get_my_id()  
#endif
#if defined (__LIBNUMA) && !defined (__HWLOC)
  core =  MOD(my_info%mp_info%myproc*my_info%nr_threads+id,ma_linux_ncores())
  CALL ma_linux_set_core(core)
  my_info%threads_info(id+1)%core = core
  my_info%threads_info(id+1)%node = ma_linux_get_mynode()
  my_info%threads_info(id+1)%id_omp = id
  my_info%threads_info(id+1)%id_real = ma_linux_get_my_id()
#endif
  END SUBROUTINE ma_linear_place

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

  INTEGER FUNCTION ma_get_proc_core() RESULT (core)
   core = -1  
#if defined (__HWLOC)  && !defined (__LIBNUMA)
   core = ma_get_proccore()
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
     WRITE(unit_num, "(T2,A,27X,I3)")      "MACHINE| Number of caches levels                 ",&
                                topology%ncaches
     WRITE(unit_num, "(T2,A,27X,I3)")      "MACHINE| Number of shared caches                 ",&
                                topology%nshared_caches
     WRITE(unit_num, "()")     
    ENDIF
  END SUBROUTINE ma_show_topology

! *****************************************************************************
!> \brief Get the core and NUMA node of the process
!> \var my_core the core where the process is running
!> \var my_node the NUMA node where the process is running
! *****************************************************************************
  SUBROUTINE ma_running_on (my_core, my_node) 
    INTEGER                                  :: my_core, my_node

     my_core = ma_get_proc_core()
     my_node = ma_get_proc_node()     
     my_info%core = my_core
     my_info%node = my_node
  END SUBROUTINE ma_running_on

! *****************************************************************************
!> \brief Update for all processes where the threads are running, also print
!> \var unit_num the output console
! *****************************************************************************
  SUBROUTINE ma_thread_running_on (unit_num) 
    INTEGER                                  :: unit_num

    INTEGER                                  :: id, p
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: thread_send
    INTEGER, ALLOCATABLE, DIMENSION(:, :, :) :: thread_recv

  IF (has_MPI) THEN

    ALLOCATE(thread_send(4,my_info%nr_threads))
    ALLOCATE(thread_recv(4,my_info%nr_threads,my_info%mp_info%numproc))

    thread_send(:,:) = 0
    thread_recv(:,:,:) = 0

    CALL ma_pack_threads(my_info,thread_send)
    CALL mp_allgather(thread_send, thread_recv, my_info%mp_info%mp_group)
    CALL ma_unpack_threads(my_info,thread_recv,thread_mapping)
  ENDIF
   
   IF(unit_num .GT. 0) THEN
   IF (my_info%nr_threads .GT. 1) THEN 
     WRITE(unit_num,"(T2,A)") "OMP | Thread placement"
     DO p=1, my_info%mp_info%numproc 
      WRITE(unit_num, "(T2,A,I4)") "Processs ", p-1
      DO id=1, my_info%nr_threads 
       IF (thread_mapping(p,id)%node .EQ. -1) THEN      
         WRITE(unit_num, "(T2,A,I4,2X,A,I3)")   "OMP | Thread ", thread_mapping(p,id)%id_omp,&
         " running on core ", thread_mapping(p,id)%core
       ELSE
        WRITE(unit_num, "(T2,A,I4,2X,A,I4,2X,A,I4)")   "OMP | Thread ", thread_mapping(p,id)%id_omp,&
        " running on NUMA node ", thread_mapping(p,id)%node, " core ", thread_mapping(p,id)%core
       ENDIF
      END DO
     END DO
     WRITE(unit_num,'()')  
   ENDIF
   ENDIF
  END SUBROUTINE ma_thread_running_on

! *****************************************************************************
!> \brief Update for all processes where the threads are running, also print
!>        To be called inside the cp2k libraries
! *****************************************************************************
  SUBROUTINE ma_thread_run ()
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: thread_send
    INTEGER, ALLOCATABLE, DIMENSION(:, :, :) :: thread_recv

  IF (has_MPI) THEN

    ALLOCATE(thread_send(4,my_info%nr_threads))
    ALLOCATE(thread_recv(4,my_info%nr_threads,my_info%mp_info%numproc))

    thread_send(:,:) = 0
    thread_recv(:,:,:) = 0

    CALL ma_pack_threads(my_info,thread_send)
    CALL mp_allgather(thread_send, thread_recv, my_info%mp_info%mp_group)
    CALL ma_unpack_threads(my_info,thread_recv,thread_mapping)
  ENDIF

  END SUBROUTINE ma_thread_run

! *****************************************************************************
!> \brief Shows threads of each MPI process
! *****************************************************************************
  SUBROUTINE ma_show_threads (unit_num)
    INTEGER                                  :: unit_num

    INTEGER                                  :: id, p

  IF ( print_thread ) THEN
   IF(unit_num .GT. 0 .AND. my_info%mp_info%myproc .EQ. 0) THEN
   IF (my_info%nr_threads .GT. 1) THEN
     WRITE(unit_num,"()")
     WRITE(unit_num,"(T2,A)") "OMP | Thread placement"
     DO p=1, my_info%mp_info%numproc
      WRITE(unit_num, "(T2,A,I4)") "Processs ", p-1
      DO id=1, my_info%nr_threads
       IF (thread_mapping(p,id)%node .EQ. -1) THEN
         WRITE(unit_num, "(T2,A,I4,2X,A,I3)")   "OMP | Thread ", thread_mapping(p,id)%id_omp,&
         " running on core ", thread_mapping(p,id)%core
       ELSE
        WRITE(unit_num, "(T2,A,I4,2X,A,I4,2X,A,I4)")   "OMP | Thread ", thread_mapping(p,id)%id_omp,&
        " running on NUMA node ", thread_mapping(p,id)%node, " core ", thread_mapping(p,id)%core
       ENDIF
      END DO
     END DO
     WRITE(unit_num,'()')
   ENDIF
   ENDIF
  ENDIF
  END SUBROUTINE ma_show_threads

! *****************************************************************************
!> \brief Print the current core and NUMA node where a thread is running
!> \var unit_num the output console
!> \var id Thread OMP id
! *****************************************************************************
  SUBROUTINE ma_current_thread_run(id, unit_num)
    INTEGER                                  :: id, unit_num

    INTEGER                                  :: core, node

     core = ma_get_mycore()
     node = ma_get_mynode()     

   IF(unit_num .GT. 0) THEN
     IF (node .EQ. -1) THEN      
         WRITE(unit_num, "(T2,A,I4,A,I4,2X,A,I3)") "OMP | Thread ", id,&
         " from process ", my_info%mp_info%myproc, " running on core ", core
     ELSE
        WRITE(unit_num, "(T2,A,I4,A,I4,2X,A,I4,2X,A,I4)")   "OMP | Thread ", id,&
        " from process ", my_info%mp_info%myproc, " running on NUMA node ", node, " core ", core
     ENDIF
   ENDIF
  END SUBROUTINE ma_current_thread_run

! *****************************************************************************
! Get the core for a thread
! *****************************************************************************
  SUBROUTINE ma_get_thread_run(id)
    INTEGER  :: id
#if defined (__HWLOC) || defined (__LIBNUMA)
  my_info%threads_info(id+1)%core = ma_get_mycore()
  my_info%threads_info(id+1)%node = ma_get_mynode()
  my_info%threads_info(id+1)%id_omp = id
#endif
#if defined (__HWLOC) && !defined (__LIBNUMA)
  my_info%threads_info(id+1)%id_real = ma_get_my_id()
#endif
#if !defined (__HWLOC) && defined (__LIBNUMA)
  my_info%threads_info(id+1)%id_real = ma_linux_get_my_id()
#endif
  END SUBROUTINE ma_get_thread_run


! *****************************************************************************
! Get/set the core for a thread
! *****************************************************************************
  SUBROUTINE ma_sched_threads(id, policy)
    INTEGER                                  :: id, policy

     SELECT CASE (policy) 
     CASE (nosched)
       CALL ma_no_place(id)
     CASE (linear) 
       CALL ma_linear_place(id)
     END SELECT
  END SUBROUTINE ma_sched_threads

! *****************************************************************************
!> \brief Verify the last and current core for a thread
!> \var unit_num the output console
!> \var id Thread OMP id
! *****************************************************************************
  SUBROUTINE ma_verify_place(id, unit_num)
    INTEGER                                  :: id, unit_num

    INTEGER                                  :: core, node, pe

     core = ma_get_mycore()
     node = ma_get_mynode()     
     pe = my_info%mp_info%myproc

     IF (thread_mapping(pe+1,id+1)%core .NE. core) THEN
      WRITE(unit_num,'(T2,A)') "SCHED | Thread placement changed"   
      IF (node .EQ. -1) THEN      
         WRITE(unit_num, "(T2,A,I4,2X,A,I3)")   "OMP | Thread ", thread_mapping(pe+1,id+1)%id_omp,&
         " WAS running on core ", thread_mapping(pe+1,id+1)%core
         WRITE(unit_num, "(T2,A,I4,2X,A,I3)")   "OMP | Thread ", id,&
         " running on core ", core        
      ELSE
        WRITE(unit_num, "(T2,A,I4,2X,A,I4,2X,A,I4)")   "OMP | Thread ", &
        thread_mapping(pe+1,id+1)%id_omp," WAS running on NUMA node ",&
        thread_mapping(pe+1,id+1)%node, " core ",thread_mapping(pe+1,id+1)%core
        WRITE(unit_num, "(T2,A,I4,2X,A,I4,2X,A,I4)")   "OMP | Thread ", id,&
        " running on NUMA node ", node, " core ", core
      ENDIF
    ENDIF
  END SUBROUTINE ma_verify_place

END MODULE machine_architecture 
