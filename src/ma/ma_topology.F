!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulationsi        !
!   Copyright (C) 2000 - 2014  CP2K developers group group!
!-----------------------------------------------------------------------------!
!> \brief   network topology support for cp2k
!> \author  Christiane Pousa Ribeiro
!> \date    2012-05-02
!> \version 1.0
!>
!> <b>Modification history:</b>
!> - Created 2012-05-02


MODULE ma_topology

  USE ISO_C_BINDING
  USE ma_errors,                       ONLY: ma_error_allocation,&
                                             ma_error_database,&
                                             ma_error_option,&
                                             ma_error_remove,&
                                             ma_error_stop
  USE ma_kinds,                        ONLY: default_string_length,&
                                             dp
  USE machine,                         ONLY: m_getpid,&
                                             m_walltime
  USE machine_architecture_types,      ONLY: cannon_graph,&
                                             complete_graph,&
                                             net_topo
  USE mathlib,                         ONLY: gcd
  USE message_passing,                 ONLY: mp_environ,&
                                             mp_isendrecv,&
                                             mp_proc_name,&
                                             mp_sum,&
                                             mp_waitall
  USE string_utilities,                ONLY: integer_to_string

IMPLICIT NONE

 PRIVATE
 
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ma_topology'

  !Maximum number of nodes in each direction
  PUBLIC :: ma_net_topology, ma_allocated_topology
  PUBLIC :: ma_compute_contention
  PUBLIC :: ma_comm_graph
  PUBLIC :: ma_2dgrid_dimensions
  PUBLIC :: ma_get_job_allocation
  PUBLIC :: ma_destroy_topology
  
 
! These are for Machine architecture internal use.
!
  PUBLIC :: net_topology
  PUBLIC :: allocated_topology

  TYPE(net_topo),DIMENSION(:),POINTER      :: net_topology
#if defined __parallel
#if defined __GEMINI || __SEASTAR
  INTEGER                                  :: max_x
  INTEGER                                  :: max_y
  INTEGER                                  :: max_z
#endif 
#endif

  INTEGER, DIMENSION(:), POINTER           :: allocated_topology
  INTEGER, DIMENSION(:,:), POINTER         :: comm_graph
  INTEGER                                  :: max_nodes

#if defined __NET
  REAL, DIMENSION(:,:), POINTER         :: dist_nodes
#endif
!
! Interfaces to deal with CRAY Interfaces
!
#if defined __parallel
#if defined __GEMINI || __SEASTAR

  INTERFACE
     FUNCTION ma_extract_topology() RESULT (istat) BIND(C, name="extract_topology")
       USE ISO_C_BINDING
    INTEGER(KIND=C_INT)                      :: istat

     END FUNCTION ma_extract_topology
  END INTERFACE

  INTERFACE
     FUNCTION ma_remove_topology() RESULT (istat) BIND(C, name="remove_topology")
       USE ISO_C_BINDING
    INTEGER(KIND=C_INT)                      :: istat

     END FUNCTION ma_remove_topology
  END INTERFACE 

  INTERFACE 
   FUNCTION ma_get_job_alloc_strategy() RESULT(istat) BIND(C,name="get_job_alloc_strategy")
       USE ISO_C_BINDING
    INTEGER(KIND=C_INT)                      :: istat

   END FUNCTION ma_get_job_alloc_strategy 
  END INTERFACE
 
#endif
#endif

CONTAINS

! *****************************************************************************
!> \brief Get the total number of nodes of the machine
!> \retval numnodes ...
! *****************************************************************************
 FUNCTION ma_num_nodes() RESULT (numnodes)
    INTEGER                                  :: numnodes

  numnodes = max_nodes

 END FUNCTION ma_num_nodes

! *****************************************************************************
!> \brief Get the job allocation strategy used on Cray machines
!> \retval strategy ...
! *****************************************************************************
 FUNCTION ma_get_job_allocation() RESULT (strategy)
    INTEGER                                  :: strategy

   strategy = 0

#if defined __parallel
#if defined __GEMINI || __SEASTAR
  strategy = ma_get_job_alloc_strategy()
#else
  strategy = -1
#endif
#endif
  END FUNCTION ma_get_job_allocation

! *****************************************************************************
!> \brief Initialize the network topology and send to all MPI processes.
!> The topology is structured differently for each type of system.
!> For instance, on Cray machines with torus it keeps the coordinates. On 
!> machines with no API support for topology information, we use definition
!> of factor (remote/local). 
!> \param mp_group ...
!> \note     2.2012 created [ Christiane Pousa ]
! *****************************************************************************
 SUBROUTINE ma_net_topology (mp_group)
    INTEGER                                :: mp_group
    
  CHARACTER(LEN=*), PARAMETER :: routineN = 'ma_net_topology', &
      routineP = moduleN//':'//routineN

#if defined(__parallel)
  INTEGER                                  :: istat
  INTEGER                                  :: source, numtask
#if defined(__SEASTAR)
  CHARACTER(LEN=default_string_length)     :: name_file, myid_str
  INTEGER                                  :: endoffile, fnid, fx, fy, fz, &
                                              myid, nfields, rst_unit
#elif defined(__GEMINI)
  CHARACTER(LEN=default_string_length)     :: name_file, myid_str
  INTEGER                                  :: endoffile, fnid, fx, fy, fz, &
                                              myid, nfields, rst_unit
  INTEGER, DIMENSION(:), POINTER           :: nid,x,y,z
#elif defined(__NET)
  INTEGER                                  :: i, j
  REAL(KIND=dp)                            :: startTime, endTime
  INTEGER, DIMENSION(2)                    :: req
  REAL, DIMENSION(:), POINTER              :: msgin, msgout
#else
  CHARACTER(LEN=default_string_length)     :: host_name
#endif
#endif

!if it is the seastar connection
#if defined __SEASTAR && __parallel
  rst_unit = -1
  nfields = 0 
 
  CALL mp_environ(numtask,source,mp_group)
  IF (source == 0) THEN
    istat =  ma_extract_topology()
    IF (istat /= 0) CALL ma_error_stop(ma_error_database)
    rst_unit = 121314  
    CALL m_getpid(myid)
    CALL integer_to_string(myid,myid_str)
    name_file = "topo_"//TRIM(myid_str)//"/topology"
    OPEN(rst_unit,FILE=name_file,ACTION="READ",STATUS="OLD",ACCESS="SEQUENTIAL")
    READ(rst_unit,*) nfields
    max_nodes = nfields
  ENDIF

  CALL mp_sum(nfields,source,mp_group)
  max_nodes = nfields
 
  ALLOCATE(net_topology(nfields),STAT=istat)
  IF (istat /= 0) THEN
        CALL ma_error_stop(ma_error_allocation)
  END IF
 
  net_topology(:)%nid = 0
  net_topology(:)%coord_x = 0
  net_topology(:)%coord_y = 0
  net_topology(:)%coord_z = 0
  max_x  = 0
  max_y  = 0
  max_z  = 0
  i = 1

  IF (source == 0 ) THEN
   DO 
     READ(rst_unit,*, IOSTAT=endoffile) fnid,fx,fy,fz
     IF (endoffile .NE. 0) EXIT
     net_topology(i)%nid = fnid
     net_topology(i)%coord_x = fx
     net_topology(i)%coord_y = fy
     net_topology(i)%coord_z = fz
     i = i + 1   
   ENDDO
   CLOSE(rst_unit)
   istat =  ma_remove_topology()
   IF (istat /= 0) CALL ma_error_stop(ma_error_remove)
 ENDIF

  CALL mp_sum(net_topology(:)%nid,mp_group)
  CALL mp_sum(net_topology(:)%coord_x,mp_group)
  CALL mp_sum(net_topology(:)%coord_y,mp_group)
  CALL mp_sum(net_topology(:)%coord_z,mp_group)

  max_x = MAXVAL(net_topology(:)%coord_x)+1
  max_y = MAXVAL(net_topology(:)%coord_y)+1
  max_z = MAXVAL(net_topology(:)%coord_z)+1

!if it is Gemini, then we have to take into account the two nodes connect to one 
!gemini connection 
#elif __GEMINI  && __parallel
  rst_unit = -1
  nfields = 0

  CALL mp_environ(numtask,source,mp_group)

#ifdef __DEBUG_HWTOPO
   IF(source==0) WRITE(*,*) 'Extracting CRAY machine topology - accessing xtdb2proc'   
#endif 

  IF (source==0) THEN
    CALL m_getpid(myid)
    CALL integer_to_string(myid,myid_str)
    name_file = "topo_"//TRIM(myid_str)//"/topology"
    istat =  ma_extract_topology()
    IF (istat /= 0) CALL ma_error_stop(ma_error_database)
    rst_unit = 121314
    OPEN(rst_unit,FILE=name_file,ACTION="READ",STATUS="OLD",ACCESS="SEQUENTIAL")
    READ(rst_unit,*) nfields
  ENDIF

  CALL mp_sum(nfields,mp_group)
  nfields = nfields + 1
  ALLOCATE(net_topology(nfields),STAT=istat)
  IF (istat /= 0) THEN
        CALL ma_error_stop(ma_error_allocation)
  END IF

  IF (source == 0) THEN
   ALLOCATE(nid(nfields),x(nfields),y(nfields),z(nfields),STAT=istat)
    IF (istat /= 0) THEN
         CALL ma_error_stop(ma_error_allocation)
    END IF

    nid(:) = 0
    x(:) = 0
    y(:) = 0
    z(:) = 0 
    endoffile = 0
    
    DO 
     READ(rst_unit,*,IOSTAT=endoffile) fnid,fx,fy,fz
     IF (endoffile .NE. 0) EXIT
     nid(fnid+1) = fnid
     x(fnid+1) = fx
     z(fnid+1) = fz
     y(fnid+1) = 2*fy+MOD(fnid+1,2) 
    ENDDO
    CLOSE(rst_unit)

    istat =  ma_remove_topology()
    IF (istat /= 0) CALL ma_error_stop(ma_error_remove)

    net_topology(:)%nid = nid(:)
    net_topology(:)%coord_x = x(:)
    net_topology(:)%coord_y = y(:)
    net_topology(:)%coord_z = z(:)   

    IF(ASSOCIATED(nid))DEALLOCATE(nid,x,y,z,STAT=istat)
    IF (istat /= 0) THEN
        CALL ma_error_stop(ma_error_allocation)
    END IF
  ENDIF

   max_x  = 0
   max_y  = 0
   max_z  = 0

   IF ( source .NE. 0) THEN
     net_topology(:)%nid     = 0
     net_topology(:)%coord_x = 0
     net_topology(:)%coord_y = 0
     net_topology(:)%coord_z = 0
   ENDIF

   CALL mp_sum(net_topology(:)%nid,mp_group)
   CALL mp_sum(net_topology(:)%coord_x,mp_group)
   CALL mp_sum(net_topology(:)%coord_y,mp_group)
   CALL mp_sum(net_topology(:)%coord_z,mp_group)

   max_x = MAXVAL(net_topology(:)%coord_x)+1
   max_y = MAXVAL(net_topology(:)%coord_y)+1
   max_z = MAXVAL(net_topology(:)%coord_z)+1

#elif defined __NET && __parallel

 CALL mp_environ(numtask,source,mp_group)
 ALLOCATE (dist_nodes(numtask,numtask),STAT=istat)
 IF (istat /= 0) THEN
        CALL ma_error_stop(ma_error_allocation)
 END IF
 dist_nodes(:,:) = 0.0
 
 ALLOCATE (msgin(1024*1024*10),STAT=istat)
 IF (istat /= 0) THEN
        CALL ma_error_stop(ma_error_allocation)
 END IF
 msgin = 0.0

 ALLOCATE (msgout(1024*1024*10),STAT=istat)
 IF (istat /= 0) THEN
        CALL ma_error_stop(ma_error_allocation)
 END IF
 msgout = 0.0 
 
! here we collect the timings to compute the distance between
 ! the nodes
 DO j=0, numtask-1
  startTime = m_walltime ( )
  CALL mp_isendrecv(msgin,source,msgout,j,mp_group,req(1),req(2))
  CALL mp_waitall(req)
  endTime = m_walltime ( ) - startTime
  dist_nodes(source+1,j+1) = endTime
  dist_nodes(j+1,source+1) = endTime
 ENDDO   
 
 !And now update all processes with the comm costs 
 CALL mp_sum(dist_nodes,mp_group) 
 
 IF(ASSOCIATED(msgin))DEALLOCATE (msgin,STAT=istat)
 IF (istat /= 0) THEN
        CALL ma_error_stop(ma_error_allocation)
 END IF

 IF(ASSOCIATED(msgout))DEALLOCATE (msgout,STAT=istat)
 IF (istat /= 0) THEN
        CALL ma_error_stop(ma_error_allocation)
 END IF

#elif defined __parallel 

  CALL mp_environ(numtask,source,mp_group)

  ALLOCATE(net_topology(numtask),STAT=istat)
  IF (istat /= 0) THEN
        CALL ma_error_stop(ma_error_allocation)
  END IF
  
  net_topology(:)%coord_x = 0
  net_topology(:)%coord_y = 0
  net_topology(:)%coord_z = 0
 
 CALL mp_proc_name(host_name)
 READ(host_name(4:LEN_TRIM(host_name)),*) net_topology(source)

#endif

 END SUBROUTINE ma_net_topology

! *****************************************************************************
!> \brief Initialize the allocated topology and update to all MPI processes
!> \param mp_group the communicator 
! *****************************************************************************
 SUBROUTINE ma_allocated_topology (mp_group)
    INTEGER                                  :: mp_group

    INTEGER                                  :: istat, numtask, &
                                                source
#if defined __GEMINI || __SEASTAR
    CHARACTER(LEN=default_string_length)     :: host_name
#elif defined __NET
    CHARACTER(LEN=default_string_length)     :: host_name
    INTEGER                                  :: host_code, i
#endif

  CALL mp_environ(numtask,source,mp_group) 
  ALLOCATE (allocated_topology(numtask),STAT=istat)
  IF (istat /= 0) THEN
        CALL ma_error_stop(ma_error_allocation)
  END IF
  
  allocated_topology(:) = 0
  
#if defined __GEMINI || __SEASTAR

#ifdef __DEBUG_HWTOPO
   IF(source==0) WRITE(*,*) 'Obtaining allocated nodes on the machine'  
#endif

 CALL mp_proc_name(host_name)
 READ(host_name(4:LEN_TRIM(host_name)),*) allocated_topology(source+1) 
 
 CALL mp_sum(allocated_topology,mp_group) 

#ifdef __DEBUG_HWTOPO
   IF(source==0) WRITE(*,*) 'Allocated nodes on the machine: ', allocated_topology              
#endif

#elif defined __NET

 CALL mp_proc_name(host_name)
 host_code = 0
 DO i=1,LEN(host_name)
  host_code = host_code + ICHAR(host_name(i:i))
 ENDDO 

 allocated_topology(source+1) = host_code
 CALL mp_sum(allocated_topology,mp_group) 

#else
  allocated_topology(:) = 0
#endif

 END SUBROUTINE ma_allocated_topology

! *****************************************************************************
!> \brief Get the allocated topology
!> \param alloc_topo  the allocated topology 
! *****************************************************************************
 SUBROUTINE ma_get_allocated_topology (alloc_topo)
    INTEGER, DIMENSION(:), POINTER           :: alloc_topo

  alloc_topo(:) = allocated_topology(:)  
 
 END SUBROUTINE ma_get_allocated_topology

! *****************************************************************************
!> \brief Set the allocated topology
!> \param alloc_topo  the allocated topology 
! *****************************************************************************
 SUBROUTINE ma_set_allocated_topology (alloc_topo)
    INTEGER, DIMENSION(:), POINTER           :: alloc_topo

  allocated_topology(:) = alloc_topo(:)  

 END SUBROUTINE ma_set_allocated_topology

! *****************************************************************************
!> \brief Deallocate the topologies structures 
! *****************************************************************************
 SUBROUTINE ma_destroy_topology ()
    INTEGER                                  :: istat

#if defined __NET
  IF (ASSOCIATED(dist_nodes)) DEALLOCATE(dist_nodes,STAT=istat)
  IF (istat /= 0) THEN
        CALL ma_error_stop(ma_error_allocation)
  END IF
#endif

  IF(ASSOCIATED(net_topology))DEALLOCATE(net_topology,STAT=istat)
  IF (istat /= 0) THEN
        CALL ma_error_stop(ma_error_allocation)
  END IF

  IF(ASSOCIATED(comm_graph))DEALLOCATE(comm_graph,STAT=istat)
  IF (istat /= 0) THEN
        CALL ma_error_stop(ma_error_allocation)
  END IF

  IF(ASSOCIATED(allocated_topology))DEALLOCATE(allocated_topology,STAT=istat)
  IF (istat /= 0) THEN
        CALL ma_error_stop(ma_error_allocation)
  END IF


 END SUBROUTINE ma_destroy_topology

! *****************************************************************************
!> \brief Get the maximum number of nodes for a 3D Torus in X direction
!> \retval maxx ...
! *****************************************************************************
 FUNCTION ma_nodes_x() RESULT (maxx)
  INTEGER               ::         maxx

#if defined __GEMINI || __SEASTAR
  maxx = max_x 
#else
  maxx = 0
#endif

 END FUNCTION ma_nodes_x 

! *****************************************************************************
!> \brief Get the maximum number of nodes for a 3D Torus in Y direction
!> \retval maxy ...
! *****************************************************************************
 FUNCTION ma_nodes_y() RESULT (maxy)
  INTEGER  ::         maxy

#if defined __GEMINI || __SEASTAR
  maxy = max_y
#else
  maxy = 0
#endif

 END FUNCTION ma_nodes_y

! *****************************************************************************
!> \brief Get the maximum number of nodes for a 3D Torus in Z direction
!> \retval maxz ...
! *****************************************************************************
 FUNCTION ma_nodes_z() RESULT (maxz)
  INTEGER  ::         maxz

#if defined __GEMINI || __SEASTAR
  maxz = max_z
#else
  maxz = 0
#endif

 END FUNCTION ma_nodes_z

! *****************************************************************************
!> \brief Get the host coordinates for a MPI process
!> \param nid ...
!> \param x ...
!> \param y ...
!> \param z ...
! *****************************************************************************
 SUBROUTINE ma_get_host_coordinates (nid, x,y,z)
  INTEGER, INTENT(IN)               :: nid
  INTEGER, INTENT(OUT)              :: x,y,z

#if defined __GEMINI || __SEASTAR
   x = net_topology(nid+1)%coord_x
   y = net_topology(nid+1)%coord_y
   z = net_topology(nid+1)%coord_z
#else
  x = 0 
  y = 0 
  z = 0
#endif

 END SUBROUTINE ma_get_host_coordinates

! *****************************************************************************
!> \brief Get the number of hops between two host
!> \param src   the id of the source node 
!> \param dest  the id of the dest node 
!> \retval nhops ...
! *****************************************************************************
  FUNCTION ma_hops (src,dest) RESULT (nhops)
   INTEGER, INTENT(IN)          :: src, dest

   INTEGER                      :: nhops

#if defined __GEMINI || __SEASTAR  
   INTEGER                      :: src_x, src_y, src_z,&
                                   dest_x, dest_y, dest_z,&
                                   x, y, z
   REAL                         :: tmp

   CALL ma_get_host_coordinates(src,src_x, src_y, src_z)
   CALL ma_get_host_coordinates(dest,dest_x, dest_y, dest_z)   

   tmp = REAL(src_x-dest_x)/REAL(ma_nodes_x())
   x=NINT((tmp-NINT(tmp))*REAL(ma_nodes_x()))
   tmp = REAL(src_y-dest_y)/REAL(ma_nodes_y())
   y=NINT((tmp-NINT(tmp))*REAL(ma_nodes_y()))
   tmp = REAL(src_z-dest_z)/REAL(ma_nodes_z())
   z=NINT((tmp-NINT(tmp))*REAL(ma_nodes_z()))
  
   nhops = x+y+z 
#else
   nhops = 0
#endif

  END FUNCTION ma_hops

! *****************************************************************************
!> \brief Compute the hop-bytes metric for the system
!> \param mp_group the communicator 
!> \retval hopsbyte ...
! *****************************************************************************
 FUNCTION ma_hop_bytes (mp_group) RESULT (hopsbyte)
  REAL                         :: hopsbyte
  INTEGER, INTENT(IN)          :: mp_group

#if defined __GEMINI || __SEASTAR     
  INTEGER                                  :: numtask,source, &
                                              ipe,jpe,x,y,z, src_x, src_y,&
                                              src_z,dest_x, dest_y, dest_z, &
                                              src,dest,sum_hb, sum_hb_ngb
  REAL                                     :: tmp

  CALL mp_environ(numtask,source,mp_group)

  tmp = 0
  sum_hb = 0
  sum_hb_ngb = 0

  DO ipe=1,numtask
    src = allocated_topology(ipe)
   DO jpe=1,SIZE(comm_graph,2)
    dest = allocated_topology(comm_graph(ipe-1,jpe-1)+1)
    IF (src .NE. dest) THEN
     CALL ma_get_host_coordinates(src,src_x, src_y, src_z)
     CALL ma_get_host_coordinates(dest,dest_x, dest_y, dest_z)
!     if(source == 0) write(*,*) 'src ', src_x,i src_y, src_z,' dest ', dest_x, dest_y, dest_z

     tmp = REAL(src_x-dest_x)/REAL(ma_nodes_x())
     x=NINT((tmp-NINT(tmp))*REAL(ma_nodes_x()))
     tmp = REAL(src_y-dest_y)/REAL(ma_nodes_y())
     y=NINT((tmp-NINT(tmp))*REAL(ma_nodes_y()))
     tmp = REAL(src_z-dest_z)/REAL(ma_nodes_z())
     z=NINT((tmp-NINT(tmp))*REAL(ma_nodes_z()))

 !    if(source == 0) write(*,*) ABS(x),ABS(y),ABS(z)

     sum_hb_ngb = sum_hb_ngb + ABS(x) + ABS(y) + ABS(z)
    ENDIF
   ENDDO
    sum_hb = sum_hb + sum_hb_ngb
    sum_hb_ngb = 0
  ENDDO
 
  hopsbyte = sum_hb/2

  IF(hopsbyte .LT. 0) hopsbyte = 0
#else
  hopsbyte = 0  
#endif 

 END FUNCTION ma_hop_bytes 

! *****************************************************************************
!> \brief Compute the hopes-byte metric for the system
!> \param mp_group the communicator 
!> \retval hopsbyte ...
! *****************************************************************************
 FUNCTION ma_hops_byte (mp_group) RESULT (hopsbyte)
  REAL                         :: hopsbyte
  INTEGER,INTENT(IN)           :: mp_group

#if defined __GEMINI || __SEASTAR     
  INTEGER                                  :: numtask,source, &
                                              ipe,jpe,x,y,z, src_x, src_y,&
                                              src_z,dest_x, dest_y, dest_z, &
                                              src,dest,sum_hb,sum_comm
  REAL                                     :: tmp

  CALL mp_environ(numtask,source,mp_group)

  sum_hb = 1 
  sum_comm = 1

  DO ipe=1,numtask
    src = allocated_topology(ipe)
   DO jpe=1,SIZE(comm_graph,2)
    dest = allocated_topology(comm_graph(ipe-1,jpe-1)+1)
    IF ( src .NE. dest ) THEN
     CALL ma_get_host_coordinates(src,src_x, src_y, src_z)
     CALL ma_get_host_coordinates(dest,dest_x, dest_y, dest_z)
     tmp = REAL(src_x-dest_x)/REAL(ma_nodes_x())
     x=NINT((tmp-NINT(tmp))*REAL(ma_nodes_x()))
     tmp = REAL(src_y-dest_y)/REAL(ma_nodes_y())
     y=NINT((tmp-NINT(tmp))*REAL(ma_nodes_y()))
     tmp = REAL(src_z-dest_z)/REAL(ma_nodes_z())
     z=NINT((tmp-NINT(tmp))*REAL(ma_nodes_z()))
     sum_hb = sum_hb + ABS(x) + ABS(y*6) + ABS(z*4)
     sum_comm = sum_comm + 1
    ENDIF
   ENDDO
  ENDDO

  hopsbyte = sum_hb/sum_comm
#else
  hopsbyte = 0.0
#endif 
 END FUNCTION ma_hops_byte

! *****************************************************************************
!> \brief ...
!> \param mp_group ...
!> \retval distbytes ...
! *****************************************************************************
 FUNCTION ma_dist_bytes (mp_group) RESULT (distbytes)
  REAL                          ::   distbytes
  INTEGER,INTENT(IN)            ::   mp_group 

#if defined __NET
  INTEGER           ::   numtask,source, jpe, &
                         ipe, src, dest, sum_comm
  REAL              ::   sum_db

  sum_db = 1
  sum_comm = 1
  CALL mp_environ(numtask,source,mp_group)

  DO ipe=1,numtask
    src = allocated_topology(ipe)
   DO jpe=1,SIZE(comm_graph,2)
    dest = allocated_topology(comm_graph(ipe-1,jpe-1)+1)
    IF ( src .NE. dest ) THEN
     sum_db = sum_db + dist_nodes(ipe,comm_graph(ipe-1,jpe-1)+1) 
     sum_comm = sum_comm + 1       
    ENDIF
   ENDDO
  ENDDO    

  distbytes = sum_db/sum_comm

#else
  distbytes = 0
#endif  

 END FUNCTION ma_dist_bytes 


! *****************************************************************************
!> \brief Generates the comm graph dim  
!> \param npcol ...
!> \param nprow ...
!> \param numtask ...
! *****************************************************************************
SUBROUTINE ma_2dgrid_dimensions(npcol,nprow,numtask)
    INTEGER, INTENT(OUT)                     :: npcol, nprow
    INTEGER, INTENT(IN)                      :: numtask

    INTEGER                                  :: gcd_max, ipe, jpe

  npcol = 0
  nprow = 0 
  gcd_max=-1

 DO ipe=1,CEILING(SQRT(REAL(numtask,dp)))
    jpe=numtask/ipe
    IF(ipe*jpe.NE.numtask) CYCLE
  IF (gcd(ipe,jpe)>=gcd_max) THEN
    nprow=ipe
    npcol=jpe
    gcd_max=gcd(ipe,jpe)
  ENDIF
 ENDDO

END SUBROUTINE ma_2dgrid_dimensions


! *****************************************************************************
!> \brief Generates the communication graph for alltoall pattern
!> \param mp_group the communicator 
! *****************************************************************************
SUBROUTINE ma_comm_pattern_alltoall (mp_group)
    INTEGER                                  :: mp_group

#if defined(__parallel)
    INTEGER                                  :: istat, numtask, source, &
                                                ipe, jpe
 CALL mp_environ(numtask,source,mp_group)
 ALLOCATE (comm_graph(0:numtask-1,0:numtask-1),STAT=istat)
  IF (istat /= 0) THEN
        CALL ma_error_stop(ma_error_allocation)
  END IF

  DO ipe=0,numtask-1 
   DO jpe=0,numtask-1
    comm_graph(ipe,jpe) = jpe
   ENDDO
  ENDDO
#endif

END SUBROUTINE ma_comm_pattern_alltoall 


! *****************************************************************************
!> \brief Generates the communication graph for cannon algorithm
!> \param mp_group the communicator 
! *****************************************************************************
 SUBROUTINE ma_comm_pattern_cannon (mp_group)
    INTEGER                                  :: mp_group

#if defined(__parallel)
    INTEGER                                  :: down, gcd_max, i, ipe, istat, &
                                                jpe, left, mycol, myrow, &
                                                npcol, nprow, numtask, right, &
                                                source, up
 CALL mp_environ(numtask,source,mp_group)

 IF (.NOT. ASSOCIATED(comm_graph))THEN

#ifdef __DEBUG_HWTOPO
   IF(source==0) WRITE(*,*) 'Generating the communication graph'
#endif

 !Cannon algorithm - each MPI has only for neighbors 
 ALLOCATE (comm_graph(0:numtask-1,0:3),STAT=istat)
  IF (istat /= 0) THEN
        CALL ma_error_stop(ma_error_allocation)
  END IF
  comm_graph(0:numtask-1,0:3) = 0

 gcd_max=-1
 DO ipe=1,CEILING(SQRT(REAL(numtask,dp)))
    jpe=numtask/ipe
    IF(ipe*jpe.NE.numtask) CYCLE
  IF (gcd(ipe,jpe)>=gcd_max) THEN
    nprow=ipe
    npcol=jpe
    gcd_max=gcd(ipe,jpe)
  ENDIF
 ENDDO

 DO i=0,numtask-1
  myrow = MOD(i,nprow)
  mycol = i/npcol

  left  = mycol*myrow + MOD((mycol-1+npcol),npcol)
  right = mycol*myrow + MOD((mycol+1+npcol),npcol)
  up    = MOD((myrow+1+nprow),nprow)*npcol+mycol
  down  = MOD((myrow-1+nprow),nprow)*npcol+mycol

  !save the neighbors
  comm_graph(i,0) = left
  comm_graph(i,1) = right
  comm_graph(i,2) = up
  comm_graph(i,3) = down
 ENDDO

ENDIF
#endif

 END SUBROUTINE ma_comm_pattern_cannon

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Get the direction to for each communication
! *****************************************************************************
!> \brief ...
!> \param src ...
!> \param dest ...
!> \param xdir ...
!> \param ydir ...
!> \param zdir ...
! *****************************************************************************
 SUBROUTINE get_direction(src,dest,xdir,ydir,zdir)
   INTEGER                        :: src, dest, xdir, ydir, zdir

#if defined __GEMINI || __SEASTAR   
   INTEGER                        :: x,y,z

   x = MOD((net_topology(dest)%coord_x - net_topology(src)%coord_x),ma_nodes_x())
   IF(x .LT. ma_nodes_x()/2) THEN
     xdir = 1 
   ELSE
     xdir = -1
   ENDIF 

   y = MOD((net_topology(dest)%coord_y - net_topology(src)%coord_y),ma_nodes_y())
   IF(y .LT. ma_nodes_y()/2) THEN
     ydir = 1
   ELSE
     ydir = -1
   ENDIF

   z = MOD((net_topology(dest)%coord_z - net_topology(src)%coord_z),ma_nodes_z())
   IF(z .LT. ma_nodes_z()/2) THEN
     zdir = 1
   ELSE
     zdir = -1
   ENDIF
#endif

 END SUBROUTINE get_direction

! *****************************************************************************
!> \brief Get the contention for a given comm pattern, use the hops-byte metric
!> \param contention The total contention on each direction of the 3D torus 
!> \param mp_group   The MPI group for communication 
! *****************************************************************************
 SUBROUTINE ma_compute_contention (contention,mp_group)
  REAL, DIMENSION(:)          ::  contention
  INTEGER                        ::  mp_group
#if defined __GEMINI || __SEASTAR  
  contention(:) = 0
  contention(1) = ma_hop_bytes(mp_group)
#elif defined __NET
  contention(:) = 0
  contention(1) = ma_dist_bytes(mp_group)
#else
 contention(:) = 0
#endif
 END SUBROUTINE ma_compute_contention

! *****************************************************************************
!> \brief ...
!> \param mp_comm ...
!> \param kernel ...
! *****************************************************************************
 SUBROUTINE ma_comm_graph (mp_comm,kernel)
  INTEGER, INTENT(IN)      :: mp_comm,kernel

#if defined(__parallel)  
  SELECT CASE (kernel)
    CASE (cannon_graph)
      CALL ma_comm_pattern_cannon(mp_comm) 
    CASE (complete_graph)
      CALL ma_comm_pattern_alltoall(mp_comm)      
    CASE DEFAULT
      CALL ma_error_stop(ma_error_option)
  END SELECT 
#endif
   
 END SUBROUTINE 

END MODULE ma_topology   
