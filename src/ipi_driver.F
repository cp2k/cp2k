!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Driver mode - To communicate with i-PI Python wrapper
!> \par History
!>      none
!> \author Michele Ceriotti 03.2012
! *****************************************************************************
MODULE ipi_driver
  USE ISO_C_BINDING,                   ONLY: C_CHAR,&
                                             C_DOUBLE,&
                                             C_INT,&
                                             C_LOC,&
                                             C_NULL_CHAR,&
                                             C_PTR
  USE bibliography,                    ONLY: Ceriotti2014,&
                                             Kapil2016,&
                                             cite_reference
  USE cell_types,                      ONLY: cell_create,&
                                             cell_type,&
                                             init_cell
  USE cp_external_control,             ONLY: external_control
  USE cp_log_handling,                 ONLY: cp_logger_get_default_io_unit
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_set,&
                                             cp_subsys_type
  USE force_env_methods,               ONLY: force_env_calc_energy_force
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE global_types,                    ONLY: global_environment_type
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp,&
                                             int_4
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_irecv,&
                                             mp_send,&
                                             mp_sync,&
                                             mp_testany
  USE virial_types,                    ONLY: virial_type
#include "./base/base_uses.f90"

  IMPLICIT NONE
  
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ipi_driver'

  PUBLIC :: run_driver

#ifndef __NO_IPI_DRIVER
  INTERFACE writebuffer
      MODULE PROCEDURE writebuffer_s, &
                       writebuffer_d, writebuffer_dv, &
                       writebuffer_i
                       
  END INTERFACE 

  INTERFACE readbuffer
      MODULE PROCEDURE readbuffer_s, &
                       readbuffer_dv, readbuffer_d, &
                       readbuffer_i
                       
  END INTERFACE 

  
  INTERFACE
      SUBROUTINE uwait(sec) BIND(C,NAME="uwait")
         USE ISO_C_BINDING, ONLY: C_DOUBLE
    REAL(C_DOUBLE)                           :: sec

      END SUBROUTINE
  END INTERFACE
  
  INTERFACE
    SUBROUTINE open_socket(psockfd, inet, port, host) BIND(C)
      IMPORT
      INTEGER(KIND=C_INT)                      :: psockfd, inet, port
      CHARACTER(KIND=C_CHAR), DIMENSION(*)     :: host

    END SUBROUTINE open_socket

    
    SUBROUTINE writebuffer_csocket(psockfd, pdata, plen) BIND(C, name="writebuffer")
      IMPORT
      INTEGER(KIND=C_INT)                      :: psockfd
      TYPE(C_PTR), VALUE                       :: pdata
      INTEGER(KIND=C_INT)                      :: plen

    END SUBROUTINE writebuffer_csocket       

    SUBROUTINE readbuffer_csocket(psockfd, pdata, plen) BIND(C, name="readbuffer")
      IMPORT
      INTEGER(KIND=C_INT)                      :: psockfd
      TYPE(C_PTR), VALUE                       :: pdata
      INTEGER(KIND=C_INT)                      :: plen

    END SUBROUTINE readbuffer_csocket   
  END INTERFACE
#endif

  CONTAINS
  
#ifndef __NO_IPI_DRIVER
! *****************************************************************************
!> \brief ...
!> \param psockfd ...
!> \param fdata ...
! *****************************************************************************
  SUBROUTINE writebuffer_d (psockfd, fdata)
    INTEGER, INTENT(IN)                      :: psockfd
    REAL(KIND=dp), INTENT(IN)                :: fdata

    CHARACTER(len=*), PARAMETER :: routineN = 'writebuffer_d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    REAL(KIND=C_DOUBLE), TARGET              :: cdata

    CALL timeset(routineN,handle)

    cdata = fdata
    CALL writebuffer_csocket(psockfd, c_loc(cdata), 8)

    CALL timestop(handle)
  END SUBROUTINE

! *****************************************************************************
!> \brief ...
!> \param psockfd ...
!> \param fdata ...
! *****************************************************************************
  SUBROUTINE writebuffer_i (psockfd, fdata)
    INTEGER, INTENT(IN)                      :: psockfd, fdata

    CHARACTER(len=*), PARAMETER :: routineN = 'writebuffer_i', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    INTEGER(KIND=C_INT), TARGET              :: cdata

    CALL timeset(routineN,handle)

    cdata = fdata
    CALL writebuffer_csocket(psockfd, c_loc(cdata), 4)

    CALL timestop(handle)
  END SUBROUTINE

! *****************************************************************************
!> \brief ...
!> \param psockfd ...
!> \param fstring ...
!> \param plen ...
! *****************************************************************************
  SUBROUTINE writebuffer_s (psockfd, fstring, plen)
    INTEGER, INTENT(IN)                      :: psockfd
    CHARACTER(LEN=*), INTENT(IN)             :: fstring
    INTEGER, INTENT(IN)                      :: plen

    CHARACTER(len=*), PARAMETER :: routineN = 'writebuffer_s', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i
    CHARACTER(LEN=1, KIND=C_CHAR), TARGET    :: cstring(plen)

    CALL timeset(routineN,handle)

    DO i = 1,plen
       cstring(i) = fstring(i:i)
    ENDDO
    CALL writebuffer_csocket(psockfd, c_loc(cstring(1)), plen)

    CALL timestop(handle)

  END SUBROUTINE

! *****************************************************************************
!> \brief ...
!> \param psockfd ...
!> \param fdata ...
!> \param plen ...
! *****************************************************************************
  SUBROUTINE writebuffer_dv(psockfd, fdata, plen)
    INTEGER, INTENT(IN)                      :: psockfd, plen
    REAL(KIND=dp), INTENT(IN), TARGET        :: fdata(plen)

    CHARACTER(len=*), PARAMETER :: routineN = 'writebuffer_dv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN,handle)

    CALL writebuffer_csocket(psockfd, c_loc(fdata(1)), 8*plen)

    CALL timestop(handle)
  END SUBROUTINE

! *****************************************************************************
!> \brief ...
!> \param psockfd ...
!> \param fdata ...
! *****************************************************************************
  SUBROUTINE readbuffer_d (psockfd, fdata)
    INTEGER, INTENT(IN)                      :: psockfd
    REAL(KIND=dp), INTENT(OUT)               :: fdata

    CHARACTER(len=*), PARAMETER :: routineN = 'readbuffer_d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    REAL(KIND=C_DOUBLE), TARGET              :: cdata

    CALL timeset(routineN,handle)

    CALL readbuffer_csocket(psockfd, c_loc(cdata), 8)
    fdata=cdata

    CALL timestop(handle)
  END SUBROUTINE

! *****************************************************************************
!> \brief ...
!> \param psockfd ...
!> \param fdata ...
! *****************************************************************************
  SUBROUTINE readbuffer_i (psockfd, fdata)
    INTEGER, INTENT(IN)                      :: psockfd
    INTEGER, INTENT(OUT)                     :: fdata

    CHARACTER(len=*), PARAMETER :: routineN = 'readbuffer_i', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    INTEGER(KIND=C_INT), TARGET              :: cdata

    CALL timeset(routineN,handle)

    CALL readbuffer_csocket(psockfd, c_loc(cdata), 4)
    fdata = cdata

    CALL timestop(handle)
  END SUBROUTINE

! *****************************************************************************
!> \brief ...
!> \param psockfd ...
!> \param fstring ...
!> \param plen ...
! *****************************************************************************
  SUBROUTINE readbuffer_s (psockfd, fstring, plen)
    INTEGER, INTENT(IN)                      :: psockfd
    CHARACTER(LEN=*), INTENT(OUT)            :: fstring
    INTEGER, INTENT(IN)                      :: plen

    CHARACTER(len=*), PARAMETER :: routineN = 'readbuffer_s', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i
    CHARACTER(LEN=1, KIND=C_CHAR), TARGET    :: cstring(plen)

    CALL timeset(routineN,handle)

    CALL readbuffer_csocket(psockfd, c_loc(cstring(1)), plen)
    fstring=""   
    DO i = 1,plen
       fstring(i:i) = cstring(i)
    ENDDO

    CALL timestop(handle)

  END SUBROUTINE

! *****************************************************************************
!> \brief ...
!> \param psockfd ...
!> \param fdata ...
!> \param plen ...
! *****************************************************************************
  SUBROUTINE readbuffer_dv(psockfd, fdata, plen)
    INTEGER, INTENT(IN)                      :: psockfd, plen
    REAL(KIND=dp), INTENT(OUT), TARGET       :: fdata(plen)

    CHARACTER(len=*), PARAMETER :: routineN = 'readbuffer_dv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN,handle)

    CALL readbuffer_csocket(psockfd, c_loc(fdata(1)), 8*plen)

    CALL timestop(handle)

  END SUBROUTINE
#endif

! *****************************************************************************
!> \brief ...
!> \param force_env ...
!> \param globenv ...
!> \par History
!>       12.2013 included in repository 
!> \author Ceriotti
! *****************************************************************************
 
  SUBROUTINE run_driver ( force_env, globenv)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(global_environment_type), POINTER   :: globenv

    CHARACTER(len=*), PARAMETER :: routineN = 'run_driver', &
      routineP = moduleN//':'//routineN

#ifdef __NO_IPI_DRIVER
    INTEGER                                  :: handle
    CALL timeset(routineN,handle)
    CPABORT("CP2K was compiled with the __NO_IPI_DRIVER option!")
    MARK_USED(globenv)
    MARK_USED(force_env)
#else
    INTEGER, PARAMETER                       :: MSGLEN = 12

    CHARACTER(len=default_path_length)       :: c_hostname, drv_hostname
    CHARACTER(LEN=default_string_length)     :: header
    INTEGER                                  :: drv_port, handle, i_drv_unix, &
                                                idir, ii, inet, ip, iwait, &
                                                nat, output_unit, socket, &
                                                wait_req(2)
    INTEGER(KIND=int_4), POINTER             :: wait_msg(:)
    LOGICAL                                  :: drv_unix, fwait, hasdata, &
                                                ionode, should_stop
    REAL(KIND=dp)                            :: cellh(3,3), cellih(3,3), &
                                                mxmat(9), pot, vir(3,3)
    REAL(KIND=dp), ALLOCATABLE               :: combuf(:)
    TYPE(cell_type), POINTER                 :: cpcell
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(section_vals_type), POINTER         :: drv_section, motion_section
    TYPE(virial_type), POINTER               :: virial

    CALL timeset(routineN,handle)

    CALL cite_reference(Ceriotti2014)
    CALL cite_reference(Kapil2016)

! server address parsing
! buffers and temporaries for communication
! access cp2k structures

    CPASSERT(ASSOCIATED(force_env))
    CALL force_env_get(force_env,para_env=para_env)

    hasdata = .FALSE.
    ionode=(para_env%source==para_env%mepos)

    output_unit= cp_logger_get_default_io_unit()

    ! reads driver parameters from input
    motion_section => section_vals_get_subs_vals(force_env%root_section,"MOTION")
    drv_section     => section_vals_get_subs_vals(motion_section,"DRIVER")

    CALL section_vals_val_get(drv_section,"HOST",c_val=drv_hostname)
    CALL section_vals_val_get(drv_section,"PORT",i_val=drv_port)
    CALL section_vals_val_get(drv_section,"UNIX",l_val=drv_unix)

    ! opens the socket
    socket=0    
    inet=1
    i_drv_unix = 1   ! a bit convoluted. socket.c uses a different convention...
    IF (drv_unix) i_drv_unix = 0 
    IF (output_unit>0) THEN
       WRITE(output_unit,*) "@ i-PI DRIVER BEING LOADED"
       WRITE(output_unit,*) "@ INPUT DATA: ", TRIM(drv_hostname), drv_port, drv_unix                          
    ENDIF    

    c_hostname=TRIM(drv_hostname)//C_NULL_CHAR
    IF (ionode) CALL open_socket(socket,i_drv_unix, drv_port, c_hostname) 
    
    NULLIFY(wait_msg)
    ALLOCATE(wait_msg(1))
    !now we have a socket, so we can initialize the CP2K environments. 
    NULLIFY(cpcell)    
    CALL cell_create(cpcell)
    driver_loop: DO
      ! do communication on master node only...
      header = ""

      CALL mp_sync(para_env%group)

      ! non-blocking sync to avoid useless CPU consumption
      IF (ionode) THEN
         CALL readbuffer(socket, header, MSGLEN)
         wait_msg=0
         DO iwait=0,para_env%num_pe-1         
           IF (iwait /= para_env%source) THEN   
              CALL mp_send(msg=wait_msg,dest=iwait,gid=para_env%group,tag=666) 
           ENDIF 
         ENDDO
      ELSE
         CALL mp_irecv(msgout=wait_msg,source=para_env%source,comm=para_env%group,&
                       tag=666, request=wait_req(2))         
         CALL mp_testany(wait_req(2:), flag=fwait)
         DO WHILE(.NOT. fwait)
           CALL mp_testany(wait_req(2:), flag=fwait)           
           CALL uwait(0.01d0)
         ENDDO
      ENDIF

       
      CALL mp_sync(para_env%group)

      CALL mp_bcast(header,para_env%source, para_env%group)
      
      IF (output_unit>0) WRITE(output_unit,*) " @ DRIVER MODE: Message from server: ", TRIM(header)
      IF (TRIM(header) == "STATUS") THEN          

         CALL mp_sync(para_env%group)
         IF (ionode) THEN  ! does not  need init (well, maybe it should, just to check atom numbers and the like... )
            IF (hasdata) THEN
               CALL writebuffer(socket,"HAVEDATA    ",MSGLEN)
            ELSE
               CALL writebuffer(socket,"READY       ",MSGLEN)
            ENDIF
         ENDIF
         CALL mp_sync(para_env%group)
      ELSE IF (TRIM(header) == "POSDATA") THEN              
         IF (ionode) THEN        
            CALL readbuffer(socket, mxmat, 9)
            cellh = RESHAPE(mxmat, (/3,3/) )
            CALL readbuffer(socket, mxmat, 9)
            cellih = RESHAPE(mxmat, (/3,3/) )
            CALL readbuffer(socket, nat)
            cellh=TRANSPOSE(cellh)
            cellih=TRANSPOSE(cellih)
         ENDIF
         CALL mp_bcast(cellh,para_env%source, para_env%group)
         CALL mp_bcast(cellih,para_env%source, para_env%group)
         CALL mp_bcast(nat,para_env%source, para_env%group)
         IF (.not.ALLOCATED(combuf)) ALLOCATE(combuf(3*nat))
         IF (ionode) CALL readbuffer(socket, combuf, nat*3)
         CALL mp_bcast(combuf,para_env%source, para_env%group)
         
         CALL force_env_get(force_env,subsys=subsys)
         IF (nat/=subsys%particles%n_els) &
                CPABORT("@DRIVER MODE: Uh-oh! Particle number mismatch between i-PI and cp2k input!")
         ii=0
         DO ip=1,subsys%particles%n_els
          DO idir=1,3
             ii=ii+1
             subsys%particles%els(ip)%r(idir)=combuf(ii)
          END DO
         END DO
         CALL init_cell(cpcell, hmat=cellh)
         CALL cp_subsys_set(subsys, cell=cpcell)
         
         CALL force_env_calc_energy_force(force_env,calc_force=.TRUE.)
    
         IF (output_unit>0) WRITE(output_unit,*) " @ DRIVER MODE: Received positions "
         
         combuf=0
         ii=0
         DO ip=1,subsys%particles%n_els
          DO idir=1,3
             ii=ii+1
             combuf(ii)=subsys%particles%els(ip)%f(idir)
          END DO
         END DO         
         CALL force_env_get(force_env, potential_energy=pot)
         CALL force_env_get(force_env,cell=cpcell)
         CALL cp_subsys_get(subsys, virial=virial)
         vir = TRANSPOSE(virial%pv_virial)
         
         CALL external_control(should_stop,"IPI",globenv=globenv)
         IF (should_stop) EXIT
       
         hasdata=.TRUE.
      ELSE IF (TRIM(header)=="GETFORCE") THEN
         IF (output_unit>0) WRITE(output_unit,*) " @ DRIVER MODE: Returning v,forces,stress "
         IF (ionode) THEN     
            CALL writebuffer(socket,"FORCEREADY  ",MSGLEN)            
            CALL writebuffer(socket,pot)
            CALL writebuffer(socket,nat)
            CALL writebuffer(socket,combuf,3*nat)
            CALL writebuffer(socket,RESHAPE(vir, (/9/) ),9)

            ! i-pi can also receive an arbitrary string, that will be printed out to the "extra" 
            ! trajectory file. this is useful if you want to return additional information, e.g.
            ! atomic charges, wannier centres, etc. one must return the number of characters, then
            ! the string. here we just send back zero characters.            
            nat=0
            CALL writebuffer(socket,nat)  ! writes out zero for the length of the "extra" field (not implemented yet!)
         ENDIF
         hasdata=.FALSE.
      ELSE 
         IF (output_unit>0) WRITE(output_unit,*) " @DRIVER MODE:  Socket disconnected, time to exit. "
         EXIT
      ENDIF
    ENDDO driver_loop  
#endif

    CALL timestop(handle)
    
  END SUBROUTINE run_driver
END MODULE ipi_driver
