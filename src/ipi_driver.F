! *****************************************************************************
!> \brief Driver mode - To communicate with i-PI Python wrapper
!> \par History
!>      none
!> \author Michele Ceriotti 03.2012
! *****************************************************************************
MODULE ipi_driver
  USE cell_types,                      ONLY: cell_create,&
                                             cell_type,&
                                             init_cell
  USE cp_external_control,             ONLY: external_control
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_set,&
                                             cp_subsys_type
  USE f77_interface,                   ONLY: default_para_env
  USE force_env_methods,               ONLY: force_env_calc_energy_force
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE global_types,                    ONLY: global_environment_type
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE iso_c_binding
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_sync
  USE termination,                     ONLY: stop_program
  USE virial_types,                    ONLY: virial_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ipi_driver'

  PUBLIC :: run_driver

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
    SUBROUTINE open_socket(psockfd, inet, port, host) BIND(C)
      USE ISO_C_BINDING
    INTEGER(KIND=C_INT)                      :: psockfd, inet, port
    CHARACTER(KIND=C_CHAR), DIMENSION(*)     :: host

    END SUBROUTINE open_socket

    
    SUBROUTINE writebuffer_csocket(psockfd, pdata, plen) BIND(C, name="writebuffer")
      USE ISO_C_BINDING
    INTEGER(KIND=C_INT)                      :: psockfd
    TYPE(C_PTR), VALUE                       :: pdata
    INTEGER(KIND=C_INT)                      :: plen

    END SUBROUTINE writebuffer_csocket       

    SUBROUTINE readbuffer_csocket(psockfd, pdata, plen) BIND(C, name="readbuffer")
      USE ISO_C_BINDING
    INTEGER(KIND=C_INT)                      :: psockfd
    TYPE(C_PTR), VALUE                       :: pdata
    INTEGER(KIND=C_INT)                      :: plen

    END SUBROUTINE readbuffer_csocket   
  END INTERFACE

  CONTAINS

! *****************************************************************************
!> \brief ...
!> \param psockfd ...
!> \param fdata ...
! *****************************************************************************
  SUBROUTINE writebuffer_d (psockfd, fdata)
      USE ISO_C_BINDING
    INTEGER, INTENT(IN)                      :: psockfd
    REAL(KIND=dp), INTENT(IN)                :: fdata

    REAL(KIND=C_DOUBLE), TARGET              :: cdata

      cdata = fdata
      CALL writebuffer_csocket(psockfd, c_loc(cdata), 8)
  END SUBROUTINE

! *****************************************************************************
!> \brief ...
!> \param psockfd ...
!> \param fdata ...
! *****************************************************************************
  SUBROUTINE writebuffer_i (psockfd, fdata)
      USE ISO_C_BINDING
    INTEGER, INTENT(IN)                      :: psockfd, fdata

    INTEGER(KIND=C_INT), TARGET              :: cdata

      cdata = fdata
      CALL writebuffer_csocket(psockfd, c_loc(cdata), 4)
  END SUBROUTINE

! *****************************************************************************
!> \brief ...
!> \param psockfd ...
!> \param fstring ...
!> \param plen ...
! *****************************************************************************
  SUBROUTINE writebuffer_s (psockfd, fstring, plen)
      USE ISO_C_BINDING
    INTEGER, INTENT(IN)                      :: psockfd
    CHARACTER(LEN=*), INTENT(IN)             :: fstring
    INTEGER, INTENT(IN)                      :: plen

    INTEGER                                  :: i
    CHARACTER(LEN=1, KIND=C_CHAR), TARGET    :: cstring(plen)

      DO i = 1,plen
         cstring(i) = fstring(i:i)
      ENDDO
      CALL writebuffer_csocket(psockfd, c_loc(cstring(1)), plen)
  END SUBROUTINE

! *****************************************************************************
!> \brief ...
!> \param psockfd ...
!> \param fdata ...
!> \param plen ...
! *****************************************************************************
  SUBROUTINE writebuffer_dv(psockfd, fdata, plen)
      USE ISO_C_BINDING  
    INTEGER, INTENT(IN)                      :: psockfd, plen
    REAL(KIND=dp), INTENT(IN), TARGET        :: fdata(plen)

      CALL writebuffer_csocket(psockfd, c_loc(fdata(1)), 8*plen)
  END SUBROUTINE

! *****************************************************************************
!> \brief ...
!> \param psockfd ...
!> \param fdata ...
! *****************************************************************************
  SUBROUTINE readbuffer_d (psockfd, fdata)
      USE ISO_C_BINDING
    INTEGER, INTENT(IN)                      :: psockfd
    REAL(KIND=dp), INTENT(OUT)               :: fdata

    REAL(KIND=C_DOUBLE), TARGET              :: cdata

      CALL readbuffer_csocket(psockfd, c_loc(cdata), 8)
      fdata=cdata
  END SUBROUTINE

! *****************************************************************************
!> \brief ...
!> \param psockfd ...
!> \param fdata ...
! *****************************************************************************
  SUBROUTINE readbuffer_i (psockfd, fdata)
      USE ISO_C_BINDING
    INTEGER, INTENT(IN)                      :: psockfd
    INTEGER, INTENT(OUT)                     :: fdata

    INTEGER(KIND=C_INT), TARGET              :: cdata

      CALL readbuffer_csocket(psockfd, c_loc(cdata), 4)
      fdata = cdata
  END SUBROUTINE

! *****************************************************************************
!> \brief ...
!> \param psockfd ...
!> \param fstring ...
!> \param plen ...
! *****************************************************************************
  SUBROUTINE readbuffer_s (psockfd, fstring, plen)
      USE ISO_C_BINDING
    INTEGER, INTENT(IN)                      :: psockfd
    CHARACTER(LEN=*), INTENT(OUT)            :: fstring
    INTEGER, INTENT(IN)                      :: plen

    INTEGER                                  :: i
    CHARACTER(LEN=1, KIND=C_CHAR), TARGET    :: cstring(plen)

      CALL readbuffer_csocket(psockfd, c_loc(cstring(1)), plen)
      fstring=""   
      WRITE(*,*) "READS", plen, cstring
      DO i = 1,plen
         fstring(i:i) = cstring(i)
      ENDDO
  END SUBROUTINE

! *****************************************************************************
!> \brief ...
!> \param psockfd ...
!> \param fdata ...
!> \param plen ...
! *****************************************************************************
  SUBROUTINE readbuffer_dv(psockfd, fdata, plen)
      USE ISO_C_BINDING  
    INTEGER, INTENT(IN)                      :: psockfd, plen
    REAL(KIND=dp), INTENT(OUT), TARGET       :: fdata(plen)

      CALL readbuffer_csocket(psockfd, c_loc(fdata(1)), 8*plen)
  END SUBROUTINE

! *****************************************************************************
!> \brief ...
!> \param force_env ...
!> \param globenv ...
!> \param error ...
!> \par History
!>       12.2013 included in repository 
!> \author Ceriotti
! *****************************************************************************
 
  SUBROUTINE run_driver ( force_env, globenv, error )
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'run_driver', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: MSGLEN = 12

    CHARACTER(len=default_path_length)       :: c_hostname, drv_hostname
    CHARACTER(LEN=default_string_length)     :: header
    INTEGER                                  :: drv_port, i_drv_unix, idir, &
                                                ii, inet, ip, nat, socket
    LOGICAL                                  :: drv_unix, hasdata = .FALSE., &
                                                ionode = .FALSE., should_stop
    REAL(KIND=dp)                            :: cellh(3,3), cellih(3,3), &
                                                mxmat(9), pot, vir(3,3)
    REAL(KIND=dp), ALLOCATABLE               :: combuf(:)
    TYPE(cell_type), POINTER                 :: cpcell
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(section_vals_type), POINTER         :: drv_section, motion_section
    TYPE(virial_type), POINTER               :: virial

! server address parsing
! buffers and temporaries for communication
! access cp2k structures

    ionode=(default_para_env%source==default_para_env%mepos)

    ! reads driver parameters from input
    motion_section => section_vals_get_subs_vals(force_env%root_section,"MOTION",error=error)
    drv_section     => section_vals_get_subs_vals(motion_section,"DRIVER",error=error)

    CALL section_vals_val_get(drv_section,"HOST",c_val=drv_hostname,error=error)
    CALL section_vals_val_get(drv_section,"PORT",i_val=drv_port,error=error)
    CALL section_vals_val_get(drv_section,"UNIX",l_val=drv_unix,error=error)

#ifdef __NO_IPI_DRIVER
    CALL stop_program(routineN,moduleN,__LINE__,"CP2K was compiled with the __NO_IPI_DRIVER option!")
#else
    
    ! opens the socket
    socket=0    
    inet=1
    i_drv_unix = 1   ! a bit convoluted. socket.c uses a different convention...
    IF (drv_unix) i_drv_unix = 0 
    IF (ionode) THEN
       WRITE(*,*) "@ i-PI DRIVER BEING LOADED"
       WRITE(*,*) "@ INPUT DATA: ", TRIM(drv_hostname), drv_port, drv_unix                          
       c_hostname=TRIM(drv_hostname)//ACHAR(0)
       CALL open_socket(socket,i_drv_unix, drv_port, c_hostname) 
    ENDIF    
    
    !now we have a socket, so we can initialize the CP2K environments. 
    NULLIFY(cpcell)
    CALL cell_create(cpcell,error=error)
    driver_loop: DO
      ! do communication on master node only...
      header = ""

      CALL mp_sync(default_para_env%group)

      IF (ionode) CALL readbuffer(socket, header, MSGLEN)
      
      CALL mp_sync(default_para_env%group)

      CALL mp_bcast(header,default_para_env%source, default_para_env%group)
      
      IF (ionode) WRITE(*,*) " @ DRIVER MODE: Message from server: ", TRIM(header)
      IF (TRIM(header) == "STATUS") THEN          

         CALL mp_sync(default_para_env%group)
         IF (ionode) THEN  ! does not  need init (well, maybe it should, just to check atom numbers and the like... )
            IF (hasdata) THEN
               CALL writebuffer(socket,"HAVEDATA    ",MSGLEN)
            ELSE
               CALL writebuffer(socket,"READY       ",MSGLEN)
            ENDIF
         ENDIF
         CALL mp_sync(default_para_env%group)
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
         CALL mp_bcast(cellh,default_para_env%source, default_para_env%group)
         CALL mp_bcast(cellih,default_para_env%source, default_para_env%group)
         CALL mp_bcast(nat,default_para_env%source, default_para_env%group)
         IF (.not.ALLOCATED(combuf)) ALLOCATE(combuf(3*nat))
         IF (ionode) CALL readbuffer(socket, combuf, nat*3)
         CALL mp_bcast(combuf,default_para_env%source, default_para_env%group)
         
         CALL force_env_get(force_env,subsys=subsys,error=error)
         IF (nat/=subsys%particles%n_els) WRITE(*,*) &
                " @DRIVER MODE: Uh-oh! Particle number mismatch between i-PI and cp2k input!"
         ii=0
         DO ip=1,subsys%particles%n_els
          DO idir=1,3
             ii=ii+1
             subsys%particles%els(ip)%r(idir)=combuf(ii)
          END DO
         END DO
         CALL init_cell(cpcell, hmat=cellh)
         CALL cp_subsys_set(subsys, cell=cpcell, error=error)
         
         CALL force_env_calc_energy_force(force_env,calc_force=.TRUE. ,error=error)
    
         IF (ionode) WRITE(*,*) " @ DRIVER MODE: Received positions "
         
         combuf=0
         ii=0
         DO ip=1,subsys%particles%n_els
          DO idir=1,3
             ii=ii+1
             combuf(ii)=subsys%particles%els(ip)%f(idir)
          END DO
         END DO         
         CALL force_env_get(force_env, potential_energy=pot, error=error)
         CALL force_env_get(force_env,cell=cpcell, error=error)
         CALL cp_subsys_get(subsys, virial=virial, error=error)
         vir = TRANSPOSE(virial%pv_virial)
         
         CALL external_control(should_stop,"IPI",globenv=globenv,error=error)
         IF (should_stop) EXIT
       
         hasdata=.TRUE.
      ELSE IF (TRIM(header)=="GETFORCE") THEN
         IF (ionode) WRITE(*,*) " @ DRIVER MODE: Returning v,forces,stress "
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
         IF (ionode) WRITE(*,*) " @DRIVER MODE:  Socket disconnected, time to exit. "
         EXIT
      ENDIF
    ENDDO driver_loop  
#endif
    
  END SUBROUTINE run_driver
END MODULE ipi_driver
