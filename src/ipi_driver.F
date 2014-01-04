! *****************************************************************************
!> \brief Driver mode - To communicate with i-PI Python wrapper
!> \par History
!>      none
!> \author Michele Ceriotti 03.2012
! *****************************************************************************
MODULE ipi_driver

  USE cell_types,                      ONLY: cell_clone,&
                                             cell_create,&
                                             cell_release,&
                                             cell_type,&
                                             compare_cells,&
                                             init_cell,&
                                             real_to_scaled,&
                                             scaled_to_real                                             
  USE cp_external_control,             ONLY: external_control
  USE cp_subsys_types,                 ONLY: cp_subsys_release,&
                                             cp_subsys_type
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp                                             
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE iso_c_binding 
  USE global_types,                    ONLY: global_environment_type
  USE f77_interface,                   ONLY: default_para_env
  USE force_env_methods,               ONLY: force_env_calc_energy_force
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type, &
                                             force_env_set_cell
  USE message_passing,                 ONLY: mp_bcast, mp_sync
  USE virial_types,                    ONLY: virial_type 

#include "cp_common_uses.h"                                               
                                             
  IMPLICIT NONE
  
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ipi_driver'

  PUBLIC :: run_driver

  interface
  subroutine usleep(useconds) bind(C)
  ! integer(c_int) function usleep(useconds) bind(C)
  use iso_c_binding
  implicit none
  integer(c_int32_t), value :: useconds
  ! end function
  end subroutine
  end interface

  CONTAINS
! *****************************************************************************
!> \brief
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


    CHARACTER(len=default_path_length)       :: drv_hostname, c_hostname
    INTEGER                                  :: drv_port, ip, ii, idir, nat
    LOGICAL                                  :: drv_unix, should_stop
    
    !MIK DRIVER
    ! server address parsing
    CHARACTER*1024  :: serveraddr, host
    INTEGER         :: socket, port, inet, nread, readbuffer, slock, swait, uwait
    ! buffers and temporaries for communication
    INTEGER, PARAMETER :: MSGLEN=12
    LOGICAL         :: isinit=.true., hasdata=.false., ionode=.false.
    CHARACTER(LEN=default_string_length)     :: header
    CHARACTER*1024  :: parbuffer
    REAL (KIND=dp)  :: cellh(3,3), cellih(3,3), vir(3,3), pot
    REAL(KIND=dp) , ALLOCATABLE :: combuf(:)
    REAL(KIND=dp)   :: sigma(3,3)
    TYPE(section_vals_type), POINTER         :: drv_section, motion_section
    TYPE(virial_type), POINTER               :: virial
    ! access cp2k structures
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(cell_type), POINTER                 :: cpcell


    ionode=(default_para_env%source==default_para_env%mepos)

    ! reads driver parameters from input
    motion_section => section_vals_get_subs_vals(force_env%root_section,"MOTION",error=error)
    drv_section     => section_vals_get_subs_vals(motion_section,"DRIVER",error=error)

    CALL section_vals_val_get(drv_section,"HOST",c_val=drv_hostname,error=error)
    CALL section_vals_val_get(drv_section,"PORT",i_val=drv_port,error=error)
    CALL section_vals_val_get(drv_section,"UNIX",l_val=drv_unix,error=error)

    
    ! opens the socket
    socket=0    
    inet=1
    if (ionode) then
       write(*,*) "@ i-PI DRIVER BEING LOADED"
       write(*,*) "@ INPUT DATA: ", TRIM(drv_hostname), drv_port, drv_unix                          
       c_hostname=TRIM(drv_hostname)//achar(0)
       CALL open_socket(socket, .not. drv_unix, drv_port, c_hostname) 
    endif    
    
    !now we have a socket, so we can initialize the CP2K environments. 
    NULLIFY(cpcell)
    call cell_create(cpcell,error=error)
    uwait=10000  ! number of MICROseconds to be waited in filesystem lock
    driver_loop: DO
      ! do communication on master node only...
      header = ""

      ! this syncs the processes, possibly (see sockets.c) without calling MPI_Barrier, 
      ! which is nice as MPI_barrier eats up a lot of CPU for nothing
      inet=slock(default_para_env%source, default_para_env%mepos) 
      CALL mp_sync(default_para_env%group)

      if (ionode) nread=readbuffer(socket, header, MSGLEN)
      if (ionode)  write(0,*) "returned from readbuffer" 

      inet=swait(uwait, default_para_env%source, default_para_env%mepos)
      CALL mp_sync(default_para_env%group)

      call mp_bcast(nread,default_para_env%source, default_para_env%group)
      if (nread .eq. 0) then
        if (ionode) write(*,*) " @ DRIVER MODE: Could not read from socket, exiting now."
        exit
      endif

      call mp_bcast(header,default_para_env%source, default_para_env%group)
      
      if (ionode) write(*,*) " @ DRIVER MODE: Message from server: ", trim(header)
      if (trim(header) == "STATUS") then          

         inet=slock(default_para_env%source, default_para_env%mepos) 
         CALL mp_sync(default_para_env%group)
         if (ionode) then  ! does not  need init (well, maybe it should, just to check atom numbers and the like... )
            if (hasdata) then
               call writebuffer(socket,"HAVEDATA    ",MSGLEN)
            else
               call writebuffer(socket,"READY       ",MSGLEN)
            endif
         endif
         inet=swait(uwait,default_para_env%source, default_para_env%mepos)
         CALL mp_sync(default_para_env%group)
      else if (trim(header) == "POSDATA") then              
         if (ionode) then        
            nread=readbuffer(socket, cellh, 9*8)
            nread=readbuffer(socket, cellih, 9*8)
            nread=readbuffer(socket, nat, 4)
            cellh=transpose(cellh)
            cellih=transpose(cellih)
         endif
         call mp_bcast(cellh,default_para_env%source, default_para_env%group)
         call mp_bcast(cellih,default_para_env%source, default_para_env%group)
         call mp_bcast(nat,default_para_env%source, default_para_env%group)
         if (.not.allocated(combuf)) allocate(combuf(3*nat))
         if (ionode) nread=readbuffer(socket, combuf, nat*3*8)
         call mp_bcast(combuf,default_para_env%source, default_para_env%group)
         
         CALL force_env_get(force_env,subsys=subsys,error=error)
         if (nat/=subsys%particles%n_els) WRITE(*,*) &
                " @DRIVER MODE: Uh-oh! Particle number mismatch between i-pi and cp2k input!"
         ii=0
         DO ip=1,subsys%particles%n_els
          DO idir=1,3
             ii=ii+1
             subsys%particles%els(ip)%r(idir)=combuf(ii)
          END DO
         END DO
         CALL init_cell(cpcell, hmat=cellh)
         CALL force_env_set_cell(force_env,cell=cpcell,error=error)
         
         CALL force_env_calc_energy_force(force_env,calc_force=.TRUE. ,error=error)
    
         if (ionode) write(*,*) " @ DRIVER MODE: Received positions "
         
         combuf=0
         ii=0
         DO ip=1,subsys%particles%n_els
          DO idir=1,3
             ii=ii+1
             combuf(ii)=subsys%particles%els(ip)%f(idir)
          END DO
         END DO         
         CALL force_env_get(force_env, potential_energy=pot, error=error)
         CALL force_env_get(force_env,cell=cpcell, virial=virial, error=error)
         vir = transpose(virial%pv_virial)
         
         CALL external_control(should_stop,"DPI",globenv=globenv,error=error)
         IF (should_stop) EXIT
       
         hasdata=.true.
      else if (trim(header)=="GETFORCE") then
         if (ionode) write(*,*) " @ DRIVER MODE: Returning v,forces,stress "
         if (ionode) then     
            call writebuffer(socket,"FORCEREADY  ",MSGLEN)            
            call writebuffer(socket,pot,8)
            call writebuffer(socket,nat,4)            
            call writebuffer(socket,combuf,3*nat*8)
            call writebuffer(socket,vir,9*8)

            ! i-pi can also receive an arbitrary string, that will be printed out to the "extra" 
            ! trajectory file. this is useful if you want to return additional information, e.g.
            ! atomic charges, wannier centres, etc. one must return the number of characters, then
            ! the string. here we just send back zero characters.            
            nat=0
            call writebuffer(socket,nat,4)  ! writes out zero for the length of the "extra" field (not implemented yet!)
         endif
         hasdata=.false.
      else 
         if (ionode) write(*,*) " @DRIVER MODE:  Socket disconnected, time to exit. "
         exit
      endif
    ENDDO driver_loop  
    
  END SUBROUTINE run_driver
  
END MODULE ipi_driver
