!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

MODULE io_parameters
  
  IMPLICIT NONE
  
! Defaults for some global variables
  
! *****************************************************************************
  
! scr: Logical unit number of the screen
! crd: Logical unit number of the coordinate file
! vel: Logical unit number of the velocity file
! ptn: Logical unit number of the pressure tensor file
! ene: Logical unit number of the energy file
! tem: Logical unit number of the temperature file
! dmp: Logical unit number of the dump file
! icrd : frequency of writing coordinate file
! ivel : frequency of writing velocity file
! iptens : frequency of writing pressure tensor file
! itemp : frequency of writing temperature file
! idump : frequency of writing dump file
! iscreen : frequency of writing screen file
  
! *****************************************************************************
  
  PRIVATE
  PUBLIC :: crd_file_name, vel_file_name, ptens_file_name, ener_file_name, &
            temp_file_name, dump_file_name

  PUBLIC :: crd, vel, ptn, ene, tem, scr, dmp

  PUBLIC :: icrd, ivel, iptens, iener, itemp, iscreen, idump
  
  CHARACTER ( LEN = 50 ) :: crd_file_name, &
                            vel_file_name, &
                            ptens_file_name,&
                            ener_file_name, &
                           temp_file_name, &
                           dump_file_name 
  INTEGER, PARAMETER     :: crd = 1, &
                            vel = 2, &
                            ptn = 3, &
                            ene = 4, &
                            tem = 5, & 
                            scr = 6, &
                            dmp = 7
  INTEGER                :: icrd, ivel, iptens,  &
                            iener, itemp, idump, &
                            iscreen
  
END MODULE io_parameters
