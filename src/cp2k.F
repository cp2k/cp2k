!-----------------------------------------------------------------------------!
!                                                                             !
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!                                                                             !
!   This program is free software; you can redistribute it and/or modify      !
!   it under the terms of the GNU General Public License as published by      !
!   the Free Software Foundation; either version 2 of the License, or         !
!   (at your option) any later version.                                       !
!                                                                             !
!   This program is distributed in the hope that it will be useful,           !
!   but WITHOUT ANY WARRANTY; without even the implied warranty of            !
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             !
!   GNU General Public License for more details.                              !
!                                                                             !
!   You should have received a copy of the GNU General Public License         !
!   along with this program; if not, write to the Free Software               !
!   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                 !
!                                                                             !
!-----------------------------------------------------------------------------!
! 
!!  Main Program
!!       - performs only basic initialisations
!!       - calls major modules

PROGRAM cp2k
  
  USE cp2k_input, ONLY : read_cp2k_section
  USE environment, ONLY : initialisation, trailer
  USE fist, ONLY : fist_main
  USE fist_global, ONLY : set_fist_global
  USE global_types, ONLY : global_environment_type
  USE kinds, ONLY : dbl, print_kind_info
  USE parallel, ONLY : start_parallel, end_parallel
  USE physcon, ONLY : print_physcon
  USE tbmd, ONLY : tbmd_main
  USE tbmd_global, ONLY : set_tbmd_global
  
  IMPLICIT NONE
  
  TYPE ( global_environment_type ) :: globenv
  
!------------------------------------------------------------------------------
  
  CALL start_parallel ( globenv )
  
  CALL initialisation ( globenv )
  
  CALL read_cp2k_section ( globenv )
  
  IF ( globenv % ionode .AND. globenv % print_level > 9 ) &
       CALL print_kind_info ( globenv % scr )
  IF ( globenv % ionode .AND. globenv % print_level > 4 ) &
       CALL print_physcon ( globenv % scr )
  
  IF ( globenv % program_name == 'FIST' ) THEN
     CALL set_fist_global ( globenv )
     CALL fist_main()
     
  ELSE IF ( globenv % program_name == 'TBMD' ) THEN
     CALL set_tbmd_global ( globenv )
     CALL tbmd_main()
  END IF
  
  CALL trailer ( globenv )
  
  CALL end_parallel()
  
END PROGRAM cp2k
