!!****** cp2k/cp2k [1.0] *
!!
!!   NAME
!!     cp2k, The main program
!!
!!   COPYRIGHT
!!I---------------------------------------------------------------------------I
!!I                                                                           I
!!I CP2K: A general program to perform molecular dynamics simulations         I
!!I Copyright (C) 2000,2001,2002,2003,2004  CP2K developers group             I
!!I                                                                           I
!!I This program is free software; you can redistribute it and/or modify      I
!!I it under the terms of the GNU General Public License as published by      I
!!I the Free Software Foundation; either version 2 of the License, or         I
!!I (at your option) any later version.                                       I
!!I                                                                           I
!!I This program is distributed in the hope that it will be useful,           I
!!I but WITHOUT ANY WARRANTY; without even the implied warranty of            I
!!I MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             I
!!I GNU General Public License for more details.                              I
!!I                                                                           I
!!I You should have received a copy of the GNU General Public License         I
!!I along with this program; if not, write to the Free Software               I
!!I Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                 I
!!I                                                                           I
!!I---------------------------------------------------------------------------I
!!
!!   FUNCTION
!!      starts the program
!!
!!   AUTHOR
!!      Joost VandeVondele 
!!
!!   NOTES
!!      should be kept as lean as possible.
!!      see cp2k_run for more comments
!!
!------------------------------------------------------------------------------
PROGRAM cp2k
  USE cp2k_runs,         ONLY : cp2k_run
  USE machine,           ONLY : m_iargc,m_getarg
  USE message_passing,   ONLY : mp_world_init, &
                                mp_world_finalize
  USE kinds,             ONLY : default_string_length, &
                                default_output_unit
  IMPLICIT NONE

  CHARACTER(LEN=default_string_length) :: input_file_name
  INTEGER                              :: output_unit
  INTEGER                              :: mpi_comm_default
  LOGICAL                              :: cp2k_run_OK

  ! output goes to the screen by default
  output_unit=default_output_unit

  ! get the default system wide communicator
  CALL mp_world_init(mpi_comm_default)

  ! get the input file name, 
  IF (m_iargc()==0) THEN
      input_file_name = "input"
  ELSE
      CALL m_getarg(1,input_file_name)
  ENDIF

  ! do a cp2k run with these arguments
  cp2k_run_OK=cp2k_run(input_file_name,output_unit,mpi_comm_default)

  ! and finalize the mpi communicator
  CALL mp_world_finalize()

END PROGRAM
