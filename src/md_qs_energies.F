!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2003  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/md_qs_energies [1.0] *
!!
!!   NAME
!!     md_qs_energies
!!
!!   FUNCTION
!!     prints all energy info per timestep to the screen or to
!!     user defined output files
!!
!!   AUTHOR
!!     Joost VandeVondele (copy from md_fist_energies)
!!
!!   MODIFICATIONS
!!     - New MD data are appended to the old data (15.09.2003,MK)
!!
!******************************************************************************

MODULE md_qs_energies

  USE global_types,                    ONLY: global_environment_type
  USE integrator,                      ONLY: pressure
  USE kinds,                           ONLY: wp => dp
  USE linklist_types,                  ONLY: linklist_internal_data_type,&
                                             get_linklist_internal_data
  USE machine,                         ONLY: m_flush
  USE mathconstants,                   ONLY: zero
  USE md,                              ONLY: simulation_parameters_type,&
                                             mdio_parameters_type
  USE md_conserved_quantities,         ONLY: compute_conserved_quantity,&
                                             compute_conserved_quantity_pimd
  USE md_environment_types,            ONLY: md_environment_type,&
                                             get_md_env
  USE message_passing,                 ONLY: mp_bcast,mp_sum, mp_comm_dup
  USE qs_parser,                       ONLY: open_file
  USE simulation_cell,                 ONLY: get_cell_param,&
                                             cell_type
  USE structure_types,                 ONLY: structure_type
  USE particle_types,                  ONLY: particle_type, &
                                             write_particle_coordinates
  USE physcon,                         ONLY: angstrom,&
                                             femtoseconds,&
                                             kelvin
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: module_name = "md_qs_energies"

  REAL(wp) :: econs_avg,t_old

  CHARACTER (LEN=50) :: crd_file_name, vel_file_name, ptens_file_name, &
                        ener_file_name, temp_file_name, dump_file_name
  INTEGER :: crd, vel, ptn, ene, tem, scr, coe, pol, group, itimes_old, source
  INTEGER :: icrd, ivel, iptens, iener, itemp, idump, iscreen
  LOGICAL :: ionode

  PUBLIC :: md_energy_qs,&
            set_qs_energies
!!***
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** md_qs_energies/set_qs_energies [1.0] *
!!
!!   NAME
!!     set_qs_energies
!!
!!   FUNCTION
!!     setup of parameters for qs_energies
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!   NOTES
!!
!!*** *************************************************************************

SUBROUTINE set_qs_energies ( globenv, mdio )

  IMPLICIT NONE

! Arguments
  TYPE ( mdio_parameters_type ), INTENT ( IN ) :: mdio
  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv

!------------------------------------------------------------------------------

  group = globenv % group
  source = globenv % source
  ionode = globenv % ionode
  scr = globenv % scr

  crd_file_name = mdio % crd_file_name
  vel_file_name = mdio % vel_file_name
  ptens_file_name = mdio % ptens_file_name
  ener_file_name = mdio % ener_file_name
  temp_file_name = mdio % temp_file_name
  dump_file_name = mdio % dump_file_name

  icrd = mdio % icrd
  ivel = mdio % ivel
  iptens = mdio % iptens
  iener = mdio % iener
  itemp = mdio % itemp
  idump = mdio % idump
  iscreen = mdio % iscreen

END SUBROUTINE set_qs_energies

!******************************************************************************
!!****** md_qs_energies/md_energy_qs [1.0] *
!!
!!   NAME
!!     energy_qs
!!
!!   FUNCTION
!!     This routine computes the conserved quantity, temperature
!!     and things like that and prints them out
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*****************************************************************************

SUBROUTINE md_energy_qs ( md_env )

  IMPLICIT NONE

! Arguments
  TYPE (md_environment_type), INTENT ( INOUT ) :: md_env

! Locals
  TYPE ( simulation_parameters_type ), POINTER :: simpar
  TYPE ( structure_type ), POINTER :: struc (:)
  INTEGER, POINTER :: itimes
  REAL (wp), POINTER :: constant
  CHARACTER(LEN=200) :: title

  INTEGER i, j, id, n, l, nnodes, beads, iw, handle
  REAL ( wp ) :: alpha, epsilon0, qsq, lup, aup
  REAL ( wp ) :: cons, temp, rcons, ekin
  REAL ( wp ), DIMENSION(3) :: cell, ang
  REAL ( wp ) :: pv_scalar,t
  TYPE(particle_type), DIMENSION(:), POINTER :: particle_set

  LOGICAL :: file_exists

!------------------------------------------------------------------------------

  CALL timeset ( 'energy', 'E', ' ', handle )

  CALL get_md_env ( md_env=md_env,&
                    itimes=itimes, &
                    constant=constant, &
                    simpar=simpar,&
                    struc=struc)

  particle_set => struc(1) % part

  IF (itimes == 0) THEN
    IF (ionode) THEN
!MK   CALL open_md_file(temp_file_name,tem)
      INQUIRE (FILE=TRIM(ener_file_name),EXIST=file_exists)
      CALL open_md_file(ener_file_name,ene)
      IF (file_exists) THEN
        BACKSPACE (UNIT=ene)
        READ (UNIT=ene,FMT=*) itimes_old,t_old
      ELSE
        itimes_old = 0
        t_old = 0.0_wp
      END IF
      CALL open_md_file(crd_file_name,crd)
      CALL open_md_file(vel_file_name,vel)
    END IF
    CALL mp_bcast(itimes_old,source,group)
    CALL mp_bcast(t_old,source,group)
    itimes = itimes_old
    t = t_old
  ELSE
    t = t_old + REAL(itimes-itimes_old,wp)*md_env%simpar%dt*femtoseconds
  END IF

  IF (itimes > itimes_old) THEN

  iw = scr

  CALL compute_conserved_quantity (md_env, md_env % force_env % qs_env % energy % total, cons,&
                                   cons_rel=rcons,&
                                   ekin=ekin,&
                                   temp=temp)

  temp = temp*kelvin

  write(title,'(A4,I8,A4,F12.3,X,A6,F20.9,X,A6,F9.3,X,A1)') " i= ",itimes," t= ",t,"fs E= ", &
                          md_env % force_env % qs_env % energy % total,"au T= ",temp,"K"

!MK  IF (ionode .AND. mod(itimes,itemp)==0) THEN
!MK     WRITE (tem, '(I8,F20.10,F20.10)') itimes,t, temp
!MK     CALL m_flush(tem)
!MK  END IF

  IF (ionode .AND. mod(itimes,iscreen)==0) THEN
     WRITE (scr, '(A1)') " "
     WRITE (scr, '(A6,I8,F12.3,F20.9,F20.9,F20.9,F20.9,F20.12)') " -MD- ",itimes, &
                         t,ekin,temp,md_env % force_env % qs_env % energy % total, cons, cons-constant
     WRITE (scr, '(A1)') " "
  END IF

  IF (ionode .AND. mod(itimes,iener)==0) THEN
     WRITE (ene, '(I8,F12.3,F20.9,F20.9,F20.9,F20.9,F20.12)') itimes, &
                         t, ekin,temp,md_env % force_env % qs_env % energy % total, cons, cons-constant
     CALL m_flush(ene)
  END IF

  IF (ionode .AND. mod(itimes,icrd)==0) THEN
     CALL write_particle_coordinates(particle_set,crd,"XMOL","POS",title)
     CALL m_flush(crd)
  END IF

  IF (ionode .AND. mod(itimes,ivel)==0) THEN
     CALL write_particle_coordinates(particle_set,vel,"XMOL","POS_VEL_FORCE",title)
     CALL m_flush(vel)
  END IF

  END IF

  CALL timestop ( zero, handle )

END SUBROUTINE md_energy_qs

! *****************************************************************************

  SUBROUTINE open_md_file(file_name,unit_number)

    CHARACTER(LEN=*), INTENT(IN) :: file_name
    INTEGER, INTENT(OUT)         :: unit_number

!   *** Local variables ***

    INTEGER :: istat
    LOGICAL :: file_exists

!   ---------------------------------------------------------------------------

    INQUIRE (FILE=TRIM(file_name),EXIST=file_exists)

    IF (file_exists) THEN
      CALL open_file(file_name=TRIM(file_name),&
                     file_action="READWRITE",&
                     file_form="FORMATTED",&
                     file_position="APPEND",&
                     file_status="OLD",&
                     unit_number=unit_number)
    ELSE
      CALL open_file(file_name=TRIM(file_name),&
                     file_action="WRITE",&
                     file_form="FORMATTED",&
                     file_position="REWIND",&
                     file_status="REPLACE",&
                     unit_number=unit_number)
    END IF

  END SUBROUTINE open_md_file

! *****************************************************************************

END MODULE md_qs_energies
