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
!******************************************************************************

MODULE md_qs_energies

  USE global_types,                    ONLY: global_environment_type
  USE integrator,                      ONLY : pressure
  USE kinds,                           ONLY : dbl
  USE linklist_types,                  ONLY : linklist_internal_data_type,  &
                                              get_linklist_internal_data
  USE machine,                         ONLY : m_flush
  USE mathconstants,                   ONLY : zero
  USE md,                              ONLY : simulation_parameters_type, &
                                              mdio_parameters_type
  USE md_conserved_quantities,         ONLY : compute_conserved_quantity,&
                                              compute_conserved_quantity_pimd
  USE md_environment_types,            ONLY : md_environment_type,&
                                              get_md_env
  USE message_passing,                 ONLY : mp_sum, mp_comm_dup
  USE simulation_cell,                 ONLY : get_cell_param,&
                                              cell_type
  USE structure_types,                 ONLY : structure_type
  USE particle_types,                  ONLY : particle_type, &
                                              write_particle_coordinates
  USE termination,                     ONLY : stop_memory,&
                                              stop_program
  USE timings,                         ONLY : timeset,&
                                              timestop
  USE util,                            ONLY : get_unit
  USE physcon,                         ONLY : angstrom , a_bohr , e_mass , joule

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: set_qs_energies, set_qs_energies_param, md_energy_qs

  CHARACTER ( LEN = 10 ) :: pres_label, e_label, temp_label
  CHARACTER ( LEN = 10 ) :: l_label, vol_label, angl_label
  REAL ( dbl ) :: econs_avg
  REAL ( dbl ) :: p_conv, temp_conv

  LOGICAL :: ionode
  CHARACTER (LEN=50) :: crd_file_name, vel_file_name, ptens_file_name, &
                        ener_file_name, temp_file_name, dump_file_name
  INTEGER :: crd, vel, ptn, ene, tem, scr, coe, pol
  INTEGER :: icrd, ivel, iptens, iener, itemp, idump, iscreen


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
!!****** md_qs_energies/set_qs_energies_params [1.0] *
!!
!!   NAME
!!     set_qs_energies_param
!!
!!   FUNCTION
!!     provides a label for the units, and a conversion factor ... rather ugly
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** *************************************************************************

SUBROUTINE set_qs_energies_param ( pconv, tempconv, l, vol, e, pres, temp, angl )

  IMPLICIT NONE

! Arguments
  REAL ( dbl ), INTENT ( IN ), OPTIONAL :: pconv, tempconv
  CHARACTER ( LEN = 10 ), INTENT ( IN ), OPTIONAL :: l, vol, e
  CHARACTER ( LEN = 10 ), INTENT ( IN ), OPTIONAL :: pres, temp, angl

!------------------------------------------------------------------------------

  IF ( PRESENT ( pconv ) ) p_conv = pconv
  IF ( PRESENT ( tempconv ) ) temp_conv = tempconv
  IF ( PRESENT ( l ) ) l_label = l
  IF ( PRESENT ( vol ) ) vol_label = vol
  IF ( PRESENT ( e ) ) e_label = e
  IF ( PRESENT ( pres ) ) pres_label = pres
  IF ( PRESENT ( temp ) ) temp_label = temp
  IF ( PRESENT ( angl ) ) angl_label = angl

END SUBROUTINE set_qs_energies_param

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
!!*** *************************************************************************


SUBROUTINE md_energy_qs ( md_env )

  IMPLICIT NONE

! Arguments
  TYPE (md_environment_type), INTENT ( INOUT ) :: md_env

! Locals
  TYPE ( simulation_parameters_type ), POINTER :: simpar
  TYPE ( structure_type ), POINTER :: struc (:)
  INTEGER, POINTER :: itimes
  REAL (dbl), POINTER :: constant
  CHARACTER(LEN=200) :: title

  INTEGER i, j, id, n, l, nnodes, beads, iw, handle
  REAL ( dbl ) :: alpha, epsilon0, qsq, lup, aup
  REAL ( dbl ) :: cons, temp, rcons, ekin
  REAL ( dbl ), DIMENSION(3) :: cell, ang
  REAL ( dbl ) :: pv_scalar,t,au2fs
  TYPE(particle_type), DIMENSION(:), POINTER :: particle_set

!------------------------------------------------------------------------------

  CALL timeset ( 'energy', 'E', ' ', handle )

  CALL get_md_env ( md_env=md_env,&
                    itimes=itimes, &
                    constant=constant, &
                    simpar=simpar,&
                    struc=struc)

  particle_set => struc(1) % part

! GT : to be changed?
  au2fs = 1.0 / ( 1.0E-15_dbl / ( a_bohr * SQRT ( e_mass / joule ) ) )
  t = itimes * md_env % simpar % dt * au2fs

  IF ( ionode .AND. itimes == 0 ) THEN
     tem = get_unit()
     OPEN ( UNIT = tem, FILE = temp_file_name )
     ene = get_unit()
     OPEN ( UNIT = ene, FILE = ener_file_name )
     crd = get_unit()
     OPEN ( UNIT = crd, FILE = crd_file_name )
     vel = get_unit()
     OPEN ( UNIT = vel, FILE = vel_file_name )
  END IF

  iw = scr

  CALL compute_conserved_quantity (md_env, md_env % force_env % qs_env % energy % total, cons,&
                                   cons_rel=rcons, &
                                   ekin=ekin, &
                                   temp=temp)

! GT : to be changed?
  temp=temp/temp_conv
! GT : to be changed?

  write(title,'(A4,I8,A4,F12.3,X,A6,F20.9,X,A6,F9.3,X,A1)') " i= ",itimes," t= ",t,"fs E= ", &
                          md_env % force_env % qs_env % energy % total,"au T= ",temp,"K"

  IF (ionode .AND. mod(itimes,itemp)==0) THEN
     WRITE (tem, '(I8,F20.10,F20.10)') itimes,t, temp
     CALL m_flush(tem)
  END IF

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

  CALL timestop ( zero, handle )

END SUBROUTINE md_energy_qs

!******************************************************************************
END MODULE md_qs_energies
!******************************************************************************
