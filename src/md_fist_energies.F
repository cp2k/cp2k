!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/md_fist_energies [1.0] *
!!
!!   NAME
!!     md_fist_energies
!!
!!   FUNCTION
!!     prints all energy info per timestep to the screen or to
!!     user defined output files
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!   gt : conserved qantity is now calculated by md_conserved_quantity
!!
!!   SOURCE
!******************************************************************************

MODULE md_fist_energies
  
!  USE cp2k_units,                      ONLY : unit_convert_type, &
!                                              get_cp2k_units
!  USE fist_environment_types,          ONLY : fist_thermodynamic_type
!  USE global_types,                    ONLY: global_environment_type
!  USE integrator,                      ONLY : pressure
!  USE kinds,                           ONLY : dbl
!  USE linklist_types,                  ONLY : linklist_internal_data_type,  &
!                                              get_linklist_internal_data
!  USE mathconstants,                   ONLY : zero
!  USE md,                              ONLY : simulation_parameters_type, &
!                                              mdio_parameters_type
!  USE md_conserved_quantities,         ONLY : compute_conserved_quantity,&
!                                              compute_conserved_quantity_pimd
!  USE md_environment_types,            ONLY : md_environment_type,&
!                                              get_md_env
!  USE message_passing,                 ONLY : mp_sum, mp_comm_dup
!  USE simulation_cell,                 ONLY : get_cell_param,&
!                                              cell_type
!  USE structure_types,                 ONLY : structure_type
!  USE termination,                     ONLY : stop_memory,&
!                                              stop_program
!  USE timings,                         ONLY : timeset,&
!                                              timestop
!  USE util,                            ONLY : get_unit
!
!  IMPLICIT NONE
!
!  PRIVATE
!  PUBLIC :: set_fist_energies, md_energy_fist
!
!  TYPE average_quantities_type
!     REAL ( dbl ) :: avetemp, avepot, avekin, avevol, aveca, avecb, avecc
!     REAL ( dbl ) :: aveal, avebe, avega, avepress, avekinc, avetempc
!  END TYPE average_quantities_type
!
!  REAL ( dbl ) :: econs_avg
!  TYPE (average_quantities_type ) :: averages
!
!  LOGICAL :: ionode
!  INTEGER :: int_group
!  CHARACTER (LEN=50) :: crd_file_name, vel_file_name, ptens_file_name, &
!                        ener_file_name, temp_file_name, dump_file_name
!  INTEGER :: crd, vel, ptn, ene, tem, scr, coe, pol
!  INTEGER :: icrd, ivel, iptens, iener, itemp, idump, iscreen
!
!  TYPE ( global_environment_type ) :: glob_env
!  TYPE ( unit_convert_type ), POINTER :: unit
!
!!!***
!!******************************************************************************
!
!CONTAINS
!
!!******************************************************************************
!!!****** md_fist_energies/set_fist_energies [1.0] *
!!!
!!!   NAME
!!!     set_fist_energies
!!!
!!!   FUNCTION
!!!     setup of parameters for fist_energies
!!!
!!!   AUTHOR
!!!     JGH
!!!
!!!   MODIFICATION HISTORY
!!!     JGH (29-Mar-2001) : duplicate communicator, instead of copy
!!!
!!!   NOTES
!!!     Call this subroutine before the first call to energy or velocity_verlet
!!!     or if you want to change ionode and/or output files
!!!
!!!*** *************************************************************************
!
!SUBROUTINE set_fist_energies ( globenv, mdio )
!
!  IMPLICIT NONE
!
!! Arguments
!  TYPE ( mdio_parameters_type ), INTENT ( IN ) :: mdio
!  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv
!
!!------------------------------------------------------------------------------
!
!  glob_env = globenv
!  CALL mp_comm_dup ( globenv % group, glob_env % group )
!
!  ionode = glob_env % ionode
!  int_group = glob_env % group
!  scr = glob_env % scr
!
!  crd_file_name = mdio % crd_file_name
!  vel_file_name = mdio % vel_file_name
!  ptens_file_name = mdio % ptens_file_name
!  ener_file_name = mdio % ener_file_name
!  temp_file_name = mdio % temp_file_name
!  dump_file_name = mdio % dump_file_name
!
!  icrd = mdio % icrd
!  ivel = mdio % ivel
!  iptens = mdio % iptens
!  iener = mdio % iener
!  itemp = mdio % itemp
!  idump = mdio % idump
!  iscreen = mdio % iscreen
!
!  CALL get_cp2k_units ( unit )
!
!END SUBROUTINE set_fist_energies
!
!!******************************************************************************
!!!****** md_fist_energies/md_energy_fist [1.0] *
!!!
!!!   NAME
!!!     energy_fist
!!!
!!!   FUNCTION
!!!     This routine computes the conserved quantity, temperature
!!!     and things like that and prints them out
!!!
!!!   AUTHOR
!!!     CJM
!!!
!!!   MODIFICATION HISTORY
!!!     none
!!!
!!!*** *************************************************************************
!
!
!SUBROUTINE md_energy_fist ( md_env, box, ll_data, thermo, pimd )
!
!  IMPLICIT NONE
!
!! Arguments
!  TYPE (md_environment_type), INTENT ( INOUT ) :: md_env
!  TYPE (cell_type), INTENT ( IN ):: box
!  TYPE (linklist_internal_data_type), INTENT (IN) :: ll_data
!  TYPE (fist_thermodynamic_type ), INTENT ( INOUT ) :: thermo
!  LOGICAL, INTENT (in) :: pimd
!
!! Locals
!  TYPE ( simulation_parameters_type ), POINTER :: simpar
!  TYPE ( structure_type ), POINTER :: struc (:)
!  INTEGER, POINTER :: itimes
!  REAL (dbl), POINTER :: constant
!
!  INTEGER i, j, id, n, l, nnodes, beads, iw, handle
!  REAL ( dbl ) :: alpha, epsilon0, qsq, lup, aup
!  REAL ( dbl ) :: cons, temp, econs
!  REAL ( dbl ), DIMENSION(3) :: cell, ang
!  REAL ( dbl ) :: pv_scalar
!
!!------------------------------------------------------------------------------
!
!  CALL timeset ( 'ENERGY', 'E', ' ', handle )
!
!  CALL get_md_env ( md_env=md_env,&
!                    itimes=itimes, &
!                    constant=constant, &
!                    simpar=simpar,&
!                    struc=struc )
!
!  nnodes = SIZE ( struc ( 1 ) % pnode )
!
!  beads = SIZE ( struc )
!
!  IF ( ionode .AND. itimes == 0 ) THEN
!     tem = get_unit()
!     OPEN ( UNIT = tem, FILE = temp_file_name )
!     ene = get_unit()
!     OPEN ( UNIT = ene, FILE = ener_file_name )
!     crd = get_unit()
!     OPEN ( UNIT = crd, FILE = crd_file_name )
!     vel = get_unit()
!     OPEN ( UNIT = vel, FILE = vel_file_name )
!     ptn = get_unit()
!     OPEN ( UNIT = ptn, FILE = ptens_file_name )
!  END IF
!
!  iw = scr
!
!  IF (pimd) THEN
!     CALL compute_conserved_quantity_pimd (md_env, thermo%pot, cons,&
!                                   econs, thermo%kin, temp)
!  ELSE
!     CALL compute_conserved_quantity (md_env, thermo%pot, cons,&
!                                      cons_rel=econs, &
!                                      ekin=thermo%kin, &
!                                      temp=temp)
!  END IF
!
!  temp = temp / unit % kelvin
!
!  IF (ionode .AND. mod(itimes,itemp)==0) THEN
!     WRITE (tem, *) itimes, temp
!  END IF
!
!  IF (ionode .AND. mod(itimes,iptens)==0) THEN
!     DO i=1,3
!        WRITE (ptn,*) itimes, (thermo%virial%pv_total(i,j), j=1,3)
!     END DO
!  END IF
!
!  IF (ionode .AND. mod(itimes,icrd)==0) THEN
!     DO j=1, beads
!        DO i=1, nnodes
!           WRITE (crd,'(i6,3(2x,f12.5) )' ) itimes, struc ( j )%pnode(i) %p%r
!        END DO
!     END DO
!  END IF
!
!  IF (ionode .AND. mod(itimes,ivel)==0) THEN
!     DO j=1, beads
!        DO i=1, nnodes
!           WRITE (vel,'(i6,3f10.5 )' ) itimes, struc ( j )%pnode(i) %p%v
!        END DO
!     END DO
!  END IF
!
!! cell parameters:
!  CALL get_cell_param ( box, cell, ang )
!
!! initializing pv_const=0
!  IF (itimes==0) THEN
!     thermo%virial%pv_constraint = 0._dbl
!     CALL pressure(struc (1) %pnode,thermo%virial)
!  END IF
!
!! computing the scalar pressure
!  pv_scalar = 0._dbl
!  DO i = 1, 3
!     pv_scalar = pv_scalar + thermo%virial%pv_total(i,i)
!  END DO
!  pv_scalar = pv_scalar/3._dbl/box%deth
!  pv_scalar = pv_scalar * unit % bar
!
!! average quantities:
!  SELECT CASE (itimes)
!  CASE (0)
!     IF (ionode) THEN
!        WRITE (iw,*)
!        WRITE ( iw, '( A )' ) ' ENERGY| Initialization proceeding'
!        WRITE (iw,*)
!     END IF
!  CASE (1)
!     averages%avepress = pv_scalar
!     averages%avetemp = temp
!     averages%avepot = thermo%pot
!     averages%avekin = thermo%kin
!     averages%avevol = box % deth
!     averages%aveca = cell(1)
!     averages%avecb = cell(2)
!     averages%avecc = cell(3)
!     averages%aveal = ang(3)
!     averages%avebe = ang(2)
!     averages%avega = ang(1)
!     econs_avg =  0._dbl
!  CASE DEFAULT
!     averages%avepress = (averages%avepress*real(itimes-1)+pv_scalar)/ &
!          REAL ( itimes, dbl )
!     averages%avetemp = (averages%avetemp*real(itimes-1)+temp)/ &
!          real(itimes)
!     averages%avepot = (averages%avepot*real(itimes-1)+thermo%pot)/ &
!          real(itimes)
!     averages%avekin = (averages%avekin*real(itimes-1)+thermo%kin)/ &
!          real(itimes)
!     averages%avevol = (averages%avevol*real(itimes-1)+box%deth)/ &
!          real(itimes)
!     averages%aveca = (averages%aveca*real(itimes-1)+cell(1))/real(itimes)
!     averages%avecb = (averages%avecb*real(itimes-1)+cell(2))/real(itimes)
!     averages%avecc = (averages%avecc*real(itimes-1)+cell(3))/real(itimes)
!     averages%aveal = (averages%aveal*real(itimes-1)+ang(3))/ &
!          real(itimes)
!     averages%avebe = (averages%avebe*real(itimes-1)+ang(2))/ &
!          real(itimes)
!     averages%avega = (averages%avega*real(itimes-1)+ang(1))/ &
!          real(itimes)
!     econs_avg = econs_avg + abs((cons-constant)/constant)
!  END SELECT
!
!  IF (ionode .AND. (mod(itimes,iener)==0) .AND. (itimes /= 0 )) THEN
!     WRITE (ene, *) itimes, econs, econs_avg/real(itimes)
!  END IF
!
!  IF (itimes==0 .AND. ionode) THEN
!     WRITE ( iw, '( )' )
!     WRITE ( iw, '( A,A )' ) ' ******************************** ', &
!          'GO CP2K GO! **********************************'
!     WRITE ( iw, '( A,A,T40,A,T65,1(1X,E15.7) )' ) ' INITIAL POTENTIAL ENERGY', &
!          adjustl( unit % e_label), '= ', thermo%pot
!     WRITE ( iw, '( A,A,T40,A,T65,1(1X,E15.7) )' ) ' INITIAL KINETIC ENERGY', &
!          adjustl( unit % e_label), '= ', thermo%kin
!     WRITE ( iw, '( A,A,T40,A,T65,1(1X,E15.7) )' ) ' INITIAL TEMPERATURE', &
!          adjustl( unit % temp_label ), '= ', temp
!     WRITE ( iw, '( A,A,T40,A,T65,1(1X,E15.7) )' ) ' INITIAL PRESSURE', &
!          adjustl( unit % pres_label), '= ', pv_scalar
!     WRITE ( iw, '( A,A,T40,A,T65,1(1X,E15.7) )' ) ' INITIAL VOLUME', &
!          adjustl( unit % vol_label), '= ', box%deth
!     WRITE ( iw, '( A,A,T29,A,T33,3(1X,E15.7) )' ) ' INITIAL CELL LNTHS', &
!          adjustl(unit % l_label), '= ', cell(1), cell(2), cell(3)
!     WRITE ( iw, '( A,A,T29,A,T33,3(1X,E15.7) )' ) ' INITIAL CELL ANGLS', &
!          adjustl(unit % angl_label), '= ', ang(3), ang(2), ang(1)
!     WRITE ( iw, '( A,A )' ) ' ******************************** ', &
!          'GO CP2K GO! **********************************'
!  END IF
!
!  IF ( ( itimes /= 0 ) .AND. ( MOD ( itimes, iscreen ) == 0 ) &
!       .AND. ionode ) THEN
!
!     WRITE ( iw, '( )' )
!     WRITE ( iw, '( A,A )' ) ' **************************************', &
!          '*****************************************'
!     WRITE ( iw, '( A,T61,A )' ) ' ENSEMBLE TYPE                = ', &
!          ADJUSTR ( simpar%ensemble)
!     WRITE ( iw, '( A,T71,I10 )' ) ' STEP NUMBER                  = ', itimes
!     WRITE ( iw, '( A,T65,1(1X,E15.7) )' ) ' CONSERVED QNTY               = ', &
!          cons
!     WRITE ( iw, '( )' )
!     WRITE ( iw, '( T52,A )' ) 'INSTANTANEOUS        AVERAGES'
!     WRITE ( iw, '( A,T49,2(1X,E15.7) )' ) ' {E-E0}/E0                    = ', &
!          econs, econs_avg/real(itimes)
!     WRITE ( iw, '( A,A,T31,A,T49,2(1X,E15.7) )' ) ' POTENTIAL ENERGY', &
!          adjustl(unit % e_label), '= ', thermo%pot, averages%avepot
!     WRITE ( iw, '( A,A,T31,A,T49,2(1X,E15.7) )' ) ' KINETIC ENERGY', &
!          adjustl(unit % e_label), '= ', thermo%kin, averages%avekin
!     WRITE ( iw, '( A,A,T31,A,T49,2(1X,E15.7) )' ) ' TEMPERATURE', &
!          adjustl(unit % temp_label), '= ', temp, averages%avetemp
!     WRITE ( iw, '( A,A,T31,A,T49,2(1X,E15.7) )' ) ' PRESSURE', &
!          adjustl(unit % pres_label), '= ', pv_scalar, averages%avepress
!     CALL get_linklist_internal_data ( ll_data, lup, aup )
!     IF (simpar%ensemble(1:3)=='NPT') THEN
!        WRITE ( iw, '( A,A,T31,A,T49,2(1X,E15.7) )' ) ' VOLUME', &
!             adjustl(unit % vol_label), '= ', box%deth, averages%avevol
!        WRITE ( iw, '( A,A,T31,A,T33,3(1X,E15.7) )' ) ' CELL LNTHS', &
!             adjustl(unit % l_label), '= ', cell(1), cell(2), cell(3)
!        WRITE ( iw, '( A,A,T31,A,T33,3(1X,E15.7) )' ) ' AVE. CELL LNTHS', &
!             adjustl(unit % l_label), '= ', averages%aveca, averages%avecb, &
!             averages%avecc
!     END IF
!     IF (simpar%ensemble=='NPT_F') THEN
!        WRITE ( iw, '( A,A,T31,A,T33,3(1X,E15.7) )' ) ' CELL ANGLS', &
!             adjustl(unit % angl_label), '= ', ang(3), ang(2), ang(1)
!        WRITE ( iw, '( A,A,T31,A,T33,3(1X,E15.7) )' ) ' AVE. CELL ANGLS', &
!             adjustl(unit % angl_label), '= ', averages%aveal, averages%avebe, &
!             averages%avega
!     END IF
!     WRITE ( iw, '( A,A,T31,A,T49,2(1X,F15.2) )' ) &
!          ' LIST UPDATES', '[steps]', '= ', lup, aup
!     WRITE ( iw, '( A,A )' ) ' **************************************', &
!          '*****************************************'
!     WRITE ( iw, '( )' )
!  END IF
!
!  CALL timestop ( zero, handle )
!
!END SUBROUTINE md_energy_fist
!
!******************************************************************************
END MODULE md_fist_energies
!******************************************************************************
