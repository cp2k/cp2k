!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/md_kg_energies [1.0] *
!!
!!   NAME
!!     md_kg_energies
!!
!!   FUNCTION
!!     prints all energy info per timestep to the screen or to
!!     user defined output files
!!
!!   AUTHOR
!!     gt, 20.10.2002
!!
!!   MODIFICATION HISTORY
!!   6Feb03  calcualtion of conserv quantity moved to md_conserved_quantities
!!
!!   SOURCE
!******************************************************************************

MODULE md_kg_energies
  
  USE cp2k_units,                        ONLY : unit_convert_type, &
                                                get_cp2k_units
  USE coefficient_types,                 ONLY : coeff_type
  USE extended_system_types,          ONLY : lnhc_parameters_type
  USE global_types,                      ONLY : global_environment_type
!  USE kg_environment_types,              ONLY : kg_thermodynamic_type
  USE kinds,                             ONLY : dbl
  USE mathconstants,                     ONLY : zero
  USE md,                                ONLY : simulation_parameters_type,&
                                                mdio_parameters_type
  USE md_conserved_quantities,           ONLY : compute_conserved_quantity
  USE md_environment_types,              ONLY : md_environment_type, &
                                                get_md_env
  USE message_passing,                   ONLY : mp_sum,&
                                                mp_comm_dup
  USE simulation_cell,                   ONLY : get_cell_param, &
                                                cell_type
  USE termination,                       ONLY : stop_memory, &
                                                stop_program
  USE timings,                           ONLY : timeset,&
                                                timestop
  PRIVATE
  PUBLIC :: set_kg_energies, md_energy_kg

  TYPE average_quantities_type
     REAL ( dbl ) :: avetemp, avepot, avekin, avecore, avexc, avehartree, avebc
     REAL ( dbl ) :: avepress, avekinc, avetempc, aveppl, aveppl2
  END TYPE average_quantities_type

  REAL ( dbl ) :: econs_avg
  REAL ( dbl ) :: ecore_avg, ehartree_avg, exc_avg
  TYPE (average_quantities_type ) :: averages

  LOGICAL :: ionode
  INTEGER :: int_group
  CHARACTER (LEN=50) :: crd_file_name, vel_file_name, ptens_file_name, &
                        ener_file_name, temp_file_name, dump_file_name, &
                        pot_file_name, coef_file_name
  INTEGER :: crd, vel, pot, ene, tem, scr, coe, pol
  INTEGER :: icrd, ivel, ipot, iener, itemp, idump, icoef, iscreen

  TYPE ( global_environment_type ) :: glob_env
  TYPE ( unit_convert_type ), POINTER :: unit

!!******************************************************************************

CONTAINS

!!******************************************************************************
!!!****** md_kg_energies/set_kg_energies [1.0] *
!!!
!!!   NAME
!!!     set_kg_energies
!!!
!!!   FUNCTION
!!!     setup of parameters for kg_energies
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
!!!   SOURCE
!!******************************************************************************
!
SUBROUTINE set_kg_energies ( globenv, mdio )

  IMPLICIT NONE

! Arguments
  TYPE ( mdio_parameters_type ), INTENT ( IN ) :: mdio
  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv
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
!  pot_file_name = mdio % pot_file_name
!  ener_file_name = mdio % ener_file_name
!  temp_file_name = mdio % temp_file_name
!  dump_file_name = mdio % dump_file_name
!  coef_file_name = mdio % coef_file_name
!
!  icrd = mdio % icrd
!  icoef = mdio % icoef
!  ivel = mdio % ivel
!  ipot = mdio % ipot
!  iener = mdio % iener
!  itemp = mdio % itemp
!  idump = mdio % idump
!  iscreen = mdio % iscreen
!
!  CALL get_cp2k_units ( unit )
!
END SUBROUTINE set_kg_energies
!
!******************************************************************************
!!****** md_kg_energies/md_energy_kg [1.0] *
!!
!!   NAME
!!     energy_kg
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
!!   SOURCE
!******************************************************************************
!
!
SUBROUTINE md_energy_kg ( md_env )
!
  IMPLICIT NONE

! Arguments
  TYPE ( md_environment_type ), INTENT ( INOUT ) :: md_env
!  TYPE (kg_thermodynamic_type ), INTENT ( INOUT ) :: thermo
!
! Locals
!  
!  INTEGER, POINTER :: itimes
!  REAL (dbl), POINTER :: constant
!  TYPE ( cell_type ), POINTER :: box
!  TYPE ( simulation_parameters_type ), POINTER :: simpar
!  TYPE ( structure_type ), POINTER :: struc 
!
!  INTEGER i, j, id, n, l, nnodes, ncoef, iw, handle, inuse
!  REAL ( dbl ) :: ekin_coefs, temp_coefs
!  REAL ( dbl ) :: temp, econs, cons
!  REAL ( dbl ), DIMENSION(3) :: cell, ang
!
!!------------------------------------------------------------------------------
!
!  CALL timeset ( 'ENERGY', 'E', ' ', handle )
!  CALL get_md_env ( md_env=md_env, itimes=itimes, constant=constant, &
!                            cell=box, simpar=simpar )                  
!
!  struc => md_env % struc ( 1 )
!
!  nnodes = SIZE ( struc % pnode )
! 
!  IF (ASSOCIATED (struc % coef_vel ) ) THEN 
!    inuse = struc % coef_pos ( 1 ) % in_use
!    SELECT CASE ( inuse )
!    CASE ( 402 )
!      ncoef = SIZE ( struc % coef_vel ( 1 ) % ao % cr )
!    CASE ( 401 )
!      ncoef = SIZE ( struc % coef_vel ( 1 ) % pw % cc )
!    END SELECT
!   END IF
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
!     pot = get_unit()
!     OPEN ( UNIT = pot, FILE = pot_file_name )
!     IF (ASSOCIATED (struc % coef_vel ) ) THEN 
!       coe = get_unit()
!       OPEN ( UNIT = coe, FILE = coef_file_name )
!     END IF
!  END IF
!
!  iw = scr
!
!  IF (ASSOCIATED (struc%coef_vel)) THEN
!    CALL compute_conserved_quantity ( md_env, thermo%pot, cons, &
!                                      cons_rel=econs, &
!                                      ekin=thermo%kin, &
!                                      temp=temp,&
!                                      ekin_coefs=ekin_coefs, &
!                                      temp_coefs=temp_coefs)
!  ELSE
!    CALL compute_conserved_quantity ( md_env, thermo%pot, cons, &
!                                      cons_rel=econs , &
!                                      ekin=thermo%kin, temp=temp )                          
!  END IF
!
!  temp=temp / unit % kelvin
!
!  IF (ionode .AND. mod(itimes,itemp)==0) THEN
!     WRITE (tem, *) itimes, temp
!  END IF
!
!  IF (ionode .AND. mod(itimes,icrd)==0) THEN
!     DO i=1, nnodes
!        WRITE (crd,'(i6,3(2x,f12.5) )' ) itimes, struc%pnode(i) %p%r
!     END DO
!  END IF
!
!  IF (ionode .AND. mod(itimes,ivel)==0) THEN
!     DO i=1, nnodes
!        WRITE (vel,'(i6,3f10.5 )' ) itimes, struc%pnode(i) %p%v
!     END DO
!  END IF
!
!  IF ( ASSOCIATED ( struc % coef_vel  ) ) THEN
!    temp_coefs = temp_coefs / unit % kelvin
!    IF (ionode .AND. mod(itimes,icoef)==0) THEN
!      SELECT CASE ( struc % coef_pos ( 1 ) % in_use )
!      CASE ( 402 )
!         DO i=1, ncoef
!            WRITE (coe,'(i6,3f20.12 )' ) itimes, struc%coef_pos(1) %ao%cr (i), &
!                   struc%coef_vel(1) %ao%cr (i),  struc%coef_force(1) %ao%cr (i)
!         END DO
!      CASE ( 401 )
!         DO i=1, ncoef
!            WRITE (coe, *) itimes, struc%coef_pos(1)%pw%cc (i), &
!                   struc%coef_vel(1) %pw%cc (i),  struc%coef_force(1) %pw%cc (i)
!         END DO
!      END SELECT
!    END IF
!  END IF
!
!! cell parameters:
!  CALL get_cell_param ( box, cell, ang )
!
!  SELECT CASE (itimes)
!  CASE (0)
!     IF (ionode) THEN
!        WRITE (iw,*)
!        WRITE ( iw, '( A )' ) ' ENERGY| Initialization proceeding'
!        WRITE (iw,*)
!     END IF
!  CASE (1)
!     averages%avetemp = temp
!     averages%avecore = thermo%e_rspace_core
!     averages%avehartree = thermo%e_hartree
!     averages%avebc = thermo%e_bond_corr
!     averages%avexc = thermo%e_xc
!     averages%aveppl = thermo%e_ppl
!     averages%aveppl2 = thermo%e_ppl2
!     averages%avepot = thermo%pot
!     averages%avekin = thermo%kin
!     IF  ( ASSOCIATED ( struc % coef_vel ) ) THEN 
!       averages%avekinc = ekin_coefs   
!       averages%avetempc = temp_coefs
!     END IF
!     econs_avg =  0._dbl
!  CASE DEFAULT
!     averages%avetemp = (averages%avetemp*real(itimes-1)+temp)/ &
!          real(itimes)
!     averages%avecore = (averages%avecore*real(itimes-1)+thermo%e_rspace_core)/ &
!          real(itimes)
!     averages%avehartree = (averages%avehartree*real(itimes-1)+thermo%e_hartree)/ &
!          real(itimes)
!     averages%avebc = (averages%avebc*real(itimes-1)+thermo%e_bond_corr)/ &
!          real(itimes)
!     averages%avexc = (averages%avexc*real(itimes-1)+thermo%e_xc)/ &
!          real(itimes)
!     averages%aveppl = (averages%aveppl*real(itimes-1)+thermo%e_ppl)/ &
!          real(itimes)
!     averages%avepot = (averages%avepot*real(itimes-1)+thermo%pot)/ &
!          real(itimes)
!     averages%avekin = (averages%avekin*real(itimes-1)+thermo%kin)/ &
!          real(itimes)
!     IF  ( ASSOCIATED ( struc % coef_vel ) ) THEN 
!       averages%avekinc = (averages%avekinc*real(itimes-1)+ ekin_coefs)/ &
!            real(itimes)
!       averages%avetempc = (averages%avetempc*real(itimes-1)+temp_coefs)/ &
!            real(itimes)
!       averages%aveppl2 = (averages%aveppl2*real(itimes-1)+thermo%e_ppl2)/ &
!            real(itimes)
!     END IF
!     econs_avg = econs_avg + abs((cons-constant)/constant)
!  END SELECT
!
!  IF (ionode .AND. (mod(itimes,iener)==0) .AND. (itimes /= 0 )) THEN
!     WRITE (ene, *) itimes, econs, econs_avg/real(itimes)
!  END IF
!
!  IF (ionode .AND. (mod(itimes,ipot)==0) .AND. (itimes /= 0 )) THEN
!    WRITE (pot, *) itimes, thermo%pot, thermo%e_xc, thermo%e_hartree, &
!    thermo%e_rspace_core, thermo%e_ppl, thermo%e_bond_corr
!  END IF
!
!
!  IF (itimes==0 .AND. ionode) THEN
!     WRITE ( iw, '( )' )
!     WRITE ( iw, '( A,A )' ) ' ******************************** ', &
!          'GO CP2K GO! **********************************'
!     WRITE ( iw, '( A,A,T40,A,T65,1(1X,E15.7) )' ) ' INITIAL TOTAL POT ENERGY', &
!          adjustl(unit % e_label), '= ', thermo%pot
!     WRITE ( iw, '( A,A,T40,A,T65,1(1X,E15.7) )' ) ' INITIAL BOND CORRECT ENERGY', &
!          adjustl(unit % e_label), '= ', thermo%e_bond_corr
!     WRITE ( iw, '( A,A,T40,A,T65,1(1X,E15.7) )' ) ' INITIAL K+XC ENERGY', &
!          adjustl(unit % e_label), '= ', thermo%e_xc
!     WRITE ( iw, '( A,A,T40,A,T65,1(1X,E15.7) )' ) ' INITIAL HARTREE ENERGY', &
!          adjustl(unit % e_label), '= ', thermo%e_hartree
!     WRITE ( iw, '( A,A,T40,A,T65,1(1X,E15.7) )' ) ' INITIAL CORE OVLP ENERGY', &
!          adjustl(unit % e_label), '= ', thermo%e_rspace_core
!     WRITE ( iw, '( A,A,T40,A,T65,1(1X,E15.7) )' ) ' INITIAL PPL ENERGY', &
!           adjustl(unit % e_label), '= ', thermo%e_ppl
!     WRITE ( iw, '( A,A,T40,A,T65,1(1X,E15.7) )' ) ' SELF-ENERGY', &
!          adjustl(unit % e_label), '= ', thermo%e_self
!     WRITE ( iw, '( A,A,T40,A,T65,1(1X,E15.7) )' ) ' INITIAL KINETIC ENERGY', &
!          adjustl(unit % e_label), '= ', thermo%kin
!     WRITE ( iw, '( A,A,T40,A,T65,1(1X,E15.7) )' ) ' INITIAL TEMPERATURE', &
!          adjustl(unit % temp_label), '= ', temp
!     IF  ( ASSOCIATED ( struc % coef_vel  ) ) THEN 
!       WRITE ( iw, '( A,A,T40,A,T65,1(1X,E15.7) )' ) ' INITIAL COEF KINETIC ENERGY', &
!            adjustl(unit % e_label), '= ', ekin_coefs
!       WRITE ( iw, '( A,A,T40,A,T65,1(1X,E15.7) )' ) ' INITIAL COEF TEMPERATURE', &
!            adjustl(unit % temp_label), '= ', temp_coefs
!       WRITE ( iw, '( A,A,T40,A,T65,1(1X,E15.7) )' ) ' INITIAL POL PPL ENERGY  ', &
!            adjustl(unit % temp_label), '= ', thermo%e_ppl2
!     END IF 
!     WRITE ( iw, '( A,A,T40,A,T65,1(1X,E15.7) )' ) ' INITIAL VOLUME', &
!          adjustl(unit % vol_label), '= ', box%deth
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
!     WRITE ( iw, '( A,T64,1(1X,E16.8) )' ) ' CONSERVED QNTY               = ', &
!          cons
!     WRITE ( iw, '( )' )
!     WRITE ( iw, '( T52,A )' ) 'INSTANTANEOUS        AVERAGES'
!     WRITE ( iw, '( A,T49,2(1X,E15.7) )' ) ' {E-E0}/E0                    = ', &
!          econs, econs_avg/real(itimes)
!     WRITE ( iw, '( A,A,T31,A,T49,2(1X,E15.7) )' ) ' TOTAL POT ENERGY', &
!          adjustl(unit % e_label), '= ', thermo%pot, averages%avepot
!     WRITE ( iw, '( A,A,T31,A,T49,2(1X,E15.7) )' ) ' BOND CORRECT ENERGY', &
!          adjustl(unit % e_label), '= ', thermo%e_bond_corr, averages%avebc
!     WRITE ( iw, '( A,A,T31,A,T49,2(1X,E15.7) )' ) ' K+XC ENERGY', &
!          adjustl(unit % e_label), '= ', thermo%e_xc, averages%avexc
!     WRITE ( iw, '( A,A,T31,A,T49,2(1X,E15.7) )' ) ' CORE OVLP ENERGY', &
!          adjustl(unit % e_label), '= ', thermo%e_rspace_core, averages%avecore
!     WRITE ( iw, '( A,A,T31,A,T49,2(1X,E15.7) )' ) ' HARTREE   ENERGY', &
!          adjustl(unit % e_label), '= ', thermo%e_hartree, averages%avehartree
!     WRITE ( iw, '( A,A,T31,A,T49,2(1X,E15.7) )' ) ' LOCAL PP  ENERGY', &
!          adjustl(unit % e_label), '= ', thermo%e_ppl, averages%aveppl
!     WRITE ( iw, '( A,A,T31,A,T49,2(1X,E15.7) )' ) ' KINETIC ENERGY', &
!          adjustl(unit % e_label), '= ', thermo%kin, averages%avekin
!     WRITE ( iw, '( A,A,T31,A,T49,2(1X,E15.7) )' ) ' TEMPERATURE', &
!          adjustl(unit % temp_label), '= ', temp, averages%avetemp
!    IF  ( ASSOCIATED ( struc % coef_vel ) ) THEN 
!       WRITE ( iw, '( A,A,T31,A,T49,2(1X,E15.7) )' ) ' COEF KINETIC ENERGY', &
!            adjustl(unit % e_label), '= ', ekin_coefs,  averages%avekinc
!       WRITE ( iw, '( A,A,T31,A,T49,2(1X,E15.7) )' ) ' COEF TEMPERATURE', &
!            adjustl(unit % temp_label), '= ', temp_coefs, averages%avetempc
!       WRITE ( iw, '( A,A,T31,A,T49,2(1X,E15.7) )' ) ' POL PPL ENERGY  ', &
!            adjustl(unit%e_label), '= ', thermo%e_ppl2, averages%aveppl2
!     END IF
!     WRITE ( iw, '( A,A )' ) ' **************************************', &
!          '*****************************************'
!     WRITE ( iw, '( )' )
!  END IF
!
!  CALL timestop ( zero, handle )
!
END SUBROUTINE md_energy_kg
!
!******************************************************************************
END MODULE md_kg_energies
!******************************************************************************
