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
!!
!!   SOURCE
!******************************************************************************

MODULE md_fist_energies
  
!  USE coefficient_types, ONLY : coeff_type
  USE fist_environment_types, ONLY : fist_thermodynamic_type
  USE global_types, ONLY: global_environment_type
  USE integrator, ONLY : pressure
  USE kinds, ONLY : dbl
  USE linklist_control, ONLY : get_ll_parm
  USE linklist_types, ONLY : linklist_internal_data_type
  USE mathconstants, ONLY : zero
  USE md, ONLY : simulation_parameters_type, mdio_parameters_type
  USE message_passing, ONLY : mp_sum, mp_comm_dup
  USE nose, ONLY : lnhc_parameters_type, npt_info_type
  USE simulation_cell, ONLY : get_cell_param, cell_type
  USE structure_types, ONLY : structure_type
  USE termination, ONLY : stop_memory, stop_program
  USE timings, ONLY : timeset, timestop
  USE util, ONLY : get_unit

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: set_fist_energies, set_fist_energies_param, md_energy_fist, &
            md_energy_fist_pimd

  TYPE average_quantities_type
     REAL ( dbl ) :: avetemp, avepot, avekin, avevol, aveca, avecb, avecc
     REAL ( dbl ) :: aveal, avebe, avega, avepress, avekinc, avetempc
  END TYPE average_quantities_type

  CHARACTER ( LEN = 10 ) :: pres_label, e_label, temp_label
  CHARACTER ( LEN = 10 ) :: l_label, vol_label, angl_label
  REAL ( dbl ) :: econs_avg
  REAL ( dbl ) :: p_conv, temp_conv
  TYPE (average_quantities_type ) :: averages

  LOGICAL :: ionode
  INTEGER :: int_group
  CHARACTER (LEN=50) :: crd_file_name, vel_file_name, ptens_file_name, &
                        ener_file_name, temp_file_name, dump_file_name
  INTEGER :: crd, vel, ptn, ene, tem, scr, coe, pol
  INTEGER :: icrd, ivel, iptens, iener, itemp, idump, iscreen

  TYPE ( global_environment_type ) :: glob_env

!******************************************************************************

CONTAINS

!******************************************************************************
!!****** md_fist_energies/set_fist_energies [1.0] *
!!
!!   NAME
!!     set_fist_energies
!!
!!   FUNCTION
!!     setup of parameters for fist_energies
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     JGH (29-Mar-2001) : duplicate communicator, instead of copy
!!
!!   NOTES
!!     Call this subroutine before the first call to energy or velocity_verlet
!!     or if you want to change ionode and/or output files
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE set_fist_energies ( globenv, mdio )

  IMPLICIT NONE

! Arguments
  TYPE ( mdio_parameters_type ), INTENT ( IN ) :: mdio
  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv

!------------------------------------------------------------------------------

  glob_env = globenv
  CALL mp_comm_dup ( globenv % group, glob_env % group )

  ionode = glob_env % ionode
  int_group = glob_env % group
  scr = glob_env % scr

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

END SUBROUTINE set_fist_energies

!******************************************************************************
!!****** md_fist_energies/set_fist_energies_params [1.0] *
!!
!!   NAME
!!     set_fist_energies_param
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE set_fist_energies_param ( pconv, tempconv, l, vol, e, pres, temp, angl )

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

END SUBROUTINE set_fist_energies_param

!******************************************************************************
!!****** md_fist_energies/md_energy_fist [1.0] *
!!
!!   NAME
!!     energy_fist
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


SUBROUTINE md_energy_fist ( itimes, constant, box, ll_data, nhc_part, nhc_baro, &
                            npt, simpar, struc, thermo )
  IMPLICIT NONE

! Arguments
  TYPE (cell_type), INTENT ( IN ):: box
  TYPE (linklist_internal_data_type), INTENT (IN) :: ll_data
  TYPE (fist_thermodynamic_type ), INTENT ( INOUT ) :: thermo
  TYPE (simulation_parameters_type ), INTENT ( IN )  :: simpar
  TYPE (structure_type ), INTENT ( IN ) :: struc
  TYPE (lnhc_parameters_type ), INTENT ( IN ) :: nhc_part, nhc_baro
  TYPE (npt_info_type ), INTENT ( IN ) :: npt ( :, : )
  INTEGER, INTENT ( IN ) :: itimes
  REAL ( dbl ), INTENT ( INOUT ):: constant

! Locals
  INTEGER i, j, id, n, l, nnodes, iw, handle
  REAL ( dbl ) :: nhc_part_pot, nhc_part_kin, baro_kin, baro_pot
  REAL ( dbl ) :: nhc_baro_pot, nhc_baro_kin
  REAL ( dbl ) :: alpha, epsilon0, qsq, lup, aup
  REAL ( dbl ) :: cons, temp, econs
  REAL ( dbl ), DIMENSION(3) :: cell, ang
  REAL ( dbl ) :: pv_scalar

!------------------------------------------------------------------------------

  CALL timeset ( 'ENERGY', 'E', ' ', handle )
 
  nnodes = SIZE ( struc % pnode )

  IF ( ionode .AND. itimes == 0 ) THEN
     tem = get_unit()
     OPEN ( UNIT = tem, FILE = temp_file_name )
     ene = get_unit()
     OPEN ( UNIT = ene, FILE = ener_file_name )
     crd = get_unit()
     OPEN ( UNIT = crd, FILE = crd_file_name )
     vel = get_unit()
     OPEN ( UNIT = vel, FILE = vel_file_name )
     ptn = get_unit()
     OPEN ( UNIT = ptn, FILE = ptens_file_name )
  END IF

  iw = scr

  CALL get_part_ke ( struc, thermo )
  
  CALL get_part_temp ( thermo, simpar, temp )


  IF (ionode .AND. mod(itimes,itemp)==0) THEN
     WRITE (tem, *) itimes, temp
  END IF

  IF (ionode .AND. mod(itimes,iptens)==0) THEN
     DO i=1,3
        WRITE (ptn,*) itimes, (thermo%virial%pv_total(i,j), j=1,3)
     END DO
  END IF

  IF (ionode .AND. mod(itimes,icrd)==0) THEN
     DO i=1, nnodes
        WRITE (crd,'(i6,3(2x,f12.5) )' ) itimes, struc%pnode(i) %p%r
     END DO
  END IF

  IF (ionode .AND. mod(itimes,ivel)==0) THEN
     DO i=1, nnodes
        WRITE (vel,'(i6,3f10.5 )' ) itimes, struc%pnode(i) %p%v
     END DO
  END IF
  IF ( simpar % ensemble ( 1 : 3 )/= 'NVE' ) &
  CALL get_nhc_energies ( nhc_part, nhc_part_kin, nhc_part_pot )

  IF ( simpar % ensemble ( 1 : 3 ) == 'NPT' ) &
  CALL get_nhc_energies (  nhc_baro, nhc_baro_kin, nhc_baro_pot )

  IF ( simpar % ensemble ( 1 : 3 ) == 'NPT' ) &
  CALL get_baro_energies ( box, simpar, nhc_baro, npt, baro_kin, baro_pot )



  SELECT CASE ( simpar % ensemble ( 1 : 3 ) )
  CASE ( 'NVE')
     cons = thermo%kin + thermo%pot
  CASE ( 'NVT')
     cons = thermo%kin + thermo%pot + nhc_part_kin + nhc_part_pot
  CASE ( 'NPT')
     cons = thermo%kin + thermo%pot + nhc_part_kin + nhc_part_pot  &
            + nhc_baro_kin + nhc_baro_pot + baro_kin + baro_pot
  END SELECT

  IF (itimes==1) THEN
     constant = cons
  END IF
  IF (itimes /= 0 ) econs = abs((cons-constant)/constant)

! cell parameters:
  CALL get_cell_param ( box, cell, ang )

! initializing pv_const=0
  IF (itimes==0) THEN
     thermo%virial%pv_constraint = 0._dbl
     CALL pressure(struc%pnode,thermo%virial)
  END IF

! computing the scalar pressure
  pv_scalar = 0._dbl
  DO i = 1, 3
     pv_scalar = pv_scalar + thermo%virial%pv_total(i,i)
  END DO
  pv_scalar = pv_scalar/3._dbl/box%deth
  pv_scalar = pv_scalar*p_conv

! average quantities:
  SELECT CASE (itimes)
  CASE (0)
     IF (ionode) THEN
        WRITE (iw,*)
        WRITE ( iw, '( A )' ) ' ENERGY| Initialization proceeding'
        WRITE (iw,*)
     END IF
  CASE (1)
     averages%avepress = pv_scalar
     averages%avetemp = temp
     averages%avepot = thermo%pot
     averages%avekin = thermo%kin
     averages%avevol = box % deth
     averages%aveca = cell(1)
     averages%avecb = cell(2)
     averages%avecc = cell(3)
     averages%aveal = ang(3)
     averages%avebe = ang(2)
     averages%avega = ang(1)
     econs_avg =  0._dbl
  CASE DEFAULT
     averages%avepress = (averages%avepress*real(itimes-1)+pv_scalar)/ &
          REAL ( itimes, dbl )
     averages%avetemp = (averages%avetemp*real(itimes-1)+temp)/ &
          real(itimes)
     averages%avepot = (averages%avepot*real(itimes-1)+thermo%pot)/ &
          real(itimes)
     averages%avekin = (averages%avekin*real(itimes-1)+thermo%kin)/ &
          real(itimes)
     averages%avevol = (averages%avevol*real(itimes-1)+box%deth)/ &
          real(itimes)
     averages%aveca = (averages%aveca*real(itimes-1)+cell(1))/real(itimes)
     averages%avecb = (averages%avecb*real(itimes-1)+cell(2))/real(itimes)
     averages%avecc = (averages%avecc*real(itimes-1)+cell(3))/real(itimes)
     averages%aveal = (averages%aveal*real(itimes-1)+ang(3))/ &
          real(itimes)
     averages%avebe = (averages%avebe*real(itimes-1)+ang(2))/ &
          real(itimes)
     averages%avega = (averages%avega*real(itimes-1)+ang(1))/ &
          real(itimes)
     econs_avg = econs_avg + abs((cons-constant)/constant)
  END SELECT

  IF (ionode .AND. (mod(itimes,iener)==0) .AND. (itimes /= 0 )) THEN
     WRITE (ene, *) itimes, econs, econs_avg/real(itimes)
  END IF

  IF (itimes==0 .AND. ionode) THEN
     WRITE ( iw, '( )' )
     WRITE ( iw, '( A,A )' ) ' ******************************** ', &
          'GO CP2K GO! **********************************'
     WRITE ( iw, '( A,A,T40,A,T65,1(1X,E15.7) )' ) ' INITIAL POTENTIAL ENERGY', &
          adjustl(e_label), '= ', thermo%pot
     WRITE ( iw, '( A,A,T40,A,T65,1(1X,E15.7) )' ) ' INITIAL KINETIC ENERGY', &
          adjustl(e_label), '= ', thermo%kin
     WRITE ( iw, '( A,A,T40,A,T65,1(1X,E15.7) )' ) ' INITIAL TEMPERATURE', &
          adjustl(temp_label), '= ', temp
     WRITE ( iw, '( A,A,T40,A,T65,1(1X,E15.7) )' ) ' INITIAL PRESSURE', &
          adjustl(pres_label), '= ', pv_scalar
     WRITE ( iw, '( A,A,T40,A,T65,1(1X,E15.7) )' ) ' INITIAL VOLUME', &
          adjustl(vol_label), '= ', box%deth
     WRITE ( iw, '( A,A,T29,A,T33,3(1X,E15.7) )' ) ' INITIAL CELL LNTHS', &
          adjustl(l_label), '= ', cell(1), cell(2), cell(3)
     WRITE ( iw, '( A,A,T29,A,T33,3(1X,E15.7) )' ) ' INITIAL CELL ANGLS', &
          adjustl(angl_label), '= ', ang(3), ang(2), ang(1)
     WRITE ( iw, '( A,A )' ) ' ******************************** ', &
          'GO CP2K GO! **********************************'
  END IF

  IF ( ( itimes /= 0 ) .AND. ( MOD ( itimes, iscreen ) == 0 ) &
       .AND. ionode ) THEN

     WRITE ( iw, '( )' )
     WRITE ( iw, '( A,A )' ) ' **************************************', &
          '*****************************************'
     WRITE ( iw, '( A,T61,A )' ) ' ENSEMBLE TYPE                = ', &
          ADJUSTR ( simpar%ensemble)
     WRITE ( iw, '( A,T71,I10 )' ) ' STEP NUMBER                  = ', itimes
     WRITE ( iw, '( A,T65,1(1X,E15.7) )' ) ' CONSERVED QNTY               = ', &
          cons
     WRITE ( iw, '( )' )
     WRITE ( iw, '( T52,A )' ) 'INSTANTANEOUS        AVERAGES'
     WRITE ( iw, '( A,T49,2(1X,E15.7) )' ) ' {E-E0}/E0                    = ', &
          econs, econs_avg/real(itimes)
     WRITE ( iw, '( A,A,T31,A,T49,2(1X,E15.7) )' ) ' POTENTIAL ENERGY', &
          adjustl(e_label), '= ', thermo%pot, averages%avepot
     WRITE ( iw, '( A,A,T31,A,T49,2(1X,E15.7) )' ) ' KINETIC ENERGY', &
          adjustl(e_label), '= ', thermo%kin, averages%avekin
     WRITE ( iw, '( A,A,T31,A,T49,2(1X,E15.7) )' ) ' TEMPERATURE', &
          adjustl(temp_label), '= ', temp, averages%avetemp
     WRITE ( iw, '( A,A,T31,A,T49,2(1X,E15.7) )' ) ' PRESSURE', &
          adjustl(pres_label), '= ', pv_scalar, averages%avepress
     CALL get_ll_parm ( ll_data, lup, aup )
     IF (simpar%ensemble(1:3)=='NPT') THEN
        WRITE ( iw, '( A,A,T31,A,T49,2(1X,E15.7) )' ) ' VOLUME', &
             adjustl(vol_label), '= ', box%deth, averages%avevol
        WRITE ( iw, '( A,A,T31,A,T33,3(1X,E15.7) )' ) ' CELL LNTHS', &
             adjustl(l_label), '= ', cell(1), cell(2), cell(3)
        WRITE ( iw, '( A,A,T31,A,T33,3(1X,E15.7) )' ) ' AVE. CELL LNTHS', &
             adjustl(l_label), '= ', averages%aveca, averages%avecb, &
             averages%avecc
     END IF
     IF (simpar%ensemble=='NPT_F') THEN
        WRITE ( iw, '( A,A,T31,A,T33,3(1X,E15.7) )' ) ' CELL ANGLS', &
             adjustl(angl_label), '= ', ang(3), ang(2), ang(1)
        WRITE ( iw, '( A,A,T31,A,T33,3(1X,E15.7) )' ) ' AVE. CELL ANGLS', &
             adjustl(angl_label), '= ', averages%aveal, averages%avebe, &
             averages%avega
     END IF
     WRITE ( iw, '( A,A,T31,A,T49,2(1X,F15.2) )' ) &
          ' LIST UPDATES', '[steps]', '= ', lup, aup
     WRITE ( iw, '( A,A )' ) ' **************************************', &
          '*****************************************'
     WRITE ( iw, '( )' )
  END IF

  CALL timestop ( zero, handle )

END SUBROUTINE md_energy_fist

!******************************************************************************
!!****** md_fist_energies/md_energy_fist_pimd [1.0] *
!!
!!   NAME
!!     energy_fist_pimd
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

SUBROUTINE md_energy_fist_pimd ( itimes, constant, box, ll_data, nhc_part, nhc_baro, &
                                 npt, simpar, struc, thermo )
  IMPLICIT NONE

! Arguments
  TYPE (cell_type), INTENT (IN) :: box
  TYPE (linklist_internal_data_type), INTENT (IN) :: ll_data
  TYPE (fist_thermodynamic_type ), INTENT ( INOUT ) :: thermo
  TYPE (simulation_parameters_type ), INTENT ( IN )  :: simpar
  TYPE (structure_type ), INTENT ( IN ) :: struc ( : )
  TYPE (lnhc_parameters_type ), INTENT ( IN ) :: nhc_part ( : ), nhc_baro ( : )
  TYPE (npt_info_type ), INTENT ( IN ) :: npt ( :, : )
  INTEGER, INTENT ( IN ) :: itimes
  REAL ( dbl ), INTENT ( INOUT ):: constant

! Locals
  INTEGER i, j, id, n, l, nnodes, iw, handle
  REAL ( dbl ) :: nhc_part_pot, nhc_part_kin, baro_kin, baro_pot
  REAL ( dbl ) :: nhc_baro_pot, nhc_baro_kin
  REAL ( dbl ) :: alpha, epsilon0, qsq, lup, aup
  REAL ( dbl ) :: cons, temp, econs
  REAL ( dbl ), DIMENSION(3) :: cell, ang
  REAL ( dbl ) :: pv_scalar

!------------------------------------------------------------------------------
  CALL timeset ( 'ENERGY', 'E', ' ', handle )
  
  nnodes = SIZE ( struc ( 1 ) % pnode )

  IF ( ionode .AND. itimes == 0 ) THEN
     tem = get_unit()
     OPEN ( UNIT = tem, FILE = temp_file_name )
     ene = get_unit()
     OPEN ( UNIT = ene, FILE = ener_file_name )
     crd = get_unit()
     OPEN ( UNIT = crd, FILE = crd_file_name )
     vel = get_unit()
     OPEN ( UNIT = vel, FILE = vel_file_name )
     ptn = get_unit()
     OPEN ( UNIT = ptn, FILE = ptens_file_name )
  END IF

  iw = scr

  CALL get_part_ke ( struc ( 1 ), thermo )
  
  CALL get_part_temp ( thermo, simpar, temp )


  IF (ionode .AND. mod(itimes,itemp)==0) THEN
     WRITE (tem, *) itimes, temp
  END IF

  IF (ionode .AND. mod(itimes,iptens)==0) THEN
     DO i=1,3
        WRITE (ptn,*) itimes, (thermo%virial%pv_total(i,j), j=1,3)
     END DO
  END IF

  IF (ionode .AND. mod(itimes,icrd)==0) THEN
     DO i=1, nnodes
        WRITE (crd,'(i6,3(2x,f12.5) )' ) itimes, struc ( 1 )%pnode(i) %p%r
     END DO
  END IF

  IF (ionode .AND. mod(itimes,ivel)==0) THEN
     DO i=1, nnodes
        WRITE (vel,'(i6,3f10.5 )' ) itimes, struc ( 1 )%pnode(i) %p%v
     END DO
  END IF
  IF ( simpar % ensemble ( 1 : 3 )/= 'NVE' ) &
  CALL get_nhc_energies ( nhc_part ( 1 ), nhc_part_kin, nhc_part_pot )

  IF ( simpar % ensemble ( 1 : 3 ) == 'NPT' ) &
  CALL get_nhc_energies (  nhc_baro ( 1 ), nhc_baro_kin, nhc_baro_pot )

  IF ( simpar % ensemble ( 1 : 3 ) == 'NPT' ) &
  CALL get_baro_energies ( box, simpar, nhc_baro ( 1 ), npt,  &
                           baro_kin, baro_pot )



  SELECT CASE ( simpar % ensemble ( 1 : 3 ) )
  CASE ( 'NVE')
     cons = thermo%kin + thermo%pot
  CASE ( 'NVT')
     cons = thermo%kin + thermo%pot + nhc_part_kin + nhc_part_pot
  CASE ( 'NPT')
     cons = thermo%kin + thermo%pot + nhc_part_kin + nhc_part_pot  &
            + nhc_baro_kin + nhc_baro_pot + baro_kin + baro_pot
  END SELECT

  IF (itimes==1) THEN
     constant = cons
  END IF
  IF (itimes /= 0 ) econs = abs((cons-constant)/constant)

! cell parameters:
  CALL get_cell_param ( box, cell, ang )

! initializing pv_const=0
  IF (itimes==0) THEN
     thermo%virial%pv_constraint = 0._dbl
     CALL pressure(struc ( 1 )%pnode,thermo%virial)
  END IF

! computing the scalar pressure
  pv_scalar = 0._dbl
  DO i = 1, 3
     pv_scalar = pv_scalar + thermo%virial%pv_total(i,i)
  END DO
  pv_scalar = pv_scalar/3._dbl/box%deth
  pv_scalar = pv_scalar*p_conv

! average quantities:
  SELECT CASE (itimes)
  CASE (0)
     IF (ionode) THEN
        WRITE (iw,*)
        WRITE ( iw, '( A )' ) ' ENERGY| Initialization proceeding'
        WRITE (iw,*)
     END IF
  CASE (1)
     averages%avepress = pv_scalar
     averages%avetemp = temp
     averages%avepot = thermo%pot
     averages%avekin = thermo%kin
     averages%avevol = box%deth
     averages%aveca = cell(1)
     averages%avecb = cell(2)
     averages%avecc = cell(3)
     averages%aveal = ang(3)
     averages%avebe = ang(2)
     averages%avega = ang(1)
     econs_avg =  0._dbl
  CASE DEFAULT
     averages%avepress = (averages%avepress*real(itimes-1)+pv_scalar)/ &
          REAL ( itimes, dbl )
     averages%avetemp = (averages%avetemp*real(itimes-1)+temp)/ &
          real(itimes)
     averages%avepot = (averages%avepot*real(itimes-1)+thermo%pot)/ &
          real(itimes)
     averages%avekin = (averages%avekin*real(itimes-1)+thermo%kin)/ &
          real(itimes)
     averages%avevol = (averages%avevol*real(itimes-1)+box%deth)/ &
          real(itimes)
     averages%aveca = (averages%aveca*real(itimes-1)+cell(1))/real(itimes)
     averages%avecb = (averages%avecb*real(itimes-1)+cell(2))/real(itimes)
     averages%avecc = (averages%avecc*real(itimes-1)+cell(3))/real(itimes)
     averages%aveal = (averages%aveal*real(itimes-1)+ang(3))/ &
          real(itimes)
     averages%avebe = (averages%avebe*real(itimes-1)+ang(2))/ &
          real(itimes)
     averages%avega = (averages%avega*real(itimes-1)+ang(1))/ &
          real(itimes)
     econs_avg = econs_avg + abs((cons-constant)/constant)
  END SELECT

  IF (ionode .AND. (mod(itimes,iener)==0) .AND. (itimes /= 0 )) THEN
     WRITE (ene, *) itimes, econs, econs_avg/real(itimes)
  END IF

  IF (itimes==0 .AND. ionode) THEN
     WRITE ( iw, '( )' )
     WRITE ( iw, '( A,A )' ) ' ******************************** ', &
          'GO CP2K GO! **********************************'
     WRITE ( iw, '( A,A,T40,A,T65,1(1X,E15.7) )' ) ' INITIAL POTENTIAL ENERGY', &
          adjustl(e_label), '= ', thermo%pot
     WRITE ( iw, '( A,A,T40,A,T65,1(1X,E15.7) )' ) ' INITIAL KINETIC ENERGY', &
          adjustl(e_label), '= ', thermo%kin
     WRITE ( iw, '( A,A,T40,A,T65,1(1X,E15.7) )' ) ' INITIAL TEMPERATURE', &
          adjustl(temp_label), '= ', temp
     WRITE ( iw, '( A,A,T40,A,T65,1(1X,E15.7) )' ) ' INITIAL PRESSURE', &
          adjustl(pres_label), '= ', pv_scalar
     WRITE ( iw, '( A,A,T40,A,T65,1(1X,E15.7) )' ) ' INITIAL VOLUME', &
          adjustl(vol_label), '= ', box%deth
     WRITE ( iw, '( A,A,T29,A,T33,3(1X,E15.7) )' ) ' INITIAL CELL LNTHS', &
          adjustl(l_label), '= ', cell(1), cell(2), cell(3)
     WRITE ( iw, '( A,A,T29,A,T33,3(1X,E15.7) )' ) ' INITIAL CELL ANGLS', &
          adjustl(angl_label), '= ', ang(3), ang(2), ang(1)
     WRITE ( iw, '( A,A )' ) ' ******************************** ', &
          'GO CP2K GO! **********************************'
  END IF

  IF ( ( itimes /= 0 ) .AND. ( MOD ( itimes, iscreen ) == 0 ) &
       .AND. ionode ) THEN

     WRITE ( iw, '( )' )
     WRITE ( iw, '( A,A )' ) ' **************************************', &
          '*****************************************'
     WRITE ( iw, '( A,T61,A )' ) ' ENSEMBLE TYPE                = ', &
          ADJUSTR ( simpar%ensemble)
     WRITE ( iw, '( A,T71,I10 )' ) ' STEP NUMBER                  = ', itimes
     WRITE ( iw, '( A,T65,1(1X,E15.7) )' ) ' CONSERVED QNTY               = ', &
          cons
     WRITE ( iw, '( )' )
     WRITE ( iw, '( T52,A )' ) 'INSTANTANEOUS        AVERAGES'
     WRITE ( iw, '( A,T49,2(1X,E15.7) )' ) ' {E-E0}/E0                    = ', &
          econs, econs_avg/real(itimes)
     WRITE ( iw, '( A,A,T31,A,T49,2(1X,E15.7) )' ) ' POTENTIAL ENERGY', &
          adjustl(e_label), '= ', thermo%pot, averages%avepot
     WRITE ( iw, '( A,A,T31,A,T49,2(1X,E15.7) )' ) ' KINETIC ENERGY', &
          adjustl(e_label), '= ', thermo%kin, averages%avekin
     WRITE ( iw, '( A,A,T31,A,T49,2(1X,E15.7) )' ) ' TEMPERATURE', &
          adjustl(temp_label), '= ', temp, averages%avetemp
     WRITE ( iw, '( A,A,T31,A,T49,2(1X,E15.7) )' ) ' PRESSURE', &
          adjustl(pres_label), '= ', pv_scalar, averages%avepress
     CALL get_ll_parm ( ll_data, lup, aup )
     IF (simpar%ensemble(1:3)=='NPT') THEN
        WRITE ( iw, '( A,A,T31,A,T49,2(1X,E15.7) )' ) ' VOLUME', &
             adjustl(vol_label), '= ', box%deth, averages%avevol
        WRITE ( iw, '( A,A,T31,A,T33,3(1X,E15.7) )' ) ' CELL LNTHS', &
             adjustl(l_label), '= ', cell(1), cell(2), cell(3)
        WRITE ( iw, '( A,A,T31,A,T33,3(1X,E15.7) )' ) ' AVE. CELL LNTHS', &
             adjustl(l_label), '= ', averages%aveca, averages%avecb, &
             averages%avecc
     END IF
     IF (simpar%ensemble=='NPT_F') THEN
        WRITE ( iw, '( A,A,T31,A,T33,3(1X,E15.7) )' ) ' CELL ANGLS', &
             adjustl(angl_label), '= ', ang(3), ang(2), ang(1)
        WRITE ( iw, '( A,A,T31,A,T33,3(1X,E15.7) )' ) ' AVE. CELL ANGLS', &
             adjustl(angl_label), '= ', averages%aveal, averages%avebe, &
             averages%avega
     END IF
     WRITE ( iw, '( A,A,T31,A,T49,2(1X,F15.2) )' ) &
          ' LIST UPDATES', '[steps]', '= ', lup, aup
     WRITE ( iw, '( A,A )' ) ' **************************************', &
          '*****************************************'
     WRITE ( iw, '( )' )
  END IF

  CALL timestop ( zero, handle )

END SUBROUTINE md_energy_fist_pimd

!!*****
!******************************************************************************
!!****** md_fist_energies/get_baro_energies [1.0] *
!!
!!   NAME
!!     get_baro_energies
!!
!!   FUNCTION
!!     Calculates kinetic energy and potential of barostat
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
!!*****
SUBROUTINE get_baro_energies ( box, simpar, nhc, npt, baro_kin, baro_pot )
  IMPLICIT NONE

  TYPE ( cell_type ), INTENT ( IN ) :: box
  TYPE ( simulation_parameters_type), INTENT ( IN ) :: simpar
  TYPE ( lnhc_parameters_type), INTENT ( IN ) :: nhc
  TYPE ( npt_info_type ), INTENT ( IN ), DIMENSION ( :, : ) :: npt
  REAL ( dbl ), INTENT ( OUT ) :: baro_pot, baro_kin

  INTEGER :: i, j

  IF ( simpar%ensemble == 'NPT_I' ) THEN
     baro_pot = simpar % p_ext * box % deth
     baro_kin = 0.5_dbl * npt ( 1, 1 ) % v * npt ( 1, 1 ) % v * &
          npt ( 1, 1 ) % mass
  ELSE IF ( simpar%ensemble=='NPT_F' ) THEN
     baro_pot =  simpar % p_ext * box % deth
     baro_kin = 0._dbl
     DO i = 1, 3
        DO j = 1, 3
           baro_kin = baro_kin + 0.5_dbl * npt ( i, j ) % v *  &
                      npt ( i, j ) % v * npt ( i, j ) % mass
        END DO
     END DO
  END IF

  IF (simpar%ensemble(1:3)=='NPT') THEN
#if defined(__parallel)
     IF (nhc%dis_type=='DIS_REP') THEN
        CALL mp_sum(baro_kin,int_group)
        CALL mp_sum(baro_pot,int_group)
     END IF
#endif
  END IF
END SUBROUTINE get_baro_energies 

!******************************************************************************
!!****** md_fist_energies/get_nhc_energies [1.0] *
!!
!!   NAME
!!     get_nhc_energies
!!
!!   FUNCTION
!!     Calculates kinetic energy and potential energy
!!     of the nhc variables
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
SUBROUTINE get_nhc_energies ( nhc, nhc_pot, nhc_kin )
  IMPLICIT NONE
  TYPE ( lnhc_parameters_type ), INTENT ( IN ) :: nhc
  REAL ( dbl ), INTENT ( OUT ) :: nhc_pot, nhc_kin

  INTEGER :: n, l

   nhc_kin = 0.0_dbl
   nhc_pot = 0.0_dbl
   DO n = 1, nhc % num_nhc
      DO l = 1, nhc % nhc_len
         nhc_kin = nhc_kin +  &
         0.5_dbl * nhc % nvt ( l, n ) % mass * &
         nhc % nvt ( l, n ) % v * nhc % nvt ( l, n ) % v

         nhc_pot = nhc_pot + nhc % nvt ( l, n ) % nkt *  &
                        nhc % nvt ( l, n ) % eta
      END DO
   END DO
#if defined(__parallel)
   IF (nhc%dis_type=='DIS_REP') THEN
      CALL mp_sum(nhc_kin,int_group)
      CALL mp_sum(nhc_pot,int_group)
   END IF
#endif
END SUBROUTINE get_nhc_energies 
!******************************************************************************
!!****** md_fist_energies/get_part_temp [1.0] *
!!
!!   NAME
!!     get_part_temp
!!
!!   FUNCTION
!!     Calculates temperature of particles
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
SUBROUTINE get_part_temp ( thermo, simpar, temperature )

  IMPLICIT NONE
  REAL ( dbl ), INTENT ( OUT ) :: temperature 
  TYPE ( simulation_parameters_type ), INTENT ( IN ) :: simpar
  TYPE ( fist_thermodynamic_type ), INTENT ( INOUT ) :: thermo

  temperature = 2.0_dbl*thermo%kin/float(simpar%nfree)/temp_conv

END SUBROUTINE  get_part_temp
!!*****
!******************************************************************************
!!****** md_fist_energies/get_part_ke [1.0] *
!!
!!   NAME
!!     get_part_ke
!!
!!   FUNCTION
!!     Calculates kinetic energy of particles
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE get_part_ke ( struc, thermo )
  IMPLICIT NONE
  TYPE ( structure_type ), INTENT ( IN ) :: struc
  TYPE ( fist_thermodynamic_type ), INTENT ( INOUT ) :: thermo

  INTEGER :: i, nnodes

  thermo % kin = zero
  nnodes = SIZE ( struc % pnode )
  DO i = 1, nnodes
     thermo%kin = thermo%kin + 0.5_dbl*struc%pnode(i) %p%prop%mass*(struc &
          %pnode(i)%p%v(1)*struc%pnode(i)%p%v(1)+struc%pnode(i)%p%v(2)*struc &
          %pnode(i)%p%v(2)+struc%pnode(i)%p%v(3)*struc%pnode(i)%p%v(3))
  END DO

! sum all contributions to energy over calculated parts on all processors
#if defined(__parallel)
  CALL mp_sum(thermo%kin,int_group)
#endif
END SUBROUTINE get_part_ke 

!******************************************************************************
!
!SUBROUTINE get_coef_temp ( ekinc, inuse, ncoef, ctemp )
!  IMPLICIT NONE
!  REAL (dbl), INTENT ( IN ) :: ekinc
!  INTEGER, INTENT ( IN ) :: inuse, ncoef
!  REAL ( dbl ), INTENT ( OUT ) :: ctemp
!! locals
!
!  SELECT CASE ( inuse )
!  CASE ( 402 )
!    ctemp = 2.0_dbl* ekinc / REAL ( ncoef, dbl ) / temp_conv
!  CASE ( 401 )
!    ctemp =  ekinc / REAL ( ncoef, dbl ) / temp_conv
!  CASE DEFAULT
!     CALL stop_program ( 'get_coef_temp','integrator not implemented')
!  END SELECT
!END SUBROUTINE  get_coef_temp
!
!******************************************************************************
!!****** md_fist_energies/get_coef_ke [1.0] *
!!
!!   NAME
!!     get_coef_ke
!!
!!   FUNCTION
!!     Calculates kinetic energy of coefs
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
!SUBROUTINE get_coef_ke ( coef, inuse, ncoef, ekinc )
!  IMPLICIT NONE
!  TYPE ( coeff_type ), INTENT ( IN ) :: coef
!  INTEGER, INTENT ( IN ) :: inuse, ncoef
!  REAL ( dbl ), INTENT ( out ) :: ekinc
!
!  INTEGER :: i
!  REAL ( dbl ) :: cmass
!  ekinc = zero
!
!  SELECT CASE ( inuse )
!  CASE ( 402 )
!    DO i = 1, ncoef 
!       ekinc = ekinc + 0.5_dbl * coef % ao % mass ( i ) * & 
!               coef % ao % cr ( i ) * coef % ao % cr ( i )
!    END DO
!  CASE ( 401 )
!    DO i = 1, ncoef 
!       cmass = coef % pw % mass_cc ( i )
!       ekinc = ekinc + 0.5_dbl * cmass * &
!        ( REAL ( coef % pw % cc ( i ) ) * REAL ( coef % pw % cc ( i ) ) + &
!          AIMAG ( coef % pw % cc ( i ) )* AIMAG ( coef % pw % cc ( i ) ) )
!    END DO
!  CASE DEFAULT
!     CALL stop_program ( 'get_coef_ke','integrator not implemented')
!  END SELECT 
!
!! sum all contributions to energy over calculated parts on all processors
!#if defined(__parallel)
!  CALL mp_sum(ekinc,int_group)
!#endif
!END SUBROUTINE get_coef_ke
!******************************************************************************
END MODULE md_fist_energies
