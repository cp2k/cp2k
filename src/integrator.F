!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/integrator [1.0] *
!!
!!   NAME
!!     integrator
!!
!!   FUNCTION
!!     Provides integrator routines (velocity verlet) for all the
!!     ensemble types
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     JGH (15-Mar-2001) : Pass logical for box change to force routine
!!     Harald Forbert (Apr-2001): added path integral routine nvt_pimd
!!     CJM (15-Apr-2001) : added coef integrators and energy routines
!!
!!   SOURCE
!******************************************************************************

MODULE integrator
  
  USE brillouin, ONLY : kpoint_type
  USE constraint, ONLY: shake_control, rattle_control, pv_constraint, getold, &
      shake_roll_control, rattle_roll_control
  USE coefficient_types, ONLY : coeff_type, PLANEWAVES, ATOMICORBITALS     
  USE dft_types, ONLY : dft_control_type
  USE eigenvalueproblems, ONLY : diagonalise
  USE ewald_parameters_types, ONLY : ewald_parameters_type
  USE force_control, ONLY : force
  USE global_types, ONLY: global_environment_type
  USE kinds, ONLY : dbl
  USE linklist_control, ONLY : get_ll_parm
  USE mathconstants, ONLY : pi, zero
  USE md, ONLY : simulation_parameters_type, thermodynamic_type, &
       mdio_parameters_type
  USE molecule_types, ONLY : molecule_structure_type, particle_node_type, &
       linklist_atoms
  USE message_passing, ONLY : mp_sum, mp_sync, mp_comm_dup
  USE nose, ONLY : lnhc_parameters_type, npt_info_type, lnhc, yoshida_coef
  USE particle_types, ONLY : particle_type
  USE simulation_cell, ONLY : get_cell_param
  USE structure_types, ONLY : structure_type, interaction_type
  USE termination, ONLY : stop_memory, stop_program
  USE timings, ONLY : timeset, timestop
  USE transformations, ONLY : u2r, fr2fu
  USE util, ONLY : get_unit, matvec_3x3, matmul_3x3, transpose_3d

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: velocity_verlet, energy, set_energy_parm, set_integrator

  TYPE old_variables_type
     REAL ( dbl ), POINTER, DIMENSION ( :, : ) :: v
     REAL ( dbl ), POINTER, DIMENSION ( :, : ) :: r
     REAL ( dbl ), POINTER, DIMENSION ( :, : ) :: eps
     REAL ( dbl ), POINTER, DIMENSION ( :, : ) :: veps
     REAL ( dbl ), POINTER, DIMENSION ( :, : ) :: h
  END TYPE old_variables_type

  TYPE ( old_variables_type ) :: old

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
                        ener_file_name, temp_file_name, dump_file_name, &
                        coef_file_name, epol_file_name
  INTEGER :: crd, vel, ptn, ene, tem, scr, coe, pol
  INTEGER :: icrd, ivel, iptens, iener, itemp, idump, iscreen, icoef, iepol

  TYPE ( global_environment_type ) :: intenv

  INTERFACE energy              
    MODULE PROCEDURE energy_fist, energy_fist_pimd, energy_coef
  END INTERFACE

  INTERFACE velocity_verlet
    MODULE PROCEDURE velocity_verlet_fist, velocity_verlet_coef
  END INTERFACE


!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** integrator/velocity_verlet_coef [1.0] *
!!
!!   NAME
!!     velocity_verlet_coef
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

SUBROUTINE velocity_verlet_coef ( itimes, constant, simpar, inter, &
     thermo, struc, ewald_param, nhc_part, nhc_baro, nhc_coef, npt, &
     dft_control, kp )
  
  IMPLICIT NONE

! Arguments
  INTEGER, INTENT ( IN ) :: itimes
  REAL ( dbl ), INTENT ( INOUT ) :: constant
  TYPE ( simulation_parameters_type ), INTENT ( IN ) :: simpar
  TYPE ( thermodynamic_type ), INTENT ( INOUT ) :: thermo
  TYPE ( structure_type ),INTENT ( INOUT ) :: struc
  TYPE ( lnhc_parameters_type ), INTENT ( INOUT ) :: nhc_part, nhc_baro
  TYPE ( lnhc_parameters_type ), INTENT ( INOUT ) :: nhc_coef 
  TYPE ( npt_info_type ), INTENT ( INOUT ), DIMENSION ( :, : ) :: npt
  TYPE ( ewald_parameters_type ), INTENT ( INOUT ) :: ewald_param
  TYPE ( interaction_type ), INTENT ( IN ) :: inter
  TYPE ( dft_control_type ) :: dft_control
  TYPE ( kpoint_type ), INTENT ( IN ) :: kp
  
  INTEGER :: handle

!------------------------------------------------------------------------------
  CALL timeset ( 'VERLET', 'I', ' ', handle )
  SELECT CASE ( simpar % ensemble )
  CASE DEFAULT
     CALL stop_program ( 'velocity_verlet','integrator not implemented')
  CASE ( 'NVE')
     CALL nve_coef ( itimes, constant, simpar, inter, thermo, struc, &
                    ewald_param, nhc_coef, dft_control, kp )
  CASE ( 'NVT')
     CALL nvt_coef ( itimes, constant, simpar, inter, thermo, struc, &
                    ewald_param, nhc_part, nhc_coef, dft_control, kp )
  END SELECT

  CALL energy_coef( itimes, constant, simpar, struc, thermo, nhc_part, & 
                   nhc_baro, nhc_coef, npt )

  CALL timestop ( zero, handle )

END SUBROUTINE velocity_verlet_coef

!!*****
!******************************************************************************
!!****** integrator/velocity_verlet_fist [1.0] *
!!
!!   NAME
!!     velocity_verlet_fist
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

SUBROUTINE velocity_verlet_fist ( itimes, constant, simpar, inter, &
     thermo, struc, ewald_param, nhc_part, nhc_baro, npt, dft_control, kp )
  
  IMPLICIT NONE

! Arguments
  INTEGER, INTENT ( IN ) :: itimes
  REAL ( dbl ), INTENT ( INOUT ) :: constant
  TYPE ( simulation_parameters_type ), INTENT ( IN ) :: simpar
  TYPE ( thermodynamic_type ), INTENT ( INOUT ) :: thermo
  TYPE ( structure_type ), DIMENSION ( : ), INTENT ( INOUT ) :: struc
  TYPE ( lnhc_parameters_type ), DIMENSION ( : ), &
                              INTENT ( INOUT ) :: nhc_part, nhc_baro
  TYPE ( npt_info_type ), INTENT ( INOUT ), DIMENSION ( :, : ) :: npt
  TYPE ( ewald_parameters_type ), INTENT ( INOUT ) :: ewald_param
  TYPE ( interaction_type ), INTENT ( IN ) :: inter
  TYPE ( dft_control_type ) :: dft_control
  TYPE ( kpoint_type ), INTENT ( IN ) :: kp
  
  INTEGER :: handle

!------------------------------------------------------------------------------
  CALL timeset ( 'VERLET', 'I', ' ', handle )
  IF ( simpar % pimd_params % beads > 0 ) THEN
     SELECT CASE ( simpar % ensemble )
     CASE DEFAULT
        CALL stop_program ( 'velocity','PIMD not possible yet with this integrator' )
     CASE ( 'NVT' )
        CALL nvt_pimd ( itimes, constant, simpar, inter, thermo, struc, &
                      ewald_param, nhc_part, dft_control, kp )
     END SELECT
     CALL energy_fist_pimd ( itimes, constant, simpar, struc, thermo, nhc_part, & 
                             nhc_baro, npt )
  ELSE
     SELECT CASE ( simpar % ensemble )
     CASE DEFAULT
        CALL stop_program ( 'velocity_verlet','integrator not implemented')
     CASE ( 'NVE')
        CALL nve ( itimes, constant, simpar, inter, thermo, struc ( 1 ), ewald_param, &
             dft_control, kp )
     CASE ( 'NVT')
        CALL nvt ( itimes, constant, simpar, inter, thermo, struc ( 1 ), ewald_param, &
             nhc_part ( 1 ), dft_control, kp )
     CASE ( 'NPT_I')
        CALL npt_i( itimes, constant, simpar, inter, thermo, struc ( 1 ), ewald_param, &
                   nhc_part ( 1 ), nhc_baro ( 1 ), npt, dft_control, kp )
     CASE ( 'NPT_F')
        CALL npt_f( itimes, constant, simpar, inter, thermo, struc ( 1 ), ewald_param, &
                   nhc_part ( 1 ), nhc_baro ( 1 ), npt, dft_control, kp )
     END SELECT
     CALL energy_fist ( itimes, constant, simpar, struc ( 1 ), thermo, nhc_part ( 1 ), &
                       nhc_baro ( 1 ), npt)
  END IF

  CALL timestop ( zero, handle )

END SUBROUTINE velocity_verlet_fist

!!*****
!******************************************************************************
!!****** integrator/set_integrator [1.0] *
!!
!!   NAME
!!     set_integrator
!!
!!   FUNCTION
!!     setup of parameters for integrator
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


SUBROUTINE set_integrator ( globenv, mdio )

  IMPLICIT NONE

! Arguments
  TYPE(mdio_parameters_type), INTENT (IN) :: mdio
  TYPE ( global_environment_type ), INTENT (IN) :: globenv

!------------------------------------------------------------------------------

  intenv = globenv
  CALL mp_comm_dup ( globenv % group, intenv % group )

  ionode = intenv % ionode
  int_group = intenv % group
  scr = intenv % scr

  crd_file_name = mdio % crd_file_name
  vel_file_name = mdio % vel_file_name
  ptens_file_name = mdio % ptens_file_name
  ener_file_name = mdio % ener_file_name
  temp_file_name = mdio % temp_file_name
  dump_file_name = mdio % dump_file_name
  coef_file_name = mdio % coef_file_name
  epol_file_name = mdio % epol_file_name

  icrd = mdio % icrd
  ivel = mdio % ivel
  iptens = mdio % iptens
  iener = mdio % iener
  itemp = mdio % itemp
  idump = mdio % idump
  icoef = mdio % icoef
  iepol = mdio % iepol
  iscreen = mdio % iscreen

END SUBROUTINE set_integrator

!!*****
!******************************************************************************
!!****** integrator/set_energy_parm [1.0] *
!!
!!   NAME
!!     set_energy_parm
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

SUBROUTINE set_energy_parm ( pconv, tempconv, l, vol, e, pres, temp, angl )

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
END SUBROUTINE set_energy_parm

!!*****
!******************************************************************************
!!****** integrator/energy_coef [1.0] *
!!
!!   NAME
!!     energy_coef
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


SUBROUTINE energy_coef ( itimes, constant, simpar, struc, thermo,  &
                    nhc_part, nhc_baro, nhc_coef, npt )

  IMPLICIT NONE

! Arguments
  INTEGER, INTENT ( IN ) :: itimes
  REAL ( dbl ), INTENT ( INOUT ) :: constant
  TYPE (simulation_parameters_type ), INTENT ( IN ) :: simpar
  TYPE (structure_type ), INTENT ( INOUT ) :: struc
  TYPE (thermodynamic_type ), INTENT ( INOUT ) :: thermo
  TYPE (lnhc_parameters_type ), INTENT ( INOUT ) :: nhc_part, nhc_baro, nhc_coef
  TYPE (npt_info_type ), INTENT ( INOUT ), DIMENSION ( :, : ) :: npt

! Locals
  INTEGER i, j, id, n, l, nnodes, ncoef, iw, handle, inuse
  REAL ( dbl ) :: nhc_part_pot, nhc_part_kin, baro_kin, baro_pot
  REAL ( dbl ) :: nhc_baro_pot, nhc_baro_kin
  REAL ( dbl ) :: alpha, epsilon0, qsq, lup, aup
  REAL ( dbl ) :: cons, temp, econs
  REAL ( dbl ) :: coef_temp, coef_kin
  REAL ( dbl ) :: nhc_coef_pot, nhc_coef_kin
  REAL ( dbl ), DIMENSION(3) :: cell, ang
  REAL ( dbl ) :: pv_scalar

!------------------------------------------------------------------------------

  CALL timeset ( 'ENERGY', 'E', ' ', handle )
  
  nnodes = SIZE ( struc % pnode )
  inuse = struc % coef_pos ( 1 ) % in_use
  SELECT CASE ( inuse )
  CASE ( 402 )
    ncoef = SIZE ( struc % coef_vel ( 1 ) % ao % cr )
  CASE ( 401 )
    ncoef = SIZE ( struc % coef_vel ( 1 ) % pw % cc )
  END SELECT

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
     coe = get_unit()
     OPEN ( UNIT = coe, FILE = coef_file_name )
     pol = get_unit()
     OPEN ( UNIT = pol, FILE = epol_file_name )
  END IF

  iw = scr

  CALL get_part_ke ( struc, thermo )
  
  CALL get_part_temp ( thermo, simpar, temp )

  CALL get_coef_ke ( struc % coef_vel ( 1 ), inuse, ncoef, coef_kin )
  
  CALL get_coef_temp ( coef_kin, inuse, ncoef, coef_temp )


  IF (ionode .AND. mod(itimes,itemp)==0) THEN
     WRITE (tem, *) itimes, temp, coef_temp
  END IF

  IF (ionode .AND. mod(itimes,iptens)==0) THEN
     DO i=1,3
        WRITE (ptn,*) itimes, (thermo%ptens(i,j), j=1,3)
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

  IF (ionode .AND. mod(itimes,icoef)==0) THEN
    SELECT CASE ( struc % coef_pos ( 1 ) % in_use )
    CASE ( 402 )
       DO i=1, ncoef
          WRITE (coe,'(i6,3f20.12 )' ) itimes, struc%coef_pos(1) %ao%cr (i), &
                 struc%coef_vel(1) %ao%cr (i),  struc%coef_force(1) %ao%cr (i)
       END DO
    CASE ( 401 )
       DO i=1, ncoef
          WRITE (coe, *) itimes, struc%coef_pos(1)%pw%cc (i), &
                 struc%coef_vel(1) %pw%cc (i),  struc%coef_force(1) %pw%cc (i)
       END DO
    END SELECT
  END IF

  IF (ionode .AND. mod(itimes,iepol)==0) THEN
     WRITE (pol,*) itimes,  thermo%e_class, thermo % e_electro, &      
            thermo % e_overlap, thermo% e_thk, thermo% pot                       
  END IF

  IF ( simpar % pol_ensemble ( 1 : 3 )/= 'NVE' ) &
  CALL get_nhc_energies ( nhc_coef, nhc_coef_kin, nhc_coef_pot )

  IF ( simpar % ensemble ( 1 : 3 )/= 'NVE' ) &
  CALL get_nhc_energies ( nhc_part, nhc_part_kin, nhc_part_pot )

  IF ( simpar % ensemble ( 1 : 3 ) == 'NPT' ) &
  CALL get_nhc_energies (  nhc_baro, nhc_baro_kin, nhc_baro_pot )

  IF ( simpar % ensemble ( 1 : 3 ) == 'NPT' ) &
  CALL get_baro_energies ( struc, simpar, nhc_baro, npt, baro_kin, baro_pot )



  SELECT CASE ( simpar % ensemble ( 1 : 3 ) )
  CASE ( 'NVE')
     SELECT CASE ( simpar % pol_ensemble ( 1 : 3 ) )
     CASE ( 'NVE')
        cons = thermo%kin + thermo%pot + coef_kin
     CASE ( 'NVT')
        cons = thermo%kin + thermo%pot + nhc_coef_kin + nhc_coef_pot + coef_kin
     END SELECT
  CASE ( 'NVT')
     SELECT CASE ( simpar % pol_ensemble ( 1 : 3 ) )
     CASE ( 'NVE')
        cons = thermo%kin + thermo%pot + nhc_part_kin + nhc_part_pot + coef_kin 
     CASE ( 'NVT')
        cons = thermo%kin + thermo%pot + nhc_part_kin + nhc_part_pot + &
                nhc_coef_kin + nhc_coef_pot + coef_kin
     END SELECT 
! to check
!  CASE ( 'NPT')
!     SELECT CASE ( simpar % pol_ensemble ( 1 : 3 ) )
!     CASE ( 'NVE')
!        cons = thermo%kin + thermo%pot + nhc_part_kin + nhc_part_pot  &
!              + nhc_baro_kin + nhc_baro_pot + baro_kin + baro_pot + coef_kin
!     CASE ( 'NVT')
!        cons = thermo%kin + thermo%pot + nhc_part_kin + nhc_part_pot  &
!              + nhc_baro_kin + nhc_baro_pot + baro_kin + baro_pot + &
!                 nhc_coef_kin + nhc_coef_pot + coef_kin
!      END SELECT
  END SELECT

  IF (itimes==1) THEN
     constant = cons
  END IF
  IF (itimes /= 0 ) econs = abs((cons-constant)/constant)

! cell parameters:
  CALL get_cell_param ( struc % box, cell, ang )

! initializing pv_const=0
  IF (itimes==0) THEN
     thermo%pv_const = 0._dbl
     CALL pressure(struc%pnode,thermo)
  END IF

! computing the scalar pressure
  pv_scalar = 0._dbl
  DO i = 1, 3
     pv_scalar = pv_scalar + thermo%ptens(i,i)
  END DO
  pv_scalar = pv_scalar/3._dbl/struc%box%deth
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
     averages%avekinc = coef_kin
     averages%avetempc = coef_temp
     averages%avevol = struc%box%deth
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
     averages%avetempc = (averages%avetempc*real(itimes-1)+coef_temp)/ &
          real(itimes)
     averages%avepot = (averages%avepot*real(itimes-1)+thermo%pot)/ &
          real(itimes)
     averages%avekin = (averages%avekin*real(itimes-1)+thermo%kin)/ &
          real(itimes)
     averages%avekinc = (averages%avekinc*real(itimes-1)+ coef_kin)/ &
          real(itimes)
     averages%avevol = (averages%avevol*real(itimes-1)+struc%box%deth)/ &
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
     WRITE ( iw, '( A,A,T40,A,T65,1(1X,E15.7) )' ) ' INITIAL COEF KINETIC ENERGY', &
          adjustl(e_label), '= ', coef_kin  
     WRITE ( iw, '( A,A,T40,A,T65,1(1X,E15.7) )' ) ' INITIAL TEMPERATURE', &
          adjustl(temp_label), '= ', temp
     WRITE ( iw, '( A,A,T40,A,T65,1(1X,E15.7) )' ) ' INITIAL COEF TEMPERATURE', &
          adjustl(temp_label), '= ', coef_temp
     WRITE ( iw, '( A,A,T40,A,T65,1(1X,E15.7) )' ) ' INITIAL PRESSURE', &
          adjustl(pres_label), '= ', pv_scalar
     WRITE ( iw, '( A,A,T40,A,T65,1(1X,E15.7) )' ) ' INITIAL VOLUME', &
          adjustl(vol_label), '= ', struc%box%deth
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
     WRITE ( iw, '( A,A,T31,A,T49,2(1X,E15.7) )' ) ' COEF KINETIC ENERGY', &
          adjustl(e_label), '= ', coef_kin,  averages%avekinc
     WRITE ( iw, '( A,A,T31,A,T49,2(1X,E15.7) )' ) ' TEMPERATURE', &
          adjustl(temp_label), '= ', temp, averages%avetemp
     WRITE ( iw, '( A,A,T31,A,T49,2(1X,E15.7) )' ) ' COEF TEMPERATURE', &
          adjustl(temp_label), '= ', coef_temp, averages%avetempc
     WRITE ( iw, '( A,A,T31,A,T49,2(1X,E15.7) )' ) ' PRESSURE', &
          adjustl(pres_label), '= ', pv_scalar, averages%avepress
     CALL get_ll_parm ( struc % ll_data (1), lup, aup )
     IF (simpar%ensemble(1:3)=='NPT') THEN
        WRITE ( iw, '( A,A,T31,A,T49,2(1X,E15.7) )' ) ' VOLUME', &
             adjustl(vol_label), '= ', struc%box%deth, averages%avevol
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

END SUBROUTINE energy_coef

!!*****
!******************************************************************************
!!****** integrator/energy_fist [1.0] *
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


SUBROUTINE energy_fist ( itimes, constant, simpar, struc, thermo,  &
                    nhc_part, nhc_baro, npt )

  IMPLICIT NONE

! Arguments
  INTEGER, INTENT ( IN ) :: itimes
  REAL ( dbl ), INTENT ( INOUT ) :: constant
  TYPE (simulation_parameters_type ), INTENT ( IN ) :: simpar
  TYPE (structure_type ), INTENT ( INOUT ) :: struc
  TYPE (thermodynamic_type ), INTENT ( INOUT ) :: thermo
  TYPE (lnhc_parameters_type ), INTENT ( INOUT ) :: nhc_part, nhc_baro
  TYPE (npt_info_type ), INTENT ( INOUT ), DIMENSION ( :, : ) :: npt

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
        WRITE (ptn,*) itimes, (thermo%ptens(i,j), j=1,3)
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
  CALL get_baro_energies ( struc, simpar, nhc_baro, npt, baro_kin, baro_pot )



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
  CALL get_cell_param ( struc % box, cell, ang )

! initializing pv_const=0
  IF (itimes==0) THEN
     thermo%pv_const = 0._dbl
     CALL pressure(struc%pnode,thermo)
  END IF

! computing the scalar pressure
  pv_scalar = 0._dbl
  DO i = 1, 3
     pv_scalar = pv_scalar + thermo%ptens(i,i)
  END DO
  pv_scalar = pv_scalar/3._dbl/struc%box%deth
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
     averages%avevol = struc % box % deth
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
     averages%avevol = (averages%avevol*real(itimes-1)+struc%box%deth)/ &
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
          adjustl(vol_label), '= ', struc%box%deth
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
     CALL get_ll_parm ( struc % ll_data (1), lup, aup )
     IF (simpar%ensemble(1:3)=='NPT') THEN
        WRITE ( iw, '( A,A,T31,A,T49,2(1X,E15.7) )' ) ' VOLUME', &
             adjustl(vol_label), '= ', struc%box%deth, averages%avevol
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

END SUBROUTINE energy_fist

!!*****
!******************************************************************************
!!****** integrator/energy_fist_pimd [1.0] *
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


SUBROUTINE energy_fist_pimd ( itimes, constant, simpar, struc, thermo,  &
                    nhc_part, nhc_baro, npt )

  IMPLICIT NONE

! Arguments
  INTEGER, INTENT ( IN ) :: itimes
  REAL ( dbl ), INTENT ( INOUT ) :: constant
  TYPE (simulation_parameters_type ), INTENT ( IN ) :: simpar
  TYPE ( structure_type ), DIMENSION ( : ), INTENT ( INOUT ) :: struc
  TYPE (thermodynamic_type ), INTENT ( INOUT ) :: thermo
  TYPE (lnhc_parameters_type ), DIMENSION ( : ), INTENT ( INOUT ) :: &
                                                 nhc_part, nhc_baro
  TYPE (npt_info_type ), INTENT ( INOUT ), DIMENSION ( :, : ) :: npt

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
        WRITE (ptn,*) itimes, (thermo%ptens(i,j), j=1,3)
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
  CALL get_baro_energies ( struc ( 1 ), simpar, nhc_baro ( 1 ), npt,  &
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
  CALL get_cell_param ( struc ( 1 ) % box, cell, ang )

! initializing pv_const=0
  IF (itimes==0) THEN
     thermo%pv_const = 0._dbl
     CALL pressure(struc ( 1 )%pnode,thermo)
  END IF

! computing the scalar pressure
  pv_scalar = 0._dbl
  DO i = 1, 3
     pv_scalar = pv_scalar + thermo%ptens(i,i)
  END DO
  pv_scalar = pv_scalar/3._dbl/struc ( 1 )%box%deth
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
     averages%avevol = struc ( 1 )%box%deth
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
     averages%avevol = (averages%avevol*real(itimes-1)+struc ( 1 )%box%deth)/ &
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
          adjustl(vol_label), '= ', struc ( 1 )%box%deth
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
     CALL get_ll_parm ( struc ( 1 ) % ll_data (1), lup, aup )
     IF (simpar%ensemble(1:3)=='NPT') THEN
        WRITE ( iw, '( A,A,T31,A,T49,2(1X,E15.7) )' ) ' VOLUME', &
             adjustl(vol_label), '= ', struc ( 1 )%box%deth, averages%avevol
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

END SUBROUTINE energy_fist_pimd

!!*****
!******************************************************************************
!!****** integrator/get_baro_energies [1.0] *
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
SUBROUTINE get_baro_energies ( struc, simpar, nhc, npt, baro_kin, baro_pot )
  IMPLICIT NONE

  TYPE ( structure_type ), INTENT (IN) :: struc
  TYPE ( simulation_parameters_type), INTENT ( IN ) :: simpar
  TYPE ( lnhc_parameters_type), INTENT ( IN ) :: nhc
  TYPE ( npt_info_type ), INTENT ( INOUT ), DIMENSION ( :, : ) :: npt
  REAL ( dbl ), INTENT ( OUT ) :: baro_pot, baro_kin

  INTEGER :: i, j

  IF ( simpar%ensemble == 'NPT_I' ) THEN
     baro_pot = simpar % p_ext * struc % box % deth
     baro_kin = 0.5_dbl * npt ( 1, 1 ) % v * npt ( 1, 1 ) % v * &
          npt ( 1, 1 ) % mass
  ELSE IF ( simpar%ensemble=='NPT_F' ) THEN
     baro_pot =  simpar % p_ext * struc % box % deth
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
!!****** integrator/get_nhc_energies [1.0] *
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
!!*****
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
!!****** integrator/get_part_temp [1.0] *
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
  TYPE ( thermodynamic_type ), INTENT ( INOUT ) :: thermo

  temperature = 2.0_dbl*thermo%kin/float(simpar%nfree)/temp_conv
END SUBROUTINE  get_part_temp
!!*****
!******************************************************************************
!!****** integrator/get_part_ke [1.0] *
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
  TYPE ( thermodynamic_type ), INTENT ( INOUT ) :: thermo

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

SUBROUTINE get_coef_temp ( ekinc, inuse, ncoef, ctemp )
  IMPLICIT NONE
  REAL (dbl), INTENT ( IN ) :: ekinc
  INTEGER, INTENT ( IN ) :: inuse, ncoef
  REAL ( dbl ), INTENT ( OUT ) :: ctemp
! locals

  SELECT CASE ( inuse )
  CASE ( 402 )
    ctemp = 2.0_dbl* ekinc / REAL ( ncoef, dbl ) / temp_conv
  CASE ( 401 )
    ctemp =  ekinc / REAL ( ncoef, dbl ) / temp_conv
  CASE DEFAULT
     CALL stop_program ( 'get_coef_temp','integrator not implemented')
  END SELECT
END SUBROUTINE  get_coef_temp
!!*****
!******************************************************************************
!!****** integrator/get_coef_ke [1.0] *
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

SUBROUTINE get_coef_ke ( coef, inuse, ncoef, ekinc )
  IMPLICIT NONE
  TYPE ( coeff_type ), INTENT ( IN ) :: coef
  INTEGER, INTENT ( IN ) :: inuse, ncoef
  REAL ( dbl ), INTENT ( out ) :: ekinc

  INTEGER :: i
  REAL ( dbl ) :: cmass
  ekinc = zero

  SELECT CASE ( inuse )
  CASE ( 402 )
    DO i = 1, ncoef 
       ekinc = ekinc + 0.5_dbl * coef % ao % mass ( i ) * & 
               coef % ao % cr ( i ) * coef % ao % cr ( i )
    END DO
  CASE ( 401 )
    DO i = 1, ncoef 
       cmass = coef % pw % mass_cc ( i )
       ekinc = ekinc + 0.5_dbl * cmass * &
        ( REAL ( coef % pw % cc ( i ) ) * REAL ( coef % pw % cc ( i ) ) + &
          AIMAG ( coef % pw % cc ( i ) )* AIMAG ( coef % pw % cc ( i ) ) )
    END DO
  CASE DEFAULT
     CALL stop_program ( 'get_coef_ke','integrator not implemented')
  END SELECT 

! sum all contributions to energy over calculated parts on all processors
#if defined(__parallel)
  CALL mp_sum(ekinc,int_group)
#endif
END SUBROUTINE get_coef_ke

!!*****
!******************************************************************************
!!****** integrator/nvt_coef [1.0] *
!!
!!   NAME
!!     nvt_coef
!!
!!   FUNCTION
!!     Parses data bases on "in_use" component of coef and
!!     sends to appropriate integrator
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE nvt_coef ( itimes, constant, simpar, inter, thermo,  &
           struc, ewald_param, nhc_part, nhc_coef, dft_control, kp )

  IMPLICIT NONE

! Arguments
  INTEGER, INTENT ( IN ) :: itimes
  REAL ( dbl ), INTENT ( INOUT ) :: constant
  TYPE ( simulation_parameters_type ), INTENT ( IN ) :: simpar
  TYPE ( interaction_type ), INTENT ( IN ) :: inter
  TYPE ( thermodynamic_type ), INTENT ( INOUT ) :: thermo
  TYPE ( structure_type ), INTENT ( INOUT ) :: struc
  TYPE ( ewald_parameters_type ), INTENT ( INOUT ) :: ewald_param
  TYPE ( lnhc_parameters_type ), INTENT ( INOUT ) :: nhc_coef
  TYPE ( lnhc_parameters_type ), INTENT ( INOUT ) :: nhc_part
  TYPE ( dft_control_type ) :: dft_control
  TYPE ( kpoint_type ), INTENT ( IN ) :: kp

  SELECT CASE ( struc % coef_pos (1) % in_use )
  CASE ( 402 )
    CALL nvt_ao ( itimes, constant, simpar, inter, thermo,  &
         struc, ewald_param, nhc_part, nhc_coef, dft_control, kp )
  CASE ( 401 )
    CALL nvt_pw ( itimes, constant, simpar, inter, thermo, struc, &
                  ewald_param, nhc_part, nhc_coef, dft_control, kp )
  CASE DEFAULT
     CALL stop_program ( 'nvt_coef','integrator not implemented')
  END SELECT 
  
END SUBROUTINE nvt_coef 

!!*****
!******************************************************************************
!!****** integrator/nve_coef [1.0] *
!!
!!   NAME
!!     nve_coef
!!
!!   FUNCTION
!!     Parses data bases on "in_use" component of coef and
!!     sends to appropriate integrator
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE nve_coef ( itimes, constant, simpar, inter, thermo,  &
                      struc, ewald_param, nhc_coef, dft_control, kp )
  IMPLICIT NONE

! Arguments
  INTEGER, INTENT ( IN ) :: itimes
  REAL ( dbl ), INTENT ( INOUT ) :: constant
  TYPE (simulation_parameters_type ), INTENT ( IN ) :: simpar
  TYPE (interaction_type ), INTENT ( IN ) :: inter
  TYPE (thermodynamic_type ), INTENT ( INOUT ) :: thermo
  TYPE (structure_type ), INTENT ( INOUT ) :: struc
  TYPE (ewald_parameters_type ), INTENT ( INOUT ) :: ewald_param
  TYPE ( dft_control_type ) :: dft_control
  TYPE ( kpoint_type ), INTENT ( IN ) :: kp
  TYPE ( lnhc_parameters_type ), INTENT ( INOUT ) :: nhc_coef

  SELECT CASE ( struc % coef_pos ( 1 ) % in_use )
  CASE ( 402 )
    CALL nve_ao ( itimes, constant, simpar, inter, thermo,  &
                  struc, ewald_param, nhc_coef, dft_control, kp )
  CASE ( 401 )
    CALL nve_pw ( itimes, constant, simpar, inter, thermo,  &
                  struc, ewald_param, nhc_coef, dft_control, kp )
  CASE DEFAULT
     CALL stop_program ( 'nve_coef','integrator not implemented')
  END SELECT 
  

END SUBROUTINE nve_coef 
!!*****
!******************************************************************************
!!****** integrator/nve_pw [1.0] *
!!
!!   NAME
!!     nve_ao
!!
!!   FUNCTION
!!     nve integrator for particle positions & momenta
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE nve_pw ( itimes, constant, simpar, inter, thermo, struc, &
                    ewald_param, nhc_coef, dft_control, kp )
  
  IMPLICIT NONE
  
! Arguments
  INTEGER, INTENT ( IN ) :: itimes
  REAL ( dbl ), INTENT ( INOUT ) :: constant
  TYPE (simulation_parameters_type ), INTENT ( IN ) :: simpar
  TYPE (thermodynamic_type ), INTENT ( INOUT ) :: thermo
  TYPE (structure_type ), INTENT ( INOUT ) :: struc
  TYPE (ewald_parameters_type ), INTENT ( INOUT ) :: ewald_param
  TYPE (interaction_type ), INTENT ( IN ) :: inter
  TYPE ( dft_control_type ) :: dft_control
  TYPE ( kpoint_type ), INTENT ( IN ) :: kp
  TYPE ( lnhc_parameters_type ), INTENT ( INOUT ) :: nhc_coef
  
! Local
  INTEGER :: i, nnodes, ncoef
  REAL ( dbl ) :: dtom, dtocm, cmass
  LOGICAL, PARAMETER :: box_change = .FALSE.

!------------------------------------------------------------------------------
  IF ( simpar % pol_ensemble == "NVT" ) & 
  CALL lnhc ( nhc_coef, struc % coef_vel ( 1 ) % pw, int_group )

  IF ( simpar % constraint ) CALL getold ( struc % molecule )
!
! first half of velocity verlet
!
  nnodes = size(struc%pnode)
  DO i = 1, nnodes
     dtom = 0.5_dbl*simpar%dt/struc%pnode(i) %p%prop%mass
     struc%pnode(i) %p%v(1) = struc%pnode(i) %p%v(1) + &
          dtom*struc%pnode(i) %p%f(1)
     struc%pnode(i) %p%v(2) = struc%pnode(i) %p%v(2) + &
          dtom*struc%pnode(i) %p%f(2)
     struc%pnode(i) %p%v(3) = struc%pnode(i) %p%v(3) + &
          dtom*struc%pnode(i) %p%f(3)
     struc%pnode(i) %p%r(1) = struc%pnode(i) %p%r(1) + &
          struc%pnode(i) %p%v(1)*simpar%dt
     struc%pnode(i) %p%r(2) = struc%pnode(i) %p%r(2) + &
          struc%pnode(i) %p%v(2)*simpar%dt
     struc%pnode(i) %p%r(3) = struc%pnode(i) %p%r(3) + &
          struc%pnode(i) %p%v(3)*simpar%dt
  END DO

  IF (simpar%constraint) CALL shake_control(struc%molecule,simpar)

  ncoef = SIZE ( struc % coef_pos ( 1 ) % pw % cc )
  DO i = 1, ncoef
     cmass = REAL ( struc % coef_vel ( 1 ) % pw % mass_cc ( i ), dbl ) 
     dtocm = 0.5_dbl * simpar % dt / cmass
     struc % coef_vel ( 1 ) % pw % cc ( i ) = struc % coef_vel ( 1 ) % pw  % cc ( i ) &
                     + dtocm * struc % coef_force ( 1 ) % pw % cc ( i )
     struc % coef_pos ( 1 ) % pw % cc ( i ) = struc % coef_pos ( 1 ) % pw % cc ( i ) &
             + struc % coef_vel ( 1 ) % pw % cc ( i ) * simpar % dt
  END DO

#if defined(__parallel)
  CALL update_structure(struc,'POS')
#endif

!
! get new forces
!
  CALL force ( struc, inter, thermo, simpar, ewald_param, box_change, &
               intenv, dft_control, kp )
  
!
! second half of velocity verlet
!
  DO i = 1, nnodes
     dtom = 0.5_dbl*simpar%dt/struc%pnode(i) %p%prop%mass
     struc%pnode(i) %p%v(1) = struc%pnode(i) %p%v(1) + &
          dtom*struc%pnode(i) %p%f(1)
     struc%pnode(i) %p%v(2) = struc%pnode(i) %p%v(2) + &
          dtom*struc%pnode(i) %p%f(2)
     struc%pnode(i) %p%v(3) = struc%pnode(i) %p%v(3) + &
          dtom*struc%pnode(i) %p%f(3)
  END DO

  DO i = 1, ncoef
     cmass = REAL ( struc % coef_vel ( 1 ) % pw % mass_cc ( i ), dbl ) 
     dtocm = 0.5_dbl * simpar % dt / cmass
     struc % coef_vel ( 1 ) % pw % cc ( i ) = struc % coef_vel ( 1 ) % pw % cc ( i ) &
                      + dtocm * struc % coef_force ( 1 ) % pw % cc ( i )
  END DO

  IF ( simpar % pol_ensemble == "NVT" ) & 
  CALL lnhc ( nhc_coef, struc % coef_vel ( 1 ) % pw, int_group )

  IF (simpar%constraint) CALL rattle_control(struc%molecule,simpar)
#if defined(__parallel)
  CALL update_structure(struc,'VEL')
#endif

! updating the virial
  IF (simpar%constraint) CALL pv_constraint(struc%molecule,thermo, &
                                            'FULL',int_group )

  CALL pressure(struc%pnode,thermo)

END SUBROUTINE nve_pw

!!*****
!******************************************************************************
!!****** integrator/nve_ao [1.0] *
!!
!!   NAME
!!     nve_ao
!!
!!   FUNCTION
!!     nve integrator for particle positions & momenta
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE nve_ao ( itimes, constant, simpar, inter, thermo, struc, &
                    ewald_param, nhc_coef, dft_control, kp )
  
  IMPLICIT NONE
  
! Arguments
  INTEGER, INTENT ( IN ) :: itimes
  REAL ( dbl ), INTENT ( INOUT ) :: constant
  TYPE (simulation_parameters_type ), INTENT ( IN ) :: simpar
  TYPE (thermodynamic_type ), INTENT ( INOUT ) :: thermo
  TYPE (structure_type ), INTENT ( INOUT ) :: struc
  TYPE (ewald_parameters_type ), INTENT ( INOUT ) :: ewald_param
  TYPE (interaction_type ), INTENT ( IN ) :: inter
  TYPE ( dft_control_type ) :: dft_control
  TYPE ( kpoint_type ), INTENT ( IN ) :: kp
  TYPE ( lnhc_parameters_type ), INTENT ( INOUT ) :: nhc_coef
  
! Local
  INTEGER :: i, nnodes, ncoef
  REAL ( dbl ) :: dtom, dtocm
  LOGICAL, PARAMETER :: box_change = .FALSE.

!------------------------------------------------------------------------------
  IF ( simpar % pol_ensemble == "NVT" ) & 
  CALL lnhc ( nhc_coef, struc % coef_vel ( 1 ) % ao, int_group )

  IF ( simpar % constraint ) CALL getold ( struc % molecule )
!
! first half of velocity verlet
!
  nnodes = size(struc%pnode)
  DO i = 1, nnodes
     dtom = 0.5_dbl*simpar%dt/struc%pnode(i) %p%prop%mass
     struc%pnode(i) %p%v(1) = struc%pnode(i) %p%v(1) + &
          dtom*struc%pnode(i) %p%f(1)
     struc%pnode(i) %p%v(2) = struc%pnode(i) %p%v(2) + &
          dtom*struc%pnode(i) %p%f(2)
     struc%pnode(i) %p%v(3) = struc%pnode(i) %p%v(3) + &
          dtom*struc%pnode(i) %p%f(3)
     struc%pnode(i) %p%r(1) = struc%pnode(i) %p%r(1) + &
          struc%pnode(i) %p%v(1)*simpar%dt
     struc%pnode(i) %p%r(2) = struc%pnode(i) %p%r(2) + &
          struc%pnode(i) %p%v(2)*simpar%dt
     struc%pnode(i) %p%r(3) = struc%pnode(i) %p%r(3) + &
          struc%pnode(i) %p%v(3)*simpar%dt
  END DO

  IF (simpar%constraint) CALL shake_control(struc%molecule,simpar)

  ncoef = SIZE ( struc % coef_pos ( 1 ) % ao % cr )
  DO i = 1, ncoef
     dtocm = 0.5_dbl * simpar % dt / struc % coef_vel ( 1 ) % ao % mass ( i )
     struc % coef_vel ( 1 ) % ao % cr ( i ) = struc % coef_vel ( 1 ) % ao  % cr ( i ) &
                     + dtocm * struc % coef_force ( 1 ) % ao % cr ( i )
     struc % coef_pos ( 1 ) % ao % cr ( i ) = struc % coef_pos ( 1 ) % ao % cr ( i ) &
             + struc % coef_vel ( 1 ) % ao % cr ( i ) * simpar % dt
  END DO

#if defined(__parallel)
  CALL update_structure(struc,'POS')
#endif

!
! get new forces
!
  CALL force ( struc, inter, thermo, simpar, ewald_param, box_change, &
               intenv, dft_control, kp )
  
!
! second half of velocity verlet
!
  DO i = 1, nnodes
     dtom = 0.5_dbl*simpar%dt/struc%pnode(i) %p%prop%mass
     struc%pnode(i) %p%v(1) = struc%pnode(i) %p%v(1) + &
          dtom*struc%pnode(i) %p%f(1)
     struc%pnode(i) %p%v(2) = struc%pnode(i) %p%v(2) + &
          dtom*struc%pnode(i) %p%f(2)
     struc%pnode(i) %p%v(3) = struc%pnode(i) %p%v(3) + &
          dtom*struc%pnode(i) %p%f(3)
  END DO

  DO i = 1, ncoef
     dtocm = 0.5_dbl * simpar % dt / struc % coef_vel ( 1 ) % ao % mass ( i )
     struc % coef_vel ( 1 ) % ao % cr ( i ) = struc % coef_vel ( 1 ) % ao % cr ( i ) &
                      + dtocm * struc % coef_force ( 1 ) % ao % cr ( i )
  END DO

  IF ( simpar % pol_ensemble == "NVT" ) & 
  CALL lnhc ( nhc_coef, struc % coef_vel ( 1 ) % ao, int_group )

  IF (simpar%constraint) CALL rattle_control(struc%molecule,simpar)
#if defined(__parallel)
  CALL update_structure(struc,'VEL')
#endif

! updating the virial
  IF (simpar%constraint) CALL pv_constraint(struc%molecule,thermo, &
                                            'FULL',int_group )

  CALL pressure(struc%pnode,thermo)

END SUBROUTINE nve_ao
!!*****
!******************************************************************************
!!****** integrator/nvt_ao [1.0] *
!!
!!   NAME
!!     nvt_pw
!!
!!   FUNCTION
!!     nvt integrator for particle positions & momenta
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE nvt_pw ( itimes, constant, simpar, inter, thermo, struc, &
                    ewald_param, nhc_part, nhc_coef, dft_control, kp )

  IMPLICIT NONE
  
! Arguments
  INTEGER, INTENT ( IN ) :: itimes
  REAL ( dbl ), INTENT ( INOUT ) :: constant
  TYPE ( simulation_parameters_type ), INTENT ( IN ) :: simpar
  TYPE ( interaction_type ), INTENT ( IN ) :: inter
  TYPE ( thermodynamic_type ), INTENT ( INOUT ) :: thermo
  TYPE ( structure_type ), INTENT ( INOUT ) :: struc
  TYPE ( ewald_parameters_type ), INTENT ( INOUT ) :: ewald_param
  TYPE ( lnhc_parameters_type ), INTENT ( INOUT ) :: nhc_part
  TYPE ( lnhc_parameters_type ), INTENT ( INOUT ) :: nhc_coef
  TYPE ( dft_control_type ) :: dft_control
  TYPE ( kpoint_type ), INTENT ( IN ) :: kp
  
! Local
  INTEGER :: i, nnodes, ncoef
  REAL ( dbl ) :: dtom,dtocm,cmass
  LOGICAL, PARAMETER :: box_change = .FALSE.

!------------------------------------------------------------------------------
  IF ( simpar % pol_ensemble == "NVT" ) & 
  CALL lnhc ( nhc_coef, struc % coef_vel ( 1 ) % pw, int_group )

  CALL lnhc ( nhc_part, struc % pnode, int_group )

  IF ( simpar % constraint ) CALL getold ( struc % molecule )
!
! first half of velocity verlet
!
  nnodes = size(struc%pnode)
  DO i = 1, nnodes
     dtom = 0.5_dbl*simpar%dt/struc%pnode(i) %p%prop%mass
     struc%pnode(i) %p%v(1) = struc%pnode(i) %p%v(1) + &
          dtom*struc%pnode(i) %p%f(1)
     struc%pnode(i) %p%v(2) = struc%pnode(i) %p%v(2) + &
          dtom*struc%pnode(i) %p%f(2)
     struc%pnode(i) %p%v(3) = struc%pnode(i) %p%v(3) + &
          dtom*struc%pnode(i) %p%f(3)
     struc%pnode(i) %p%r(1) = struc%pnode(i) %p%r(1) + &
          struc%pnode(i) %p%v(1)*simpar%dt
     struc%pnode(i) %p%r(2) = struc%pnode(i) %p%r(2) + &
          struc%pnode(i) %p%v(2)*simpar%dt
     struc%pnode(i) %p%r(3) = struc%pnode(i) %p%r(3) + &
          struc%pnode(i) %p%v(3)*simpar%dt
  END DO

  IF (simpar%constraint) CALL shake_control(struc%molecule,simpar)

  ncoef = SIZE ( struc % coef_pos (1) % pw % cc )
  DO i = 1, ncoef
     cmass = REAL ( struc % coef_vel ( 1 ) % pw % mass_cc ( i ), dbl )
     dtocm = 0.5_dbl * simpar % dt / cmass
     struc % coef_vel ( 1 ) % pw % cc ( i ) = struc % coef_vel ( 1 ) % pw  % cc ( i ) &
                      + dtocm * struc % coef_force ( 1 ) % pw % cc ( i )
     struc % coef_pos ( 1 ) % pw % cc ( i ) = struc % coef_pos ( 1 ) % pw % cc ( i ) &
                      + struc % coef_vel ( 1 ) % pw % cc ( i ) * simpar % dt
  END DO

#if defined(__parallel)
  CALL update_structure(struc,'POS')
#endif

!
! get new forces
!
  CALL force ( struc, inter, thermo, simpar, ewald_param, box_change, &
               intenv, dft_control, kp )
  
!
! second half of velocity verlet
!
  DO i = 1, nnodes
     dtom = 0.5_dbl*simpar%dt/struc%pnode(i) %p%prop%mass
     struc%pnode(i) %p%v(1) = struc%pnode(i) %p%v(1) + &
          dtom*struc%pnode(i) %p%f(1)
     struc%pnode(i) %p%v(2) = struc%pnode(i) %p%v(2) + &
          dtom*struc%pnode(i) %p%f(2)
     struc%pnode(i) %p%v(3) = struc%pnode(i) %p%v(3) + &
          dtom*struc%pnode(i) %p%f(3)
  END DO

  DO i = 1, ncoef
     cmass = REAL ( struc % coef_vel ( 1 ) % pw % mass_cc ( i ), dbl )
     dtocm = 0.5_dbl * simpar % dt / cmass
     struc % coef_vel ( 1 ) % pw % cc ( i ) = struc % coef_vel ( 1 ) % pw % cc ( i ) &
                      + dtocm * struc % coef_force ( 1 ) % pw % cc ( i )
  END DO

  IF ( simpar % pol_ensemble == "NVT" ) & 
  CALL lnhc ( nhc_coef, struc % coef_vel ( 1 ) % pw, int_group )

  CALL lnhc ( nhc_part, struc % pnode, int_group )

  IF (simpar%constraint) CALL rattle_control(struc%molecule,simpar)
#if defined(__parallel)
  CALL update_structure(struc,'VEL')
#endif

! updating the virial
  IF (simpar%constraint) CALL pv_constraint(struc%molecule,thermo, &
                                            'FULL',int_group )

  CALL pressure(struc%pnode,thermo)

END SUBROUTINE nvt_pw

!!*****
!******************************************************************************
!!****** integrator/nvt_ao [1.0] *
!!
!!   NAME
!!     nvt_ao
!!
!!   FUNCTION
!!     nvt integrator for particle positions & momenta
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE nvt_ao ( itimes, constant, simpar, inter, thermo, struc, &
                    ewald_param, nhc_part, nhc_coef, dft_control, kp )

  IMPLICIT NONE
  
! Arguments
  INTEGER, INTENT ( IN ) :: itimes
  REAL ( dbl ), INTENT ( INOUT ) :: constant
  TYPE ( simulation_parameters_type ), INTENT ( IN ) :: simpar
  TYPE ( interaction_type ), INTENT ( IN ) :: inter
  TYPE ( thermodynamic_type ), INTENT ( INOUT ) :: thermo
  TYPE ( structure_type ), INTENT ( INOUT ) :: struc
  TYPE ( ewald_parameters_type ), INTENT ( INOUT ) :: ewald_param
  TYPE ( lnhc_parameters_type ), INTENT ( INOUT ) :: nhc_part
  TYPE ( lnhc_parameters_type ), INTENT ( INOUT ) :: nhc_coef
  TYPE ( dft_control_type ) :: dft_control
  TYPE ( kpoint_type ), INTENT ( IN ) :: kp
  
! Local
  INTEGER :: i, nnodes, ncoef
  REAL ( dbl ) :: dtom, dtocm
  LOGICAL, PARAMETER :: box_change = .FALSE.

!------------------------------------------------------------------------------
  IF ( simpar % pol_ensemble == "NVT" ) & 
  CALL lnhc ( nhc_coef, struc % coef_vel ( 1 ) % ao, int_group )

  CALL lnhc ( nhc_part, struc % pnode, int_group )

  IF ( simpar % constraint ) CALL getold ( struc % molecule )
!
! first half of velocity verlet
!
  nnodes = size(struc%pnode)
  DO i = 1, nnodes
     dtom = 0.5_dbl*simpar%dt/struc%pnode(i) %p%prop%mass
     struc%pnode(i) %p%v(1) = struc%pnode(i) %p%v(1) + &
          dtom*struc%pnode(i) %p%f(1)
     struc%pnode(i) %p%v(2) = struc%pnode(i) %p%v(2) + &
          dtom*struc%pnode(i) %p%f(2)
     struc%pnode(i) %p%v(3) = struc%pnode(i) %p%v(3) + &
          dtom*struc%pnode(i) %p%f(3)
     struc%pnode(i) %p%r(1) = struc%pnode(i) %p%r(1) + &
          struc%pnode(i) %p%v(1)*simpar%dt
     struc%pnode(i) %p%r(2) = struc%pnode(i) %p%r(2) + &
          struc%pnode(i) %p%v(2)*simpar%dt
     struc%pnode(i) %p%r(3) = struc%pnode(i) %p%r(3) + &
          struc%pnode(i) %p%v(3)*simpar%dt
  END DO

  IF (simpar%constraint) CALL shake_control(struc%molecule,simpar)

  ncoef = SIZE ( struc % coef_pos (1) % ao % cr )
  DO i = 1, ncoef
     dtocm = 0.5_dbl * simpar % dt / struc % coef_vel ( 1 ) % ao % mass ( i )
     struc % coef_vel ( 1 ) % ao % cr ( i ) = struc % coef_vel ( 1 ) % ao  % cr ( i ) &
                      + dtocm * struc % coef_force ( 1 ) % ao % cr ( i )
     struc % coef_pos ( 1 ) % ao % cr ( i ) = struc % coef_pos ( 1 ) % ao % cr ( i ) &
                      + struc % coef_vel ( 1 ) % ao % cr ( i ) * simpar % dt
  END DO

#if defined(__parallel)
  CALL update_structure(struc,'POS')
#endif

!
! get new forces
!
  CALL force ( struc, inter, thermo, simpar, ewald_param, box_change, &
               intenv, dft_control, kp )
  
!
! second half of velocity verlet
!
  DO i = 1, nnodes
     dtom = 0.5_dbl*simpar%dt/struc%pnode(i) %p%prop%mass
     struc%pnode(i) %p%v(1) = struc%pnode(i) %p%v(1) + &
          dtom*struc%pnode(i) %p%f(1)
     struc%pnode(i) %p%v(2) = struc%pnode(i) %p%v(2) + &
          dtom*struc%pnode(i) %p%f(2)
     struc%pnode(i) %p%v(3) = struc%pnode(i) %p%v(3) + &
          dtom*struc%pnode(i) %p%f(3)
  END DO

  DO i = 1, ncoef
     dtocm = 0.5_dbl * simpar % dt / struc % coef_vel ( 1 ) % ao % mass ( i )
     struc % coef_vel ( 1 ) % ao % cr ( i ) = struc % coef_vel ( 1 ) % ao % cr ( i ) &
                      + dtocm * struc % coef_force ( 1 ) % ao % cr ( i )
  END DO

  IF ( simpar % pol_ensemble == "NVT" ) & 
  CALL lnhc ( nhc_coef, struc % coef_vel ( 1 ) % ao, int_group )

  CALL lnhc ( nhc_part, struc % pnode, int_group )

  IF (simpar%constraint) CALL rattle_control(struc%molecule,simpar)
#if defined(__parallel)
  CALL update_structure(struc,'VEL')
#endif

! updating the virial
  IF (simpar%constraint) CALL pv_constraint(struc%molecule,thermo, &
                                            'FULL',int_group )

  CALL pressure(struc%pnode,thermo)

END SUBROUTINE nvt_ao

!!*****
!******************************************************************************
!!****** integrator/nve [1.0] *
!!
!!   NAME
!!     nve
!!
!!   FUNCTION
!!     nve integrator for particle positions & momenta
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE nve ( itimes, constant, simpar, inter, thermo, struc, ewald_param, &
     dft_control, kp )
  
  IMPLICIT NONE
  
! Arguments
  INTEGER, INTENT ( IN ) :: itimes
  REAL ( dbl ), INTENT ( INOUT ) :: constant
  TYPE (simulation_parameters_type ), INTENT ( IN ) :: simpar
  TYPE (thermodynamic_type ), INTENT ( INOUT ) :: thermo
  TYPE (structure_type ), INTENT ( INOUT ) :: struc
  TYPE (ewald_parameters_type ), INTENT ( INOUT ) :: ewald_param
  TYPE (interaction_type ), INTENT ( IN ) :: inter
  TYPE ( dft_control_type ) :: dft_control
  TYPE ( kpoint_type ), INTENT ( IN ) :: kp
  
! Local
  INTEGER :: i, nnodes
  REAL ( dbl ) :: dtom
  LOGICAL, PARAMETER :: box_change = .FALSE.

!------------------------------------------------------------------------------

  IF ( simpar % constraint ) CALL getold ( struc % molecule )

!
! first half of velocity verlet
!
  nnodes = size(struc%pnode)
  DO i = 1, nnodes
     dtom = 0.5_dbl*simpar%dt/struc%pnode(i) %p%prop%mass
     struc%pnode(i) %p%v(1) = struc%pnode(i) %p%v(1) + &
          dtom*struc%pnode(i) %p%f(1)
     struc%pnode(i) %p%v(2) = struc%pnode(i) %p%v(2) + &
          dtom*struc%pnode(i) %p%f(2)
     struc%pnode(i) %p%v(3) = struc%pnode(i) %p%v(3) + &
          dtom*struc%pnode(i) %p%f(3)
     struc%pnode(i) %p%r(1) = struc%pnode(i) %p%r(1) + &
          struc%pnode(i) %p%v(1)*simpar%dt
     struc%pnode(i) %p%r(2) = struc%pnode(i) %p%r(2) + &
          struc%pnode(i) %p%v(2)*simpar%dt
     struc%pnode(i) %p%r(3) = struc%pnode(i) %p%r(3) + &
          struc%pnode(i) %p%v(3)*simpar%dt
  END DO

  IF (simpar%constraint) CALL shake_control(struc%molecule,simpar)
#if defined(__parallel)
  CALL update_structure(struc,'POS')
#endif

!
! get new forces
!
  CALL force ( struc, inter, thermo, simpar, ewald_param, box_change, &
               intenv, dft_control, kp )
  
!
! second half of velocity verlet
!
  DO i = 1, nnodes
     dtom = 0.5_dbl*simpar%dt/struc%pnode(i) %p%prop%mass
     struc%pnode(i) %p%v(1) = struc%pnode(i) %p%v(1) + &
          dtom*struc%pnode(i) %p%f(1)
     struc%pnode(i) %p%v(2) = struc%pnode(i) %p%v(2) + &
          dtom*struc%pnode(i) %p%f(2)
     struc%pnode(i) %p%v(3) = struc%pnode(i) %p%v(3) + &
          dtom*struc%pnode(i) %p%f(3)
  END DO

  IF (simpar%constraint) CALL rattle_control(struc%molecule,simpar)
#if defined(__parallel)
  CALL update_structure(struc,'VEL')
#endif

! updating the virial
  IF (simpar%constraint) CALL pv_constraint(struc%molecule,thermo, &
                                            'FULL',int_group )

  CALL pressure(struc%pnode,thermo)

END SUBROUTINE nve

!!*****
!******************************************************************************
!!****** integrator/nvt [1.0] *
!!
!!   NAME
!!     nvt
!!
!!   FUNCTION
!!     nvt integrator for particle positions & momenta
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE nvt ( itimes, constant, simpar, inter, thermo, struc, &
     ewald_param, nhc, dft_control, kp )
  
  IMPLICIT NONE

! Arguments
  INTEGER, INTENT ( IN ) :: itimes
  REAL ( dbl ), INTENT ( INOUT ) :: constant
  TYPE (simulation_parameters_type ), INTENT ( IN ) :: simpar
  TYPE (thermodynamic_type ), INTENT ( INOUT ) :: thermo
  TYPE (structure_type ), INTENT ( INOUT ) :: struc
  TYPE (lnhc_parameters_type ), INTENT ( INOUT ) :: nhc
  TYPE (ewald_parameters_type ), INTENT ( INOUT ) :: ewald_param
  TYPE (interaction_type ), INTENT ( IN ) :: inter
  TYPE ( dft_control_type ) :: dft_control
  TYPE ( kpoint_type ), INTENT ( IN ) :: kp
  
! Locals
  INTEGER :: i, nnodes
  REAL ( dbl ) :: dtom
  LOGICAL, PARAMETER :: box_change = .FALSE.

!------------------------------------------------------------------------------

  CALL lnhc ( nhc, struc % pnode, int_group )

  IF ( simpar % constraint ) CALL getold ( struc % molecule )

!
! first half of velocity verlet
!
  nnodes = SIZE ( struc % pnode )
  DO i = 1, nnodes
     dtom = 0.5_dbl * simpar%dt / struc%pnode(i)%p%prop%mass
     struc%pnode(i) %p%v(1) = struc%pnode(i) %p%v(1) + &
          dtom*struc%pnode(i) %p%f(1)
     struc%pnode(i) %p%v(2) = struc%pnode(i) %p%v(2) + &
          dtom*struc%pnode(i) %p%f(2)
     struc%pnode(i) %p%v(3) = struc%pnode(i) %p%v(3) + &
          dtom*struc%pnode(i) %p%f(3)
     struc%pnode(i) %p%r(1) = struc%pnode(i) %p%r(1) + &
          struc%pnode(i) %p%v(1)*simpar%dt
     struc%pnode(i) %p%r(2) = struc%pnode(i) %p%r(2) + &
          struc%pnode(i) %p%v(2)*simpar%dt
     struc%pnode(i) %p%r(3) = struc%pnode(i) %p%r(3) + &
          struc%pnode(i) %p%v(3)*simpar%dt
  END DO

  IF ( simpar % constraint ) CALL shake_control ( struc % molecule, simpar )
#if defined(__parallel)
  CALL update_structure ( struc, 'POS' )
#endif

!
! get new forces
!
  CALL force ( struc, inter, thermo, simpar, ewald_param, box_change, &
               intenv, dft_control, kp )
!
! second half of velocity verlet
!
  DO i = 1, nnodes
     dtom = 0.5_dbl*simpar%dt/struc%pnode(i) %p%prop%mass
     struc%pnode(i) %p%v(1) = struc%pnode(i) %p%v(1) + &
          dtom*struc%pnode(i) %p%f(1)
     struc%pnode(i) %p%v(2) = struc%pnode(i) %p%v(2) + &
          dtom*struc%pnode(i) %p%f(2)
     struc%pnode(i) %p%v(3) = struc%pnode(i) %p%v(3) + &
          dtom*struc%pnode(i) %p%f(3)
  END DO

  IF ( simpar % constraint ) CALL rattle_control ( struc % molecule, simpar )

  CALL lnhc ( nhc, struc % pnode, int_group )

#if defined(__parallel)
  CALL update_structure(struc,'VEL')
#endif

! updating the virial
  IF (simpar%constraint) CALL pv_constraint(struc%molecule,thermo, &
                                            'FULL',int_group )
  CALL pressure(struc%pnode,thermo)

END SUBROUTINE nvt
!!*****
!******************************************************************************
!!****** integrator/npt_i [1.0] *
!!
!!   NAME
!!     npt_i
!!
!!   FUNCTION
!!     npt_i integrator for particle positions & momenta
!!     isotropic box changes
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE npt_i ( itimes, constant, simpar, inter, thermo, struc, &
     ewald_param, nhc_part, nhc_baro, npt, dft_control, kp )
  
  IMPLICIT NONE

! Arguments
  INTEGER, INTENT ( IN ) :: itimes
  REAL ( dbl ), INTENT ( INOUT ) :: constant
  TYPE (simulation_parameters_type ), INTENT ( IN ) :: simpar
  TYPE (thermodynamic_type ), INTENT ( INOUT ) :: thermo
  TYPE (structure_type ), INTENT ( INOUT ) :: struc
  TYPE (lnhc_parameters_type ), INTENT ( INOUT ) :: nhc_part, nhc_baro
  TYPE (npt_info_type ), INTENT ( INOUT ), DIMENSION ( :, : ) :: npt
  TYPE (ewald_parameters_type ), INTENT ( INOUT ) :: ewald_param
  TYPE (interaction_type ), INTENT ( IN ) :: inter
  TYPE ( dft_control_type ) :: dft_control
  TYPE ( kpoint_type ), INTENT ( IN ) :: kp
  
! Locals
  INTEGER :: i, ii,  nnodes, iroll
  REAL ( dbl ) :: dtom
  REAL ( dbl ), PARAMETER :: e2 = 1._dbl/6._dbl, e4 = e2/20._dbl, &
       e6 = e4/42._dbl, e8 = e6/72._dbl
  REAL ( dbl ), SAVE :: eps_0
  REAL ( dbl ) :: arg_r, poly_r, scale_r, roll_tol
  REAL ( dbl ) :: arg_v, poly_v, scale_v, kin, fdotr, infree
  REAL ( dbl ), DIMENSION (3) :: vector_r, vector_v
  LOGICAL, PARAMETER :: box_change = .TRUE.
  LOGICAL :: first

!------------------------------------------------------------------------------
  infree = 1.0_dbl / REAL ( simpar % nfree, dbl ) 



! initialize eps_0 the first time through
  IF ( itimes == 1 ) eps_0 = npt(1,1) % eps


! Nose-Hoover for the Barostat
     CALL lnhc ( nhc_baro, npt, int_group )
! Nose-Hoover for the Particles
     CALL lnhc ( nhc_part, struc%pnode, int_group )

! setting up for ROLL: saving old variables
  roll_tol = -1.0E+10_dbl
  IF (simpar%constraint) THEN
     iroll = 1
     CALL set ( struc, npt, 'F' )
     CALL getold ( struc % molecule )
  END IF

  SR: DO WHILE (abs(roll_tol)>=1.E-10_dbl) ! SHAKE-ROLL LOOP

     IF (simpar%constraint) THEN
        CALL set( struc, npt, 'B' )
     END IF
     CALL get_veps ( npt, simpar, struc, thermo, infree )

!
! first half of velocity verlet
!
     arg_r = ( 0.5_dbl * npt ( 1, 1 ) % v * simpar % dt ) * &
             ( 0.5_dbl * npt ( 1, 1 ) % v * simpar % dt)
     poly_r = 1._dbl + e2*arg_r + e4*arg_r*arg_r + e6*arg_r**3 + e8*arg_r**4

     arg_v = ( 0.25_dbl * npt ( 1, 1 ) % v * simpar % dt *  &
             (1._dbl + 3._dbl * infree ) )*( 0.25_dbl * npt ( 1, 1 ) % v * &
              simpar % dt * ( 1._dbl + 3._dbl * infree ) )
     poly_v = 1._dbl + e2*arg_v + e4*arg_v*arg_v + e6*arg_v**3 + e8*arg_v**4

     scale_r = EXP( 0.5_dbl * simpar % dt * npt ( 1, 1 ) %v )
     scale_v = EXP( -0.25_dbl * simpar % dt * npt ( 1, 1 ) % v * &
                  (1._dbl + 3._dbl * infree ) )
     nnodes = size(struc%pnode)
     DO i = 1, nnodes
        dtom = 0.5_dbl*simpar%dt/struc%pnode(i) %p%prop%mass

        struc%pnode(i) %p%v(1) = struc%pnode(i) %p%v(1)*scale_v*scale_v + &
             scale_v*poly_v*struc%pnode(i) %p%f(1)* dtom
        struc%pnode(i) %p%v(2) = struc%pnode(i) %p%v(2)*scale_v*scale_v + &
             scale_v*poly_v*struc%pnode(i) %p%f(2)* dtom
        struc%pnode(i) %p%v(3) = struc%pnode(i) %p%v(3)*scale_v*scale_v + &
             scale_v*poly_v*struc%pnode(i) %p%f(3)* dtom

        struc%pnode(i) %p%r(1) = struc%pnode(i) %p%r(1)*scale_r*scale_r + &
             scale_r*poly_r*struc%pnode(i) %p%v(1)*simpar%dt
        struc%pnode(i) %p%r(2) = struc%pnode(i) %p%r(2)*scale_r*scale_r + &
             scale_r*poly_r*struc%pnode(i) %p%v(2)*simpar%dt
        struc%pnode(i) %p%r(3) = struc%pnode(i) %p%r(3)*scale_r*scale_r + &
             scale_r*poly_r*struc%pnode(i) %p%v(3)*simpar%dt
     END DO
     roll_tol = 0._dbl
     vector_r ( : ) = scale_r * poly_r
     vector_v ( : ) = scale_v * poly_v

     IF ( simpar%constraint ) CALL shake_roll_control(struc%molecule,simpar, &
          thermo,roll_tol,iroll,vector_r,vector_v,int_group )
  END DO SR

! Update eps:
  npt ( :, : ) % eps =  npt ( :, : ) % eps &
                       + simpar % dt * npt ( :, : ) % v

! Update h_mat
  struc % box % hmat ( :, : ) = struc % box % hmat ( :, : ) &
       * EXP ( npt ( 1, 1 ) % eps - eps_0 )

  eps_0 =  npt ( 1, 1 ) % eps

#if defined(__parallel)
  CALL update_structure(struc,'POS')
#endif

!
! get new forces
!
  CALL force ( struc, inter, thermo, simpar, ewald_param, box_change, &
               intenv, dft_control, kp )
  
!
! second half of velocity verlet
!
  DO i = 1, nnodes
     dtom = 0.5_dbl*simpar%dt/struc%pnode(i) %p%prop%mass
     struc%pnode(i) %p%v(1) = struc%pnode(i) %p%v(1)*scale_v*scale_v + &
          scale_v*poly_v*struc%pnode(i) %p%f(1)* dtom
     struc%pnode(i) %p%v(2) = struc%pnode(i) %p%v(2)*scale_v*scale_v + &
          scale_v*poly_v*struc%pnode(i) %p%f(2)* dtom
     struc%pnode(i) %p%v(3) = struc%pnode(i) %p%v(3)*scale_v*scale_v + &
          scale_v*poly_v*struc%pnode(i) %p%f(3)* dtom
  END DO


  roll_tol = -1.E10_dbl
  IF (simpar%constraint) THEN
     first = .TRUE.
     iroll = 1
     CALL set( struc, npt, 'F' )
  END IF

  RR: DO WHILE (abs(roll_tol)>=1.E-10_dbl) ! RATTLE-ROLL LOOP
     roll_tol = 0._dbl
     IF (simpar%constraint) CALL rattle_roll_setup ( npt, struc, &
            thermo, simpar, vector_v, roll_tol, iroll, infree, first )

     CALL get_veps ( npt, simpar, struc, thermo, infree )


  END DO RR


! Nose-Hoover for the Particles
  CALL lnhc ( nhc_part, struc%pnode, int_group )
! Nose-Hoover for the Barostat
  CALL lnhc ( nhc_baro, npt, int_group )


#if defined(__parallel)
  CALL update_structure(struc,'VEL')
#endif

  CALL pressure(struc%pnode,thermo)

END SUBROUTINE npt_i

!!*****
!******************************************************************************
!!****** integrator/pressure [1.0] *
!!
!!   NAME
!!     pressure
!!
!!   FUNCTION
!!     Computes the kinetic part of the pressure tensor and updates
!!     the full VIRIAL (PV)
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE pressure ( pnode, thermo )

  IMPLICIT NONE

! Arguments
  TYPE ( particle_node_type ), DIMENSION ( : ), INTENT ( IN ) :: pnode
  TYPE ( thermodynamic_type ), INTENT ( INOUT ) :: thermo

! Locals
  INTEGER :: i, j, iatom, nnodes, handle
  REAL ( dbl ) :: mfl

!------------------------------------------------------------------------------

  CALL timeset ( 'PRESSURE', 'E', 'Mflops', handle )

  thermo % pv_kin = 0._dbl
  nnodes = size(pnode)
  DO i = 1, 3
     DO j = 1, 3
        DO iatom = 1, nnodes
           thermo % pv_kin ( i, j ) = thermo % pv_kin ( i, j ) + &
                pnode ( iatom ) % p % prop % mass * &
                pnode ( iatom ) % p % v ( i )* &
                pnode ( iatom ) % p % v ( j )
        END DO
     END DO
  END DO
  mfl = REAL( 9 * nnodes, dbl ) * 2._dbl * 1.e-6_dbl

#if defined(__parallel)
  CALL mp_sum(thermo%pv_kin,int_group)
#endif

! total virial
  thermo%ptens = thermo%pv + thermo%pv_kin + thermo%pv_const

  CALL timestop ( mfl, handle )

END SUBROUTINE pressure

!!*****
!******************************************************************************
!!****** integrator/update_structure [1.0] *
!!
!!   NAME
!!     update_structure
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

SUBROUTINE update_structure(struc,task)

  IMPLICIT NONE

! Arguments
  TYPE (structure_type ), INTENT ( INOUT ) :: struc
  CHARACTER ( LEN = * ), INTENT ( IN ) :: task

! Locals
  REAL ( dbl ), ALLOCATABLE, DIMENSION ( :, : ) :: atot
  INTEGER :: natoms, ios, imol, atombase, iat, i, handle
  TYPE (linklist_atoms), POINTER :: current_atom

!------------------------------------------------------------------------------

  CALL timeset ( 'UPDATE_STRUCTURE', 'E', ' ', handle )

  natoms = SIZE ( struc % part )
  ALLOCATE (atot(3,natoms),STAT=ios)
  IF ( ios /= 0 ) CALL stop_memory ( 'update_structure', 'atot', 3 * natoms )
  atot = zero
  DO imol = 1, size(struc%molecule)
     current_atom => struc%molecule(imol) %ll_atoms
     DO i = 1, struc%molecule(imol) %natoms_mol
        iat = struc%molecule(imol) %atombase + i - 1
        IF (task=='POS') THEN
           atot(1,iat) = current_atom%part%p%r(1)
           atot(2,iat) = current_atom%part%p%r(2)
           atot(3,iat) = current_atom%part%p%r(3)
        ELSE IF (task=='VEL') THEN
           atot(1,iat) = current_atom%part%p%v(1)
           atot(2,iat) = current_atom%part%p%v(2)
           atot(3,iat) = current_atom%part%p%v(3)
        END IF
        current_atom => current_atom%next
     END DO
  END DO

#ifdef __parallel
  CALL mp_sum(atot,int_group)
#endif

  IF (task=='POS') THEN
     DO iat = 1, natoms
        struc%part(iat) %r(1) = atot(1,iat)
        struc%part(iat) %r(2) = atot(2,iat)
        struc%part(iat) %r(3) = atot(3,iat)
     END DO
  ELSE IF (task=='VEL') THEN
     DO iat = 1, natoms
        struc%part(iat) %v(1) = atot(1,iat)
        struc%part(iat) %v(2) = atot(2,iat)
        struc%part(iat) %v(3) = atot(3,iat)
     END DO
  END IF

! DEALLOCATE ( atot, STAT = ios )
! IF ( ios /= 0 ) CALL stop_memory ( 'update_structure', 'atot' )

  CALL timestop ( zero, handle )

END SUBROUTINE update_structure

!!*****
!******************************************************************************
!!****** integrator/set [1.0] *
!!
!!   NAME
!!     set
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

SUBROUTINE set(struc,npt,char)

  IMPLICIT NONE

! Arguments
  TYPE (structure_type ), INTENT ( INOUT ) :: struc
  TYPE (npt_info_type ), DIMENSION ( :, : ), INTENT ( INOUT ) :: npt
  CHARACTER ( LEN = * ), INTENT ( IN ) :: char

! Locals
  INTEGER :: idim
  INTEGER :: isos

!------------------------------------------------------------------------------

  isos = 0
  IF ( .NOT. ASSOCIATED ( old % v ) ) &
    ALLOCATE ( old % v ( SIZE ( struc % part ), 3 ), STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'set', 'old%v', 0 )
  IF ( .NOT. ASSOCIATED ( old % r ) ) &
    ALLOCATE ( old % r ( SIZE ( struc % part ), 3 ), STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'set', 'old%r', 0 )
  IF ( .NOT. ASSOCIATED ( old % eps ) )  &
    ALLOCATE ( old % eps ( SIZE ( npt, 1), SIZE ( npt, 2 ) ),STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'set', 'old%eps', 0 )
  IF ( .NOT. ASSOCIATED ( old % veps ) ) &
    ALLOCATE ( old % veps ( SIZE ( npt, 1), SIZE ( npt, 2 ) ),STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'set', 'old%veps', 0 )
  IF ( .NOT. ASSOCIATED ( old % h ) ) &
    ALLOCATE (old%h(SIZE(struc%box%hmat,1), SIZE(struc%box%hmat,2)),STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'set', 'old%h', 0 )

  SELECT CASE (char)
  CASE ( 'F') ! forward assigning the old
     DO idim = 1, 3
        old % v ( :, idim ) = struc % part ( : ) % v ( idim )
        old % r ( :, idim ) = struc % part ( : ) % r ( idim )
     END DO
     old % eps ( :, : ) = npt ( :, : ) % eps
     old % veps ( :, : ) = npt ( :, : ) % v
     old % h ( :, : ) = struc%box%hmat ( :, : )
  CASE ( 'B') ! back assigning the original variables
     DO idim = 1, 3
        struc % part ( : ) % r ( idim ) = old % r ( :, idim )
        struc % part ( : ) % v ( idim ) = old % v ( :, idim )
     END DO
     npt ( :, : ) % eps = old % eps ( :, : )
     npt ( :, : ) % v = old % veps ( :, : )
     struc % box % hmat ( :, : )=old % h ( :, : )
  END SELECT

END SUBROUTINE set

!!*****
!******************************************************************************
!!****** integrator/npt_f [1.0] *
!!
!!   NAME
!!     npt_f
!!
!!   FUNCTION
!!     Velocity Verlet integrator for the NPT ensemble with fully flexible cell
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE npt_f ( itimes, constant, simpar, inter, thermo, struc, &
     ewald_param, nhc_part, nhc_baro, npt, dft_control, kp )
  
  IMPLICIT NONE

! Arguments
  INTEGER, INTENT ( IN ) :: itimes
  REAL ( dbl ), INTENT ( INOUT ) :: constant
  TYPE ( simulation_parameters_type ), INTENT ( IN ) :: simpar
  TYPE ( thermodynamic_type ), INTENT ( INOUT ) :: thermo
  TYPE ( structure_type ), INTENT ( INOUT ) :: struc
  TYPE ( lnhc_parameters_type ), INTENT ( INOUT ) :: nhc_part, nhc_baro
  TYPE ( npt_info_type ), INTENT ( INOUT ), DIMENSION ( :, : ) :: npt
  TYPE ( ewald_parameters_type ), INTENT ( INOUT ) :: ewald_param
  TYPE ( interaction_type ), INTENT ( IN ) :: inter
  TYPE ( dft_control_type ) :: dft_control
  TYPE ( kpoint_type ), INTENT ( IN ) :: kp
  
! Locals
  INTEGER :: i, j, nnodes, iroll
  REAL ( dbl ) :: dtom
  REAL ( dbl ), PARAMETER :: e2 = 1._dbl/6._dbl, e4 = e2/20._dbl, &
       e6 = e4/42._dbl, e8 = e6/72._dbl
  REAL ( dbl ), SAVE :: eps_0
  REAL ( dbl ) :: roll_tol, infree, trvg
  REAL ( dbl ), DIMENSION (3,3) :: pv_kinetic, u, uh
  REAL ( dbl ), DIMENSION (3) :: arg_r, poly_r, scale_r, e_val, ur, uv, uf
  REAL ( dbl ), DIMENSION (3) :: arg_v, poly_v, scale_v
  LOGICAL, PARAMETER :: box_change = .TRUE.
  LOGICAL :: first

!------------------------------------------------------------------------------
  infree = 1.0_dbl / REAL ( simpar % nfree, dbl ) 
  


! Nose-Hoover for the Barostat
     CALL lnhc ( nhc_baro, npt, int_group )
! Nose-Hoover for the Particles
     CALL lnhc ( nhc_part, struc%pnode, int_group )

! setting up for ROLL: saving old variables
  roll_tol = -1.0E+10_dbl
  IF (simpar%constraint) THEN
     iroll = 1
     CALL set(struc, npt, 'F' )
     CALL getold ( struc % molecule )
  END IF
  SR: DO WHILE (abs(roll_tol)>=1.E-10_dbl) ! SHAKE-ROLL LOOP

     IF (simpar%constraint) THEN
        CALL set ( struc, npt, 'B' )
     END IF

     CALL get_veps ( npt, simpar, struc, thermo, infree )

     trvg = npt ( 1, 1 ) % v  + npt ( 2, 2 ) % v + npt ( 3, 3 ) % v 
     
! find eigenvalues and eigenvectors of npt ( :, : ) % v
!

     CALL diagonalise ( matrix =  npt % v, mysize = 3, &
          storageform = "UPPER", eigenvalues = e_val, eigenvectors = u )
!
! first half of velocity verlet
!
     arg_r ( : ) =  0.5_dbl * e_val ( : ) * simpar % dt * &
                    0.5_dbl * e_val ( : ) * simpar % dt
     poly_r = 1._dbl + e2*arg_r + e4*arg_r*arg_r + e6*arg_r**3 + e8*arg_r**4
     scale_r ( : ) = EXP ( 0.5_dbl * simpar % dt * e_val ( : ) )

     arg_v ( : ) =  0.25_dbl * simpar % dt * ( e_val ( : ) + trvg * infree ) *&
                    0.25_dbl * simpar % dt * ( e_val ( : ) + trvg * infree )
     poly_v = 1._dbl + e2*arg_v + e4*arg_v*arg_v + e6*arg_v**3 + e8*arg_v**4
     scale_v ( : ) = EXP ( -0.25_dbl * simpar % dt * ( &
                           e_val ( : ) + trvg * infree ) )

     nnodes = size(struc%pnode)
     DO i = 1, nnodes
        dtom = 0.5_dbl*simpar%dt/struc%pnode(i) %p%prop%mass
! Transform positions and velocities and forces
        ur = MATVEC_3X3 ( TRANSPOSE_3D ( u ), struc % pnode ( i ) % p % r )
        uv = MATVEC_3X3 ( TRANSPOSE_3D ( u ), struc % pnode ( i ) % p % v )
        uf = MATVEC_3X3 ( TRANSPOSE_3D ( u ), struc % pnode ( i ) % p % f )

        uv(1) = uv(1)*scale_v(1)*scale_v(1) + uf(1)*scale_v(1)*poly_v(1)*dtom
        uv(2) = uv(2)*scale_v(2)*scale_v(2) + uf(2)*scale_v(2)*poly_v(2)*dtom
        uv(3) = uv(3)*scale_v(3)*scale_v(3) + uf(3)*scale_v(3)*poly_v(3)*dtom

        ur(1) = ur(1)*scale_r(1)*scale_r(1) +  &
                uv(1)*scale_r(1)*poly_r(1)*simpar%dt
        ur(2) = ur(2)*scale_r(2)*scale_r(2) +  &
                uv(2)*scale_r(2)*poly_r(2)*simpar%dt
        ur(3) = ur(3)*scale_r(3)*scale_r(3) +  &
                uv(3)*scale_r(3)*poly_r(3)*simpar%dt

        struc % pnode ( i ) % p % r ( : ) = MATVEC_3x3 ( u, ur )
        struc % pnode ( i ) % p % v ( : ) = MATVEC_3x3 ( u, uv )
     END DO

     roll_tol = 0.0_dbl
     IF ( simpar % constraint ) CALL shake_roll_control ( struc % molecule, &
          simpar, thermo,roll_tol,iroll,scale_r*poly_r,scale_v*poly_v,  &
          int_group ,u)
  END DO SR

! Update h_mat

  uh = MATMUL_3X3 ( TRANSPOSE_3D ( u ), struc % box % hmat )

  DO i = 1, 3
     DO j = 1, 3
       uh ( i, j ) = uh ( i, j ) * scale_r ( i ) * scale_r ( i )
     END DO
  END DO

  struc % box % hmat = matmul_3x3 ( u, uh )

#if defined(__parallel)
  CALL update_structure(struc,'POS')
#endif

!
! get new forces
!
  CALL force ( struc, inter, thermo, simpar, ewald_param, box_change, &
               intenv, dft_control, kp )
  
!
! second half of velocity verlet
!
  DO i = 1, nnodes
     dtom = 0.5_dbl*simpar%dt/struc%pnode(i) %p%prop%mass
     uv = MATVEC_3X3 ( TRANSPOSE_3D ( u ), struc % pnode ( i ) % p % v )
     uf = MATVEC_3x3 ( TRANSPOSE_3D ( u ), struc % pnode ( i ) % p % f )

     uv(1) = uv(1)*scale_v(1)*scale_v(1) + uf(1)*scale_v(1)*poly_v(1)*dtom
     uv(2) = uv(2)*scale_v(2)*scale_v(2) + uf(2)*scale_v(2)*poly_v(2)*dtom
     uv(3) = uv(3)*scale_v(3)*scale_v(3) + uf(3)*scale_v(3)*poly_v(3)*dtom

     struc%pnode(i) %p%v ( : ) = MATVEC_3x3(u,uv)
  END DO

  roll_tol = -1.0E+10_dbl
  IF (simpar%constraint) THEN
     first = .TRUE.
     iroll = 1
     CALL set ( struc, npt, 'F' )
  END IF

  RR: DO WHILE ( ABS ( roll_tol ) >= 1.0E-10_dbl ) ! RATTLE-ROLL LOOP
     roll_tol = 0.0_dbl
     IF ( simpar%constraint ) CALL rattle_roll_setup ( npt, struc, &
        thermo, simpar, scale_v * poly_v, roll_tol, iroll, infree, first )

     CALL get_veps ( npt, simpar, struc, thermo, infree )

  END DO RR

#if defined(__parallel)
  CALL update_structure(struc,'VEL')
#endif

! Nose-Hoover for the Particles
  CALL lnhc ( nhc_part, struc%pnode, int_group )
! Nose-Hoover for the Barostat
  CALL lnhc ( nhc_baro, npt, int_group )

  CALL pressure(struc%pnode,thermo)

END SUBROUTINE npt_f

!!*****
!******************************************************************************
!!****** integrator/rattle_roll_setup [1.0] *
!!
!!   NAME
!!     rattle_roll_setup
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

SUBROUTINE rattle_roll_setup ( npt, struc, thermo, simpar, vector_v,  &
                               roll_tol, iroll, infree, first )

  IMPLICIT NONE

! Arguments
  TYPE ( simulation_parameters_type ), INTENT ( IN ) :: simpar
  TYPE ( structure_type ), INTENT ( INOUT ) :: struc
  TYPE ( npt_info_type ), INTENT ( INOUT ), DIMENSION ( :, : ) :: npt
  TYPE ( thermodynamic_type ), INTENT ( INOUT ) :: thermo
  LOGICAL, INTENT ( INOUT ) :: first
  REAL ( dbl ), INTENT ( OUT ) :: roll_tol
  REAL ( dbl ), INTENT ( IN ) :: infree
  REAL ( dbl ), INTENT ( IN ), DIMENSION ( : ) :: vector_v
  INTEGER, INTENT ( INOUT ) :: iroll

! Locals
  TYPE ( npt_info_type ), DIMENSION ( 3, 3 ) :: npt_loc

!------------------------------------------------------------------------------

  IF ( first ) THEN
    CALL get_veps ( npt, simpar, struc, thermo, infree )
  END IF
  first = .FALSE.

! assigning local variable
  SELECT CASE ( simpar % ensemble )
  CASE ( 'NPT_I' )
    npt_loc ( :, : ) % v = 0._dbl
    npt_loc ( :, : ) % mass = 0._dbl
    npt_loc ( 1, 1 ) % v = npt ( 1, 1 ) % v
    npt_loc ( 2, 2 ) % v = npt ( 1, 1 ) % v
    npt_loc ( 3, 3 ) % v = npt ( 1, 1 ) % v
    npt_loc ( 1, 1 ) % mass = npt ( 1, 1 ) % mass
    npt_loc ( 2, 2 ) % mass = npt ( 1, 1 ) % mass
    npt_loc ( 3, 3 ) % mass = npt ( 1, 1 ) % mass
  CASE ( 'NPT_F' )
    npt_loc = npt
  END SELECT

! resetting

  CALL set ( struc, npt, 'B' )
  CALL rattle_roll_control( struc%molecule, simpar, thermo, npt_loc % v, &
                            npt_loc % mass, roll_tol, iroll, vector_v,  &
                            int_group )

END SUBROUTINE rattle_roll_setup

!!*****
!******************************************************************************
!!****** integrator/get_veps [1.0] *
!!
!!   NAME
!!     get_veps
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

SUBROUTINE get_veps ( npt, simpar, struc, thermo, infree )

! Arguments
  TYPE (simulation_parameters_type ), INTENT ( IN ) :: simpar
  TYPE (thermodynamic_type ), INTENT ( INOUT ) :: thermo
  TYPE (structure_type ), INTENT ( IN ) :: struc
  TYPE (npt_info_type ), INTENT ( INOUT ), DIMENSION ( :, : ) :: npt
  REAL ( dbl ), INTENT ( IN ) :: infree

! locals
  INTEGER :: i, j, iatom, ii, jj, nnodes
  REAL ( dbl ) :: kin, fdotr
  REAL ( dbl ), DIMENSION ( 3, 3 ) :: pv_kin, unit
! dbg ISOTROPIC LIMIT
!  REAL ( dbl ) :: trace
! dbg

! initializing locals
  unit = 0.0_dbl
  pv_kin = 0.0_dbl
  unit(1,1) = 1.0_dbl
  unit(2,2) = 1.0_dbl
  unit(3,3) = 1.0_dbl

  nnodes = SIZE ( struc % pnode )
  DO i = 1, 3
     DO j = 1, 3
        DO iatom = 1, nnodes
           pv_kin ( i, j ) = pv_kin ( i, j ) + &
                struc % pnode ( iatom ) % p % prop % mass * &
                struc % pnode ( iatom ) % p % v ( i )* &
                struc % pnode ( iatom ) % p % v ( j )
        END DO
     END DO
  END DO

  kin = 0.0_dbl
  DO ii = 1, 3
     kin = kin +  pv_kin ( ii, ii )
  END DO

#if defined(__parallel)
    CALL mp_sum( pv_kin, int_group )
#endif


! updating the constraint virial
  IF ( simpar % constraint )  &
    CALL pv_constraint ( struc % molecule, thermo, 'FULL', int_group )

  IF ( simpar % ensemble == 'NPT_I' ) THEN
! get force on barostat
     fdotr = 0.0_dbl
     DO ii = 1, 3
       fdotr = fdotr + thermo % pv ( ii, ii ) +  &
                       thermo % pv_const ( ii, ii )
     END DO

     npt ( :, : ) % f = ( 1.0_dbl + ( 3.0_dbl * infree ) ) * kin + fdotr - &
       3.0_dbl * simpar % p_ext * struc % box % deth
  ELSEIF ( simpar % ensemble == 'NPT_F' ) THEN
     npt ( :, : ) % f = thermo % pv ( :, : ) +                   &
       pv_kin ( :, : ) +  thermo % pv_const ( :, : ) -  &
       unit ( :, : ) * simpar % p_ext * struc % box % deth +       &
       infree * kin * unit ( :, : )
! dbg ISOTROPIC LIMIT
!    trace = npt ( 1, 1 ) % f + npt ( 2, 2 ) % f + npt ( 3, 3 ) % f
!    trace = trace / 3.0_dbl
!    npt ( :, : ) % f = trace * unit ( :, : )
! dbg
  ENDIF

! update barostat velocities
     npt ( :, : ) % v = npt ( :, : ) % v + &
     0.5_dbl * simpar % dt * npt ( :, : ) % f / npt ( :, : ) % mass

END SUBROUTINE get_veps
!!*****
!******************************************************************************
!!****** integrator/nvt_fist_pimd [1.0] *
!!
!!   NAME
!!     nvt_pimd
!!
!!   FUNCTION
!!     nvt_pimd integrator for particle positions & momenta
!!
!!   AUTHOR
!!     HAF
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE nvt_pimd(itimes,constant,simpar,inter,thermo,struc, &
     ewald_param,nhc,dft_control,kp)
  
  IMPLICIT NONE
  
! Arguments
  INTEGER, INTENT ( IN ) :: itimes
  REAL ( dbl ), INTENT ( INOUT ) :: constant
  TYPE (simulation_parameters_type ), INTENT ( IN ) :: simpar
  TYPE (thermodynamic_type ), INTENT ( INOUT ) :: thermo
  TYPE (structure_type ), DIMENSION ( : ), INTENT ( INOUT ) :: struc
  TYPE (lnhc_parameters_type ), DIMENSION ( : ), INTENT ( INOUT ) :: nhc
  TYPE (ewald_parameters_type ), INTENT ( INOUT ) :: ewald_param
  TYPE (interaction_type ), INTENT ( IN ) :: inter
  TYPE ( dft_control_type ) :: dft_control
  TYPE ( kpoint_type ), INTENT ( IN ) :: kp
  
! Locals
  INTEGER :: i, j, k, nnodes, beads
  REAL ( dbl ) :: dtom, transdt, old_u, fac
  LOGICAL, PARAMETER :: box_change=.false.
  
!------------------------------------------------------------------------------

  beads = simpar % pimd_params % beads
  
  DO i = 1, beads
     CALL lnhc ( nhc ( i ), struc ( i ) % pnode, int_group )
  END DO
  
!  IF ( simpar % constraint ) CALL getold ( struc % molecule )
  IF ( simpar % constraint ) CALL stop_program('pimd_nvt', &
                                     'contraints not supported atm with PIMD')
#if defined(__parallel)
  CALL stop_program('pimd_nvt','parallel computation not supported with PIMD atm')
#endif

!
! first half of velocity verlet
!
  fac = 0.5_dbl/beads
  nnodes = size(struc(1)%pnode)
  DO j=1, beads
     transdt = simpar%dt*struc(j)%pimd_data%ismeff
     DO i = 1, nnodes
        dtom = fac*transdt/struc(j)%pnode(i) %p%prop%mass
        struc(j)%pnode(i) %p%v(1) = struc(j)%pnode(i) %p%v(1) + &
             dtom*struc(j)%pnode(i) %p%f(1)
        struc(j)%pnode(i) %p%v(2) = struc(j)%pnode(i) %p%v(2) + &
             dtom*struc(j)%pnode(i) %p%f(2)
        struc(j)%pnode(i) %p%v(3) = struc(j)%pnode(i) %p%v(3) + &
             dtom*struc(j)%pnode(i) %p%f(3)
     END DO
  END DO
  
!  IF ( simpar % constraint ) CALL shake_control ( struc % molecule, simpar )
#if defined(__parallel)
!  CALL update_structure ( struc, 'POS' )
#endif

! Harmonic bead-bead forces, assuming diagonalizing transformation for now

  DO j=1, beads
    DO i=1, nnodes
      old_u = struc(j)%pimd_data%u(1,i)
      struc(j)%pimd_data%u(1,i)= old_u * struc(j)%pimd_data%harmonic_prop(1,1)&
            + struc(j)%pnode(i)%p%v(1) * struc(j)%pimd_data%harmonic_prop(1,2)
      struc(j)%pnode(i)%p%v(1) = old_u * struc(j)%pimd_data%harmonic_prop(2,1)&
            + struc(j)%pnode(i)%p%v(1) * struc(j)%pimd_data%harmonic_prop(2,2)

      old_u = struc(j)%pimd_data%u(2,i)
      struc(j)%pimd_data%u(2,i)= old_u * struc(j)%pimd_data%harmonic_prop(1,1)&
            + struc(j)%pnode(i)%p%v(2) * struc(j)%pimd_data%harmonic_prop(1,2)
      struc(j)%pnode(i)%p%v(2) = old_u * struc(j)%pimd_data%harmonic_prop(2,1)&
            + struc(j)%pnode(i)%p%v(2) * struc(j)%pimd_data%harmonic_prop(2,2)

      old_u = struc(j)%pimd_data%u(3,i)
      struc(j)%pimd_data%u(3,i)= old_u * struc(j)%pimd_data%harmonic_prop(1,1)&
            + struc(j)%pnode(i)%p%v(3) * struc(j)%pimd_data%harmonic_prop(1,2)
      struc(j)%pnode(i)%p%v(3) = old_u * struc(j)%pimd_data%harmonic_prop(2,1)&
            + struc(j)%pnode(i)%p%v(3) * struc(j)%pimd_data%harmonic_prop(2,2)
    END DO
  END DO
  IF (simpar%pimd_params%ndelta > 1) THEN
     ! double time step for thermostats
     DO i=1, beads
        CALL yoshida_coef(nhc(i),2.0_dbl*simpar % dt / &
                          simpar % pimd_params % ndelta)
     END DO
     DO k= 2, simpar%pimd_params%ndelta
        DO i = 1, beads
           CALL lnhc ( nhc(i), struc(i) % pnode, int_group )
           DO j=1, nnodes
              old_u = struc(i)%pimd_data%u(1,j)
              struc(i)%pimd_data%u(1,j)= old_u * &
                         struc(i)%pimd_data%harmonic_prop(1,1) &
                       + struc(i)%pnode(j)%p%v(1) * &
                         struc(i)%pimd_data%harmonic_prop(1,2)
              struc(i)%pnode(j)%p%v(1) = old_u * &
                         struc(i)%pimd_data%harmonic_prop(2,1) &
                       + struc(i)%pnode(j)%p%v(1) * &
                         struc(i)%pimd_data%harmonic_prop(2,2)

              old_u = struc(i)%pimd_data%u(2,j)
              struc(i)%pimd_data%u(2,j)= old_u * &
                         struc(i)%pimd_data%harmonic_prop(1,1) &
                       + struc(i)%pnode(j)%p%v(2) * &
                         struc(i)%pimd_data%harmonic_prop(1,2)
              struc(i)%pnode(j)%p%v(2) = old_u * &
                         struc(i)%pimd_data%harmonic_prop(2,1) &
                       + struc(i)%pnode(j)%p%v(2) * &
                         struc(i)%pimd_data%harmonic_prop(2,2)

              old_u = struc(i)%pimd_data%u(3,j)
              struc(i)%pimd_data%u(3,j)= old_u * &
                         struc(i)%pimd_data%harmonic_prop(1,1) &
                       + struc(i)%pnode(j)%p%v(3) * &
                         struc(i)%pimd_data%harmonic_prop(1,2)
              struc(i)%pnode(j)%p%v(3) = old_u * &
                         struc(i)%pimd_data%harmonic_prop(2,1) &
                       + struc(i)%pnode(j)%p%v(3) * &
                         struc(i)%pimd_data%harmonic_prop(2,2)
           END DO
        END DO
     END DO
     ! set time step for thermostats back to normal
     DO i=1, beads
        CALL yoshida_coef(nhc(i),simpar % dt / &
                          simpar % pimd_params % ndelta)
     END DO
  END IF

  CALL u2r (struc, simpar)
!
! get new forces
!
  CALL force ( struc, inter, thermo, simpar, ewald_param, box_change, &
                  intenv, dft_control, kp )

  CALL fr2fu (struc, simpar)
  
!
! second half of velocity verlet
!

  DO j=1, beads
     transdt = 0.5_dbl*simpar%dt*struc(j)%pimd_data%ismeff / beads
     DO i = 1, nnodes
        dtom = transdt/struc(j)%pnode(i) %p%prop%mass
        struc(j)%pnode(i) %p%v(1) = struc(j)%pnode(i) %p%v(1) + &
             dtom*struc(j)%pnode(i) %p%f(1)
        struc(j)%pnode(i) %p%v(2) = struc(j)%pnode(i) %p%v(2) + &
             dtom*struc(j)%pnode(i) %p%f(2)
        struc(j)%pnode(i) %p%v(3) = struc(j)%pnode(i) %p%v(3) + &
             dtom*struc(j)%pnode(i) %p%f(3)
     END DO
  END DO

!  IF ( simpar % constraint ) CALL rattle_control ( struc % molecule, simpar )
  
  DO i = 1, beads
     CALL lnhc ( nhc(i), struc(i) % pnode, int_group )
  END DO
  
#if defined(__parallel)
!  CALL update_structure(struc,'VEL')
#endif
  
! updating the virial
!  IF (simpar%constraint) CALL pv_constraint(struc%molecule,thermo, &
!                                            'FULL',int_group )
!+++ that is not right yet either:
  CALL pressure ( struc(1)%pnode,thermo)
  
END SUBROUTINE nvt_pimd

!!*****
!******************************************************************************

END MODULE integrator

!******************************************************************************
