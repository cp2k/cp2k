!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

MODULE integrator
  
!  USE constraint, ONLY : shake_control, rattle_control, getold, &
!       pv_constraint, initialize_roll, rattle_roll_control, &
!       shake_roll_control
  USE constraint
  USE global_types, ONLY: global_environment_type
  USE eigenvalueproblems, ONLY : diagonalise
  USE ewald_parameters_types, ONLY : ewald_parameters_type
  USE force_control, ONLY : force
  USE kinds, ONLY : dbl
!*apsi      USE lapack, ONLY : rs
  USE linklist_control, ONLY : get_ll_parm
  USE mathconstants, ONLY : pi, zero
  USE md, ONLY : simulation_parameters_type, thermodynamic_type, &
       mdio_parameters_type
  USE molecule_types, ONLY : molecule_structure_type, particle_node_type, &
       linklist_atoms
  USE message_passing, ONLY : mp_sum
  USE nose, ONLY : extended_parameters_type, lnhc, lnhcp, lnhcpf
  USE particle_types, ONLY : particle_type
  USE simulation_cell, ONLY : get_cell_param
  USE stop_program, ONLY : stop_prg, stop_memory
  USE structure_types, ONLY : structure_type, interaction_type
  USE timings, ONLY : timeset, timestop
  USE util, ONLY : get_unit
  
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: velocity_verlet, force, energy, set_energy_parm, set_integrator
  
  TYPE old_variables_type
     REAL ( dbl ), POINTER, DIMENSION ( :, : ) :: v
     REAL ( dbl ), POINTER, DIMENSION ( :, : ) :: r
     REAL ( dbl ), POINTER, DIMENSION ( :, : ) :: eps
     REAL ( dbl ), POINTER, DIMENSION ( :, : ) :: veps
     REAL ( dbl ), POINTER, DIMENSION ( :, : ) :: eta
     REAL ( dbl ), POINTER, DIMENSION ( :, : ) :: veta
     REAL ( dbl ), POINTER, DIMENSION ( :, : ) :: h
  END TYPE old_variables_type
  
  TYPE ( old_variables_type ) :: old
  
  TYPE average_quantities_type
     REAL ( dbl ) :: avetemp, avepot, avekin, avevol, aveca, avecb, avecc
     REAL ( dbl ) :: aveal, avebe, avega, avepress, avekincoef
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
  INTEGER :: crd, vel, ptn, ene, tem, scr
  INTEGER :: icrd, ivel, iptens, iener, itemp, idump, iscreen
  
  TYPE ( global_environment_type ) :: intenv

!******************************************************************************

CONTAINS

!******************************************************************************

SUBROUTINE velocity_verlet ( itimes, constant, simpar, inter, &
     thermo, struc, ewald_param, nhcp )
  
  IMPLICIT NONE
  
! Arguments
  INTEGER, INTENT ( IN ) :: itimes
  REAL ( dbl ), INTENT ( INOUT ) :: constant
  TYPE ( simulation_parameters_type ), INTENT ( IN ) :: simpar
  TYPE ( thermodynamic_type ), INTENT ( INOUT ) :: thermo
  TYPE ( structure_type ), INTENT ( INOUT ) :: struc
  TYPE ( extended_parameters_type ), INTENT ( INOUT ) :: nhcp
  TYPE ( ewald_parameters_type ), INTENT ( INOUT ) :: ewald_param
  TYPE ( interaction_type ), INTENT ( IN ) :: inter
  
  INTEGER :: handle

!------------------------------------------------------------------------------
  CALL timeset ( 'VERLET', 'I', ' ', handle )
  SELECT CASE ( simpar % ensemble )
  CASE DEFAULT
     CALL stop_prg ( 'velocity_verlet','integrator not implemented')
  CASE ( 'NVE')
     CALL nve(itimes,constant,simpar,inter,thermo,struc,ewald_param)
     CALL energy(itimes,constant,simpar,struc,thermo,nhcp)
  CASE ( 'NVT')
     CALL nvt(itimes,constant,simpar,inter,thermo,struc,ewald_param,nhcp)
     CALL energy(itimes,constant,simpar,struc,thermo,nhcp)
  CASE ( 'NPT_I')
     CALL npt_i(itimes,constant,simpar,inter,thermo,struc,ewald_param,nhcp)
     CALL energy(itimes,constant,simpar,struc,thermo,nhcp)
  CASE ( 'NPT_F')
     CALL npt_f(itimes,constant,simpar,inter,thermo,struc,ewald_param,nhcp)
     CALL energy(itimes,constant,simpar,struc,thermo,nhcp)
  END SELECT

  CALL timestop ( zero, handle )
  
END SUBROUTINE velocity_verlet

!******************************************************************************

! call this subroutine before the first call to energy or velocity_verlet
!    or if you want to change ionode and/or output files

SUBROUTINE set_integrator ( globenv, mdio )
  IMPLICIT NONE
  
! Arguments
  TYPE(mdio_parameters_type), INTENT (IN) :: mdio
  TYPE ( global_environment_type ), INTENT (IN) :: globenv
  
!------------------------------------------------------------------------------
  
  intenv = globenv 

  ionode = globenv % ionode
  int_group = globenv % group
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
  
END SUBROUTINE set_integrator

!******************************************************************************

SUBROUTINE set_energy_parm ( pconv, tempconv, l, vol, e, pres, temp, angl )
  
  IMPLICIT NONE
  
! Arguments
  REAL ( dbl ), INTENT ( IN ), OPTIONAL :: pconv, tempconv
  CHARACTER ( LEN = 10 ), INTENT ( IN ), OPTIONAL :: l, vol, e
  CHARACTER ( LEN = 10 ), INTENT ( IN ), OPTIONAL :: pres, temp, angl
  
!------------------------------------------------------------------------------
  
  IF ( PRESENT ( pconv)) p_conv = pconv
  IF ( PRESENT ( tempconv)) temp_conv = tempconv
  IF ( PRESENT ( l)) l_label = l
  IF ( PRESENT ( vol)) vol_label = vol
  IF ( PRESENT ( e)) e_label = e
  IF ( PRESENT ( pres)) pres_label = pres
  IF ( PRESENT ( temp)) temp_label = temp
  IF ( PRESENT ( angl)) angl_label = angl
END SUBROUTINE set_energy_parm

!******************************************************************************
! this routine computes the conserved quantity, temperature
!    and things like that and prints them out

SUBROUTINE energy ( itimes, constant, simpar, struc, thermo, nhcp )
  
  IMPLICIT NONE
  
! Arguments
  INTEGER, INTENT ( IN ) :: itimes
  REAL ( dbl ), INTENT ( INOUT ) :: constant
  TYPE (simulation_parameters_type ), INTENT ( IN ) :: simpar
  TYPE (structure_type ), INTENT ( INOUT ) :: struc
  TYPE (thermodynamic_type ), INTENT ( INOUT ) :: thermo
  TYPE (extended_parameters_type ), INTENT ( INOUT ) :: nhcp
  
! Locals
  INTEGER i, j, id, n, l, nnodes, iw, handle
  REAL ( dbl ) :: nhc_pot, nhc_kin, baro_kin, baro_pot
  REAL ( dbl ) :: alpha, epsilon0, qsq, lup, aup
  REAL ( dbl ) :: cons, temp, econs
  REAL ( dbl ), DIMENSION(3) :: cell, ang
  REAL ( dbl ) :: pv_scalar
  
!------------------------------------------------------------------------------
  
  CALL timeset ( 'ENERGY', 'E', ' ', handle )
  
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
  
  temp = 2.0_dbl*thermo%kin/float(simpar%nfree)/temp_conv
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
        WRITE (crd,'(i6,3(x,f12.5) )' ) itimes, struc%pnode(i) %p%r
     END DO
  END IF
  
IF (ionode .AND. mod(itimes,ivel)==0) THEN
     DO i=1, nnodes
        WRITE (vel,'(i6,3f10.5 )' ) itimes, struc%pnode(i) %p%v
     END DO
  END IF
  
  IF (simpar%ensemble/='NVE') THEN
     nhc_kin = 0.0_dbl
     nhc_pot = 0.0_dbl
     DO n = 1, nhcp%num_nhc
        DO l = 1, nhcp%nhc_len
           nhc_kin = nhc_kin + 0.5_dbl*nhcp%nvt(l,n) %mass*nhcp%nvt(l,n) %v &
                **2
           nhc_pot = nhc_pot + nhcp%nvt(l,n) %nkt*nhcp%nvt(l,n) %eta
        END DO
     END DO
#if defined(__parallel)
     IF (nhcp%dis_type=='DIS_REP') THEN
        CALL mp_sum(nhc_kin,int_group)
        CALL mp_sum(nhc_pot,int_group)
     END IF
#endif
  END IF
  
  IF (simpar%ensemble=='NPT_I') THEN
     baro_pot = simpar%p_ext*struc%box%deth
     baro_kin = 0.5_dbl*nhcp%npt(1,1) %v*nhcp%npt(1,1) %v* &
          nhcp%npt(1,1) %mass
  ELSE IF (simpar%ensemble=='NPT_F') THEN
     baro_pot = simpar%p_ext*struc%box%deth
     baro_kin = 0._dbl
     DO i = 1, 3
        DO j = 1, 3
           baro_kin = baro_kin + 0.5_dbl*nhcp%npt(i,j) %v*nhcp%npt(i,j) %v* &
                nhcp%npt(i,j) %mass
        END DO
     END DO
  END IF
  
  IF (simpar%ensemble(1:3)=='NPT') THEN
#if defined(__parallel)
     IF (nhcp%dis_type=='DIS_REP') THEN
        CALL mp_sum(baro_kin,int_group)
        CALL mp_sum(baro_pot,int_group)
     END IF
#endif
  END IF
  
  SELECT CASE (simpar%ensemble(1:3))
  CASE ( 'NVE')
     cons = thermo%kin + thermo%pot
  CASE ( 'NVT')
     cons = thermo%kin + thermo%pot + nhc_kin + nhc_pot
  CASE ( 'NPT')
     cons = thermo%kin + thermo%pot + nhc_kin + nhc_pot + baro_kin + &
          baro_pot
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
     CALL get_ll_parm(lup,aup)
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
  
END SUBROUTINE energy

!******************************************************************************

! nve integrator for particle positions & momenta
!    velocity Verlet

SUBROUTINE nve(itimes,constant,simpar,inter,thermo,struc,ewald_param)
  IMPLICIT NONE
  INTEGER, INTENT ( IN ) :: itimes
  REAL ( dbl ), INTENT ( INOUT ) :: constant
  TYPE (simulation_parameters_type ), INTENT ( IN ) :: simpar
  TYPE (thermodynamic_type ), INTENT ( INOUT ) :: thermo
  TYPE (structure_type ), INTENT ( INOUT ) :: struc
  TYPE (ewald_parameters_type ), INTENT ( INOUT ) :: ewald_param
  TYPE (interaction_type ), INTENT ( IN ) :: inter
  
  INTEGER :: i, nnodes
  REAL ( dbl ) :: dtom
  
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
  CALL force ( struc, inter, thermo, simpar, ewald_param, intenv )
  
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
       'FULL', int_group )
  
  CALL pressure(struc%pnode,thermo)
  
END SUBROUTINE nve

!******************************************************************************

SUBROUTINE nvt(itimes,constant,simpar,inter,thermo,struc, &
     ewald_param,nhcp)
  
  IMPLICIT NONE
  
! Arguments
  INTEGER, INTENT ( IN ) :: itimes
  REAL ( dbl ), INTENT ( INOUT ) :: constant
  TYPE (simulation_parameters_type ), INTENT ( IN ) :: simpar
  TYPE (thermodynamic_type ), INTENT ( INOUT ) :: thermo
  TYPE (structure_type ), INTENT ( INOUT ) :: struc
  TYPE (extended_parameters_type ), INTENT ( INOUT ) :: nhcp
  TYPE (ewald_parameters_type ), INTENT ( INOUT ) :: ewald_param
  TYPE (interaction_type ), INTENT ( IN ) :: inter
  
! Locals
  INTEGER :: i, nnodes
  REAL ( dbl ) :: dtom
  
!------------------------------------------------------------------------------
  
  CALL lnhc ( nhcp, struc % pnode, int_group )
  
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
  
  IF ( simpar % constraint ) CALL shake_control ( struc % molecule, simpar )
#if defined(__parallel)
  CALL update_structure ( struc, 'POS' )
#endif
  
!
! get new forces
!
  CALL force ( struc, inter, thermo, simpar, ewald_param, intenv )
  
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
  
  CALL lnhc ( nhcp, struc % pnode, int_group )
  
#if defined(__parallel)
  CALL update_structure(struc,'VEL')
#endif
  
! updating the virial
  IF (simpar%constraint) CALL pv_constraint(struc%molecule,thermo, &
       'FULL', int_group)
  CALL pressure(struc%pnode,thermo)
  
END SUBROUTINE nvt

!******************************************************************************

SUBROUTINE npt_i ( itimes, constant, simpar, inter, thermo, struc, &
     ewald_param, nhcp )
  
  IMPLICIT NONE
  
! Arguments
  INTEGER, INTENT ( IN ) :: itimes
  REAL ( dbl ), INTENT ( INOUT ) :: constant
  TYPE (simulation_parameters_type ), INTENT ( IN ) :: simpar
  TYPE (thermodynamic_type ), INTENT ( INOUT ) :: thermo
  TYPE (structure_type ), INTENT ( INOUT ) :: struc
  TYPE (extended_parameters_type ), INTENT ( INOUT ) :: nhcp
  TYPE (ewald_parameters_type ), INTENT ( INOUT ) :: ewald_param
  TYPE (interaction_type ), INTENT ( IN ) :: inter
  
! Locals
  INTEGER :: i, nnodes, iroll
  REAL ( dbl ) :: dtom
  REAL ( dbl ), PARAMETER :: e2 = 1._dbl/6._dbl, e4 = e2/20._dbl, &
       e6 = e4/42._dbl, e8 = e6/72._dbl
  REAL ( dbl ), SAVE :: eps_0
  REAL ( dbl ) :: arg, poly, scale, roll_tol
  REAL ( dbl ), DIMENSION (3,3) :: pv_kinetic
  REAL ( dbl ), DIMENSION (3) :: vector
  LOGICAL :: first
  
!------------------------------------------------------------------------------
  
! setting up for ROLL: saving old variables
  roll_tol = -1.0E+10_dbl
  IF (simpar%constraint) THEN
     iroll = 1
     CALL set(struc,nhcp,'F')
     CALL getold(struc%molecule)
  END IF
  
! initialize eps_0 the first time through
  IF ( itimes == 1 ) eps_0 = nhcp%npt(1,1) % eps
  
  SR: DO WHILE (abs(roll_tol)>=1.E-10_dbl) ! SHAKE-ROLL LOOP
     
     IF (simpar%constraint) THEN
        CALL set(struc,nhcp,'B')
        CALL initialize_roll(struc%molecule)
     END IF
     CALL lnhcp(nhcp,struc%molecule,struc%pnode,thermo%pv,thermo%pv_kin, &
          thermo%pv_const,simpar%p_ext,struc%box%deth,simpar%nfree, &
          simpar%constraint, int_group)
     
!
! first half of velocity verlet
!
     arg = ( 0.5_dbl*nhcp%npt(1,1)%v*simpar%dt)*(0.5_dbl*nhcp%npt(1,1)%v* &
          simpar%dt)
     poly = 1._dbl + e2*arg + e4*arg*arg + e6*arg**3 + e8*arg**4
     
     scale = EXP(0.5_dbl*simpar%dt*nhcp%npt(1,1)%v)
     nnodes = size(struc%pnode)
     DO i = 1, nnodes
        dtom = 0.5_dbl*simpar%dt/struc%pnode(i) %p%prop%mass
        struc%pnode(i) %p%v(1) = struc%pnode(i) %p%v(1) + &
             dtom*struc%pnode(i) %p%f(1)
        struc%pnode(i) %p%v(2) = struc%pnode(i) %p%v(2) + &
             dtom*struc%pnode(i) %p%f(2)
        struc%pnode(i) %p%v(3) = struc%pnode(i) %p%v(3) + &
             dtom*struc%pnode(i) %p%f(3)
        struc%pnode(i) %p%r(1) = struc%pnode(i) %p%r(1)*scale*scale + &
             scale*poly*struc%pnode(i) %p%v(1)*simpar%dt
        struc%pnode(i) %p%r(2) = struc%pnode(i) %p%r(2)*scale*scale + &
             scale*poly*struc%pnode(i) %p%v(2)*simpar%dt
        struc%pnode(i) %p%r(3) = struc%pnode(i) %p%r(3)*scale*scale + &
             scale*poly*struc%pnode(i) %p%v(3)*simpar%dt
     END DO
     roll_tol = 0._dbl
     vector ( : ) = scale*poly
     
     IF ( simpar%constraint ) CALL shake_roll_control(struc%molecule,simpar, &
          thermo,roll_tol,iroll,vector, int_group )
! Update eps:
     nhcp%npt ( :, : ) %eps = nhcp%npt ( :, : ) %eps &
          + simpar%dt*nhcp%npt ( :, : ) %v
     
! Update h_mat
     struc%box%hmat ( :, : ) = struc%box%hmat ( :, : ) &
          * EXP ( nhcp%npt(1,1)%eps - eps_0 )
  END DO SR
  
  eps_0 = nhcp%npt(1,1) % eps
  
#if defined(__parallel)
  CALL update_structure(struc,'POS')
#endif
  
!
! get new forces
!
  CALL force(struc,inter,thermo,simpar,ewald_param,intenv)
  
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
  
  roll_tol = -1.E10_dbl
  IF (simpar%constraint) THEN
     first = .TRUE.
     iroll = 1
     CALL set(struc,nhcp,'F')
  END IF
  
  RR: DO WHILE (abs(roll_tol)>=1.E-10_dbl) ! RATTLE-ROLL LOOP
     roll_tol = 0._dbl
     IF (simpar%constraint) CALL rattle_init(nhcp,struc,thermo,simpar, &
          roll_tol,iroll,first)
     CALL lnhcp(nhcp,struc%molecule,struc%pnode,thermo%pv,thermo%pv_kin, &
          thermo%pv_const,simpar%p_ext,struc%box%deth,simpar%nfree, &
          simpar%constraint, int_group )
  END DO RR
  
#if defined(__parallel)
  CALL update_structure(struc,'VEL')
#endif
  
! updating the virial
  IF (simpar%constraint) CALL pv_constraint(struc%molecule,thermo, &
       'FULL', int_group )
  
  CALL pressure(struc%pnode,thermo)
  
END SUBROUTINE npt_i

!******************************************************************************

! computes the kinetic part of the pressure tensor and updates
!    the full VIRIAL (PV)

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

  thermo%pv_kin = 0._dbl
  nnodes = size(pnode)
  DO i = 1, 3
     DO j = 1, 3
        DO iatom = 1, nnodes
           thermo%pv_kin(i,j) = thermo%pv_kin(i,j) + &
                pnode(iatom) %p%prop%mass*pnode(iatom) %p%v(i)* &
                pnode(iatom) %p%v(j)
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

  natoms = size(struc%part)
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
  
  DEALLOCATE ( atot, STAT = ios )
  IF ( ios /= 0 ) CALL stop_memory ( 'update_structure', 'atot' )

  CALL timestop ( zero, handle )
  
END SUBROUTINE update_structure

!******************************************************************************

SUBROUTINE set(struc,nhcp,char)
  
  IMPLICIT NONE
  
! Arguments
  TYPE (structure_type ), INTENT ( INOUT ) :: struc
  TYPE (extended_parameters_type ), INTENT ( INOUT ) :: nhcp
  CHARACTER ( LEN = * ), INTENT ( IN ) :: char
  
! Locals
  INTEGER :: idim
  INTEGER :: isos
  
!------------------------------------------------------------------------------
  
  isos = 0
  IF ( .NOT. associated(old%v)) ALLOCATE (old%v(size(struc%part),3), &
       STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'set', 'old%v', 0 )
  IF ( .NOT. associated(old%r)) ALLOCATE (old%r(size(struc%part),3), &
       STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'set', 'old%r', 0 )
  IF ( .NOT. associated(old%eps)) ALLOCATE (old%eps(size(nhcp%npt, &
       1),size(nhcp%npt,2)),STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'set', 'old%eps', 0 )
  IF ( .NOT. associated(old%veps)) ALLOCATE (old%veps(size(nhcp%npt, &
       1),size(nhcp%npt,2)),STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'set', 'old%veps', 0 )
  IF ( .NOT. associated(old%eta)) ALLOCATE (old%eta(size(nhcp%nvt, &
       1),size(nhcp%nvt,2)),STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'set', 'old%eta', 0 )
  IF ( .NOT. associated(old%veta)) ALLOCATE (old%veta(size(nhcp%nvt, &
       1),size(nhcp%nvt,2)),STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'set', 'old%veta', 0 )
  IF ( .NOT. associated(old%h)) ALLOCATE (old%h(size(struc%box%hmat,1), &
       size(struc%box%hmat,2)),STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'set', 'old%h', 0 )
  
  SELECT CASE (char)
  CASE ( 'F') ! forward assigning the old
     DO idim = 1, 3
        old%v(:,idim) = struc%part ( : ) %v(idim)
        old%r(:,idim) = struc%part ( : ) %r(idim)
     END DO
     old%eta ( :, : ) = nhcp%nvt ( :, : ) %eta
     old%veta ( :, : ) = nhcp%nvt ( :, : ) %v
     old%eps ( :, : ) = nhcp%npt ( :, : ) %eps
     old%veps ( :, : ) = nhcp%npt ( :, : ) %v
     old%h ( :, : ) = struc%box%hmat ( :, : )
  CASE ( 'B') ! back assigning the original variables
     DO idim = 1, 3
        struc%part ( : ) %r(idim) = old%r(:,idim)
        struc%part ( : ) %v(idim) = old%v(:,idim)
     END DO
     nhcp%nvt ( :, : ) %eta = old%eta ( :, : )
     nhcp%nvt ( :, : ) %v = old%veta ( :, : )
     nhcp%npt ( :, : ) %eps = old%eps ( :, : )
     nhcp%npt ( :, : ) %v = old%veps ( :, : )
     struc%box%hmat ( :, : )=old%h ( :, : )
  END SELECT
  
END SUBROUTINE set

!******************************************************************************

SUBROUTINE npt_f ( itimes, constant, simpar, inter, thermo, struc, &
     ewald_param, nhcp )
  
  IMPLICIT NONE
  
! Arguments
  INTEGER, INTENT ( IN ) :: itimes
  REAL ( dbl ), INTENT ( INOUT ) :: constant
  TYPE (simulation_parameters_type ), INTENT ( IN ) :: simpar
  TYPE (thermodynamic_type ), INTENT ( INOUT ) :: thermo
  TYPE (structure_type ), INTENT ( INOUT ) :: struc
  TYPE (extended_parameters_type ), INTENT ( INOUT ) :: nhcp
  TYPE (ewald_parameters_type ), INTENT ( INOUT ) :: ewald_param
  TYPE (interaction_type ), INTENT ( IN ) :: inter
  
! Locals
  INTEGER :: i, j, nnodes, iroll
  REAL ( dbl ) :: dtom
  REAL ( dbl ), PARAMETER :: e2 = 1._dbl/6._dbl, e4 = e2/20._dbl, &
       e6 = e4/42._dbl, e8 = e6/72._dbl
  REAL ( dbl ), SAVE :: eps_0
  REAL ( dbl ) :: roll_tol
  REAL ( dbl ), DIMENSION (3,3) :: pv_kinetic, u, uh
  REAL ( dbl ), DIMENSION (3) :: arg, poly, scale, e_val, ur, uv
  LOGICAL :: first
  
!------------------------------------------------------------------------------
  
! setting up for ROLL: saving old variables
  IF (simpar%constraint) THEN
     iroll = 1
     roll_tol = -1.0E+10_dbl
     CALL set(struc,nhcp,'F')
     CALL getold(struc%molecule)
  END IF
  
  SR: DO WHILE (abs(roll_tol)>=1.E-10_dbl) ! SHAKE-ROLL LOOP
     
     IF (simpar%constraint) THEN
        CALL set(struc,nhcp,'B')
        CALL initialize_roll(struc%molecule)
     END IF
     CALL lnhcpf(nhcp,struc%molecule,struc%pnode,thermo%pv,thermo%pv_kin, &
          thermo%pv_const,simpar%p_ext,struc%box%deth,simpar%nfree, &
          simpar%constraint, int_group )
     
! find eigenvalues and eigenvectors of nhcp%npt ( :, : )%v
!
!*apsi          CALL rs(nhcp%npt%v,e_val,u)
     CALL diagonalise ( matrix = nhcp % npt % v, mysize = 3, &
          storageform = "UPPER", eigenvalues = e_val, eigenvectors = u )
     
!
! first half of velocity verlet
!
     arg ( : ) = (0.5_dbl*e_val ( : )*simpar%dt)*(0.5_dbl*e_val ( : )*simpar%dt)
     poly = 1._dbl + e2*arg + e4*arg*arg + e6*arg**3 + e8*arg**4
     scale ( : ) = EXP(0.5_dbl*simpar%dt*e_val ( : ))
     
     nnodes = size(struc%pnode)
     DO i = 1, nnodes
        dtom = 0.5_dbl*simpar%dt/struc%pnode(i) %p%prop%mass
        struc%pnode(i) %p%v(1) = struc%pnode(i) %p%v(1) + &
             dtom*struc%pnode(i) %p%f(1)
        struc%pnode(i) %p%v(2) = struc%pnode(i) %p%v(2) + &
             dtom*struc%pnode(i) %p%f(2)
        struc%pnode(i) %p%v(3) = struc%pnode(i) %p%v(3) + &
             dtom*struc%pnode(i) %p%f(3)
        
! Transform positions and velocities
        ur = matmul(transpose(u),struc%pnode(i)%p%r)
        uv = matmul(transpose(u),struc%pnode(i)%p%v)
        ur(1) = ur(1)*scale(1)*scale(1) + uv(1)*scale(1)*poly(1)*simpar%dt
        ur(2) = ur(2)*scale(2)*scale(2) + uv(2)*scale(2)*poly(2)*simpar%dt
        ur(3) = ur(3)*scale(3)*scale(1) + uv(3)*scale(3)*poly(3)*simpar%dt
        struc%pnode(i) %p%r ( : ) = matmul(u,ur)
     END DO
     
     roll_tol = 0.0_dbl
     IF ( simpar % constraint ) CALL shake_roll_control ( struc % molecule, &
          simpar, thermo,roll_tol,iroll,scale*poly, int_group ,u)
  END DO SR
  
! Update h_mat
  uh = matmul(transpose(u),struc%box%hmat)
  DO i = 1, 3
     DO j = 1, 3
        uh(j,i) = uh(j,i)*scale(j)*scale(j)
     END DO
  END DO
  struc%box%hmat = matmul(u,uh)
  
#if defined(__parallel)
  CALL update_structure(struc,'POS')
#endif
  
!
! get new forces
!
  CALL force ( struc, inter, thermo, simpar, ewald_param, intenv )
  
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
  
  IF (simpar%constraint) THEN
     roll_tol = -1.0E+10_dbl
     first = .TRUE.
     iroll = 1
     CALL set ( struc, nhcp, 'F' )
  END IF
  
  RR: DO WHILE ( ABS ( roll_tol ) >= 1.0E-10_dbl ) ! RATTLE-ROLL LOOP
     roll_tol = 0.0_dbl
     IF ( simpar%constraint ) CALL rattle_init(nhcp,struc,thermo,simpar, &
          roll_tol,iroll,first)
     CALL lnhcpf(nhcp,struc%molecule,struc%pnode,thermo%pv,thermo%pv_kin, &
          thermo%pv_const,simpar%p_ext,struc%box%deth,simpar%nfree, &
          simpar%constraint, int_group )
     
  END DO RR
  
#if defined(__parallel)
  CALL update_structure(struc,'VEL')
#endif
  
! updating the virial
  IF (simpar%constraint) CALL pv_constraint(struc%molecule,thermo, &
       'FULL', int_group )
  CALL pressure(struc%pnode,thermo)
  
END SUBROUTINE npt_f

!******************************************************************************

SUBROUTINE rattle_init ( nhcp, struc, thermo, simpar, roll_tol, iroll, first )
  
  IMPLICIT NONE
  
! Arguments
  TYPE (simulation_parameters_type ), INTENT ( IN ) :: simpar
  TYPE (structure_type ), INTENT ( INOUT ) :: struc
  TYPE (extended_parameters_type ), INTENT ( INOUT ) :: nhcp
  TYPE (thermodynamic_type ), INTENT ( INOUT ) :: thermo
  LOGICAL, INTENT ( INOUT ) :: first
  REAL ( dbl ), INTENT ( OUT ) :: roll_tol
  INTEGER, INTENT ( INOUT ) :: iroll
  
! Locals
  REAL ( dbl ), DIMENSION (3,3) :: veps
  
!------------------------------------------------------------------------------
  
  IF ( first ) THEN
     CALL initialize_roll ( struc%molecule )
     CALL lnhcp ( nhcp,struc%molecule,struc%pnode,thermo%pv,thermo%pv_kin, &
          thermo%pv_const,simpar%p_ext,struc%box%deth,simpar%nfree, &
          simpar%constraint, int_group )
  END IF
  first = .FALSE.
  
! assigning local variable
  veps ( :, : ) =  0._dbl
  veps(1,1) = nhcp%npt(1,1) %v
  veps(2,2) = veps(1,1)
  veps(3,3) = veps(1,1)
  
! resetting
  CALL set(struc,nhcp,'B')
  CALL rattle_roll_control(struc%molecule,simpar,thermo,veps, &
       roll_tol,iroll, int_group )
  CALL initialize_roll(struc%molecule)
  
END SUBROUTINE rattle_init

!******************************************************************************

END MODULE integrator
