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
!!
!!   SOURCE
!******************************************************************************

MODULE md_kg_energies
  
  USE coefficient_types,                 ONLY : coeff_type
  USE global_types,                      ONLY : global_environment_type
  USE kg_environment_types,              ONLY : kg_thermodynamic_type
  USE kinds,                             ONLY : dbl
  USE mathconstants,                     ONLY : zero
  USE md,                                ONLY : simulation_parameters_type,&
                                                mdio_parameters_type
  USE md_environment_types,              ONLY : md_environment_type, &
                                                get_md_env
  USE message_passing,                   ONLY : mp_sum,&
                                                mp_comm_dup
  USE nose,                              ONLY : lnhc_parameters_type
  USE simulation_cell,                   ONLY : get_cell_param, &
                                                cell_type
  USE structure_types,                   ONLY : structure_type
  USE termination,                       ONLY : stop_memory, &
                                                stop_program
  USE timings,                           ONLY : timeset,&
                                                timestop
  USE util,                              ONLY : get_unit

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: set_kg_energies, set_kg_energies_param, md_energy_kg

  TYPE average_quantities_type
     REAL ( dbl ) :: avetemp, avepot, avekin, avecore, avexc, avehartree, avebc
     REAL ( dbl ) :: avepress, avekinc, avetempc
  END TYPE average_quantities_type

  CHARACTER ( LEN = 10 ) :: pres_label, e_label, temp_label
  CHARACTER ( LEN = 10 ) :: l_label, vol_label, angl_label
  REAL ( dbl ) :: econs_avg
  REAL ( dbl ) :: ecore_avg, ehartree_avg, exc_avg
  REAL ( dbl ) :: p_conv, temp_conv
  TYPE (average_quantities_type ) :: averages

  LOGICAL :: ionode
  INTEGER :: int_group
  CHARACTER (LEN=50) :: crd_file_name, vel_file_name, ptens_file_name, &
                        ener_file_name, temp_file_name, dump_file_name, &
                        pot_file_name, coef_file_name
  INTEGER :: crd, vel, pot, ene, tem, scr, coe, pol
  INTEGER :: icrd, ivel, ipot, iener, itemp, idump, icoef, iscreen

  TYPE ( global_environment_type ) :: glob_env

!******************************************************************************

CONTAINS

!******************************************************************************
!!****** md_kg_energies/set_kg_energies [1.0] *
!!
!!   NAME
!!     set_kg_energies
!!
!!   FUNCTION
!!     setup of parameters for kg_energies
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

SUBROUTINE set_kg_energies ( globenv, mdio )

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
  pot_file_name = mdio % pot_file_name
  ener_file_name = mdio % ener_file_name
  temp_file_name = mdio % temp_file_name
  dump_file_name = mdio % dump_file_name
  coef_file_name = mdio % coef_file_name

  icrd = mdio % icrd
  icoef = mdio % icoef
  ivel = mdio % ivel
  ipot = mdio % ipot
  iener = mdio % iener
  itemp = mdio % itemp
  idump = mdio % idump
  iscreen = mdio % iscreen

END SUBROUTINE set_kg_energies

!******************************************************************************
!!****** md_kg_energies/set_kg_energies_params [1.0] *
!!
!!   NAME
!!     set_kg_energies_param
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

SUBROUTINE set_kg_energies_param ( pconv, tempconv, l, vol, e, pres, temp, angl )

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

END SUBROUTINE set_kg_energies_param

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


SUBROUTINE md_energy_kg ( md_env, thermo )

  IMPLICIT NONE

! Arguments
  TYPE ( md_environment_type ), INTENT ( INOUT ) :: md_env
  TYPE (kg_thermodynamic_type ), INTENT ( INOUT ) :: thermo

! Locals
  
  INTEGER, POINTER :: itimes
  REAL (dbl), POINTER :: constant
  TYPE ( cell_type ), POINTER :: box
  TYPE ( lnhc_parameters_type ), POINTER :: nhc_part 
  TYPE ( lnhc_parameters_type ), POINTER :: nhc_coef 
  TYPE ( simulation_parameters_type ), POINTER :: simpar
  TYPE ( structure_type ), POINTER :: struc 

  INTEGER i, j, id, n, l, nnodes, ncoef, iw, handle, inuse
  REAL ( dbl ) :: nhc_part_pot, nhc_part_kin
  REAL ( dbl ) :: coef_temp, coef_kin
  REAL ( dbl ) :: nhc_coef_pot, nhc_coef_kin
  REAL ( dbl ) :: alpha, epsilon0, qsq, lup, aup
  REAL ( dbl ) :: cons, temp, econs
  REAL ( dbl ), DIMENSION(3) :: cell, ang
  REAL ( dbl ) :: pv_scalar

!------------------------------------------------------------------------------

  CALL timeset ( 'ENERGY', 'E', ' ', handle )
  CALL get_md_env ( md_env=md_env, itimes=itimes, constant=constant, &
                            cell=box, simpar=simpar )                  
  nhc_part => md_env % nhc_part ( 1 )
  struc => md_env % struc ( 1 )

  nnodes = SIZE ( struc % pnode )
 
  IF (ASSOCIATED (struc % coef_vel ) ) THEN 
    nhc_coef => md_env % nhc_coef         
    inuse = struc % coef_pos ( 1 ) % in_use
    SELECT CASE ( inuse )
    CASE ( 402 )
      ncoef = SIZE ( struc % coef_vel ( 1 ) % ao % cr )
    CASE ( 401 )
      ncoef = SIZE ( struc % coef_vel ( 1 ) % pw % cc )
    END SELECT
   END IF

  IF ( ionode .AND. itimes == 0 ) THEN
     tem = get_unit()
     OPEN ( UNIT = tem, FILE = temp_file_name )
     ene = get_unit()
     OPEN ( UNIT = ene, FILE = ener_file_name )
     crd = get_unit()
     OPEN ( UNIT = crd, FILE = crd_file_name )
     vel = get_unit()
     OPEN ( UNIT = vel, FILE = vel_file_name )
     pot = get_unit()
     OPEN ( UNIT = pot, FILE = pot_file_name )
     IF (ASSOCIATED (struc % coef_vel ) ) THEN 
       coe = get_unit()
       OPEN ( UNIT = coe, FILE = coef_file_name )
     END IF
  END IF

  iw = scr

  CALL get_part_ke ( struc, thermo )
  
  CALL get_part_temp ( thermo, simpar, temp )

  IF ( ASSOCIATED ( struc % coef_vel  ) ) THEN
  
    CALL get_coef_ke ( struc % coef_vel ( 1 ), inuse, ncoef, coef_kin )

    CALL get_coef_temp ( coef_kin, inuse, ncoef, coef_temp )

  ENDIF 

  IF (ionode .AND. mod(itimes,itemp)==0) THEN
     WRITE (tem, *) itimes, temp
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

  IF ( ASSOCIATED ( struc % coef_vel  ) ) THEN
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
  END IF

  IF ( simpar % ensemble ( 1 : 3 )/= 'NVE' ) &
  CALL get_nhc_energies ( nhc_part, nhc_part_kin, nhc_part_pot )

  IF ( simpar % ensemble_coef ( 1 : 3 )/= 'NVE' ) &
  CALL get_nhc_energies ( nhc_coef, nhc_coef_kin, nhc_coef_pot )

  IF ( simpar % ensemble ( 1 : 3 ) == 'NPT' ) &
  CALL stop_program ( "md_kg_energies", "npt not implemented " )

  SELECT CASE ( simpar % ensemble ( 1 : 3 ) )
  CASE ( 'NVE')
     cons = thermo%kin + thermo%pot
  CASE ( 'NVT')
     cons = thermo%kin + thermo%pot + nhc_part_kin + nhc_part_pot
  END SELECT

  IF  ( ASSOCIATED ( struc % coef_vel  ) ) THEN 
     SELECT CASE ( simpar % ensemble_coef ( 1 : 3 ) )
     CASE ( 'NVE')
        cons = cons + coef_kin
     CASE ( 'NVT')
        cons = cons + nhc_coef_kin + nhc_coef_pot + coef_kin
     END SELECT
  END IF

  IF (itimes==1) THEN
     constant = cons
  END IF
  IF (itimes /= 0 ) econs = abs((cons-constant)/constant)

! cell parameters:
  CALL get_cell_param ( box, cell, ang )

  SELECT CASE (itimes)
  CASE (0)
     IF (ionode) THEN
        WRITE (iw,*)
        WRITE ( iw, '( A )' ) ' ENERGY| Initialization proceeding'
        WRITE (iw,*)
     END IF
  CASE (1)
     averages%avetemp = temp
     averages%avecore = thermo%e_rspace_core
     averages%avehartree = thermo%e_hartree
     averages%avebc = thermo%e_bond_corr
     averages%avexc = thermo%e_xc
     averages%avepot = thermo%pot
     averages%avekin = thermo%kin
     IF  ( ASSOCIATED ( struc % coef_vel ) ) THEN 
       averages%avekinc = coef_kin   
       averages%avetempc = coef_temp
     END IF
     econs_avg =  0._dbl
  CASE DEFAULT
     averages%avetemp = (averages%avetemp*real(itimes-1)+temp)/ &
          real(itimes)
     averages%avecore = (averages%avecore*real(itimes-1)+thermo%e_rspace_core)/ &
          real(itimes)
     averages%avehartree = (averages%avehartree*real(itimes-1)+thermo%e_hartree)/ &
          real(itimes)
     averages%avebc = (averages%avebc*real(itimes-1)+thermo%e_bond_corr)/ &
          real(itimes)
     averages%avexc = (averages%avexc*real(itimes-1)+thermo%e_xc)/ &
          real(itimes)
     averages%avepot = (averages%avepot*real(itimes-1)+thermo%pot)/ &
          real(itimes)
     averages%avekin = (averages%avekin*real(itimes-1)+thermo%kin)/ &
          real(itimes)
     IF  ( ASSOCIATED ( struc % coef_vel ) ) THEN 
       averages%avekinc = (averages%avekinc*real(itimes-1)+ coef_kin)/ &
            real(itimes)
       averages%avetempc = (averages%avetempc*real(itimes-1)+coef_temp)/ &
            real(itimes)
     END IF
     econs_avg = econs_avg + abs((cons-constant)/constant)
  END SELECT

  IF (ionode .AND. (mod(itimes,iener)==0) .AND. (itimes /= 0 )) THEN
     WRITE (ene, *) itimes, econs, econs_avg/real(itimes)
  END IF

 IF (ionode .AND. (mod(itimes,ipot)==0) .AND. (itimes /= 0 )) THEN
     WRITE (pot, *) itimes, thermo%e_bond_corr, thermo%e_xc, &
     thermo%e_hartree, thermo%e_rspace_core, thermo%pot
  END IF


  IF (itimes==0 .AND. ionode) THEN
     WRITE ( iw, '( )' )
     WRITE ( iw, '( A,A )' ) ' ******************************** ', &
          'GO CP2K GO! **********************************'
     WRITE ( iw, '( A,A,T40,A,T65,1(1X,E15.7) )' ) ' INITIAL TOTAL POT ENERGY', &
          adjustl(e_label), '= ', thermo%pot
     WRITE ( iw, '( A,A,T40,A,T65,1(1X,E15.7) )' ) ' INITIAL BOND CORRECT ENERGY', &
          adjustl(e_label), '= ', thermo%e_bond_corr
     WRITE ( iw, '( A,A,T40,A,T65,1(1X,E15.7) )' ) ' INITIAL K+XC ENERGY', &
          adjustl(e_label), '= ', thermo%e_xc
     WRITE ( iw, '( A,A,T40,A,T65,1(1X,E15.7) )' ) ' INITIAL HARTREE ENERGY', &
          adjustl(e_label), '= ', thermo%e_hartree
     WRITE ( iw, '( A,A,T40,A,T65,1(1X,E15.7) )' ) ' INITIAL CORE OVLP ENERGY', &
          adjustl(e_label), '= ', thermo%e_rspace_core
     WRITE ( iw, '( A,A,T40,A,T65,1(1X,E15.7) )' ) ' SELF-ENERGY', &
          adjustl(e_label), '= ', thermo%e_self
     WRITE ( iw, '( A,A,T40,A,T65,1(1X,E15.7) )' ) ' INITIAL KINETIC ENERGY', &
          adjustl(e_label), '= ', thermo%kin
     WRITE ( iw, '( A,A,T40,A,T65,1(1X,E15.7) )' ) ' INITIAL TEMPERATURE', &
          adjustl(temp_label), '= ', temp
     IF  ( ASSOCIATED ( struc % coef_vel  ) ) THEN 
       WRITE ( iw, '( A,A,T40,A,T65,1(1X,E15.7) )' ) ' INITIAL COEF KINETIC ENERGY', &
            adjustl(e_label), '= ', coef_kin
       WRITE ( iw, '( A,A,T40,A,T65,1(1X,E15.7) )' ) ' INITIAL COEF TEMPERATURE', &
            adjustl(temp_label), '= ', coef_temp
     END IF 
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
     WRITE ( iw, '( A,A,T31,A,T49,2(1X,E15.7) )' ) ' TOTAL POT ENERGY', &
          adjustl(e_label), '= ', thermo%pot, averages%avepot
     WRITE ( iw, '( A,A,T31,A,T49,2(1X,E15.7) )' ) ' BOND CORRECT ENERGY', &
          adjustl(e_label), '= ', thermo%e_bond_corr, averages%avebc
     WRITE ( iw, '( A,A,T31,A,T49,2(1X,E15.7) )' ) ' K+XC ENERGY', &
          adjustl(e_label), '= ', thermo%e_xc, averages%avexc
     WRITE ( iw, '( A,A,T31,A,T49,2(1X,E15.7) )' ) ' CORE OVLP ENERGY', &
          adjustl(e_label), '= ', thermo%e_rspace_core, averages%avecore
     WRITE ( iw, '( A,A,T31,A,T49,2(1X,E15.7) )' ) ' HARTREE   ENERGY', &
          adjustl(e_label), '= ', thermo%e_hartree, averages%avehartree
     WRITE ( iw, '( A,A,T31,A,T49,2(1X,E15.7) )' ) ' KINETIC ENERGY', &
          adjustl(e_label), '= ', thermo%kin, averages%avekin
     WRITE ( iw, '( A,A,T31,A,T49,2(1X,E15.7) )' ) ' TEMPERATURE', &
          adjustl(temp_label), '= ', temp, averages%avetemp
     IF  ( ASSOCIATED ( struc % coef_vel ) ) THEN 
       WRITE ( iw, '( A,A,T31,A,T49,2(1X,E15.7) )' ) ' COEF KINETIC ENERGY', &
            adjustl(e_label), '= ', coef_kin,  averages%avekinc
       WRITE ( iw, '( A,A,T31,A,T49,2(1X,E15.7) )' ) ' COEF TEMPERATURE', &
            adjustl(temp_label), '= ', coef_temp, averages%avetempc
     END IF
     WRITE ( iw, '( A,A )' ) ' **************************************', &
          '*****************************************'
     WRITE ( iw, '( )' )
  END IF

  CALL timestop ( zero, handle )

END SUBROUTINE md_energy_kg

!******************************************************************************
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
  TYPE ( kg_thermodynamic_type ), INTENT ( INOUT ) :: thermo

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
  TYPE ( kg_thermodynamic_type ), INTENT ( INOUT ) :: thermo

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
 SUBROUTINE get_coef_temp ( ekinc, inuse, ncoef, ctemp )
   IMPLICIT NONE
   REAL (dbl), INTENT ( IN ) :: ekinc
   INTEGER, INTENT ( IN ) :: inuse, ncoef
   REAL ( dbl ), INTENT ( OUT ) :: ctemp
!! locals
!
   SELECT CASE ( inuse )
   CASE ( 402 )
     ctemp = 2.0_dbl* ekinc / REAL ( ncoef, dbl ) / temp_conv
   CASE ( 401 )
     ctemp =  ekinc / REAL ( ncoef, dbl ) / temp_conv
   CASE DEFAULT
      CALL stop_program ( 'get_coef_temp','integrator not implemented')
   END SELECT
 END SUBROUTINE  get_coef_temp
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
 SUBROUTINE get_coef_ke ( coef, inuse, ncoef, ekinc )
   IMPLICIT NONE
   TYPE ( coeff_type ), INTENT ( IN ) :: coef
   INTEGER, INTENT ( IN ) :: inuse, ncoef
   REAL ( dbl ), INTENT ( out ) :: ekinc
!
   INTEGER :: i
   REAL ( dbl ) :: cmass
   ekinc = zero
!
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
 
!! sum all contributions to energy over calculated parts on all processors
!#if defined(__parallel)
!  CALL mp_sum(ekinc,int_group)
!#endif
 END SUBROUTINE get_coef_ke
!******************************************************************************
END MODULE md_kg_energies
!******************************************************************************
