!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2003  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/md_energies [1.0] *
!!
!!   NAME
!!     md_energies
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

MODULE md_energies

  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_logger_type,&
                                             cp_to_string
  USE cp_subsystem_types,              ONLY: cp_subsys_get,&
                                             cp_subsystem_p_type
  USE cp2k_units,                      ONLY: unit_convert_type, &
                                             get_cp2k_units
  USE dynamical_coeff_types,           ONLY: dyn_coeff_set_p_type
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_p_type, &
                                             use_fist_force, &
                                             use_qs_force,& 
                                             use_kg_force,& 
                                             use_pol_force,& 
                                             use_qmmm, &
                                             use_list_force
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp, dbl
  USE machine,                         ONLY: m_flush
  USE mathconstants,                   ONLY: zero
  USE md,                              ONLY: mdio_parameters_type,&
                                             simulation_parameters_type
  USE md_conserved_quantities,         ONLY: compute_conserved_quantity

  USE md_environment_types,            ONLY: get_md_env,&
                                             md_environment_type,&
                                             set_md_env
  USE message_passing,                 ONLY: mp_bcast
  USE output_utilities,                ONLY: print_warning
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: write_particle_coordinates,&
                                             write_qs_particle_coordinates
  USE physcon,                         ONLY: femtoseconds,&
                                             kelvin
  USE qs_parser,                       ONLY: open_file
  USE simulation_cell,                 ONLY: cell_type, get_cell_param
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE virial_types,                    ONLY: virial_type

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: module_name = "md_energies"

  TYPE average_quantities_type
     REAL ( dbl ) :: avetemp, avepot, avekin, avevol, aveca, avecb, avecc
     REAL ( dbl ) :: aveal, avebe, avega, avepress, avekinc, avetempc
  END TYPE average_quantities_type

 
  TYPE ( unit_convert_type ), POINTER :: unit

  TYPE (average_quantities_type ) :: averages
  REAL(KIND = dp) :: econs_avg,t_old

  CHARACTER (LEN=50) :: crd_file_name, vel_file_name, ptens_file_name, &
                        ener_file_name, temp_file_name, dump_file_name
  INTEGER :: crd, vel, ptn, ene, tem, scr, coe, pol, group, itimes_old, source
  INTEGER :: icrd, ivel, iptens, iener, itemp, idump, iscreen
  LOGICAL :: ionode

  PUBLIC :: md_energy,&
            set_energies
!!***
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** md_energies/set_energies [1.0] *
!!
!!   NAME
!!     set_energies
!!
!!   FUNCTION
!!     setup of parameters for energies
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!   NOTES
!!
!!*** *************************************************************************

SUBROUTINE set_energies ( globenv, mdio )


    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv
    TYPE(mdio_parameters_type), INTENT(IN)   :: mdio

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
 
  CALL get_cp2k_units ( unit )

END SUBROUTINE set_energies

!******************************************************************************
!!****** md_energies/md_energy [1.0] *
!!
!!   NAME
!!     md_energy
!!
!!   FUNCTION
!!     This routine computes the conserved quantity, temperature
!!     and things like that and prints them out
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     - New MD data are appended to the old data (15.09.2003,MK)
!!
!!*****************************************************************************

  SUBROUTINE md_energy(md_env, error)

    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER              :: routine_name = "md_energy"

    CHARACTER(LEN=200)                       :: title
    CHARACTER(len=50)                        :: filename
    INTEGER                                  :: handle, istat, isubsys, iw, delta_t
    INTEGER, POINTER                         :: itimes
    LOGICAL                                  :: file_exists, replace_file
    REAL(KIND=dp)                            :: cons, e_diff, e_old, ekin, &
                                                ekin_coefs, epot, rcons, t, &
                                                temp, temp_coefs, x
    REAL(KIND=dp), POINTER                   :: constant, used_time
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    TYPE(dyn_coeff_set_p_type), &
      DIMENSION(:), POINTER                  :: dyn_coeff_set
    TYPE(force_env_p_type), POINTER          :: force_env_p ( : )
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar
    TYPE ( virial_type ), POINTER            :: virial
    REAL ( dbl ) :: cell_length ( 3 ), cell_angle ( 3 )
    REAL ( dbl ) :: pv_scalar, econs
    INTEGER :: i

!   ---------------------------------------------------------------------------

    CALL timeset(routine_name,"E","",handle)

    iw = scr

    NULLIFY(dyn_coeff_set)
    CALL get_md_env(md_env=md_env,&
                    itimes=itimes,&
                    constant=constant,&
                    used_time = used_time, &
                    simpar=simpar,&
                    force_env_p=force_env_p )
                    

    CALL force_env_get(force_env=force_env_p(1)%force_env,&
                       cell=cell,&
                       subsys=subsys,&
                       dyn_coeff_set=dyn_coeff_set,&
                       virial = virial, &
                       error=error)

! computing the scalar pressure
      IF ( virial % pv_availability ) THEN
        pv_scalar = 0._dbl
        DO i = 1, 3
          pv_scalar = pv_scalar + virial%pv_total(i,i)
        END DO
        pv_scalar = pv_scalar/3._dbl/cell%deth
        pv_scalar = pv_scalar * unit % bar
      ENDIF


    CALL force_env_get(force_env_p ( 1 ) % force_env, potential_energy=epot, error=error)
    IF (ASSOCIATED (dyn_coeff_set)) THEN
      CALL compute_conserved_quantity(md_env,&
                                      epot,&
                                      cons,&
                                      ekin=ekin,&
                                      temp=temp,&
                                      ekin_coefs=ekin_coefs,&
                                      temp_coefs=temp_coefs)
     temp_coefs=temp_coefs*kelvin
     temp = temp*kelvin

   ELSE
     CALL compute_conserved_quantity(md_env,&
                                     epot,&
                                     cons,&
                                     ekin=ekin,&
                                     temp=temp)
     temp = temp*kelvin

   END IF

  CALL get_cell_param ( cell, cell_length, cell_angle )

    replace_file = .FALSE.

    IF (itimes == 0) THEN
      IF (ionode) THEN
        WRITE (iw,*)
        WRITE ( iw, '( A )' ) ' MD_ENERGIES| Initialization proceeding'
        WRITE (iw,*)
        replace_file = (simpar%read_type == "INIT")
        filename=ener_file_name
        INQUIRE (FILE=filename,EXIST=file_exists)
        CALL open_md_file(ener_file_name,ene,replace_file)
        IF (.NOT.replace_file.AND.file_exists) THEN
          BACKSPACE (UNIT=ene,IOSTAT=istat)
          IF (istat /= 0) THEN
            CALL stop_program(routine_name,module_name,__LINE__,&
                              "BACKSPACE command failed for file <"//&
                              TRIM(ener_file_name)//">")
          END IF
          READ (UNIT=ene,FMT=*,IOSTAT=istat) itimes_old,t_old,x,x,x, &
                                             e_old,e_diff
          IF (istat /= 0) THEN
            CALL print_warning(routine_name,&
                               "Could not read the last line of the file <"//&
                               TRIM(ener_file_name)//">. MD step number "//&
                               "and time are set to zero.")
            itimes_old = 0
            t_old = 0.0_dp
            e_old = 0.0_dp
            e_diff = 0.0_dp
          END IF
        ELSE
          itimes_old = 0
          t_old = 0.0_dp
          e_old = 0.0_dp
          e_diff = 0.0_dp
        END IF
        CALL open_md_file(crd_file_name,crd,replace_file)
        CALL open_md_file(vel_file_name,vel,replace_file)
      END IF
      CALL mp_bcast(itimes_old,source,group)
      CALL mp_bcast(t_old,source,group)
      CALL mp_bcast(e_old,source,group)
      CALL mp_bcast(e_diff,source,group)
      itimes = itimes_old
      t = t_old
      constant = cons
      CALL set_md_env(md_env=md_env,constant=constant)
      IF (ionode) THEN
        WRITE ( iw, '( )' )
        WRITE ( iw, '( A,A )' ) ' ******************************** ', &
           'GO CP2K GO! **********************************'
        WRITE ( iw, '( A,A,T40,A,T60,1(1X,E20.12) )' ) ' INITIAL POTENTIAL ENERGY', &
           adjustl( unit % e_label), '= ', epot
        WRITE ( iw, '( A,A,T40,A,T60,1(1X,E20.12) )' ) ' INITIAL KINETIC ENERGY', &
           adjustl( unit % e_label), '= ', ekin
        WRITE ( iw, '( A,A,T40,A,T60,1(1X,E20.12) )' ) ' INITIAL TEMPERATURE', &
           adjustl( unit % temp_label ), '= ', temp
        IF ( virial % pv_availability ) &
        WRITE ( iw, '( A,A,T40,A,T60,1(1X,E20.12) )' ) ' INITIAL PRESSURE', &
           adjustl( unit % pres_label), '= ', pv_scalar
        WRITE ( iw, '( A,A,T40,A,T60,1(1X,E20.12) )' ) ' INITIAL VOLUME', &
           adjustl( unit % vol_label), '= ', cell%deth
        WRITE ( iw, '( A,A,T29,A,T33,3(1X,E15.7) )' ) ' INITIAL CELL LNTHS', &
           adjustl(unit % l_label), '= ', cell_length(1), cell_length(2), cell_length(3)
        WRITE ( iw, '( A,A,T29,A,T33,3(1X,E15.7) )' ) ' INITIAL CELL ANGLS', &
           adjustl(unit % angl_label), '= ', cell_angle(3), cell_angle(2), cell_angle(1)
        WRITE ( iw, '( A,A )' ) ' ******************************** ', &
           'GO CP2K GO! **********************************'
      ENDIF
    ELSE
      t = t_old + REAL(itimes-itimes_old,dp)*simpar%dt*femtoseconds
      econs = ABS((cons-constant)/constant)
    END IF


    IF (itimes > itimes_old) THEN

      IF (ionode) THEN

        IF (simpar%output_format == "XMOL") THEN
!         *** This information can be digested by Molden ***
          WRITE (UNIT=title,FMT="(A,I8,A,F20.10)")&
            " i =",itimes,", E =",epot
        ELSE
          title = ""
        END IF

! average quantities:
       delta_t = itimes - itimes_old
       SELECT CASE (delta_t)
       CASE (1)
          IF ( virial % pv_availability ) averages%avepress = pv_scalar
          averages%avetemp = temp
          averages%avepot = epot
          averages%avekin = ekin
          averages%avevol = cell % deth
          averages%aveca = cell_length ( 1 )
          averages%avecb = cell_length ( 2 )
          averages%avecc = cell_length ( 3 )
          averages%aveal = cell_angle ( 3 )
          averages%avebe = cell_angle ( 2 )
          averages%avega = cell_angle ( 1 )
          econs_avg =  0._dbl
       CASE DEFAULT
          IF ( virial % pv_availability ) &
          averages%avepress = (averages%avepress*REAL(delta_t-1,dbl)+pv_scalar)/ &
          REAL ( delta_t, dbl )
          averages%avetemp = (averages%avetemp*REAL(delta_t-1,dbl)+temp)/ &
               REAL(delta_t,dbl)
          averages%avepot = (averages%avepot*REAL(delta_t-1,dbl)+epot)/ &
               REAL(delta_t,dbl)
          averages%avekin = (averages%avekin*REAL(delta_t-1,dbl)+ekin)/ &
               REAL(delta_t,dbl)
          averages%avevol = (averages%avevol*REAL(delta_t-1,dbl)+cell%deth)/ &
               REAL(delta_t,dbl)
          averages%aveca = (averages%aveca*REAL(delta_t-1,dbl)+cell_length(1))/REAL(delta_t,dbl)
          averages%avecb = (averages%avecb*REAL(delta_t-1,dbl)+cell_length(2))/REAL(delta_t,dbl)
          averages%avecc = (averages%avecc*REAL(delta_t-1,dbl)+cell_length(3))/REAL(delta_t,dbl)
          averages%aveal = (averages%aveal*REAL(delta_t-1,dbl)+cell_angle(3))/ &
               REAL(delta_t,dbl)
          averages%avebe = (averages%avebe*REAL(delta_t-1,dbl)+cell_angle(2))/ &
               REAL(delta_t,dbl)
          averages%avega = (averages%avega*REAL(delta_t-1,dbl)+cell_angle(1))/ &
          REAL(delta_t,dbl)
          econs_avg = econs_avg + abs((cons-constant)/constant)
       END SELECT
!

       IF ( ( delta_t /= 0 ) .AND. ( MOD ( itimes, iscreen ) == 0 ) ) THEN

          WRITE ( iw, '( )' )
          WRITE ( iw, '( A,A )' ) ' **************************************', &
               '*****************************************'
          WRITE ( iw, '( A,T61,A )' ) ' ENSEMBLE TYPE                = ', &
               ADJUSTR ( simpar%ensemble)
          WRITE ( iw, '( A,T71,I10 )' ) ' STEP NUMBER                  = ', itimes
          WRITE ( iw, '( A,T60,1(1X,E20.12) )' ) ' TIME [FS]                    = ', t
          WRITE ( iw, '( A,T60,1(1X,E20.12) )' ) ' CONSERVED QNTY               = ', &
               cons
          WRITE ( iw, '( )' )
          WRITE ( iw, '( T42,A )' ) 'INSTANTANEOUS        AVERAGES'
          WRITE ( iw, '( A,T39,2(1X,E20.12) )' ) ' {E-E0}/E0                    = ', &
                econs, econs_avg/REAL( delta_t,dbl)
          WRITE ( iw, '( A,A,T31,A,T39,2(1X,E20.12) )' ) ' POTENTIAL ENERGY', &
               adjustl(unit % e_label), '= ', epot, averages%avepot
          WRITE ( iw, '( A,A,T31,A,T39,2(1X,E20.12) )' ) ' KINETIC ENERGY', &
               adjustl(unit % e_label), '= ', ekin, averages%avekin
          WRITE ( iw, '( A,A,T31,A,T39,2(1X,E20.12) )' ) ' TEMPERATURE', &
               adjustl(unit % temp_label), '= ', temp, averages%avetemp
          IF ( virial % pv_availability ) &
          WRITE ( iw, '( A,A,T31,A,T39,2(1X,E20.12) )' ) ' PRESSURE', &
               adjustl(unit % pres_label), '= ', pv_scalar, averages%avepress
          IF (simpar%ensemble(1:3)=='NPT') THEN
             WRITE ( iw, '( A,A,T31,A,T39,2(1X,E20.12) )' ) ' VOLUME', &
                  adjustl(unit % vol_label), '= ', cell%deth, averages%avevol
             WRITE ( iw, '( A,A,T31,A,T33,3(1X,E15.7) )' ) ' CELL LNTHS', &
                  adjustl(unit % l_label), '= ', cell_length(1), cell_length(2), cell_length(3)
             WRITE ( iw, '( A,A,T31,A,T33,3(1X,E15.7) )' ) ' AVE. CELL LNTHS', &
                  adjustl(unit % l_label), '= ', averages%aveca, averages%avecb, &
                  averages%avecc
          END IF
          IF (simpar%ensemble=='NPT_F') THEN
             WRITE ( iw, '( A,A,T31,A,T33,3(1X,E15.7) )' ) ' CELL ANGLS', &
                  adjustl(unit % angl_label), '= ', cell_angle(3), cell_angle(2), cell_angle(1)
             WRITE ( iw, '( A,A,T31,A,T33,3(1X,E15.7) )' ) ' AVE. CELL ANGLS', &
                  adjustl(unit % angl_label), '= ', averages%aveal, averages%avebe, &
                  averages%avega
          END IF
          WRITE ( iw, '( A,A )' ) ' **************************************', &
               '*****************************************'
          WRITE ( iw, '( )' )
       END IF

!dbg
!        IF (MODULO(itimes,iscreen) == 0) THEN
!          IF (ASSOCIATED (dyn_coeff_set)) THEN
!            WRITE (scr, "(/,T2,A6,I8,F12.3,F16.9,F16.9,F16.9,F16.9,F16.9,F16.9,F20.12,F12.6/)")&
!              "-MD-",itimes,t,ekin,ekin_coefs,temp,temp_coefs,epot,&
!              cons,cons-constant, used_time
!          ELSE
!            WRITE (scr, "(/,T2,A6,I8,F12.3,F20.9,F20.9,F20.9,F20.9,F20.12,F20.6/)")&
!              "-MD-",itimes,t,ekin,temp,epot,&
!              cons,cons-constant, used_time
!          END IF
!        END IF
!dbg
!MK
!        IF (MODULO(itimes,iscreen) == 0) THEN
!          DO isubsys=1,SIZE(subsys)
!            CALL cp_subsys_get(subsys(isubsys)%subsys,&
!                               particles=particles,&
!                               error=error)
!            CALL write_qs_particle_coordinates(particle_set=particles%els,&
!                                               cell=cell,&
!                                               error=error)
!          END DO
!        END IF
!MK
        IF (MODULO(itimes,iener) == 0) THEN
          WRITE (ene,"(I8,F12.3,F20.9,F20.9,F20.9,F20.9,F20.12)")&
            itimes,t,ekin,temp,epot,&
            econs, econs_avg/REAL(delta_t,dbl)
          CALL m_flush(ene)
        END IF

        IF (MODULO(itimes,icrd) == 0) THEN
          DO isubsys=1,SIZE(subsys)
            CALL cp_subsys_get(subsys(isubsys)%subsys,&
                              particles=particles,&
                              error=error)
            CALL write_particle_coordinates(particles%els,cell,crd,&
                                            simpar%output_format,"POS",title)
          END DO
          CALL m_flush(crd)
        END IF

        IF (MODULO(itimes,ivel) == 0) THEN
          DO isubsys=1,SIZE(subsys)
            CALL cp_subsys_get(subsys(isubsys)%subsys,&
                              particles=particles,&
                              error=error)
            CALL write_particle_coordinates(particles%els,cell,vel,&
                                            simpar%output_format,"POS_VEL_FORCE",title)
          END DO
          CALL m_flush(vel)
        END IF

      END IF

    END IF

    CALL timestop(zero,handle)

  END SUBROUTINE md_energy

! *****************************************************************************

  SUBROUTINE open_md_file(file_name,unit_number,replace_file)

    CHARACTER(LEN=*), INTENT(IN)             :: file_name
    INTEGER, INTENT(OUT)                     :: unit_number
    LOGICAL, INTENT(IN), OPTIONAL            :: replace_file

    LOGICAL                                  :: file_exists, keep_file

!   ---------------------------------------------------------------------------

    INQUIRE (FILE=file_name,EXIST=file_exists)

    IF (PRESENT(replace_file)) THEN
      keep_file = (.NOT.replace_file.AND.file_exists)
    ELSE
      keep_file = file_exists
    END IF

    IF (keep_file) THEN
      CALL open_file(file_name=TRIM(file_name),&
                     file_action="READWRITE",&
                     file_form="FORMATTED",&
                     file_position="APPEND",&
                     file_status="OLD",&
                     unit_number=unit_number)
    ELSE
      CALL open_file(file_name=TRIM(file_name),&
                     file_action="READWRITE",&
                     file_form="FORMATTED",&
                     file_position="REWIND",&
                     file_status="REPLACE",&
                     unit_number=unit_number)
    END IF

  END SUBROUTINE open_md_file

! *****************************************************************************

END MODULE md_energies
