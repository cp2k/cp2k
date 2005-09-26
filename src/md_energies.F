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
  USE kinds,                           ONLY: dp
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
  USE cp_files,                        ONLY: open_file
  USE simulation_cell,                 ONLY: cell_type, get_cell
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE virial_types,                    ONLY: virial_type
  USE input_constants,                 ONLY: nve_ensemble,&
                                             nvt_ensemble,&
                                             npt_i_ensemble,&
                                             npt_f_ensemble,&
                                             nph_ensemble,&
                                             nph_uniaxial_ensemble,&
                                             nph_uniaxial_damped_ensemble,&
                                             isokin_ensemble, &
                                             ensemble_name_types,&
                                             md_init_start,&
                                             md_pos_restart,&
                                             md_pos_vel_restart,&
                                             dump_xmol
  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: moduleN = "md_energies"

  TYPE average_quantities_type
     REAL ( dp ) :: avetemp, avepot, avekin, avevol, aveca, avecb, avecc
     REAL ( dp ) :: avetemp_baro, avehugoniot, avecpu
     REAL ( dp ) :: aveal, avebe, avega, avepress, avekinc, avetempc, avepxx
  END TYPE average_quantities_type


  PUBLIC :: md_energy

!!***
!******************************************************************************

CONTAINS

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

  SUBROUTINE md_energy(md_env, globenv, error)

    TYPE ( md_environment_type ), POINTER       :: md_env
    TYPE ( global_environment_type ), POINTER   :: globenv
    TYPE ( cp_error_type ), INTENT ( inout ), &
      OPTIONAL                                  :: error

! locals
    CHARACTER(LEN=*), PARAMETER              :: routineN = "md_energy"

    CHARACTER (LEN=200)                      :: title
    CHARACTER (len=50)                       :: filename
    CHARACTER (LEN=50)                       :: crd_file_name, vel_file_name, ptens_file_name, &
                                                ener_file_name, temp_file_name, dump_file_name
    INTEGER, SAVE                             :: crd, vel, ene 
    INTEGER                                  :: handle, istat, isubsys, iw, delta_t, i
    INTEGER                                  :: group, source
    INTEGER                                  :: icrd, ivel, iptens, iener, itemp, idump, iscreen, nat
    INTEGER, POINTER                         :: itimes
    LOGICAL                                  :: file_exists, replace_file
    LOGICAL                                  :: ionode
    REAL(KIND = dp), SAVE                    :: econs_avg, t_old, itimes_old
    REAL ( dp )                              :: abc ( 3 ), cell_angle ( 3 )
    REAL ( dp )                              :: pv_scalar, econs, pv_xx
    REAL(KIND=dp)                            :: cons, ekin, &
                                                ekin_coefs, epot, rcons, t, &
                                                temp, temp_coefs, temp_baro, &
                                                hugoniot
    REAL(KIND=dp), POINTER                   :: constant, used_time
    TYPE ( unit_convert_type ), POINTER      :: unit
    TYPE(cell_type), POINTER                 :: cell, cell_ref
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    TYPE(dyn_coeff_set_p_type), &
      DIMENSION(:), POINTER                  :: dyn_coeff_set
    TYPE(force_env_p_type), POINTER          :: force_env_p ( : )
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar
    TYPE(mdio_parameters_type), &
      POINTER                                :: mdio
    TYPE ( virial_type ), POINTER            :: virial
    TYPE (average_quantities_type ), SAVE    :: averages

!   ---------------------------------------------------------------------------

    CALL timeset(routineN,"E","",handle)

    CALL get_cp2k_units ( unit )

! Zeroing
    hugoniot = 0._dp


    NULLIFY(dyn_coeff_set)
    CALL get_md_env(md_env=md_env,&
                    itimes=itimes,&
                    constant=constant,&
                    used_time = used_time, &
                    simpar=simpar,         &
                    mdio = mdio,           &
                    force_env_p=force_env_p )
                    

!------------------------------------------------------------------------------

    group = globenv % group
    source = globenv  % source
    ionode = globenv % ionode
    iw = globenv % scr

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

    CALL force_env_get(force_env=force_env_p(1)%force_env,&
                       cell=cell,&
                       cell_ref=cell_ref,&
                       subsys=subsys,&
                       dyn_coeff_set=dyn_coeff_set,&
                       virial = virial, &
                       error=error)

    nat=0
    DO isubsys=1,SIZE(subsys)
       CALL cp_subsys_get(subsys(isubsys)%subsys,&
            particles=particles,&
            error=error)
       nat=nat+particles%n_els
    END DO

! computing the scalar pressure
      IF ( virial % pv_availability ) THEN
        pv_scalar = 0._dp
        DO i = 1, 3
          pv_scalar = pv_scalar + virial%pv_total(i,i)
        END DO
        pv_scalar = pv_scalar/3._dp/cell%deth
        pv_scalar = pv_scalar * unit % bar
        pv_xx = virial%pv_total(1,1)/cell%deth*unit%bar
      ENDIF

     IF (itimes == 0) THEN
       IF ( simpar % ensemble == nph_ensemble .OR.&
            simpar % ensemble == nph_uniaxial_ensemble .OR.&
            simpar % ensemble == nph_uniaxial_damped_ensemble ) THEN
         IF ( simpar % v0 == 0._dp ) simpar % v0 = cell_ref % deth
       ENDIF
     ENDIF

    CALL force_env_get(force_env_p ( 1 ) % force_env, potential_energy=epot, error=error)
    IF (ASSOCIATED (dyn_coeff_set)) THEN
      CALL compute_conserved_quantity(md_env,&
                                      epot,&
                                      cons,&
                                      ekin=ekin,&
                                      temp_part=temp,&
                                      temp_baro=temp_baro,&
                                      ekin_coefs=ekin_coefs,&
                                      temp_coefs=temp_coefs)
     temp_coefs=temp_coefs*kelvin
     temp = temp*kelvin
     temp_baro = temp_baro*kelvin

   ELSE
     CALL compute_conserved_quantity(md_env,&
                                     epot,&
                                     cons,&
                                     ekin=ekin,&
                                     temp_part=temp,&
                                     temp_baro=temp_baro )
     temp = temp*kelvin
     temp_baro = temp_baro*kelvin

   END IF
 
  CALL get_cell ( cell, abc = abc,  &
                  alpha = cell_angle ( 3 ), &
                  beta = cell_angle ( 2 ), &
                  gamma = cell_angle ( 1 ) )

  replace_file = .FALSE.

  IF (itimes == 0) THEN
     IF ( simpar % ensemble == nph_ensemble .OR.&
          simpar % ensemble == nph_uniaxial_ensemble .OR.&
          simpar % ensemble == nph_uniaxial_damped_ensemble ) THEN
      IF ( simpar % e0 == 0._dp ) simpar % e0 = epot + ekin
    ENDIF
    IF (ionode) THEN
        WRITE (iw,*)
        WRITE ( iw, '( A )' ) ' MD_ENERGIES| Initialization proceeding'
        WRITE (iw,*)
        replace_file = (simpar%read_type == md_init_start)
        filename=ener_file_name
        INQUIRE (FILE=filename,EXIST=file_exists)
        CALL open_md_file(ener_file_name,ene,replace_file)
        WRITE (iw,*) "Format of the energy file"
        WRITE (iw,*) "Step Nr.    Time[fs]  Kin.[a.u.]     Temp[K]  Pot.[a.u.]    Cons Qty UsedTime[s]"
        IF ((.NOT.replace_file).AND.file_exists) THEN
          BACKSPACE (UNIT=ene,IOSTAT=istat)
          IF (istat /= 0) THEN
            CALL stop_program(routineN,moduleN,__LINE__,&
                              "BACKSPACE command failed for file <"//&
                              TRIM(ener_file_name)//">")
          END IF
          READ (UNIT=ene,FMT=*,IOSTAT=istat) itimes_old,t_old
                                             
          IF (istat /= 0) THEN
            CALL print_warning(routineN,&
                               "Could not read the last line of the file <"//&
                               TRIM(ener_file_name)//">. MD step number "//&
                               "and time are set to zero.")
            itimes_old = 0
            t_old = 0.0_dp
          END IF
        ELSE
          itimes_old = 0
          t_old = 0.0_dp
        END IF
        CALL open_md_file(crd_file_name,crd,replace_file)
        CALL open_md_file(vel_file_name,vel,replace_file)
    END IF
    CALL mp_bcast(itimes_old,source,group)
    CALL mp_bcast(t_old,source,group)
    itimes = itimes_old
    t = t_old
    constant = cons
    CALL set_md_env(md_env=md_env,constant=constant)
    IF ( simpar % ensemble == nph_ensemble .OR.&
         simpar % ensemble == nph_uniaxial_ensemble .OR.&
         simpar % ensemble == nph_uniaxial_damped_ensemble ) THEN
        hugoniot = epot + ekin - simpar % e0 -  &
                  .5_dp * ( pv_xx/unit%bar + simpar % p0 ) * &
                  ( simpar % v0 - cell % deth )
    ENDIF
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
        IF ( simpar % ensemble == nph_ensemble .OR.&
             simpar % ensemble == nph_uniaxial_ensemble .OR.&
             simpar % ensemble == nph_uniaxial_damped_ensemble.OR.&
             simpar % ensemble == npt_i_ensemble.OR.&
             simpar % ensemble == npt_f_ensemble ) &
        WRITE ( iw, '( A,A,T40,A,T60,1(1X,E20.12) )' ) ' INITIAL BAROSTAT TEMP', &
           adjustl( unit % temp_label), '= ', temp_baro
        IF ( virial % pv_availability ) &
        WRITE ( iw, '( A,A,T40,A,T60,1(1X,E20.12) )' ) ' INITIAL PRESSURE', &
           adjustl( unit % pres_label), '= ', pv_scalar
        IF ( simpar % ensemble == nph_ensemble .OR.&
             simpar % ensemble == nph_uniaxial_ensemble .OR.&
             simpar % ensemble == nph_uniaxial_damped_ensemble)&
        WRITE ( iw, '( A,A,T40,A,T60,1(1X,E20.12) )' ) ' INITIAL HUGONIOT CONSTRAINT', &
           adjustl( unit % e_label), '= ', hugoniot
        WRITE ( iw, '( A,A,T40,A,T60,1(1X,E20.12) )' ) ' INITIAL VOLUME', &
           adjustl( unit % vol_label), '= ', cell%deth
        WRITE ( iw, '( A,A,T29,A,T33,3(1X,E15.7) )' ) ' INITIAL CELL LNTHS', &
           adjustl(unit % l_label), '= ', abc(1), abc(2), abc(3)
        WRITE ( iw, '( A,A,T29,A,T33,3(1X,E15.7) )' ) ' INITIAL CELL ANGLS', &
           adjustl(unit % angl_label), '= ', cell_angle(3), cell_angle(2), cell_angle(1)
        WRITE ( iw, '( A,A )' ) ' ******************************** ', &
           'GO CP2K GO! **********************************'
    ENDIF
  ELSE
     IF ( simpar % ensemble == nph_ensemble .OR.&
          simpar % ensemble == nph_uniaxial_ensemble .OR.&
          simpar % ensemble == nph_uniaxial_damped_ensemble)&
          hugoniot = epot + ekin - simpar % e0 -  &
              .5_dp * ( pv_xx/unit%bar + simpar % p0 ) * &
              ( simpar % v0 - cell % deth )
     t = t_old + REAL(itimes-itimes_old,dp)*simpar%dt*femtoseconds
     econs = ABS((cons-constant)/nat*kelvin)
  END IF


    IF (itimes > itimes_old) THEN

      IF (ionode) THEN

        IF (simpar%output_format == dump_xmol) THEN
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
          IF ( virial % pv_availability ) averages%avepxx = pv_xx
          averages%avecpu = used_time
          averages%avehugoniot = hugoniot
          averages%avetemp = temp
          averages%avetemp_baro = temp_baro
          averages%avepot = epot
          averages%avekin = ekin
          averages%avevol = cell % deth
          averages%aveca = abc ( 1 )
          averages%avecb = abc ( 2 )
          averages%avecc = abc ( 3 )
          averages%aveal = cell_angle ( 3 )
          averages%avebe = cell_angle ( 2 )
          averages%avega = cell_angle ( 1 )
          econs_avg =  0._dp
       CASE DEFAULT
          IF ( virial % pv_availability ) &
          averages%avepress = (averages%avepress*REAL(delta_t-1,dp)+pv_scalar)/ &
          REAL ( delta_t, dp )
          IF ( virial % pv_availability ) &
          averages%avepxx = (averages%avepxx*REAL(delta_t-1,dp)+pv_xx)/ &
          REAL ( delta_t, dp )
          averages%avetemp = (averages%avetemp*REAL(delta_t-1,dp)+temp)/ &
               REAL(delta_t,dp)
          averages%avecpu = (averages%avecpu*REAL(delta_t-1,dp)+used_time)/ &
               REAL(delta_t,dp)
          averages%avehugoniot = (averages%avehugoniot*REAL(delta_t-1,dp)+hugoniot)/ &
               REAL(delta_t,dp)
          averages%avetemp_baro= (averages%avetemp_baro*REAL(delta_t-1,dp)+temp_baro)/ &
               REAL(delta_t,dp)
          averages%avepot = (averages%avepot*REAL(delta_t-1,dp)+epot)/ &
               REAL(delta_t,dp)
          averages%avekin = (averages%avekin*REAL(delta_t-1,dp)+ekin)/ &
               REAL(delta_t,dp)
          averages%avevol = (averages%avevol*REAL(delta_t-1,dp)+cell%deth)/ &
               REAL(delta_t,dp)
          averages%aveca = (averages%aveca*REAL(delta_t-1,dp)+abc(1))/REAL(delta_t,dp)
          averages%avecb = (averages%avecb*REAL(delta_t-1,dp)+abc(2))/REAL(delta_t,dp)
          averages%avecc = (averages%avecc*REAL(delta_t-1,dp)+abc(3))/REAL(delta_t,dp)
          averages%aveal = (averages%aveal*REAL(delta_t-1,dp)+cell_angle(3))/ &
               REAL(delta_t,dp)
          averages%avebe = (averages%avebe*REAL(delta_t-1,dp)+cell_angle(2))/ &
               REAL(delta_t,dp)
          averages%avega = (averages%avega*REAL(delta_t-1,dp)+cell_angle(1))/ &
          REAL(delta_t,dp)
          econs_avg = econs_avg + (cons-constant)/nat*kelvin
       END SELECT
!

       IF ( ( delta_t /= 0 ) .AND. ( MOD ( itimes, iscreen ) == 0 ) ) THEN

          WRITE ( iw, '( )' )
          WRITE ( iw, '( A,A )' ) ' **************************************', &
               '*****************************************'
          WRITE ( iw, '( A,T61,A )' ) ' ENSEMBLE TYPE                = ', &
               ADJUSTR ( ensemble_name_types(simpar%ensemble)(1:20))
          WRITE ( iw, '( A,T71,I10 )' ) ' STEP NUMBER                  = ', itimes
          WRITE ( iw, '( A,T60,1(1X,F20.12) )' ) ' TIME [FS]                    = ', t
          WRITE ( iw, '( A,T60,1(1X,E20.12) )' ) ' CONSERVED QNTY               = ', &
               cons
          WRITE ( iw, '( )' )
          WRITE ( iw, '( T42,A )' ) 'INSTANTANEOUS        AVERAGES'
          WRITE ( iw, '( A,T39,2(1X,F20.2) )' ) ' CPU [S]                      = ', &
                used_time, averages%avecpu
          WRITE ( iw, '( A,T39,2(1X,E20.12) )' ) ' {E-E0}/{k_b*N_at}            = ', &
                econs, econs_avg/REAL( delta_t,dp)
          WRITE ( iw, '( A,A,T31,A,T39,2(1X,E20.12) )' ) ' POTENTIAL ENERGY', &
               adjustl(unit % e_label), '= ', epot, averages%avepot
          WRITE ( iw, '( A,A,T31,A,T39,2(1X,E20.12) )' ) ' KINETIC ENERGY', &
               adjustl(unit % e_label), '= ', ekin, averages%avekin
          WRITE ( iw, '( A,A,T31,A,T39,2(1X,E20.12) )' ) ' TEMPERATURE', &
               adjustl(unit % temp_label), '= ', temp, averages%avetemp
          IF ( virial % pv_availability ) &
          WRITE ( iw, '( A,A,T31,A,T39,2(1X,E20.12) )' ) ' PRESSURE', &
               adjustl(unit % pres_label), '= ', pv_scalar, averages%avepress
          IF ( simpar % ensemble == nph_ensemble .OR.&
               simpar % ensemble == nph_uniaxial_ensemble .OR.&
               simpar % ensemble == nph_uniaxial_damped_ensemble) THEN
            WRITE ( iw, '( A,A,T31,A,T39,2(1X,E20.12) )' ) ' P_XX', &
               adjustl(unit % pres_label), '= ', pv_xx, averages%avepxx
            WRITE ( iw, '( A,A,T31,A,T39,2(1X,E20.12) )' ) ' HUGONIOT', &
               adjustl(unit % e_label), '= ', hugoniot, averages%avehugoniot
          END IF
        IF ( simpar % ensemble == nph_ensemble .OR.&
             simpar % ensemble == nph_uniaxial_ensemble .OR.&
             simpar % ensemble == nph_uniaxial_damped_ensemble.OR.&
             simpar % ensemble == npt_i_ensemble.OR.&
             simpar % ensemble == npt_f_ensemble ) THEN
             WRITE ( iw, '( A,A,T31,A,T39,2(1X,E20.12) )' ) ' BAROSTAT TEMP', &
               adjustl(unit % temp_label), '= ', temp_baro, averages%avetemp_baro
             WRITE ( iw, '( A,A,T31,A,T39,2(1X,E20.12) )' ) ' VOLUME', &
                  adjustl(unit % vol_label), '= ', cell%deth, averages%avevol
             WRITE ( iw, '( A,A,T31,A,T33,3(1X,E15.7) )' ) ' CELL LNTHS', &
                  adjustl(unit % l_label), '= ', abc(1), abc(2), abc(3)
             WRITE ( iw, '( A,A,T31,A,T33,3(1X,E15.7) )' ) ' AVE. CELL LNTHS', &
                  adjustl(unit % l_label), '= ', averages%aveca, averages%avecb, &
                  averages%avecc
          END IF
          IF (simpar%ensemble==npt_f_ensemble) THEN
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

        IF (MODULO(itimes,iener) == 0) THEN
          ! please change also the corresponding format explaination above
          ! keep the constant of motion the true constant of motion !
          ! "Step Nr.    Time[fs]  Kin.[a.u.]     Temp[K]  Pot.[a.u.]    Cons Qty UsedTime[s]"
          !
          WRITE (ene,"(I8,F12.3,F20.9,F20.9,F20.9,F20.9,F20.9)")&
            itimes,t,ekin,temp,epot,cons,used_time
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

    ! XXX are these files ever closed ???
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
