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
  USE dynamical_coeff_types,           ONLY: dyn_coeff_set_p_type
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp
  USE machine,                         ONLY: m_flush
  USE mathconstants,                   ONLY: zero
  USE md,                              ONLY: mdio_parameters_type,&
                                             simulation_parameters_type
  USE md_conserved_quantities,         ONLY: compute_conserved_quantity
!compute_conserved_quantity_pimd

  USE md_environment_types,            ONLY: get_md_env,&
                                             md_environment_type,&
                                             set_md_env
  USE message_passing,                 ONLY: mp_bcast
  USE output_utilities,                ONLY: print_warning
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: write_particle_coordinates
  USE physcon,                         ONLY: femtoseconds,&
                                             kelvin
  USE qs_parser,                       ONLY: open_file
  USE simulation_cell,                 ONLY: cell_type
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: module_name = "md_energies"

  REAL(KIND = dp) :: econs_avg,t_old

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
!!****** md_energies/set_qs_energies [1.0] *
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

END SUBROUTINE set_qs_energies

!******************************************************************************
!!****** md_energies/md_energy_qs [1.0] *
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
!!     - New MD data are appended to the old data (15.09.2003,MK)
!!
!!*****************************************************************************

  SUBROUTINE md_energy_qs(md_env, error)

    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER              :: routine_name = "md_energy_qs"

    CHARACTER(LEN=200)                       :: title
    INTEGER                                  :: handle, istat, isubsys, iw
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
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar

!   ---------------------------------------------------------------------------

    CALL timeset(routine_name,"E","",handle)

    NULLIFY(dyn_coeff_set)
    CALL get_md_env(md_env=md_env,&
                    itimes=itimes,&
                    constant=constant,&
                    used_time = used_time, &
                    simpar=simpar,&
                    force_env=force_env)

    CALL force_env_get(force_env=force_env,&
                       cell=cell,&
                       subsys=subsys,&
                       dyn_coeff_set=dyn_coeff_set,&
                       error=error)

    replace_file = .FALSE.

    IF (itimes == 0) THEN
      IF (ionode) THEN
        replace_file = (simpar%read_type == "INIT")
        INQUIRE (FILE=TRIM(ener_file_name),EXIST=file_exists)
        CALL open_md_file(ener_file_name,ene,replace_file)
        IF (.NOT.replace_file.AND.file_exists) THEN
          BACKSPACE (UNIT=ene,IOSTAT=istat)
          IF (istat /= 0) THEN
            CALL stop_program(routine_name,module_name,__LINE__,&
                              "BACKSPACE command failed for file <"//&
                              TRIM(ener_file_name)//">")
          END IF
          READ (UNIT=ene,FMT=*,IOSTAT=istat) itimes_old,t_old,x,x,x,e_old,e_diff
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
      constant = e_old - e_diff
      CALL set_md_env(md_env=md_env,constant=constant)
    ELSE
      t = t_old + REAL(itimes-itimes_old,dp)*simpar%dt*femtoseconds
    END IF

    IF (itimes > itimes_old) THEN

       iw = scr
       CALL force_env_get(force_env,&
                          potential_energy=epot,&
                          error=error)
       IF (ASSOCIATED (dyn_coeff_set)) THEN
         CALL compute_conserved_quantity(md_env,&
                                         epot,&
                                         cons,&
                                         cons_rel=rcons,&
                                         ekin=ekin,&
                                         temp=temp,&
                                         ekin_coefs=ekin_coefs,&
                                         temp_coefs=temp_coefs)
        temp_coefs=temp_coefs*kelvin
        temp = temp*kelvin

        WRITE (title,"(A4,I8,A4,F12.3,1X,A6,F20.9,1X,A6,F9.3,1X,A6,F9.3,1X,A1)")&
          " i =",itimes," t =",t,"fs E =",epot,&
          "au T= ",temp,"K Tc=",temp_coefs,"K"
      ELSE
        CALL compute_conserved_quantity(md_env,&
                                        epot,&
                                        cons,&
                                        cons_rel=rcons,&
                                        ekin=ekin,&
                                        temp=temp)
        temp = temp*kelvin

        WRITE (title,"(A4,I8,A4,F12.3,1X,A6,F20.9,1X,A6,F9.3,1X,A1)")&
          " i =",itimes," t =",t,"fs E =",epot,&
          "au T= ",temp,"K"
      END IF


      IF (ionode) THEN

        IF (MODULO(itimes,iscreen) == 0) THEN
          IF (ASSOCIATED (dyn_coeff_set)) THEN
            WRITE (scr, "(/,T2,A6,I8,F12.3,F16.9,F16.9,F16.9,F16.9,F16.9,F16.9,F20.12,F12.6/)")&
              "-MD-",itimes,t,ekin,ekin_coefs,temp,temp_coefs,epot,&
              cons,cons-constant, used_time
          ELSE
            WRITE (scr, "(/,T2,A6,I8,F12.3,F20.9,F20.9,F20.9,F20.9,F20.12,F20.6/)")&
              "-MD-",itimes,t,ekin,temp,epot,&
              cons,cons-constant, used_time
          END IF
        END IF

        IF (MODULO(itimes,iener) == 0) THEN
          WRITE (ene,"(I8,F12.3,F20.9,F20.9,F20.9,F20.9,F20.12)")&
            itimes,t,ekin,temp,epot,&
            cons,cons-constant
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

  END SUBROUTINE md_energy_qs

! *****************************************************************************

  SUBROUTINE open_md_file(file_name,unit_number,replace_file)

    CHARACTER(LEN=*), INTENT(IN)             :: file_name
    INTEGER, INTENT(OUT)                     :: unit_number
    LOGICAL, INTENT(IN), OPTIONAL            :: replace_file

    LOGICAL                                  :: file_exists, keep_file

!   ---------------------------------------------------------------------------

    INQUIRE (FILE=TRIM(file_name),EXIST=file_exists)

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
