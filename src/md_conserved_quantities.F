!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000-2003  CP2K developers group                            !
!-----------------------------------------------------------------------------!
!!****** cp2k/md_conserved_quantities [1.0] *
!!
!!   NAME
!!     conserved_quantities
!!
!!   FUNCTION
!!     computes the conserved quantities for a given md ensemble 
!!     and also kinetic energies, thermo/barostat stuff
!!
!!   AUTHOR
!!     gtb, 05.02.2003 
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

MODULE md_conserved_quantities
  
  USE atomic_kind_list_types,          ONLY: atomic_kind_list_type
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsystem_types,              ONLY: cp_subsys_get,&
                                             cp_subsystem_p_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE dynamical_coeff_types,           ONLY: dyn_coeff_mass_global,&
                                             dyn_coeff_set_p_type,&
                                             dyn_coeff_set_type,&
                                             dyn_coeff_type
  USE extended_system_types,           ONLY: lnhc_parameters_type,&
                                             npt_info_type
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_p_type
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: zero
  USE md,                              ONLY: simulation_parameters_type
  USE md_environment_types,            ONLY: get_md_env,&
                                             md_environment_type,&
                                             set_md_env
  USE message_passing,                 ONLY: mp_sum
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: particle_type
  USE simulation_cell,                 ONLY: cell_type
  USE termination,                     ONLY: stop_program

  IMPLICIT NONE

  PRIVATE

  TYPE md_ener_type
  REAL(KIND=dp) :: epot  
  REAL(KIND=dp) :: ekin  
  REAL(KIND=dp) :: temp_part
  REAL(KIND=dp) :: temp_nose  
  REAL(KIND=dp) :: temp_baro  
  REAL(KIND=dp) :: ekin_coefs
  REAL(KIND=dp) :: temp_coefs  
  END TYPE md_ener_type


  PUBLIC :: compute_conserved_quantity
           ! compute_conserved_quantity_pimd                                    


!******************************************************************************
CONTAINS

!!*****************************************************************************
!!******
!!
!!   NAME
!!     compute_conserved_quantity
!!
!!   FUNCTION
!!     calculates conserved quantity.                                  
!!
!!   INPUT ARGUMENTS    
!!    md_env is the md_environment 
!!    epot is the total potential energy
!!   OUTPUT ARGUMENTS
!!    cons is the conserved quantity
!!   OUTPUT OPTIONAL ARGUMENTS 
!!    cons_rel : relative cons. quantity (to the first md step)
!!    ekin : kinetic energy of particles
!!    temp : temperature
!!    ekin_coefs : ke of the coefficients (KG code)
!!    temp_coefs : temperature of the coefficients (KG code)
!!
!!   AUTHOR
!!    gloria
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** *************************************************************************

SUBROUTINE compute_conserved_quantity (md_env, epot, cons, cons_rel, ekin,&
                                       temp_part, temp_baro, temp_nose, ekin_coefs, &
                                       temp_coefs )                            


    TYPE(md_environment_type), POINTER       :: md_env
    REAL(KIND=dp), INTENT(IN)                    :: epot
    REAL(KIND=dp), INTENT(OUT)                   :: cons
    REAL(KIND=dp), INTENT(OUT), OPTIONAL         :: cons_rel, ekin, temp_part, &
                                                ekin_coefs, temp_coefs, temp_baro, &
                                                temp_nose

    INTEGER                                  :: int_grp
    INTEGER, POINTER                         :: itimes
    REAL(KIND=dp), POINTER                       :: constant
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(md_ener_type)                       :: md_ener
    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar
    CHARACTER(LEN=*), PARAMETER              :: RoutineN="compute_conserved_quantity"

!------------------------------------------------------------------------------
    md_ener%ekin        = 0.0_dp 
    md_ener%epot        = 0.0_dp
    md_ener%temp_part   = 0.0_dp
    md_ener%temp_baro   = 0.0_dp
    md_ener%ekin_coefs  = 0.0_dp
    md_ener%temp_coefs  = 0.0_dp


  md_ener%epot=epot 

  CALL get_md_env (md_env=md_env,&
                   constant=constant,&
                   itimes=itimes,&
                   simpar=simpar,&
                   para_env = para_env )

  int_grp =  para_env % group  
 
  CALL get_part_ke ( md_env, md_ener % ekin, int_grp )

  md_ener % temp_part = 2.0_dp * md_ener % ekin / REAL ( simpar%nfree,KIND=dp)
  
  SELECT CASE ( TRIM(simpar % ensemble)  )
  CASE DEFAULT
    ! hence, should be implemented 
    CALL stop_program( RoutineN, 'Unknown Ensemble')
  CASE ( 'ISOKIN' )
    cons = md_ener % ekin
  CASE ( 'NVE')
    CALL get_econs_nve (md_env, md_ener, cons, int_grp)
  CASE ( 'NVT')
    CALL get_econs_nvt (md_env, md_ener, cons, int_grp)
  CASE ( 'NPT_I', 'NPT_F' )
    CALL get_econs_npt (md_env, md_ener, cons, int_grp)
  CASE ( 'NPH_UNIAXIAL' )
    CALL get_econs_nph_uniaxial (md_env, md_ener, cons, int_grp)
  CASE ( 'NPH_UNIAXIAL_DAMPED' )
    CALL get_econs_nph_damped (md_env, md_ener, cons, int_grp)
  END SELECT

  IF (itimes==1) THEN
     constant = cons
     CALL set_md_env (md_env=md_env, constant=constant)
  END IF

  IF ((PRESENT(cons_rel)).AND.(itimes /= 0 )) THEN
    cons_rel = ABS((cons-constant)/constant)
  END IF
  IF (PRESENT(ekin)) ekin=md_ener%ekin 
  IF (PRESENT(temp_part)) temp_part=md_ener%temp_part
  IF (PRESENT(temp_baro)) temp_baro=md_ener%temp_baro
  IF (PRESENT(ekin_coefs)) ekin_coefs=md_ener%ekin_coefs 
  IF (PRESENT(temp_coefs)) temp_coefs=md_ener%temp_coefs 

END SUBROUTINE compute_conserved_quantity

!*****************************************************************************

!SUBROUTINE compute_conserved_quantity_pimd (md_env, epot, cons,&
!                                            cons_rel, ekin, temp)

!  IMPLICIT NONE

! Arguments
!  TYPE (md_environment_type), POINTER :: md_env
!  REAL(KIND=dp), INTENT (INOUT) :: epot
!  REAL(KIND=dp), INTENT (OUT):: cons
!  REAL(KIND=dp), INTENT (OUT):: cons_rel
!  REAL(KIND=dp), INTENT (OUT):: ekin
!  REAL(KIND=dp), INTENT (OUT):: temp

!Locals
!  TYPE (simulation_parameters_type), POINTER :: simpar
!  TYPE (structure_type), POINTER :: struc ( : )

!  TYPE (cell_type), POINTER :: box
!  TYPE (lnhc_parameters_type), POINTER :: nhc_baro ( : )
!  TYPE (lnhc_parameters_type), POINTER :: nhc_part ( : )
!  TYPE (npt_info_type), POINTER :: npt ( :, : )

!  INTEGER, POINTER :: itimes
!  REAL(KIND=dp), POINTER :: constant
!  TYPE (atomic_kind_type), POINTER :: atomic_kind

!  INTEGER i,j,nnodes,beads                             
!  REAL(KIND=dp) :: ke, t, atmp, btmp, mass 
!  REAL(KIND=dp) :: nhc_part_kin, nhc_part_pot, nhc_baro_kin 
!  REAL(KIND=dp) :: nhc_baro_pot, baro_kin, baro_pot 
!  TYPE(cp_para_env_type), POINTER :: para_env

!------------------------------------------------------------------------------

!  CALL get_md_env (md_env=md_env,&
!                   cell=box,&
!                   constant=constant,&
!                   itimes=itimes,&
!                   nhc_baro=nhc_baro,&
!                   nhc_part=nhc_part,&
!                   npt=npt,&
!                   simpar=simpar,&
!                   struc=struc, &
!                   para_env = para_env)
!
!  int_grp = para_env % group  

!  nnodes = SIZE ( struc ( 1 ) % pnode )
!  beads = SIZE ( struc )
!  atmp = zero
!  cons = zero

!  DO i = 1, beads

!     CALL get_part_ke ( struc ( i ), ke )
!     CALL get_part_temp ( ke, simpar, t )
!     atmp = atmp + t
!     cons = cons + ke            
!
!  END DO

!  temp = atmp / beads

!  IF ( simpar % ensemble ( 1 : 3 )/= 'NVE' ) THEN
!     nhc_part_kin = zero
!     nhc_part_pot = zero
!     DO i = 1, beads
!        CALL get_nhc_energies ( nhc_part ( i ), atmp, btmp )
!        nhc_part_kin = nhc_part_kin + atmp
!        nhc_part_pot = nhc_part_pot + btmp
!     END DO
!  END IF

!  IF ( simpar % ensemble ( 1 : 3 ) == 'NPT' ) &
!  CALL get_nhc_energies (  nhc_baro ( 1 ), nhc_baro_kin, nhc_baro_pot )

!  IF ( simpar % ensemble ( 1 : 3 ) == 'NPT' ) &
!  CALL get_baro_energies ( box, simpar, nhc_baro ( 1 ), npt,  &
!                           baro_kin, baro_pot )
!  btmp = zero
!  DO i = 1, beads
!     atmp = zero
!     DO j = 1, nnodes
!        atomic_kind => struc (i)%pnode(j) %p % atomic_kind
!        CALL get_atomic_kind (atomic_kind=atomic_kind, mass=mass)
!        atmp = atmp + mass * (struc(i)%pimd_env%u(1,j)**2 + &
!           struc(i)%pimd_env%u(2,j)**2 + struc(i)%pimd_env%u(3,j)**2 )  
!     END DO
!     btmp = btmp + struc(i)%pimd_env%lbd * atmp
!  END DO
!  btmp = 0.5_dp*simpar%pimd_params%wp**2*btmp
!  cons = cons + btmp

!  ekin = 1.5_dp*beads*nnodes/simpar%pimd_params%beta - btmp
!  epot = epot/beads

!  SELECT CASE ( simpar % ensemble ( 1 : 3 ) )
!  CASE ( 'NVE')
!     cons = cons + epot
!  CASE ( 'NVT')
!     cons = cons + epot + nhc_part_kin + nhc_part_pot
!  CASE ( 'NPT')
!     cons = cons + epot + nhc_part_kin + nhc_part_pot  &
!                 + nhc_baro_kin + nhc_baro_pot + baro_kin + baro_pot
!  END SELECT

!  IF (itimes==1) THEN
!     constant = cons
!     CALL set_md_env (md_env=md_env, constant=constant)
!  END IF
!  IF (itimes /= 0 ) cons_rel = abs((cons-constant)/constant)

!END SUBROUTINE compute_conserved_quantity_pimd

!******************************************************************************
!!******
!!
!!   NAME
!!     get_econs_nve
!!
!!   FUNCTION
!!     calculates conserved quantity for nvt ensemble
!!
!!   AUTHOR
!!    gloria
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** *************************************************************************
SUBROUTINE get_econs_nve ( md_env, md_ener, cons, int_grp )


    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(md_ener_type), INTENT(inout)        :: md_ener
    REAL(KIND=dp), INTENT(OUT)                   :: cons
    INTEGER, INTENT(IN)                      :: int_grp

    REAL(KIND=dp)                                :: nhc_coef_kin, nhc_coef_pot
    TYPE(dyn_coeff_set_p_type), &
      DIMENSION(:), POINTER                  :: dyn_coeff_set
    TYPE(force_env_p_type), DIMENSION(:), POINTER :: force_env_p
    TYPE(lnhc_parameters_type), &
      DIMENSION(:), POINTER                  :: nhc_coef

    INTEGER :: iforce_env

!----------------------------------------------------------------------------

  NULLIFY(dyn_coeff_set, force_env_p, nhc_coef)

  CALL get_md_env (md_env, force_env_p=force_env_p, nhc_coef=nhc_coef)

  cons = md_ener % ekin + md_ener % epot 
  
  DO iforce_env = 1, SIZE(force_env_p)  

     CALL force_env_get (force_env = force_env_p(iforce_env)%force_env, &
                         dyn_coeff_set=dyn_coeff_set)

     IF (ASSOCIATED (dyn_coeff_set)) THEN

        CALL get_coef_ke (dyn_coeff_set(1)%dyn_coeff_set, md_ener%ekin_coefs, &
                      md_ener%temp_coefs, int_grp)
        cons = cons + md_ener%ekin_coefs
        IF (ASSOCIATED(nhc_coef)) THEN
           CALL get_nhc_energies ( nhc_coef(1), nhc_coef_pot, nhc_coef_kin,&
                                   int_grp )
           cons = cons + nhc_coef_kin + nhc_coef_pot
        END IF
     END IF
  END DO

END SUBROUTINE get_econs_nve

!******************************************************************************
!!******
!!
!!   NAME
!!     get_econs_nvt    
!!
!!   FUNCTION
!!     calculates conserved quantity for nvt ensemble           
!!
!!   AUTHOR
!!    gloria
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** *************************************************************************

SUBROUTINE get_econs_nvt ( md_env, md_ener, cons, int_grp)                            


    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(md_ener_type), INTENT(inout)        :: md_ener
    REAL(KIND=dp), INTENT(OUT)                   :: cons
    INTEGER, INTENT(IN)                      :: int_grp

    REAL(KIND=dp)                                :: nhc_coef_kin, nhc_coef_pot, &
                                                nhc_part_kin, nhc_part_pot
    TYPE(dyn_coeff_set_p_type), &
      DIMENSION(:), POINTER                  :: dyn_coeff_set
    TYPE(force_env_p_type), DIMENSION(:), POINTER  :: force_env_p
    TYPE(lnhc_parameters_type), &
      DIMENSION(:), POINTER                  :: nhc_coef, nhc_part

    INTEGER :: iforce_env

!----------------------------------------------------------------------------

  NULLIFY(dyn_coeff_set, force_env_p, nhc_part, nhc_coef)

  CALL get_md_env (md_env, force_env_p=force_env_p, nhc_part=nhc_part,&
                   nhc_coef=nhc_coef)
  CALL get_nhc_energies ( nhc_part ( 1 ), nhc_part_pot, nhc_part_kin, int_grp )
  cons = md_ener% ekin + md_ener%epot + nhc_part_kin + nhc_part_pot
  DO iforce_env = 1, SIZE(force_env_p)
     CALL force_env_get (force_env=force_env_p(iforce_env)%force_env, &
                         dyn_coeff_set=dyn_coeff_set)
     IF (ASSOCIATED (dyn_coeff_set)) THEN
        CALL get_coef_ke (dyn_coeff_set(1)%dyn_coeff_set, md_ener%ekin_coefs, &
                          md_ener%temp_coefs, int_grp)
        cons = cons + md_ener%ekin_coefs
        IF (ASSOCIATED(nhc_coef)) THEN
           CALL get_nhc_energies ( nhc_coef(1), nhc_coef_pot, nhc_coef_kin,&
                                   int_grp )
           cons = cons + nhc_coef_kin + nhc_coef_pot 
        END IF
     END IF
  END DO

END SUBROUTINE get_econs_nvt

!******************************************************************************
!!******
!!
!!   NAME
!!     get_econs_npt
!!
!!   FUNCTION
!!     calculates conserved quantity for npt ensemble
!!
!!   AUTHOR
!!     gloria
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*****************************************************************************

SUBROUTINE get_econs_npt ( md_env, md_ener, cons, int_grp)


    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(md_ener_type), INTENT(inout)        :: md_ener
    REAL(KIND=dp), INTENT(OUT)                   :: cons
    INTEGER, INTENT(IN)                      :: int_grp

    REAL(KIND=dp)                                :: baro_kin, baro_pot, &
                                                nhc_baro_kin, nhc_baro_pot,&
                                                nhc_part_kin, nhc_part_pot
    TYPE(cell_type), POINTER                 :: box
    TYPE(lnhc_parameters_type), POINTER      :: nhc_baro( : ), nhc_part( : )
    TYPE(npt_info_type), POINTER             :: npt( :, : )
    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar

!----------------------------------------------------------------------------

  CALL get_md_env ( md_env, nhc_part=nhc_part, nhc_baro=nhc_baro, &
                    simpar = simpar, npt = npt, cell = box )


  CALL get_nhc_energies ( nhc_part ( 1 ), nhc_part_pot, nhc_part_kin, int_grp )

  CALL get_nhc_energies ( nhc_baro ( 1 ), nhc_baro_pot, nhc_baro_kin, int_grp )

  CALL get_baro_energies ( box, simpar, npt, baro_kin, baro_pot, int_grp )

  cons = md_ener%ekin + md_ener%epot +  nhc_part_kin + nhc_part_pot  &
       + nhc_baro_kin + nhc_baro_pot + baro_kin + baro_pot


END SUBROUTINE get_econs_npt

!******************************************************************************
!!******
!!
!!   NAME
!!     get_econs_nph_uniaxial
!!
!!   FUNCTION
!!     calculates conserved quantity for nph_uniaxial
!!
!!   AUTHOR
!!     cjm
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*****************************************************************************

SUBROUTINE get_econs_nph_uniaxial ( md_env, md_ener, cons, int_grp)


    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(md_ener_type), INTENT(inout)        :: md_ener
    REAL(KIND=dp), INTENT(OUT)                   :: cons
    INTEGER, INTENT(IN)                      :: int_grp

    REAL(KIND=dp)                                :: baro_kin, baro_pot
    TYPE(cell_type), POINTER                 :: box
    TYPE(npt_info_type), POINTER             :: npt( :, : )
    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar

!----------------------------------------------------------------------------

  CALL get_md_env ( md_env, simpar = simpar, npt = npt, cell = box )

  CALL get_baro_energies ( box, simpar, npt, baro_kin, baro_pot, int_grp )

  md_ener % temp_baro = 2.0_dp * baro_kin 

  cons = md_ener%ekin + md_ener%epot + baro_kin + baro_pot


END SUBROUTINE get_econs_nph_uniaxial
!******************************************************************************
!!******
!!
!!   NAME
!!     get_econs_nph_damped
!!
!!   FUNCTION
!!     calculates conserved quantity for nph_uniaxial_damped
!!
!!   AUTHOR
!!     cjm
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*****************************************************************************

SUBROUTINE get_econs_nph_damped ( md_env, md_ener, cons, int_grp)


    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(md_ener_type), INTENT(inout)        :: md_ener
    REAL(KIND=dp), INTENT(OUT)                   :: cons
    INTEGER, INTENT(IN)                      :: int_grp

    REAL(KIND=dp)                                :: baro_kin, baro_pot, &
                                                    nhc_baro_kin, nhc_baro_pot
    TYPE(cell_type), POINTER                 :: box
    TYPE(lnhc_parameters_type), POINTER      :: nhc_baro( : )
    TYPE(npt_info_type), POINTER             :: npt( :, : )
    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar

!----------------------------------------------------------------------------

  CALL get_md_env ( md_env, simpar = simpar, npt = npt, nhc_baro = nhc_baro, &
                    cell = box )

  CALL get_nhc_energies ( nhc_baro ( 1 ), nhc_baro_pot, nhc_baro_kin, int_grp )

  CALL get_baro_energies ( box, simpar, npt, baro_kin, baro_pot, int_grp )

  md_ener % temp_baro = 2.0_dp * baro_kin 

  cons = md_ener%ekin + md_ener%epot + baro_kin + baro_pot + nhc_baro_kin + nhc_baro_pot


END SUBROUTINE get_econs_nph_damped
!******************************************************************************
!!****** md_conserved_quantity/get_baro_energies [1.0] *
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
!!*** *************************************************************************
SUBROUTINE get_baro_energies ( box, simpar, npt, baro_kin, baro_pot, int_grp )

    TYPE(cell_type), INTENT(IN)                        :: box
    TYPE(simulation_parameters_type), &
      INTENT(IN)                                       :: simpar
    TYPE(npt_info_type), DIMENSION(:, :), &
      INTENT(IN)                                       :: npt
    REAL(KIND=dp), INTENT(OUT)                         :: baro_kin, baro_pot
    INTEGER, INTENT(IN)                                :: int_grp
! locals
    INTEGER                                  :: i, j
    REAL ( dp )                             :: v0, iv0, v_shock 

  IF ( simpar%ensemble == 'NPT_I' ) THEN
     baro_pot = simpar % p_ext * box % deth
     baro_kin = 0.5_dp * npt ( 1, 1 ) % v * npt ( 1, 1 ) % v * &
          npt ( 1, 1 ) % mass
  ELSE IF ( simpar%ensemble=='NPT_F' ) THEN
     baro_pot =  simpar % p_ext * box % deth
     baro_kin = 0.0_dp
     DO i = 1, 3
        DO j = 1, 3
           baro_kin = baro_kin + 0.5_dp * npt ( i, j ) % v *  &
                      npt ( i, j ) % v * npt ( i, j ) % mass
        END DO
     END DO
  ELSEIF ( simpar%ensemble (1:3) == 'NPH' ) THEN
     v0 = simpar % v0
     iv0 = 1._dp / v0
     v_shock = simpar % v_shock

!!! Valid only for orthorhombic box !!!
     baro_pot = -0.5_dp * v_shock * v_shock * ( 1._dp - box % deth * iv0 )**2 - &
                 simpar % p0 * ( v0 - box % deth ) 
! dbg UNIAXIAL LIMIT
!     baro_pot =  simpar % p0 * box % deth
! dbg UNIAXIAL LIMIT
!!! Valid only for orthorhombic box !!!

     baro_kin = 0.5_dp * npt ( 1, 1 ) % v * npt ( 1, 1 ) % v * &
          npt ( 1, 1 ) % mass
  END IF

END SUBROUTINE get_baro_energies 

!******************************************************************************
!!****** md_conserved_quantity/get_nhc_energies [1.0] *
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
!!***8*************************************************************************
SUBROUTINE get_nhc_energies ( nhc, nhc_pot, nhc_kin, int_grp )
    TYPE(lnhc_parameters_type), INTENT(IN)   :: nhc
    REAL(KIND=dp), INTENT(OUT)               :: nhc_pot, nhc_kin
    INTEGER                                  :: int_grp

    INTEGER                                  :: l, n

   nhc_kin = 0.0_dp
   nhc_pot = 0.0_dp
   DO n = 1, nhc % num_nhc
      DO l = 1, nhc % nhc_len
         nhc_kin = nhc_kin +  &
         0.5_dp * nhc % nvt ( l, n ) % mass * &
         nhc % nvt ( l, n ) % v * nhc % nvt ( l, n ) % v

         nhc_pot = nhc_pot + nhc % nvt ( l, n ) % nkt *  &
                        nhc % nvt ( l, n ) % eta
      END DO
   END DO
   IF (nhc%dis_type=='DIS_REP') THEN
      CALL mp_sum(nhc_kin,int_grp)
      CALL mp_sum(nhc_pot,int_grp)
   END IF
END SUBROUTINE get_nhc_energies 
!******************************************************************************
!!****** md_conserved_quantity/get_part_ke [1.0] *
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
!!*** *************************************************************************

SUBROUTINE get_part_ke ( md_env, ekin, int_grp )


    TYPE(md_environment_type), POINTER       :: md_env
    REAL(KIND=dp), INTENT(OUT)                   :: ekin
    INTEGER, INTENT(IN)                      :: int_grp

    INTEGER                                  :: iparticle, iparticle_kind, &
                                                iparticle_local, &
                                                nparticle_kind, &
                                                nparticle_local
    REAL(KIND=dp)                                :: mass
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(force_env_p_type), DIMENSION(:), POINTER :: force_env_p
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

  ekin  = zero

  CALL get_md_env ( md_env, force_env_p = force_env_p )
  CALL force_env_get ( force_env = force_env_p(1)%force_env, subsys = subsys )
!HAF: subsys(1)?? ARGHHHHH
!HAF: ok, in that case I can just as well use force_env_p(1)%force_env 
  CALL cp_subsys_get(subsys=subsys(1)%subsys,&
                     atomic_kinds=atomic_kinds,&
                     local_particles=local_particles,&
                     particles=particles)

  nparticle_kind = atomic_kinds%n_els
  atomic_kind_set => atomic_kinds%els

  particle_set => particles%els

  DO iparticle_kind=1,nparticle_kind
    atomic_kind => atomic_kind_set(iparticle_kind)
    CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
    nparticle_local = local_particles%n_el(iparticle_kind)
    DO iparticle_local=1,nparticle_local
      iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
      ekin = ekin +  0.5_dp * mass *  &
                     ( particle_set(iparticle)%v(1)*particle_set(iparticle)%v(1) &
                    + particle_set(iparticle)%v(2)* particle_set(iparticle)%v(2) &
                    + particle_set(iparticle)%v(3)* particle_set(iparticle)%v(3) )
    END DO
  END DO


! sum all contributions to energy over calculated parts on all processors
  CALL mp_sum ( ekin, int_grp )
END SUBROUTINE get_part_ke 

!******************************************************************************
!!****** md_conserved_quantity/get_coef_ke [1.0] *
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
!!     gt 06-02-04 adapted to new dynamical coeffs type
!!
!!*****************************************************************************
!
 SUBROUTINE get_coef_ke ( dyn_coeff_set, ekinc, ctemp, int_group )

    TYPE(dyn_coeff_set_type), POINTER        :: dyn_coeff_set
    REAL(KIND=dp), INTENT(out)                   :: ekinc, ctemp
    INTEGER                                  :: int_group

    INTEGER                                  :: icoef, iel, ikind, ncoefs
    LOGICAL                                  :: gmass
    TYPE(dyn_coeff_type), POINTER            :: coeffs

!----------------------------------------------------------------------------

   NULLIFY(coeffs)
   ekinc = zero
   ncoefs = 0
   gmass=(dyn_coeff_set%mass_storage==dyn_coeff_mass_global)

   DO ikind=1,SIZE(dyn_coeff_set%coeffs_of_kind)
     IF (ASSOCIATED(dyn_coeff_set%coeffs_of_kind(ikind)%coeffs)) THEN
       coeffs=>dyn_coeff_set%coeffs_of_kind(ikind)%coeffs
       DO iel = 1, coeffs%n_els
         DO icoef = 1, coeffs%ncoef_atom
           ncoefs=ncoefs+1
           IF(gmass) THEN
             ekinc = ekinc + 0.5_dp*dyn_coeff_set%global_mass* &
                     coeffs%vel(iel,icoef)*coeffs%vel(iel,icoef)
           ELSE
             ekinc = ekinc + 0.5_dp*coeffs%masses(iel,icoef)* &
                     coeffs%vel(iel,icoef)*coeffs%vel(iel,icoef)
           END IF
         END DO
       END DO
     END IF
   END DO
 
 ! sum all contributions to energy over calculated parts on all processors
  CALL mp_sum ( ekinc, int_group )
  CALL mp_sum ( ncoefs, int_group )

! get coefficient temperature

  ctemp=2.0_dp*ekinc/REAL(ncoefs,dp) 
  
  END SUBROUTINE get_coef_ke
!******************************************************************************

END MODULE md_conserved_quantities
!******************************************************************************
