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
  
  USE atomic_kind_list_types,        ONLY: atomic_kind_list_type
  USE atomic_kind_types,             ONLY : atomic_kind_type,&
                                            get_atomic_kind
  USE coefficient_types,             ONLY : coeff_type
  USE cp_para_types,                 ONLY : cp_para_env_type
  USE cp_subsystem_types,            ONLY : cp_subsystem_p_type,&
                                            cp_subsys_get
  USE distribution_1d_types,         ONLY : distribution_1d_type
  USE extended_system_types,         ONLY : lnhc_parameters_type, &
                                            npt_info_type
  USE kinds,                         ONLY : dbl
  USE force_env_types,               ONLY : force_env_get,&
                                            force_env_type, &
                                            force_env_write_variables
  USE mathconstants,                 ONLY : zero
  USE md,                            ONLY : simulation_parameters_type
  USE md_environment_types,          ONLY : md_environment_type,&
                                            get_md_env,&
                                            set_md_env
  USE message_passing,               ONLY : mp_sum, &
                                            mp_comm_dup
  USE particle_types,                ONLY: particle_type
  USE particle_list_types,           ONLY: particle_list_type
  USE simulation_cell,               ONLY : cell_type
  USE termination,                   ONLY : stop_memory, &
                                            stop_program
  USE util,                          ONLY : get_unit

  IMPLICIT NONE

  PRIVATE

  TYPE md_ener_type
  REAL (dbl) :: epot  
  REAL (dbl) :: ekin  
  REAL (dbl) :: temp  
  REAL (dbl) :: ekin_coefs
  REAL (dbl) :: temp_coefs  
  END TYPE md_ener_type


  PUBLIC :: compute_conserved_quantity
           ! compute_conserved_quantity_pimd                                    


!******************************************************************************
CONTAINS

!!*****************************************************************************
!!******
!!
!!   NAME
!!     compute_conserved_quantity_md
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
                                       temp, ekin_coefs, temp_coefs )                            

  IMPLICIT NONE

! Arguments
  TYPE (md_environment_type), POINTER :: md_env
  REAL (dbl), INTENT (IN) :: epot
  REAL (dbl), INTENT (OUT) :: cons
  REAL (dbl), INTENT (OUT), OPTIONAL :: cons_rel
  REAL (dbl), INTENT (OUT), OPTIONAL :: ekin
  REAL (dbl), INTENT (OUT), OPTIONAL :: temp
  REAL (dbl), INTENT (OUT), OPTIONAL :: ekin_coefs
  REAL (dbl), INTENT (OUT), OPTIONAL :: temp_coefs

!Locals
  TYPE (simulation_parameters_type), POINTER :: simpar
  INTEGER, POINTER :: itimes
  INTEGER ::  int_grp
  REAL ( dbl ), POINTER :: constant
  type ( md_ener_type ) :: md_ener
  type ( cp_para_env_type ) ,POINTER :: para_env

!------------------------------------------------------------------------------
  md_ener%epot=epot 

  CALL get_md_env (md_env=md_env,&
                   constant=constant,&
                   itimes=itimes,&
                   simpar=simpar,&
                   para_env = para_env )

  int_grp =  para_env % group  
 
  CALL get_part_ke ( md_env, md_ener%ekin, int_grp )
  
  CALL get_part_temp ( md_ener%ekin, simpar, md_ener%temp )

  SELECT CASE ( simpar % ensemble ( 1 : 3 ) )
  CASE ( 'NVE')
    CALL get_econs_nve (md_env, md_ener, cons)
  CASE ( 'NVT')
    CALL get_econs_nvt (md_env, md_ener, cons, int_grp)
  CASE ( 'NPT')
    CALL get_econs_npt (md_env, md_ener, cons, int_grp)
  END SELECT

  IF (itimes==1) THEN
     constant = cons
     CALL set_md_env (md_env=md_env, constant=constant)
  END IF

  IF ((PRESENT(cons_rel)).AND.(itimes /= 0 )) THEN
    cons_rel = abs((cons-constant)/constant)
  END IF
  IF (PRESENT(ekin)) ekin=md_ener%ekin 
  IF (PRESENT(temp)) temp=md_ener%temp 
  IF (PRESENT(ekin_coefs)) ekin_coefs=md_ener%ekin_coefs 
  IF (PRESENT(temp_coefs)) temp_coefs=md_ener%temp_coefs 

END SUBROUTINE compute_conserved_quantity

!*****************************************************************************

!SUBROUTINE compute_conserved_quantity_pimd (md_env, epot, cons,&
!                                            cons_rel, ekin, temp)

!  IMPLICIT NONE

! Arguments
!  TYPE (md_environment_type), POINTER :: md_env
!  REAL (dbl), INTENT (INOUT) :: epot
!  REAL (dbl), INTENT (OUT):: cons
!  REAL (dbl), INTENT (OUT):: cons_rel
!  REAL (dbl), INTENT (OUT):: ekin
!  REAL (dbl), INTENT (OUT):: temp

!Locals
!  TYPE (simulation_parameters_type), POINTER :: simpar
!  TYPE (structure_type), POINTER :: struc ( : )

!  TYPE (cell_type), POINTER :: box
!  TYPE (lnhc_parameters_type), POINTER :: nhc_baro ( : )
!  TYPE (lnhc_parameters_type), POINTER :: nhc_part ( : )
!  TYPE (npt_info_type), POINTER :: npt ( :, : )

!  INTEGER, POINTER :: itimes
!  REAL (dbl), POINTER :: constant
!  TYPE (atomic_kind_type), POINTER :: atomic_kind

!  INTEGER i,j,nnodes,beads                             
!  REAL(dbl) :: ke, t, atmp, btmp, mass 
!  REAL(dbl) :: nhc_part_kin, nhc_part_pot, nhc_baro_kin 
!  REAL(dbl) :: nhc_baro_pot, baro_kin, baro_pot 
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
!  btmp = 0.5_dbl*simpar%pimd_params%wp**2*btmp
!  cons = cons + btmp

!  ekin = 1.5_dbl*beads*nnodes/simpar%pimd_params%beta - btmp
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
SUBROUTINE get_econs_nve ( md_env, md_ener, cons )

  IMPLICIT NONE

! Arguments
  TYPE (md_environment_type), POINTER :: md_env
  TYPE (md_ener_type), INTENT (inout) :: md_ener
  REAL (dbl), INTENT (OUT) :: cons

  TYPE (lnhc_parameters_type), POINTER :: nhc_coef

! Locals
  REAL ( dbl ) :: nhc_coef_pot, nhc_coef_kin

!----------------------------------------------------------------------------

  cons = md_ener % ekin + md_ener % epot 


!  IF (ASSOCIATED (struc%coef_vel)) THEN

!    CALL get_coef_ke (struc%coef_vel(1), md_ener%ekin_coefs)
!
!    CALL get_coef_temp (md_ener%ekin_coefs, struc%coef_vel(1), md_ener%temp_coefs)
!
!    CALL get_md_env (md_env=md_env, nhc_coef=nhc_coef)
!
!    CALL get_nhc_energies ( nhc_coef, nhc_coef_kin, nhc_coef_pot )
!
!    cons = cons + nhc_coef_kin + nhc_coef_pot + md_ener%ekin_coefs
!
!  END IF

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

  IMPLICIT NONE

! Arguments
  TYPE (md_environment_type), POINTER :: md_env
  TYPE (md_ener_type), INTENT (inout) :: md_ener
  REAL (dbl), INTENT (OUT) :: cons
  INTEGER, INTENT ( IN ) :: int_grp

! Locals
  TYPE (lnhc_parameters_type), POINTER :: nhc_part ( : )
  TYPE (lnhc_parameters_type), POINTER :: nhc_coef       

  REAL ( dbl ) :: nhc_part_pot, nhc_part_kin
  REAL ( dbl ) :: nhc_coef_pot, nhc_coef_kin

!----------------------------------------------------------------------------



  CALL get_md_env ( md_env, nhc_part = nhc_part )
  CALL get_nhc_energies ( nhc_part ( 1 ), nhc_part_kin, nhc_part_pot, int_grp )

  cons = md_ener% ekin + md_ener%epot + nhc_part_kin + nhc_part_pot


!  IF (ASSOCIATED (struc%coef_vel)) THEN

!    CALL get_coef_ke ( struc % coef_vel (1), md_ener%ekin_coefs )

!    CALL get_coef_temp ( md_ener%ekin_coefs, struc%coef_vel(1), md_ener%temp_coefs )
 
!    CALL get_md_env (md_env=md_env, nhc_coef=nhc_coef)

!    CALL get_nhc_energies ( nhc_coef, nhc_coef_kin, nhc_coef_pot )

!    cons = cons + nhc_coef_kin + nhc_coef_pot + md_ener%ekin_coefs

!  END IF

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

  IMPLICIT NONE

! Arguments
  TYPE (md_environment_type), POINTER :: md_env
  TYPE (md_ener_type), INTENT (inout) :: md_ener
  REAL (dbl), INTENT (OUT) :: cons
  INTEGER, INTENT ( IN ) :: int_grp

  TYPE (cell_type), POINTER :: box
  TYPE (lnhc_parameters_type), POINTER :: nhc_baro ( : )
  TYPE (lnhc_parameters_type), POINTER :: nhc_part ( : )
  TYPE (npt_info_type), POINTER :: npt ( :, : )
  TYPE (simulation_parameters_type), POINTER :: simpar

! Locals
  REAL ( dbl ) :: nhc_part_pot, nhc_part_kin, nhc_baro_pot, nhc_baro_kin
  REAL ( dbl ) :: baro_pot, baro_kin

!----------------------------------------------------------------------------

  CALL get_md_env ( md_env, nhc_part=nhc_part, nhc_baro=nhc_baro, &
                    simpar = simpar, npt = npt, cell = box )


  CALL get_nhc_energies ( nhc_part ( 1 ), nhc_part_kin, nhc_part_pot, int_grp )

  CALL get_nhc_energies ( nhc_baro ( 1 ), nhc_baro_kin, nhc_baro_pot, int_grp )

  CALL get_baro_energies ( box, simpar, nhc_baro ( 1 ), npt, baro_kin, baro_pot, int_grp )

  cons = md_ener%ekin + md_ener%epot +  nhc_part_kin + nhc_part_pot  &
       + nhc_baro_kin + nhc_baro_pot + baro_kin + baro_pot


END SUBROUTINE get_econs_npt


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
SUBROUTINE get_baro_energies ( box, simpar, nhc, npt, baro_kin, baro_pot, int_grp )
  IMPLICIT NONE

  TYPE ( cell_type ), INTENT ( IN ) :: box
  TYPE ( simulation_parameters_type), INTENT ( IN ) :: simpar
  TYPE ( lnhc_parameters_type), INTENT ( IN ) :: nhc
  TYPE ( npt_info_type ), INTENT ( IN ), DIMENSION ( :, : ) :: npt
  REAL ( dbl ), INTENT ( OUT ) :: baro_pot, baro_kin
  INTEGER, INTENT ( IN ) :: int_grp

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
     IF (nhc%dis_type=='DIS_REP') THEN
        CALL mp_sum(baro_kin,int_grp)
        CALL mp_sum(baro_pot,int_grp)
     END IF
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
!!***8*************************************************************************
SUBROUTINE get_nhc_energies ( nhc, nhc_pot, nhc_kin, int_grp )
  IMPLICIT NONE
  TYPE ( lnhc_parameters_type ), INTENT ( IN ) :: nhc
  REAL ( dbl ), INTENT ( OUT ) :: nhc_pot, nhc_kin
  INTEGER :: int_grp

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
   IF (nhc%dis_type=='DIS_REP') THEN
      CALL mp_sum(nhc_kin,int_grp)
      CALL mp_sum(nhc_pot,int_grp)
   END IF
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
!!*** *************************************************************************
SUBROUTINE get_part_temp ( ekin, simpar, temperature )

  IMPLICIT NONE
  REAL ( dbl ), INTENT ( OUT ) :: temperature 
  TYPE ( simulation_parameters_type ), INTENT ( IN ) :: simpar
  REAL ( dbl ), INTENT ( in ) :: ekin         

! to be converted in K in the output routine

  temperature = 2.0_dbl*ekin/float(simpar%nfree)

END SUBROUTINE  get_part_temp
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
!!*** *************************************************************************

SUBROUTINE get_part_ke ( md_env, ekin, int_grp )

  IMPLICIT NONE

!arguments
  TYPE ( md_environment_type ), POINTER :: md_env
  REAL ( dbl ), INTENT ( OUT ) :: ekin         
  INTEGER, INTENT ( IN ) :: int_grp
!locals
  TYPE(atomic_kind_list_type), POINTER :: atomic_kinds
  TYPE(atomic_kind_type), POINTER      :: atomic_kind
  TYPE(distribution_1d_type), POINTER  :: local_particles
  TYPE(force_env_type), POINTER        :: force_env
  TYPE(particle_list_type), POINTER    :: particles
  TYPE(cp_subsystem_p_type), POINTER, DIMENSION(:) :: subsys
  TYPE(atomic_kind_type), DIMENSION(:), POINTER   :: atomic_kind_set
  TYPE(particle_type), DIMENSION(:), POINTER      :: particle_set
  INTEGER  :: iparticle,iparticle_kind,iparticle_local,&
              nparticle_kind,nparticle_local
  REAL ( dbl ) :: mass

  ekin  = zero

  CALL get_md_env ( md_env, force_env = force_env )
  CALL force_env_get ( force_env = force_env, subsys = subsys)
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
      ekin = ekin +  0.5_dbl * mass *  &
                     ( particle_set(iparticle)%v(1)*particle_set(iparticle)%v(1) &
                    + particle_set(iparticle)%v(2)* particle_set(iparticle)%v(2) &
                    + particle_set(iparticle)%v(3)* particle_set(iparticle)%v(3) )
    END DO
  END DO


! sum all contributions to energy over calculated parts on all processors
  CALL mp_sum(ekin,int_grp)
END SUBROUTINE get_part_ke 

!******************************************************************************
!
! SUBROUTINE get_coef_temp ( ekinc, coef, ctemp )

!   IMPLICIT NONE

! arguments
!   TYPE ( coeff_type ), INTENT ( IN ) :: coef
!   REAL (dbl), INTENT ( IN ) :: ekinc
!   REAL ( dbl ), INTENT ( OUT ) :: ctemp
! locals
!   INTEGER :: inuse, ncoef

! to be converted in K in the output routine
 
!   inuse = coef % in_use
!   SELECT CASE ( inuse )
!   CASE ( 402 )
!     ncoef = SIZE (coef%ao%cr)
!     ctemp = 2.0_dbl* ekinc / REAL ( ncoef, dbl ) 
!   CASE ( 401 )
!     ncoef = SIZE (coef%pw%cc)
!     ctemp =  ekinc / REAL ( ncoef, dbl ) 
!   CASE DEFAULT
!      CALL stop_program ( 'get_coef_temp','integrator not implemented')
!   END SELECT
! END SUBROUTINE  get_coef_temp
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
!!*****************************************************************************
!
! SUBROUTINE get_coef_ke ( coef, ekinc )

!   IMPLICIT NONE
! arguments
!   TYPE ( coeff_type ), INTENT ( IN ) :: coef
!   REAL ( dbl ), INTENT ( out ) :: ekinc
!locals 
!   INTEGER :: i
!   INTEGER :: inuse
!   REAL ( dbl ) :: cmass

!   ekinc = zero
!   inuse = coef % in_use
!   SELECT CASE ( inuse )
!   CASE ( 402 )
!     DO i = 1, SIZE (coef%ao%cr)
!        ekinc = ekinc + 0.5_dbl * coef % ao % mass ( i ) * & 
!                coef % ao % cr ( i ) * coef % ao % cr ( i )
!     END DO
!   CASE ( 401 )
!     DO i = 1, SIZE (coef%pw%cc)
!        cmass = coef % pw % mass_cc ( i )
!        ekinc = ekinc + 0.5_dbl * cmass * &
!         ( REAL ( coef % pw % cc ( i ) ) * REAL ( coef % pw % cc ( i ) ) + &
!           AIMAG ( coef % pw % cc ( i ) )* AIMAG ( coef % pw % cc ( i ) ) )
!     END DO
!   CASE DEFAULT
!      CALL stop_program ( 'get_coef_ke','integrator not implemented')
!   END SELECT 
 
 ! sum all contributions to energy over calculated parts on all processors
! #if defined(__parallel)
!   CALL mp_sum(ekinc,int_grp)
! #endif
! END SUBROUTINE get_coef_ke

!******************************************************************************

END MODULE md_conserved_quantities

!******************************************************************************
