!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C ) 2002 - 2003  CP2K developers group                         !
!-----------------------------------------------------------------------------!
#include "cp_prep_globals.h"
!!****** cp2k/kg_force [1.0] *
!!
!!   NAME
!!     kg_force
!!
!!   FUNCTION
!!     kg force calculation routine      
!!
!!   AUTHOR
!!     gloria,30.09.2002 
!!
!!   MODIFICATION HISTORY
!!     JGH (22-Feb-03) PW grid options added
!!
!!   SOURCE
!******************************************************************************

MODULE kg_force
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE checkpoint_handler,              ONLY: write_checkpoint_information
  USE coefficient_types,               ONLY: coeff_sumup,&
                                             coeff_transform_space,&
                                             coeff_type,&
                                             coeff_zero,&
                                             coeff_init
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_subsystem_types,              ONLY: cp_subsystem_p_type
  USE dynamical_coeff_types,           ONLY: dyn_coeff_set_type,&
                                             dyn_coeff_type
  USE global_types,                    ONLY: global_environment_type
  USE hartree,                         ONLY: calculate_hartree
  USE kg_core,                         ONLY: calculate_ecore_rspace,&
                                             calculate_eself
  USE kg_density,                      ONLY: calculate_density,&
                                             calculate_epc_density,&
                                             calculate_epc_rspace_forces,&
                                             calculate_p_density,&
                                             calculate_total_density,&
                                             calculate_v_rspace_forces,&
                                             calculate_vp_rspace_forces
  USE kg_energy_types,                 ONLY: init_kg_energy,&
                                             kg_energy_type
  USE kg_environment_types,            ONLY: get_kg_env,&
                                             kg_environment_type
  USE kg_force_types,                  ONLY: init_kg_force,&
                                             kg_force_type
  USE kg_intra,                        ONLY: calculate_ebond_corr
  USE kg_kxc,                          ONLY: calculate_kxc_derivatives
  USE kg_neighbor_lists,               ONLY: build_kg_neighbor_lists
  USE kg_ppl,                          ONLY: calculate_drho_ppl,&
                                             calculate_rho0_ppl
  USE kg_rho_types,                    ONLY: kg_rho_type, &
                                             kg_rho_get
  USE kg_rspw_types,                   ONLY: kg_rspw_type,&
                                             kg_rspw_get
  USE kinds,                           ONLY: dp
  USE machine,                         ONLY: m_flush
  USE mathconstants,                   ONLY: zero
  USE message_passing,                 ONLY: mp_sum
  USE particle_types,                  ONLY: particle_type
  USE pw_pool_types,                   ONLY: pw_pool_give_back_coeff,&
                                             pw_pool_init_coeff,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_p_type
  USE pws,                             ONLY: cleanup_pw_poisson_solver,&
                                             init_pw_poisson_solver
  USE simulation_cell,                 ONLY: cell_type
  USE termination,                     ONLY: stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xc,                              ONLY: xc_vxc_pw_create1

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: kg_force_control

!******************************************************************************

CONTAINS

!!****** kg_force/kg_force_control [1.1] *
!!
!!   NAME
!!     kg_force_control
!!
!!   SYNOPSIS
!!     Subroutine kg_force_control(kg_env, globenv, debug)
!!       Implicit None
!!       Type(kg_environment_type), Intent (INOUT):: kg_env
!!       Type(global_environment_type), Intent (IN):: globenv
!!       Type(debug_variables_type), Intent (OUT), Optional:: debug
!!     End Subroutine kg_force_control
!!
!!   FUNCTION
!!     Calculates the total potential energy, total force, and the
!!     total pressure tensor from the potentials
!!
!!   AUTHOR
!!     gt
!!
!!   MODIFICATION HISTORY
!!
!!*** *************************************************************************

SUBROUTINE kg_force_control ( kg_env, globenv,error )

    TYPE(kg_environment_type), POINTER       :: kg_env
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    LOGICAL                                  :: kgpol
    TYPE(dft_control_type), POINTER          :: dft_control

!  TYPE ( debug_variables_type ), INTENT ( OUT ), OPTIONAL :: debug
!------------------------------------------------------------------------------

  kgpol = .FALSE.
  CALL get_kg_env( kg_env=kg_env, dft_control=dft_control)
  kgpol = dft_control % qs_control % polarization
 
  IF (kgpol) THEN
    CALL kgpol_calculate_forces (kg_env, globenv)
  ELSE
    CALL kg_calculate_forces (kg_env, globenv )
  END IF
  END SUBROUTINE kg_force_control


!******************************************************************************
!!****** kg_force/kgpol_calculate_forces [1.1] *
!!
!!   NAME
!!     kgpol_calculate_forces
!!
!!   FUNCTION
!!     Calculates the total potential energy, total force
!!
!!   AUTHOR
!!     gt         
!!
!!   MODIFICATION HISTORY
!!
!!*** *************************************************************************

SUBROUTINE kgpol_calculate_forces ( kg_env, globenv, error )


    TYPE(kg_environment_type), POINTER       :: kg_env
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "SUBROUTINE kgpol_calculate_forces (MODULE kg_force)"

    INTEGER                                  :: handle, i, &
                                                iatom, ikind, stat, natoms, &
                                                nkind, output_unit
    INTEGER, ALLOCATABLE, DIMENSION(:), SAVE :: natom_of_kind
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: calculate_forces, ionode,&
                                                failure
    REAL(KIND=dp)                                :: dvol, e1_xc, e2_xc, &
                                                total_rho0_rspace, &
                                                total_rhop_rspace
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(coeff_type), POINTER                :: rho_core, rho_r, rho_g,&
                                                rhop_r, rhop_g, &  
                                                v_gspace, &
                                                v_rspace, vxc_r,&
                                                v1xc_r, v2xc_r
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: my_vxc_r, my_vxc_tau, &
                                                my_rho_r, my_rho_g,&
                                                my_tau
    TYPE(cp_subsystem_p_type), DIMENSION(1)  :: subsys
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(dyn_coeff_set_type), POINTER        :: dyn_coeff_set
    TYPE(dyn_coeff_type), POINTER            :: local_coeffs
    TYPE(kg_energy_type), POINTER            :: energy
    TYPE(kg_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(kg_rho_type), POINTER               :: rho
    TYPE(kg_rspw_type), POINTER              :: rspw
    TYPE(particle_type), POINTER             :: part( : )
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool

!  TYPE ( debug_variables_type ), INTENT ( OUT ), OPTIONAL :: debug
!  CHARACTER(LEN=8) :: force_type
!------------------------------------------------------------------------------

  CALL write_checkpoint_information("entering "//routine,globenv)
  CALL timeset ( 'FORCE','I',' ',handle )
  ionode=globenv%ionode
  output_unit=globenv%scr
  NULLIFY(force,energy,part,atomic_kind_set,atomic_kind,cell,dft_control,&
          dyn_coeff_set,local_coeffs,rho,rspw,auxbas_pw_pool)
  NULLIFY(rho_r,rho_g,rhop_r,rhop_g,rho_core,v_gspace,v_rspace,&
          vxc_r,v1xc_r,v2xc_r)
  NULLIFY(my_vxc_r,my_vxc_tau,my_rho_r,my_rho_g,my_tau,subsys(1)%subsys)

  CALL get_kg_env( kg_env=kg_env, atomic_kind_set=atomic_kind_set,&
                    dft_control=dft_control, energy=energy, cell=cell,&
                    dyn_coeff_set=dyn_coeff_set,subsys=subsys(1)%subsys,&
                    force=force,particle_set=part,rho=rho,rspw=rspw)

  CPPrecondition(ASSOCIATED(rho),cp_failure_level,routine,error,failure)
  CPPrecondition(ASSOCIATED(rspw),cp_failure_level,routine,error,failure)
  IF(.NOT.failure) THEN
    CALL kg_rho_get(kg_rho=rho, rho_r=rho_r, rho_g=rho_g, rhop_r=rhop_r,&
                    rhop_g=rhop_g,rho_core=rho_core)
    CALL kg_rspw_get(kg_rspw=rspw,auxbas_pw_pool=auxbas_pw_pool,error=error)
  END IF

  ALLOCATE(v_gspace,stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routine,error,failure)
  ALLOCATE(v_rspace,stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routine,error,failure)
  ALLOCATE(vxc_r,stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routine,error,failure)
  ALLOCATE(v1xc_r,stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routine,error,failure)
  ALLOCATE(v2xc_r,stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routine,error,failure)

  IF (.NOT.failure) THEN
    CALL pw_pool_init_coeff(auxbas_pw_pool, v_gspace, &
                            use_data = COMPLEXDATA1D,&
                            in_space = RECIPROCALSPACE, error=error)
    CALL coeff_zero(v_gspace)
    CALL pw_pool_init_coeff(auxbas_pw_pool, v_rspace,&
                            use_data = REALDATA3D,&
                            in_space = REALSPACE, error=error)
    CALL coeff_zero(v_rspace)
    CALL pw_pool_init_coeff(auxbas_pw_pool, v1xc_r,&
                            use_data = REALDATA3D,&
                            in_space = REALSPACE, error=error)
    CALL coeff_zero(v1xc_r)
    CALL pw_pool_init_coeff(auxbas_pw_pool, v2xc_r,&
                            use_data = REALDATA3D,&
                            in_space = REALSPACE, error=error)
    CALL coeff_zero(v2xc_r)
  END IF

  calculate_forces=.TRUE. 
  natoms = SIZE ( part )
  nkind = SIZE( atomic_kind_set )
  stat = 0
  dvol=v_rspace%pw%pw_grid%dvol

  IF ( .NOT. ALLOCATED ( natom_of_kind ) ) THEN
     ALLOCATE (natom_of_kind(nkind),STAT=stat)
     IF (stat /= 0) CALL stop_memory(routine,"natom_of_kind",nkind)
     CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                              natom_of_kind=natom_of_kind)
  END IF

  CALL build_kg_neighbor_lists(kg_env,globenv)

! initialize forces and energies                    
 
  CALL init_kg_force(force)
  CALL init_kg_energy(energy)

!reinitialize total forces
  DO i = 1, natoms
     part ( i ) % f ( 1 ) = 0.0_dp
     part ( i ) % f ( 2 ) = 0.0_dp
     part ( i ) % f ( 3 ) = 0.0_dp
  END DO

!reinitialize forces on coefficents
  DO ikind= 1, nkind
    local_coeffs => dyn_coeff_set%coeffs_of_kind(ikind)%coeffs
    IF(ASSOCIATED(local_coeffs)) THEN
      local_coeffs%forces(:,:)=0.0_dp
    END IF
  END DO


! self energy 
  CALL calculate_eself ( kg_env, energy % core_self )

! Calculate the rspace energy of the core charge distribution
! and its force contributions

  CALL calculate_ecore_rspace(kg_env,globenv,calculate_forces)

! calculate intramolecular bonded_correction energy  

  CALL calculate_ebond_corr(kg_env,globenv,energy%bond_corr,calculate_forces)

! Calculate the local pseudopotential energy         

  CALL calculate_rho0_ppl(kg_env,globenv,calculate_forces)
  CALL calculate_drho_ppl(kg_env,globenv,calculate_forces)

! Calculate the density generated by the core charges

  CALL calculate_epc_density(rho_core,rho%total_rho_core_rspace, kg_env)

! calculate the frozen electronic pseudo density

  CALL calculate_density (rho_r, rho_g, total_rho0_rspace, kg_env)        

  CALL calculate_p_density(rhop_r, rhop_g, total_rhop_rspace, kg_env)      

! from the frozen rho, calculate the thomas-fermi and xc potentials
!  CALL xc_calculate_pw1 (rho_r%pw, rho_g%pw, vxc_r%pw, energy%exc, &
!                         dft_control,auxbas_pw_pool)
  ALLOCATE(my_rho_r(1))
  my_rho_r(1)%pw => rho_r%pw
  ALLOCATE(my_rho_g(1))
  my_rho_g(1)%pw => rho_g%pw
  CALL xc_vxc_pw_create1(my_vxc_r,my_vxc_tau,my_rho_r,my_rho_g,my_tau,&
       energy%exc,dft_control,cell,auxbas_pw_pool,subsys,&
       calculate_forces=.FALSE.)
! ALLOCATE(vxc_r)
  CALL coeff_init(vxc_r,pw=my_vxc_r(1)%pw)
  NULLIFY(my_rho_r(1)%pw,my_rho_g(1)%pw, my_vxc_r(1)%pw)
  DEALLOCATE(my_rho_r,my_rho_g, my_vxc_r)
  NULLIFY(my_rho_r,my_rho_g, my_vxc_r)

  vxc_r%pw%cr3d(:,:,:) =dvol*vxc_r%pw%cr3d(:,:,:)

  CALL calculate_kxc_derivatives(auxbas_pw_pool, rho_r, rho_g, rhop_r, vxc_r, v1xc_r,&
                                 v2xc_r, e1_xc, e2_xc, dft_control, calculate_forces)

! calculate the total density
  rho%total_rho_rspace = total_rho0_rspace+ rho%total_rho_core_rspace

! add the core density to the frozen electronic density
  CALL coeff_sumup(rho_core,rho_g)

  rho%total_rho_gspace = calculate_total_density(rho_g)

! add the polarization density to get the total density
  CALL coeff_sumup(rhop_g,rho_g)

! calculate electrostatic potential
  CALL init_pw_poisson_solver(cell%green)
  CALL calculate_hartree(rho_g, energy % hartree,v_gspace)
  CALL cleanup_pw_poisson_solver(cell%green)

!  calclation of the forces on the ions

  CALL coeff_transform_space(v_gspace,v_rspace)
  CALL pw_pool_give_back_coeff(auxbas_pw_pool,v_gspace,error=error)
  DEALLOCATE(v_gspace,stat=stat)
  CPPostconditionNoFail(stat==0,cp_warning_level,routine,error)

  v_rspace%pw%cr3d(:,:,:) = dvol*v_rspace%pw%cr3d(:,:,:)

  CALL calculate_epc_rspace_forces(v_rspace, kg_env)         

  v_rspace%pw%cr3d(:,:,:) =-1.0_dp*(vxc_r%pw%cr3d(:,:,:) +&
               v1xc_r%pw%cr3d(:,:,:)+ v_rspace%pw%cr3d(:,:,:))

  CALL calculate_vp_rspace_forces(v_rspace,kg_env, calculate_forces)

  v_rspace%pw%cr3d(:,:,:) = -1.0_dp*v_rspace%pw%cr3d(:,:,:)+&
                            v2xc_r%pw%cr3d(:,:,:)

  CALL calculate_v_rspace_forces(v_rspace,kg_env) 

! deallocate work storage
  CALL pw_pool_give_back_coeff(auxbas_pw_pool,v_rspace,error=error)
  DEALLOCATE(v_rspace,stat=stat)
  CPPostconditionNoFail(stat==0,cp_warning_level,routine,error)
  CALL pw_pool_give_back_coeff(auxbas_pw_pool,vxc_r,error=error)
  DEALLOCATE(vxc_r,stat=stat)
  CPPostconditionNoFail(stat==0,cp_warning_level,routine,error)
  CALL pw_pool_give_back_coeff(auxbas_pw_pool,v1xc_r,error=error)
  DEALLOCATE(v1xc_r,stat=stat)
  CPPostconditionNoFail(stat==0,cp_warning_level,routine,error)
  CALL pw_pool_give_back_coeff(auxbas_pw_pool,v2xc_r,error=error)
  DEALLOCATE(v2xc_r,stat=stat)
  CPPostconditionNoFail(stat==0,cp_warning_level,routine,error)

! add up all the potential energies

  energy % total = energy % core_self + energy % core_overlap + energy % exc + &
                 e1_xc + e2_xc + energy % hartree + energy % pseudo + energy%ppseudo + &
                 energy% bond_corr
     
! print
    IF (ionode.AND.(globenv%print%level>3)) THEN
       WRITE (UNIT=output_unit,FMT="(/,(T3,A,T40,F20.10))")&
            "Total electronic density (r-space): ",&
            total_rho0_rspace,  &
            "Total core charge density (r-space):",&
            rho%total_rho_core_rspace
       WRITE (UNIT=output_unit,FMT="(T3,A,T60,F20.10)")&
            "Total charge density (r-space):     ",&
            rho%total_rho_rspace
       WRITE (UNIT=output_unit,FMT="(T3,A,T60,F20.10)")&
            "Total rho_coefs charge density (r-space):     ",&
            total_rhop_rspace
       WRITE (UNIT=output_unit,FMT="(/,(T3,A,T55,F25.14))")&
            "Overlap energy of the core charge distribution:",energy%core_overlap,&
            "Bond correction energy:                        ",energy%bond_corr,&
            "Self energy of the core charge distribution:   ",energy%core_self,&
            "Pseudopotential  energy:                       ",energy%pseudo,&
            "Pseudopotential  energy  due to coefs:         ",energy%ppseudo,&
            "Hartree energy:                                ",energy%hartree,&
            "Exchange-correlation energy:                   ",energy%exc,&
            "1st order correction to xc energy coefs        ",e1_xc,&
            "2nd order correction to xc energy coefs        ",e2_xc,&
            "Total energy:                                  ",energy%total
       CALL m_flush(output_unit)
    END IF

! add up all the forces

  DO ikind = 1, nkind         
    CALL mp_sum(force(ikind)%f_rho,globenv%group)
    CALL mp_sum(force(ikind)%f_rspace_core,globenv%group)
    CALL mp_sum(force(ikind)%f_hartree_core,globenv%group)
    CALL mp_sum(force(ikind)%f_ppl,globenv%group)
    force(ikind)%f_total(:,:) = force(ikind)%f_rho(:,:) +&
                                force(ikind)%f_rspace_core(:,:) +&
                                force(ikind)%f_hartree_core(:,:) +&
                                force(ikind)%f_ppl(:,:) 
  END DO


  DO ikind = 1, nkind         
    atomic_kind => atomic_kind_set(ikind)
    CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           atom_list=atom_list) 
    DO iatom=1, SIZE(atom_list)
      i = atom_list(iatom)
! fill in the forces in part
      part(i)%f(1) = force(ikind)%f_total(1,iatom)
      part(i)%f(2) = force(ikind)%f_total(2,iatom)
      part(i)%f(3) = force(ikind)%f_total(3,iatom)
    END DO
  END DO

! deallocating all local variables
  IF ( ALLOCATED ( natom_of_kind ) ) THEN
     DEALLOCATE (natom_of_kind ,STAT=stat)
     IF (stat /= 0) CALL stop_memory(routine,"natom_of_kind",nkind)
  END IF

 CALL timestop ( zero, handle )
 CALL write_checkpoint_information("leaving "//routine,globenv)

END SUBROUTINE kgpol_calculate_forces


!******************************************************************************
!!****** kg_force/kg_calculate_forces [1.1] *
!!
!!   NAME
!!     kg_calculate_forces
!!
!!   FUNCTION
!!     Calculates the total potential energy, total force
!!
!!   AUTHOR
!!     gt         
!!
!!   MODIFICATION HISTORY
!!
!!*** *************************************************************************

SUBROUTINE kg_calculate_forces ( kg_env, globenv,error )

    TYPE(kg_environment_type), POINTER       :: kg_env
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error


    CHARACTER(LEN=*), PARAMETER :: &
      routine = "SUBROUTINE kg_calculate_forces (MODULE kg_force)"

    INTEGER                                  :: handle, i, &
                                                iatom, ikind, stat, natoms, &
                                                nkind, output_unit
    INTEGER, ALLOCATABLE, DIMENSION(:), SAVE :: natom_of_kind
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: calculate_forces, ionode,&
                                                failure
    REAL(KIND=dp)                                :: dvol, total_rho0_rspace
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(coeff_type), POINTER                :: rho_core, rho_r, rho_g,&
                                                v_gspace, &
                                                v_rspace, vxc_r
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: my_vxc_r, my_vxc_tau, &
                                                my_rho_r, my_rho_g,&
                                                my_tau
    TYPE(cp_subsystem_p_type), DIMENSION(1)  :: subsys
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(kg_energy_type), POINTER            :: energy
    TYPE(kg_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(kg_rho_type), POINTER               :: rho
    TYPE(kg_rspw_type), POINTER              :: rspw
    TYPE(particle_type), POINTER             :: part( : )
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool

!  TYPE ( debug_variables_type ), INTENT ( OUT ), OPTIONAL :: debug
!  CHARACTER(LEN=8) :: force_type
!------------------------------------------------------------------------------

  CALL write_checkpoint_information("entering "//routine,globenv)
  CALL timeset ( 'FORCE','I',' ',handle )
  failure=.FALSE.
  ionode=globenv%ionode
  output_unit=globenv%scr
  NULLIFY(force,energy,part,atomic_kind_set,atomic_kind,cell,dft_control,&
          rho,rspw,auxbas_pw_pool)
  NULLIFY(rho_r,rho_g,rho_core,v_gspace,v_rspace,vxc_r)
  NULLIFY(my_vxc_r,my_vxc_tau,my_rho_r,my_rho_g,my_tau,subsys(1)%subsys)
  CALL get_kg_env( kg_env=kg_env, atomic_kind_set=atomic_kind_set,&
                    dft_control=dft_control, energy=energy, cell=cell,&
                    subsys=subsys(1)%subsys,&
                    force=force,particle_set=part,rho=rho,rspw=rspw)
  CPPrecondition(ASSOCIATED(rho),cp_failure_level,routine,error,failure)
  CPPrecondition(ASSOCIATED(rspw),cp_failure_level,routine,error,failure)
  IF(.NOT.failure) THEN
    CALL kg_rho_get(kg_rho=rho, rho_r=rho_r, rho_g=rho_g, rho_core=rho_core)
    CALL kg_rspw_get(kg_rspw=rspw,auxbas_pw_pool=auxbas_pw_pool,error=error)
  END IF

  ALLOCATE(v_gspace,stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routine,error,failure)
  ALLOCATE(v_rspace,stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routine,error,failure)
  ALLOCATE(vxc_r,stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routine,error,failure)
  IF (.NOT.failure) THEN
    CALL pw_pool_init_coeff(auxbas_pw_pool, v_gspace, &
                            use_data = COMPLEXDATA1D,&
                            in_space = RECIPROCALSPACE, error=error)
    CALL coeff_zero(v_gspace)
    CALL pw_pool_init_coeff(auxbas_pw_pool, v_rspace,&
                            use_data = REALDATA3D,&
                            in_space = REALSPACE, error=error)
    CALL coeff_zero(v_rspace)
  END IF
 
  calculate_forces=.TRUE. 
  natoms = SIZE ( part )
  nkind = SIZE( atomic_kind_set )
  stat = 0
  dvol=v_rspace%pw%pw_grid%dvol

  IF ( .NOT. ALLOCATED ( natom_of_kind ) ) THEN
     ALLOCATE (natom_of_kind(nkind),STAT=stat)
     IF (stat /= 0) CALL stop_memory(routine,"natom_of_kind",nkind)
     CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                              natom_of_kind=natom_of_kind)
  END IF

!  IF (.not.ASSOCIATED(kg_env%sac_ppl)) THEN
     CALL build_kg_neighbor_lists(kg_env,globenv)
!  ELSE
!     CALL rebuild_neighbor_lists(kg_env,globenv)
!  END IF

! initialize forces and energies                    
 
  CALL init_kg_force(force)
  CALL init_kg_energy(energy)

  DO i = 1, natoms
     part ( i ) % f ( 1 ) = 0.0_dp
     part ( i ) % f ( 2 ) = 0.0_dp
     part ( i ) % f ( 3 ) = 0.0_dp
  END DO

! self energy 
  CALL calculate_eself ( kg_env, energy % core_self )

! Calculate the rspace energy of the core charge distribution
! and its force contributions

  CALL calculate_ecore_rspace(kg_env,globenv,calculate_forces)

! calculate intramolecular bonded_correction energy  

  CALL calculate_ebond_corr(kg_env,globenv,energy%bond_corr,calculate_forces)

! Calculate the local pseudopotential energy         

  CALL calculate_rho0_ppl(kg_env,globenv,calculate_forces)

! Calculate the density generated by the core charges

  CALL calculate_epc_density(rho_core, rho%total_rho_core_rspace, kg_env)

! calculate the frozen electronic pseudo density

  CALL calculate_density (rho_r, rho_g, total_rho0_rspace, kg_env)

! from the frozen rho, calculate the thomas-fermi and xc potentials
  ALLOCATE(my_rho_r(1))
  my_rho_r(1)%pw => rho_r%pw
  ALLOCATE(my_rho_g(1))
  my_rho_g(1)%pw => rho_g%pw
  CALL xc_vxc_pw_create1(my_vxc_r,my_vxc_tau,my_rho_r,my_rho_g,my_tau,&
       energy%exc,dft_control,cell,auxbas_pw_pool,subsys,&
       calculate_forces=.FALSE.)
! ALLOCATE(vxc_r)
  CALL coeff_init(vxc_r,pw=my_vxc_r(1)%pw)
  NULLIFY(my_rho_r(1)%pw, my_vxc_r(1)%pw,my_rho_g(1)%pw)
  DEALLOCATE(my_rho_r, my_rho_g, my_vxc_r)
  NULLIFY(my_rho_r, my_rho_g, my_vxc_r)

  vxc_r%pw%cr3d(:,:,:) =dvol*vxc_r%pw%cr3d(:,:,:)

! calculate the total density
  rho%total_rho_rspace = total_rho0_rspace+ rho%total_rho_core_rspace

! add the core density to the frozen electronic density

  CALL coeff_sumup(rho_core,rho_g)
  rho%total_rho_gspace = calculate_total_density(rho_g)

! calculate electrostatic potential
  CALL init_pw_poisson_solver(cell%green)
  CALL calculate_hartree(rho_g, energy % hartree,v_gspace)
  CALL cleanup_pw_poisson_solver(cell%green)

!  calclation of the forces on the ions
  CALL coeff_transform_space(v_gspace,v_rspace)
  CALL pw_pool_give_back_coeff(auxbas_pw_pool,v_gspace,error=error)
  DEALLOCATE(v_gspace,stat=stat)
  CPPostconditionNoFail(stat==0,cp_warning_level,routine,error)

  v_rspace%pw%cr3d(:,:,:) = dvol*v_rspace%pw%cr3d(:,:,:)

  CALL calculate_epc_rspace_forces(v_rspace, kg_env)

  v_rspace%pw%cr3d(:,:,:) = v_rspace%pw%cr3d(:,:,:)+ vxc_r%pw%cr3d(:,:,:)

  CALL calculate_v_rspace_forces(v_rspace,kg_env) 

! deallocate work storage
  CALL pw_pool_give_back_coeff(auxbas_pw_pool,v_rspace,error=error)
  DEALLOCATE(v_rspace,stat=stat)
  CPPostconditionNoFail(stat==0,cp_warning_level,routine,error)
  CALL pw_pool_give_back_coeff(auxbas_pw_pool,vxc_r,error=error)
  DEALLOCATE(vxc_r,stat=stat)
  CPPostconditionNoFail(stat==0,cp_warning_level,routine,error)

! add up all the potential energies

  energy % total = energy % core_self + energy % core_overlap + energy % exc + &
                 energy % hartree + energy % pseudo + energy % bond_corr
     
! print
!    IF (ionode) THEN
    IF (ionode.AND.(globenv%print%level>3)) THEN
       WRITE (UNIT=output_unit,FMT="(/,(T3,A,T40,F20.10))")&
            "Total electronic density (r-space): ",&
            total_rho0_rspace,  &
            "Total core charge density (r-space):",&
            rho%total_rho_core_rspace
       WRITE (UNIT=output_unit,FMT="(T3,A,T60,F20.10)")&
            "Total charge density (r-space):     ",&
            rho%total_rho_rspace
       WRITE (UNIT=output_unit,FMT="(/,(T3,A,T55,F25.14))")&
            "Overlap energy of the core charge distribution:",energy%core_overlap,&
            "Self energy of the core charge distribution:   ",energy%core_self,&
            "Pseudopotential  energy:                       ",energy%pseudo,&
            "Hartree energy:                                ",energy%hartree,&
            "Exchange-correlation energy:                   ",energy%exc,&
            "Bonded--correction energy:                     ",energy%bond_corr,&
            "Total energy:                                  ",energy%total
       CALL m_flush(output_unit)
    END IF

! add up all the forces

  DO ikind = 1, nkind         
    CALL mp_sum(force(ikind)%f_rho,globenv%group)
    CALL mp_sum(force(ikind)%f_rspace_core,globenv%group)
    CALL mp_sum(force(ikind)%f_hartree_core,globenv%group)
    CALL mp_sum(force(ikind)%f_ppl,globenv%group)
    CALL mp_sum(force(ikind)%f_bc,globenv%group)
    force(ikind)%f_total(:,:) = force(ikind)%f_rho(:,:) +&
                                force(ikind)%f_rspace_core(:,:) +&
                                force(ikind)%f_hartree_core(:,:) +&
                                force(ikind)%f_ppl(:,:) +&
                                force(ikind)%f_bc(:,:)
  END DO


  DO ikind = 1, nkind         
    atomic_kind => atomic_kind_set(ikind)
    CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           atom_list=atom_list) 
    DO iatom=1, SIZE(atom_list)
      i = atom_list(iatom)
! fill in the forces in part
      part(i)%f(1) = force(ikind)%f_total(1,iatom)
      part(i)%f(2) = force(ikind)%f_total(2,iatom)
      part(i)%f(3) = force(ikind)%f_total(3,iatom)
    END DO
  END DO

! deallocating all local variables
  IF ( ALLOCATED ( natom_of_kind ) ) THEN
     DEALLOCATE (natom_of_kind ,STAT=stat)
     IF (stat /= 0) CALL stop_memory(routine,"natom_of_kind",nkind)
  END IF

 CALL timestop ( zero, handle )
 CALL write_checkpoint_information("leaving "//routine,globenv)

END SUBROUTINE kg_calculate_forces

!******************************************************************************

END MODULE kg_force  
