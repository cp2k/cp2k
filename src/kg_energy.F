!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C ) 2002 - 2003  CP2K developers group                         !
!-----------------------------------------------------------------------------!
#include "cp_prep_globals.h"
!!****** cp2k/kg_energy [1.0] *
!!
!!   NAME
!!     kg_energy
!!
!!   FUNCTION
!!     kg_energy minim routine      
!!
!!   AUTHOR
!!     gloria,30.12.2003 
!!
!!   MODIFICATION HISTORY
!!   SOURCE
!******************************************************************************

MODULE kg_energy 
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE checkpoint_handler,              ONLY: write_checkpoint_information
  USE coefficient_types,               ONLY: coeff_copy,&
                                             coeff_sumup,&
                                             coeff_transform_space,&
                                             coeff_type,&
                                             coeff_init,&
                                             coeff_zero
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
                                             dyn_coeff_type,&
                                             get_dyn_coeff,&
                                             read_coeff_restart,&
                                             write_coeff_restart
  USE global_types,                    ONLY: global_environment_type
  USE hartree,                         ONLY: calculate_hartree
  USE kg_core,                         ONLY: calculate_ecore_rspace,&
                                             calculate_eself
  USE kg_density,                      ONLY: calculate_density,&
                                             calculate_epc_density,&
                                             calculate_p_density,&
                                             calculate_total_density,&
                                             calculate_vp_rspace_forces
  USE kg_energy_types,                 ONLY: init_kg_energy,&
                                             kg_energy_type
  USE kg_energy_utils,                 ONLY: diis_info_create,&
                                             diis_info_destroy,&
                                             diis_info_type,&
                                             do_diis_step
  USE kg_intra,                        ONLY: calculate_ebond_corr

  USE kg_environment_types,            ONLY: get_kg_env,&
                                             kg_environment_type
  USE kg_environment_methods,          ONLY: kg_qs_env_update
  USE kg_kxc,                          ONLY: calculate_kxc_derivatives
  USE kg_neighbor_lists,               ONLY: build_kg_neighbor_lists
  USE kg_ppl,                          ONLY: calculate_drho_ppl,&
                                             calculate_rho0_ppl
  USE kg_rho_types,                    ONLY: kg_rho_type, &
                                             kg_rho_get
  USE kg_rspw_types,                   ONLY: kg_rspw_type,&
                                             kg_rspw_get
  USE kg_scf,                          ONLY: kg_qs_scf
  USE kinds,                           ONLY: dbl
  USE machine,                         ONLY: m_flush,&
                                             m_walltime
  USE mathconstants,                   ONLY: zero
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_max
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
  USE qs_core_hamiltonian,             ONLY: build_core_hamiltonian_matrix
  USE qs_environment_methods,          ONLY: qs_env_update_s_mstruct
  USE qs_environment_types,            ONLY: qs_environment_type
  USE qs_neighbor_lists,               ONLY: build_qs_neighbor_lists
  USE scf_control_types,               ONLY: scf_control_type
  USE simulation_cell,                 ONLY: cell_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xc,                              ONLY: xc_vxc_pw_create1

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: kg_energies

!******************************************************************************

CONTAINS

!******************************************************************************
SUBROUTINE kg_energies ( kg_env, globenv )

    TYPE(kg_environment_type), POINTER       :: kg_env
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    TYPE(dft_control_type), POINTER          :: dft_control

!------------------------------------------------------------------------------

    NULLIFY ( dft_control )

    CALL get_kg_env(kg_env=kg_env,&
                    dft_control=dft_control)

    IF ( dft_control % qs_control % method == "KG_GPW" ) THEN

      CALL kg_qs_energies ( kg_env, globenv )

    ELSE

      CALL kg_std_energies ( kg_env, globenv )

    END IF

END SUBROUTINE kg_energies

!******************************************************************************
!!****** kg_energy/kg_std_energies [1.1] *
!!
!!   NAME
!!     kg_std_energies
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
SUBROUTINE kg_std_energies ( kg_env, globenv, error )

    TYPE(kg_environment_type), POINTER       :: kg_env
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "SUBROUTINE kg_std_energies (MODULE kg_energy)"

    INTEGER                                  :: handle, handle2, i, stat, &
                                                iter, output_unit
    LOGICAL                                  :: ionode, should_stop, failure
    REAL(dbl)                                :: maxerr, maxgrad, t1, t2
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE (coeff_type), POINTER               :: vxc_r
    TYPE(diis_info_type), POINTER            :: diis_info
    TYPE(dyn_coeff_set_type), POINTER        :: dyn_coeff_set
    TYPE(kg_rspw_type), POINTER              :: rspw
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool

!------------------------------------------------------------------------------

  CALL write_checkpoint_information("entering "//routine,globenv)
  CALL timeset ( 'kg_std_energies','I',' ',handle )

  NULLIFY(atomic_kind_set,dyn_coeff_set,diis_info,scf_control)
  NULLIFY(rspw,auxbas_pw_pool,vxc_r)
  ionode=globenv%ionode
  output_unit=globenv%scr
  stat = 0
  iter = 0

  t1 = m_walltime()

  CALL get_kg_env(kg_env=kg_env,&
                  atomic_kind_set=atomic_kind_set,&
                  dyn_coeff_set=dyn_coeff_set,&
                  rspw=rspw,&
                  scf_control=scf_control)           

  CALL kg_rspw_get(kg_rspw=rspw,auxbas_pw_pool=auxbas_pw_pool,error=error)
  ALLOCATE(vxc_r,stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routine,error,failure)
   
  CALL build_kg_neighbor_lists(kg_env,globenv)

  CALL diis_info_create(diis_info,scf_control%max_diis,dyn_coeff_set)

  CALL init_energy_calculation(kg_env,vxc_r,globenv)

  CALL calc_c_gradients(kg_env,vxc_r,maxgrad,globenv)

  IF (ionode) THEN
     WRITE (UNIT=output_unit,&
        FMT="(/,T3,A,T16,A,T49,A,T68,A,/,T3,A)")&
        "Step","maxerr","Time","Convergence","Total energy",&
         REPEAT("-",77)
  END IF

scf_loop: DO                      

    CALL timeset("scf_iter","I","",handle2)
    IF (ionode) CALL m_flush(output_unit)
    iter = iter + 1

    CALL do_diis_step(kg_env,diis_info,maxerr)
    CALL calc_c_gradients(kg_env,vxc_r,maxgrad,globenv)

    t2 = m_walltime()

    IF (ionode) THEN
       WRITE (UNIT=output_unit,&
          FMT="(T2,I5,2X,E10.2,T32,F8.2,T40,2F20.10)")&
          iter,maxerr,t2 - t1,maxgrad,kg_env%energy%total
    END IF

! ** convergence check
    IF (kg_env%para_env%source==kg_env%para_env%mepos) &
            INQUIRE(file="EXIT",exist=should_stop)
    CALL mp_bcast(should_stop,kg_env%para_env%source,kg_env%para_env%group)
    IF (maxgrad < scf_control%eps_scf) THEN
      IF (ionode) THEN
         WRITE(UNIT=output_unit,FMT="(/,T3,A,I5,A/)")&
               "*** SCF run converged in", iter, "steps ***"
      END IF
      CALL timestop(0.0_dbl,handle2)
      EXIT scf_loop
    ELSE IF (should_stop.OR.&
        iter == scf_control%max_scf) THEN
      IF (ionode) THEN
             WRITE(UNIT=output_unit,FMT="(/,T3,A,/)")&
                  "*** SCF run NOT converged ***"
      END IF
      CALL timestop(0.0_dbl,handle2)
      EXIT 
    END IF

!   *** Write restart file ***
    IF (scf_control%write_restart_each>0.AND.&
         MODULO(iter,scf_control%write_restart_each)==0) THEN
       CALL write_coeff_restart(dyn_coeff_set,atomic_kind_set,globenv)  
    END IF

    t1 = m_walltime()

    CALL timestop(0.0_dbl,handle2)

  END DO scf_loop

  CALL pw_pool_give_back_coeff(auxbas_pw_pool,vxc_r,error=error)
  DEALLOCATE(vxc_r,stat=stat)
  CPPostconditionNoFail(stat==0,cp_warning_level,routine,error)
  
!   *** Write restart file ***
  IF (scf_control%write_restart_each>=0) THEN
    CALL write_coeff_restart(dyn_coeff_set,atomic_kind_set,globenv)  
  END IF
  
  CALL diis_info_destroy(diis_info)

  CALL timestop(0.0_dbl,handle)

  CALL write_checkpoint_information("leaving "//routine,globenv)

END SUBROUTINE kg_std_energies

!******************************************************************************
!!****** kg_energy/kg_qs_energies [1.1] *
!!
!!   NAME
!!     kg_qs_energies
!!
!!   FUNCTION
!!     Calculates the total potential energy, total force
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!
!!*** *************************************************************************
SUBROUTINE kg_qs_energies ( kg_env, globenv )

    TYPE(kg_environment_type), POINTER       :: kg_env
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    CHARACTER(LEN=*), PARAMETER :: routine = "kg_qs_energies"

    INTEGER                                  :: handle
    TYPE(qs_environment_type), POINTER       :: qs_env

!------------------------------------------------------------------------------

    CALL write_checkpoint_information("entering "//routine,globenv)
    CALL timeset(routine,"I","",handle)

    NULLIFY ( qs_env )
    CALL get_kg_env(kg_env=kg_env,&
                    sub_qs_env=qs_env)

    CALL build_qs_neighbor_lists(qs_env,globenv,.TRUE.)

    ! Calculate the overlap and the core Hamiltonian integral matrix
    CALL build_core_hamiltonian_matrix(qs_env=qs_env,&
                                       globenv=globenv,&
                                       calculate_forces=.FALSE.)

    CALL kg_qs_env_update(kg_env,qs_env)

    ! Perform a SCF run
    CALL kg_qs_scf(kg_env,globenv)

    CALL timestop(0.0_dbl,handle)
    CALL write_checkpoint_information("leaving "//routine,globenv)

END SUBROUTINE kg_qs_energies

!------------------------------------------------------------------------------

SUBROUTINE init_energy_calculation ( kg_env, vxc_r, globenv, error)
  ! vxc_r is allocated here

    TYPE(kg_environment_type), POINTER       :: kg_env
    TYPE(coeff_type),POINTER                 :: vxc_r
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "SUBROUTINE init_energy_calculation (MODULE kg_energy)"

    INTEGER                                  :: handle, i, stat, &
                                                output_unit
    LOGICAL                                  :: calculate_forces, ionode, failure
    REAL(dbl)                                :: dvol, total_rho0_rspace
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(coeff_type), POINTER                :: rho_core, rho_r, rho_g,&
                                                rhop_r, rhop_g
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: my_vxc_r, my_vxc_tau, &
                                                my_rho_r, my_rho_g,&
                                                my_tau
    TYPE(cp_subsystem_p_type), &
         DIMENSION(:), POINTER               :: my_atoms
    TYPE(cell_type), POINTER                 :: cell
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(dyn_coeff_set_type), POINTER        :: dyn_coeff_set
    TYPE(kg_energy_type), POINTER            :: energy
    TYPE(kg_rho_type), POINTER               :: rho
    TYPE(kg_rspw_type), POINTER              :: rspw
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(scf_control_type), POINTER          :: scf_control

!------------------------------------------------------------------------------

  CALL write_checkpoint_information("entering "//routine,globenv)
  CALL timeset ( 'init_energy_calculation','I',' ',handle )
  ionode=globenv%ionode
  output_unit=globenv%scr
  NULLIFY(auxbas_pw_pool,energy,atomic_kind_set,dft_control,&
          cell,rho,rspw,scf_control,dyn_coeff_set)
  NULLIFY(rho_r,rho_g,rho_core)

  CALL get_kg_env( kg_env=kg_env, atomic_kind_set=atomic_kind_set,&
                    dft_control=dft_control, energy=energy, cell=cell,&
                    rho=rho,rspw=rspw,dyn_coeff_set=dyn_coeff_set,&
                    scf_control=scf_control)

  CPPrecondition(ASSOCIATED(rho),cp_failure_level,routine,error,failure)
  CPPrecondition(ASSOCIATED(rspw),cp_failure_level,routine,error,failure)
  IF(.NOT.failure) THEN
    CALL kg_rho_get(kg_rho=rho, rho_r=rho_r, rho_g=rho_g, &
                    rho_core=rho_core)
    CALL kg_rspw_get(kg_rspw=rspw,auxbas_pw_pool=auxbas_pw_pool,error=error)
  END IF

  calculate_forces=.FALSE. 
  stat = 0
  dvol=auxbas_pw_pool%pw_grid%dvol

  CALL build_kg_neighbor_lists(kg_env,globenv)

! initialize energies                    
 
  CALL init_kg_energy(energy)

! calculate the constant parts of the energy 
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

  CALL calculate_epc_density(rho_core,rho%total_rho_core_rspace, kg_env)

! calculate the frozen electronic pseudo density

  CALL calculate_density (rho_r, rho_g, total_rho0_rspace, kg_env)

! calculate the total density
  rho%total_rho_rspace = total_rho0_rspace+ rho%total_rho_core_rspace

! from the frozen rho, calculate the thomas-fermi and xc potentials
!  CALL xc_calculate_pw1 (rho_r%pw, rho_g%pw, vxc_r%pw, energy%exc, &
!                         dft_control,auxbas_pw_pool)
  ALLOCATE(my_rho_r(1))
  my_rho_r(1)%pw => rho_r%pw
  CALL xc_vxc_pw_create1(my_vxc_r,my_vxc_tau,my_rho_r,my_rho_g,my_tau,&
       energy%exc,dft_control,cell,auxbas_pw_pool,my_atoms,&
       calculate_forces=.FALSE.)
  allocate(vxc_r)
  CALL coeff_init(vxc_r,pw=my_vxc_r(1)%pw)
  NULLIFY(my_rho_r(1)%pw, my_vxc_r(1)%pw)
  DEALLOCATE(my_rho_r, my_vxc_r)
  NULLIFY(my_rho_r, my_vxc_r)

  vxc_r%pw%cr3d(:,:,:) =dvol*vxc_r%pw%cr3d(:,:,:)
! add the frozen electronic density to the core density 
  CALL coeff_sumup(rho_g,rho_core)

!  read coef from restart file
  IF (scf_control%density_guess == "RESTART") THEN
    CALL read_coeff_restart(dyn_coeff_set,globenv)
  END IF    
  CALL timestop(0.0_dbl,handle)
  CALL write_checkpoint_information("leaving "//routine,globenv)

END SUBROUTINE init_energy_calculation 
!------------------------------------------------------------------------------

SUBROUTINE calc_c_gradients ( kg_env, vxc_r, maxgrad, globenv,error )


    TYPE(kg_environment_type), POINTER       :: kg_env
    TYPE(coeff_type), POINTER                :: vxc_r
    REAL(dbl), intent(inout)                 :: maxgrad
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "SUBROUTINE calc_c_gradients (MODULE kg_energy)"

    INTEGER                                  :: handle, i, &
                                                ikind, nkind, output_unit, stat
    LOGICAL                                  :: calculate_forces, ionode,&
                                                failure
    REAL(dbl)                                :: dvol, e1_xc, e2_xc, &
                                                total_rhop_rspace
    REAL(dbl), DIMENSION(:, :), POINTER      :: forces
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(coeff_type), POINTER                :: rho_core, rho_r, rho_g,&
                                                rhop_r, rhop_g, &
                                                v_gspace, &
                                                v_rspace, &
                                                v1xc_r, v2xc_r
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(dyn_coeff_set_type), POINTER        :: dyn_coeff_set
    TYPE(dyn_coeff_type), POINTER            :: local_coeffs
    TYPE(kg_energy_type), POINTER            :: energy
    TYPE(kg_rho_type), POINTER               :: rho
    TYPE(kg_rspw_type), POINTER              :: rspw
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool

!------------------------------------------------------------------------------

  CALL write_checkpoint_information("entering "//routine,globenv)
  CALL timeset ( 'calc_c_gradients','I',' ',handle )

  ionode=globenv%ionode
  output_unit=globenv%scr
  calculate_forces=.FALSE.

  NULLIFY(energy,atomic_kind_set,cell,dft_control,&
          dyn_coeff_set,rho,rspw,local_coeffs,forces)
  NULLIFY(rho_r,rho_g,rhop_r,rhop_g,rho_core,v_gspace,v_rspace,&
          v1xc_r,v2xc_r)

  CALL get_kg_env( kg_env=kg_env, atomic_kind_set=atomic_kind_set,&
                    dft_control=dft_control, energy=energy, cell=cell,&
                    rho=rho,rspw=rspw,dyn_coeff_set=dyn_coeff_set)
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
        
  dvol=v_rspace%pw%pw_grid%dvol           
  nkind=SIZE(dyn_coeff_set%coeffs_of_kind)
! reinitialize forces
  DO ikind= 1, nkind
    local_coeffs => dyn_coeff_set%coeffs_of_kind(ikind)%coeffs
    IF(ASSOCIATED(local_coeffs)) THEN
      CALL get_dyn_coeff(coeffs=local_coeffs,forces=forces)
      forces(:,:)=0.0_dbl
    END IF
  END DO
  CALL coeff_zero(rho_g)
  CALL coeff_copy(rho_core,rho_g)
  rho%total_rho_gspace = calculate_total_density(rho_g)

  CALL calculate_drho_ppl(kg_env,globenv,calculate_forces)

  CALL calculate_p_density(rhop_r, rhop_g, total_rhop_rspace, kg_env)

  CALL calculate_kxc_derivatives(auxbas_pw_pool, rho_r, rho_g, rhop_r, vxc_r, v1xc_r,&
                                 v2xc_r, e1_xc, e2_xc, dft_control, calculate_forces)

! add the polarization density to get the total density
  CALL coeff_sumup(rhop_g,rho_g)

! calculate electrostatic potential
  CALL init_pw_poisson_solver(cell%green)
  CALL calculate_hartree(rho_g, energy % hartree,v_gspace)        
  CALL cleanup_pw_poisson_solver(cell%green)

  CALL coeff_transform_space(v_gspace,v_rspace)
  CALL pw_pool_give_back_coeff(auxbas_pw_pool,v_gspace,error=error)
  DEALLOCATE(v_gspace,stat=stat)
  CPPostconditionNoFail(stat==0,cp_warning_level,routine,error)

  v_rspace%pw%cr3d(:,:,:) = dvol*v_rspace%pw%cr3d(:,:,:)
  v_rspace%pw%cr3d(:,:,:) =-1.0_dbl*(vxc_r%pw%cr3d(:,:,:) +&
               v1xc_r%pw%cr3d(:,:,:)+ v_rspace%pw%cr3d(:,:,:))

  CALL calculate_vp_rspace_forces(v_rspace,kg_env, calculate_forces)            

! deallocate work storage
  CALL pw_pool_give_back_coeff(auxbas_pw_pool,v_rspace,error=error)
  DEALLOCATE(v_rspace,stat=stat)
  CPPostconditionNoFail(stat==0,cp_warning_level,routine,error)
  CALL pw_pool_give_back_coeff(auxbas_pw_pool,v1xc_r,error=error)
  DEALLOCATE(v1xc_r,stat=stat)
  CPPostconditionNoFail(stat==0,cp_warning_level,routine,error)
  CALL pw_pool_give_back_coeff(auxbas_pw_pool,v2xc_r,error=error)
  DEALLOCATE(v2xc_r,stat=stat)
  CPPostconditionNoFail(stat==0,cp_warning_level,routine,error)

! add up all the potential energies

  energy % total = energy % core_self + energy % core_overlap + energy % exc + &
                 e1_xc + e2_xc+ energy % hartree + energy % pseudo + energy%ppseudo + &
                 energy%bond_corr
     
! get the maximum coefficent gradient 

  maxgrad=0._dbl
  DO ikind= 1, nkind
    local_coeffs => dyn_coeff_set%coeffs_of_kind(ikind)%coeffs
    IF(ASSOCIATED(local_coeffs)) THEN
      CALL get_dyn_coeff(coeffs=local_coeffs,forces=forces)      
      maxgrad=MAX(maxgrad,MAXVAL(ABS(forces)))
    END IF
  END DO
  CALL mp_max(maxgrad,kg_env%para_env%group)

! print
    IF (ionode) THEN
       WRITE (UNIT=output_unit,FMT="(/,(T3,A,T40,F20.10))")&
            "Total core charge density (r-space):",&
            rho%total_rho_core_rspace
       WRITE (UNIT=output_unit,FMT="(T3,A,T60,F20.10)")&
            "Total charge density (r-space):     ",&
            rho%total_rho_rspace
       WRITE (UNIT=output_unit,FMT="(T3,A,T60,F20.10)")&
            "Total rho_coefs charge density (r-space):     ",&
            total_rhop_rspace
       WRITE (UNIT=output_unit,FMT="(T3,A,T60,F20.10)")&
            "Maximum coefficent gradient         ",&
            maxgrad                  
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

 CALL timestop ( zero, handle )
 CALL write_checkpoint_information("leaving "//routine,globenv)

END SUBROUTINE calc_c_gradients

!******************************************************************************
END MODULE kg_energy  
