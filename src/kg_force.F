!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief kg force calculation routine
!> \par History
!>      JGH (22-Feb-03) PW grid options added
!> \author gloria,30.09.2002
! *****************************************************************************
MODULE kg_force
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE cell_types,                      ONLY: cell_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_operations,             ONLY: dbcsr_allocate_matrix_set,&
                                             dbcsr_deallocate_matrix_set,&
                                             dbcsr_from_sm,&
                                             dbcsr_obj_p_type,&
                                             sm_from_dbcsr
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE dynamical_coeff_types,           ONLY: dyn_coeff_set_type,&
                                             dyn_coeff_type
  USE f77_blas
  USE global_types,                    ONLY: global_environment_type
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type
  USE kg_core,                         ONLY: calculate_ecore_rspace,&
                                             calculate_eself
  USE kg_density,                      ONLY: calculate_density,&
                                             calculate_epc_density,&
                                             calculate_epc_rspace_forces,&
                                             calculate_p_density,&
                                             calculate_v_rspace_forces,&
                                             calculate_vp_rspace_forces
  USE kg_energy,                       ONLY: kg_qs_energies
  USE kg_energy_types,                 ONLY: init_kg_energy,&
                                             kg_energy_type
  USE kg_environment_types,            ONLY: get_kg_env,&
                                             kg_environment_type
  USE kg_force_types,                  ONLY: init_kg_force,&
                                             kg_force_type
  USE kg_gpw_fm_mol_types,             ONLY: kg_fm_p_type
  USE kg_gpw_fm_mol_utils,             ONLY: calculate_w_matrix_per_molecule
  USE kg_intra,                        ONLY: calculate_ebond_corr
  USE kg_kxc,                          ONLY: calculate_kxc_derivatives
  USE kg_neighbor_lists,               ONLY: build_kg_neighbor_lists
  USE kg_ppl,                          ONLY: calculate_drho_ppl,&
                                             calculate_rho0_ppl
  USE kg_rho_types,                    ONLY: kg_rho_get,&
                                             kg_rho_type
  USE kg_rspw_types,                   ONLY: kg_rspw_get,&
                                             kg_rspw_type
  USE kinds,                           ONLY: dp
  USE machine,                         ONLY: m_flush
  USE message_passing,                 ONLY: mp_sum
  USE mulliken,                        ONLY: mulliken_restraint
  USE particle_types,                  ONLY: particle_type
  USE pw_methods,                      ONLY: pw_axpy,&
                                             pw_integrate_function,&
                                             pw_transfer,&
                                             pw_zero
  USE pw_poisson_methods,              ONLY: pw_poisson_solve
  USE pw_poisson_types,                ONLY: pw_poisson_type
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_p_type
  USE qs_core_energies,                ONLY: calculate_ecore_overlap
  USE qs_core_hamiltonian,             ONLY: build_core_hamiltonian_matrix
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_force,                        ONLY: write_forces
  USE qs_force_types,                  ONLY: qs_force_type,&
                                             zero_qs_force
  USE qs_ks_methods,                   ONLY: qs_ks_update_qs_env
  USE qs_ks_types,                     ONLY: qs_ks_env_type
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE scf_control_types,               ONLY: scf_control_type
  USE sparse_matrix_output,            ONLY: write_sparse_matrix
  USE sparse_matrix_types,             ONLY: allocate_matrix_set,&
                                             deallocate_matrix_set,&
                                             real_matrix_p_type,&
                                             replicate_matrix_structure,&
                                             set_matrix
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xc,                              ONLY: xc_vxc_pw_create1
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'kg_force'

  PUBLIC :: kg_force_control

CONTAINS

! *****************************************************************************
!> \brief Calculates the total potential energy, total force, and the
!>      total pressure tensor from the potentials
!> \author gt
! *****************************************************************************
SUBROUTINE kg_force_control ( kg_env, globenv, error )

    TYPE(kg_environment_type), POINTER       :: kg_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout)       :: error

    LOGICAL                                  :: kgpol
    TYPE(dft_control_type), POINTER          :: dft_control

  kgpol = .FALSE.
  CALL get_kg_env( kg_env=kg_env, dft_control=dft_control,error=error)
  kgpol = dft_control % qs_control % polarization

  IF (kgpol) THEN
    CALL kgpol_calculate_forces ( kg_env,  error)
  ELSEIF ( dft_control % qs_control % method == "KG_GPW" ) THEN
    ! The qs_forces are calculated therefore kg_forces are not allocated here
    CALL kg_qs_calculate_forces( kg_env, globenv, error )
  ELSE
    CALL kg_calculate_forces ( kg_env, error )
  END IF
  END SUBROUTINE kg_force_control

! *****************************************************************************
!> \brief Calculates the total potential energy, total force
!> \author gt
! *****************************************************************************
SUBROUTINE kgpol_calculate_forces ( kg_env, error )

    TYPE(kg_environment_type), POINTER       :: kg_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'kgpol_calculate_forces', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, iatom, ikind, &
                                                natoms, nkind, output_unit, &
                                                stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: natom_of_kind
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: calculate_forces, failure
    REAL(KIND=dp)                            :: dvol, e1_xc, e2_xc, &
                                                total_rho0_rspace, &
                                                total_rhop_rspace
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(dyn_coeff_set_type), POINTER        :: dyn_coeff_set
    TYPE(dyn_coeff_type), POINTER            :: local_coeffs
    TYPE(kg_energy_type), POINTER            :: energy
    TYPE(kg_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(kg_rho_type), POINTER               :: rho
    TYPE(kg_rspw_type), POINTER              :: rspw
    TYPE(particle_type), POINTER             :: part( : )
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: my_rho_g, my_rho_r, my_tau, &
                                                my_vxc_r, my_vxc_tau
    TYPE(pw_p_type), POINTER                 :: rho_core, rho_g, rho_r, &
                                                rhop_g, rhop_r, v1xc_r, &
                                                v2xc_r, v_gspace, v_rspace, &
                                                vxc_r
    TYPE(pw_poisson_type), POINTER           :: poisson_env
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(section_vals_type), POINTER         :: xc_section

  CALL timeset ( routineN,handle )
  failure=.FALSE.
  para_env=>kg_env%para_env
  logger=>cp_error_get_logger(error)
  ! better to hook it to a print_key in the input
  output_unit=cp_logger_get_default_io_unit(logger)

  NULLIFY(force,energy,part,atomic_kind_set,atomic_kind,cell,dft_control,&
          dyn_coeff_set,local_coeffs,rho,rspw,auxbas_pw_pool)
  NULLIFY(rho_r,rho_g,rhop_r,rhop_g,rho_core,v_gspace,v_rspace,&
          vxc_r,v1xc_r,v2xc_r)
  NULLIFY(my_vxc_r,my_vxc_tau,my_rho_r,my_rho_g,my_tau, poisson_env)

  xc_section => section_vals_get_subs_vals(kg_env%input,&
       "DFT%XC",error=error)

  CALL get_kg_env( kg_env=kg_env, atomic_kind_set=atomic_kind_set,&
                    dft_control=dft_control, energy=energy, cell=cell,&
                    dyn_coeff_set=dyn_coeff_set, force=force,particle_set=part,rho=rho,rspw=rspw,error=error)

  CPPrecondition(ASSOCIATED(rho),cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(rspw),cp_failure_level,routineP,error,failure)
  IF(.NOT.failure) THEN
    CALL kg_rho_get(kg_rho=rho, rho_r=rho_r, rho_g=rho_g, rhop_r=rhop_r,&
                    rhop_g=rhop_g,rho_core=rho_core,error=error)
    CALL kg_rspw_get(kg_rspw=rspw,auxbas_pw_pool=auxbas_pw_pool,&
         poisson_env=poisson_env,error=error)
  END IF

  ALLOCATE(v_gspace,stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  ALLOCATE(v_rspace,stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  ALLOCATE(vxc_r,stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  ALLOCATE(v1xc_r,stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  ALLOCATE(v2xc_r,stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

  IF (.NOT.failure) THEN
    CALL pw_pool_create_pw(auxbas_pw_pool, v_gspace%pw, &
                            use_data = COMPLEXDATA1D,&
                            in_space = RECIPROCALSPACE, error=error)
    CALL pw_zero(v_gspace%pw,error=error)
    CALL pw_pool_create_pw(auxbas_pw_pool, v_rspace%pw,&
                            use_data = REALDATA3D,&
                            in_space = REALSPACE, error=error)
    CALL pw_zero(v_rspace%pw,error=error)
    CALL pw_pool_create_pw(auxbas_pw_pool, v1xc_r%pw,&
                            use_data = REALDATA3D,&
                            in_space = REALSPACE, error=error)
    CALL pw_zero(v1xc_r%pw,error=error)
    CALL pw_pool_create_pw(auxbas_pw_pool, v2xc_r%pw,&
                            use_data = REALDATA3D,&
                            in_space = REALSPACE, error=error)
    CALL pw_zero(v2xc_r%pw,error=error)
  END IF

  calculate_forces=.TRUE.
  natoms = SIZE ( part )
  nkind = SIZE( atomic_kind_set )
  stat = 0
  dvol=v_rspace%pw%pw_grid%dvol

  IF ( .NOT. ALLOCATED ( natom_of_kind ) ) THEN
     ALLOCATE (natom_of_kind(nkind),STAT=stat)
     IF (stat /= 0) CALL stop_memory(routineP,"natom_of_kind",nkind)
     CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                              natom_of_kind=natom_of_kind)
  END IF

  CALL build_kg_neighbor_lists(kg_env,para_env,error)

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
  CALL calculate_eself ( kg_env, energy % core_self ,error=error)

! Calculate the rspace energy of the core charge distribution
! and its force contributions

  CALL calculate_ecore_rspace(kg_env,para_env,calculate_forces,error)

! calculate intramolecular bonded_correction energy

  CALL calculate_ebond_corr(kg_env,energy%bond_corr,calculate_forces,error=error)

! Calculate the local pseudopotential energy

  CALL calculate_rho0_ppl(kg_env,calculate_forces,error)
  CALL calculate_drho_ppl(kg_env,calculate_forces,error)

! Calculate the density generated by the core charges

  CALL calculate_epc_density(rho_core,rho%total_rho_core_rspace, kg_env,error=error)

! calculate the frozen electronic pseudo density

  CALL calculate_density (rho_r, rho_g, total_rho0_rspace, kg_env,error=error)

  CALL calculate_p_density(rhop_r, rhop_g, total_rhop_rspace, kg_env,error=error)

  ALLOCATE(my_rho_r(1))
  my_rho_r(1)%pw => rho_r%pw
  ALLOCATE(my_rho_g(1))
  my_rho_g(1)%pw => rho_g%pw
  CALL xc_vxc_pw_create1(my_vxc_r,my_vxc_tau,my_rho_r,my_rho_g,my_tau,&
       energy%exc,xc_section,cell,auxbas_pw_pool,error=error)

  vxc_r%pw => my_vxc_r(1)%pw
  NULLIFY(my_rho_r(1)%pw,my_rho_g(1)%pw, my_vxc_r(1)%pw)
  DEALLOCATE(my_rho_r,my_rho_g, my_vxc_r)
  NULLIFY(my_rho_r,my_rho_g, my_vxc_r)

  vxc_r%pw%cr3d(:,:,:) =dvol*vxc_r%pw%cr3d(:,:,:)

  CALL calculate_kxc_derivatives(auxbas_pw_pool, rho_r, rho_g, rhop_r, vxc_r, v1xc_r,&
                                 v2xc_r, e1_xc, e2_xc, xc_section, calculate_forces,error)

! calculate the total density
  rho%total_rho_rspace = total_rho0_rspace+ rho%total_rho_core_rspace

! add the core density to the frozen electronic density
  CALL pw_axpy(rho_core%pw,rho_g%pw,error=error)

  rho%total_rho_gspace = pw_integrate_function(rho_g%pw,error=error)

! add the polarization density to get the total density
  CALL pw_axpy(rhop_g%pw,rho_g%pw,error=error)

! calculate electrostatic potential
  CALL pw_poisson_solve(poisson_env,rho_g%pw, energy % hartree,v_gspace%pw,error=error)

!  calclation of the forces on the ions

  CALL pw_transfer(v_gspace%pw,v_rspace%pw,error=error)
  CALL pw_pool_give_back_pw(auxbas_pw_pool,v_gspace%pw,error=error)
  DEALLOCATE(v_gspace,stat=stat)
  CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)

  v_rspace%pw%cr3d(:,:,:) = dvol*v_rspace%pw%cr3d(:,:,:)

  CALL calculate_epc_rspace_forces(v_rspace, kg_env,error=error)

  v_rspace%pw%cr3d(:,:,:) =-1.0_dp*(vxc_r%pw%cr3d(:,:,:) +&
               v1xc_r%pw%cr3d(:,:,:)+ v_rspace%pw%cr3d(:,:,:))

  CALL calculate_vp_rspace_forces(v_rspace,kg_env, calculate_forces,error=error)

  v_rspace%pw%cr3d(:,:,:) = -1.0_dp*v_rspace%pw%cr3d(:,:,:)+&
                            v2xc_r%pw%cr3d(:,:,:)

  CALL calculate_v_rspace_forces(v_rspace,kg_env,error=error)

! deallocate work storage
  CALL pw_pool_give_back_pw(auxbas_pw_pool,v_rspace%pw,error=error)
  DEALLOCATE(v_rspace,stat=stat)
  CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
  CALL pw_pool_give_back_pw(auxbas_pw_pool,vxc_r%pw,error=error)
  DEALLOCATE(vxc_r,stat=stat)
  CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
  CALL pw_pool_give_back_pw(auxbas_pw_pool,v1xc_r%pw,error=error)
  DEALLOCATE(v1xc_r,stat=stat)
  CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
  CALL pw_pool_give_back_pw(auxbas_pw_pool,v2xc_r%pw,error=error)
  DEALLOCATE(v2xc_r,stat=stat)
  CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)

! add up all the potential energies

  energy % total = energy % core_self + energy % core_overlap + energy % exc + &
                 e1_xc + e2_xc + energy % hartree + energy % pseudo + energy%ppseudo + &
                 energy% bond_corr

! print
    IF (output_unit>0) THEN
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
    CALL mp_sum(force(ikind)%f_rho,para_env%group)
    CALL mp_sum(force(ikind)%f_rspace_core,para_env%group)
    CALL mp_sum(force(ikind)%f_hartree_core,para_env%group)
    CALL mp_sum(force(ikind)%f_ppl,para_env%group)
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
     IF (stat /= 0) CALL stop_memory(routineP,"natom_of_kind",nkind)
  END IF

 CALL timestop(handle)

END SUBROUTINE kgpol_calculate_forces

! *****************************************************************************
!> \brief Calculates the total potential energy, total force
!> \author gt
! *****************************************************************************
SUBROUTINE kg_calculate_forces ( kg_env, error )

    TYPE(kg_environment_type), POINTER       :: kg_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'kg_calculate_forces', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, iatom, ikind, &
                                                natoms, nkind, output_unit, &
                                                stat
    INTEGER, ALLOCATABLE, DIMENSION(:), SAVE :: natom_of_kind
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: calculate_forces, failure
    REAL(KIND=dp)                            :: dvol, total_rho0_rspace
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(kg_energy_type), POINTER            :: energy
    TYPE(kg_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(kg_rho_type), POINTER               :: rho
    TYPE(kg_rspw_type), POINTER              :: rspw
    TYPE(particle_type), POINTER             :: part( : )
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: my_rho_g, my_rho_r, my_tau, &
                                                my_vxc_r, my_vxc_tau
    TYPE(pw_p_type), POINTER                 :: rho_core, rho_g, rho_r, &
                                                v_gspace, v_rspace, vxc_r
    TYPE(pw_poisson_type), POINTER           :: poisson_env
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(section_vals_type), POINTER         :: xc_section

!  TYPE ( debug_variables_type ), INTENT ( OUT ), OPTIONAL :: debug
!  CHARACTER(LEN=8) :: force_type
!------------------------------------------------------------------------------

  CALL timeset (routineN,handle )
  failure=.FALSE.
  para_env=>kg_env%para_env
  logger=>cp_error_get_logger(error)
  output_unit=cp_logger_get_default_io_unit(logger)
  NULLIFY(force,energy,part,atomic_kind_set,atomic_kind,cell,dft_control,&
          rho,rspw,auxbas_pw_pool)
  NULLIFY(rho_r,rho_g,rho_core,v_gspace,v_rspace,vxc_r)
  NULLIFY(my_vxc_r,my_vxc_tau,my_rho_r,my_rho_g,my_tau)
  xc_section => section_vals_get_subs_vals(kg_env%input,"DFT%XC", error=error)
  CALL get_kg_env( kg_env=kg_env, atomic_kind_set=atomic_kind_set,&
                    dft_control=dft_control, energy=energy, cell=cell,&
                    force=force,particle_set=part,rho=rho,rspw=rspw,error=error)
  CPPrecondition(ASSOCIATED(rho),cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(rspw),cp_failure_level,routineP,error,failure)
  IF(.NOT.failure) THEN
    CALL kg_rho_get(kg_rho=rho, rho_r=rho_r, rho_g=rho_g, rho_core=rho_core,error=error)
    CALL kg_rspw_get(kg_rspw=rspw,auxbas_pw_pool=auxbas_pw_pool,&
         poisson_env=poisson_env,error=error)
  END IF

  ALLOCATE(v_gspace,stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  ALLOCATE(v_rspace,stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  ALLOCATE(vxc_r,stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  IF (.NOT.failure) THEN
    CALL pw_pool_create_pw(auxbas_pw_pool, v_gspace%pw, &
                            use_data = COMPLEXDATA1D,&
                            in_space = RECIPROCALSPACE, error=error)
    CALL pw_zero(v_gspace%pw,error=error)
    CALL pw_pool_create_pw(auxbas_pw_pool, v_rspace%pw,&
                            use_data = REALDATA3D,&
                            in_space = REALSPACE, error=error)
    CALL pw_zero(v_rspace%pw,error=error)
  END IF

  calculate_forces=.TRUE.
  natoms = SIZE ( part )
  nkind = SIZE( atomic_kind_set )
  stat = 0
  dvol=v_rspace%pw%pw_grid%dvol

  IF ( .NOT. ALLOCATED ( natom_of_kind ) ) THEN
     ALLOCATE (natom_of_kind(nkind),STAT=stat)
     IF (stat /= 0) CALL stop_memory(routineP,"natom_of_kind",nkind)
     CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                              natom_of_kind=natom_of_kind)
  END IF

  CALL build_kg_neighbor_lists(kg_env,para_env,error)

! initialize forces and energies

  CALL init_kg_force(force)
  CALL init_kg_energy(energy)

  DO i = 1, natoms
     part ( i ) % f ( 1 ) = 0.0_dp
     part ( i ) % f ( 2 ) = 0.0_dp
     part ( i ) % f ( 3 ) = 0.0_dp
  END DO

! self energy
  CALL calculate_eself ( kg_env, energy % core_self ,error=error)

! Calculate the rspace energy of the core charge distribution
! and its force contributions

  CALL calculate_ecore_rspace(kg_env,para_env,calculate_forces,error)

! calculate intramolecular bonded_correction energy

  CALL calculate_ebond_corr(kg_env,energy%bond_corr,calculate_forces,error=error)

! Calculate the local pseudopotential energy

  CALL calculate_rho0_ppl(kg_env,calculate_forces,error)

! Calculate the density generated by the core charges

  CALL calculate_epc_density(rho_core, rho%total_rho_core_rspace, kg_env,error=error)

! calculate the frozen electronic pseudo density

  CALL calculate_density (rho_r, rho_g, total_rho0_rspace, kg_env,error=error)

! from the frozen rho, calculate the thomas-fermi and xc potentials
  ALLOCATE(my_rho_r(1))
  my_rho_r(1)%pw => rho_r%pw
  ALLOCATE(my_rho_g(1))
  my_rho_g(1)%pw => rho_g%pw
  CALL xc_vxc_pw_create1(my_vxc_r,my_vxc_tau,my_rho_r,my_rho_g,my_tau,&
       energy%exc,xc_section,cell,auxbas_pw_pool,error=error)
! ALLOCATE(vxc_r)
  vxc_r%pw => my_vxc_r(1)%pw
  NULLIFY(my_rho_r(1)%pw, my_vxc_r(1)%pw,my_rho_g(1)%pw)
  DEALLOCATE(my_rho_r, my_rho_g, my_vxc_r)
  NULLIFY(my_rho_r, my_rho_g, my_vxc_r)

  vxc_r%pw%cr3d(:,:,:) =dvol*vxc_r%pw%cr3d(:,:,:)

! calculate the total density
  rho%total_rho_rspace = total_rho0_rspace+ rho%total_rho_core_rspace

! add the core density to the frozen electronic density

  CALL pw_axpy(rho_core%pw,rho_g%pw,error=error)
  rho%total_rho_gspace = pw_integrate_function(rho_g%pw,error=error)

! calculate electrostatic potential
  CALL pw_poisson_solve(poisson_env,rho_g%pw, energy % hartree,v_gspace%pw,error=error)

!  calclation of the forces on the ions
  CALL pw_transfer(v_gspace%pw,v_rspace%pw,error=error)
  CALL pw_pool_give_back_pw(auxbas_pw_pool,v_gspace%pw,error=error)
  DEALLOCATE(v_gspace,stat=stat)
  CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)

  v_rspace%pw%cr3d(:,:,:) = dvol*v_rspace%pw%cr3d(:,:,:)

  CALL calculate_epc_rspace_forces(v_rspace, kg_env,error=error)

  v_rspace%pw%cr3d(:,:,:) = v_rspace%pw%cr3d(:,:,:)+ vxc_r%pw%cr3d(:,:,:)

  CALL calculate_v_rspace_forces(v_rspace,kg_env,error=error)

! deallocate work storage
  CALL pw_pool_give_back_pw(auxbas_pw_pool,v_rspace%pw,error=error)
  DEALLOCATE(v_rspace,stat=stat)
  CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
  CALL pw_pool_give_back_pw(auxbas_pw_pool,vxc_r%pw,error=error)
  DEALLOCATE(vxc_r,stat=stat)
  CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)

! add up all the potential energies

  energy % total = energy % core_self + energy % core_overlap + energy % exc + &
                 energy % hartree + energy % pseudo + energy % bond_corr

! print
    IF (output_unit>0) THEN
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
    CALL mp_sum(force(ikind)%f_rho,para_env%group)
    CALL mp_sum(force(ikind)%f_rspace_core,para_env%group)
    CALL mp_sum(force(ikind)%f_hartree_core,para_env%group)
    CALL mp_sum(force(ikind)%f_ppl,para_env%group)
    CALL mp_sum(force(ikind)%f_bc,para_env%group)
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
     IF (stat /= 0) CALL stop_memory(routineP,"natom_of_kind",nkind)
  END IF

 CALL timestop(handle)

END SUBROUTINE kg_calculate_forces

! *****************************************************************************
!> \brief Calculates the total potential energy, total force
!>      using the KG GPW method. This implies the calculation of the
!>      GPW KS energy  and the KS matrix via qs_ks_build_kohn_sham_matrix
!>      where the KE energy correction is added to be able to treat the
!>      molecule separatedly in the construction of the density.
!>      In the SCF iteration the optimization of the wavefunctions is
!>      done by the standard diagonalization of the diagonal blocks
!>      of the MO matrix. Each block is handled independently
!>      This is possible because the wavefunctions located
!>      on different molecules are orthogonal by construction.
!> \param kg_env Kg environment which contains the qs envirement as subenvironment
!> \param globenv global environment
!> \note
!>      Even if it is a KG method energies and forces are of the type qs_energy and qs_force
!>      because the calculation goes through the QS routines for GPW.
!>      Therefore the kg_energy and kg_force remain not allocated.
!>      An additional term has been added in either qs_energy and qs_force, which
!>      is called    kg_gpw_ekin_mol and contains the contributions from the
!>      correction calculated as by the kinetic energy functional applied to
!>      the molecular densities.
!> \par History
!>      created 11-02-05
!> \author MI
! *****************************************************************************
  SUBROUTINE kg_qs_calculate_forces( kg_env, globenv, error )

    TYPE(kg_environment_type), POINTER       :: kg_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'kg_qs_calculate_forces', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, iatom, ikind, &
                                                ispin, istat, iw, natom, &
                                                nspins, output_unit
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind, kind_of
    LOGICAL                                  :: failure
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dbcsr_obj_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s_b, matrix_w_b
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(kg_fm_p_type), POINTER              :: kg_fm
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s, matrix_w
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(section_vals_type), POINTER         :: print_section

    failure = .FALSE.
    NULLIFY(logger)
    logger => cp_error_get_logger(error)

    CALL timeset(routineN,handle)

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(kg_env),cp_failure_level,routineP,error,failure)
    para_env=>kg_env%para_env

    NULLIFY( qs_env , kg_fm )
    CALL get_kg_env(kg_env=kg_env, sub_qs_env=qs_env, kg_fm_set=kg_fm, &
                    scf_control=scf_control, error=error)

    NULLIFY (atomic_kind_set)
    NULLIFY (dft_control)
    NULLIFY (force)
    NULLIFY (ks_env)
    NULLIFY (particle_set)
    NULLIFY (matrix_s)
!!    NULLIFY (scf_control)
    NULLIFY (matrix_w)
    NULLIFY (rho)

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    dft_control=dft_control,&
                    force=force,&
                    particle_set=particle_set,&
 !!                 scf_control=scf_control,&
                    error=error)

    natom = SIZE(particle_set)

    ! zero out the forces
    DO iatom=1,natom
       particle_set(iatom)%f=0.0_dp
    END DO

    IF((dft_control%qs_control%method=="GAPW") .OR. &
        dft_control%qs_control%semi_empirical .OR. &
        dft_control%qs_control%dftb) CALL stop_program(routineP,&
         "KG_GPW is not implemented with GAPW, DFTB or semiempirical methods" )

    ALLOCATE (atom_of_kind(natom),STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE (kind_of(natom),STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             atom_of_kind=atom_of_kind,&
                             kind_of=kind_of)

    CALL zero_qs_force(force)

    ! Do the SCF calculation for KG GPW
    CALL calculate_ecore_overlap(qs_env,para_env,.TRUE.,molecular=.TRUE., error=error)
    CALL kg_qs_energies(kg_env,globenv,error=error)

    CALL get_qs_env(qs_env=qs_env,&
                    ks_env=ks_env, matrix_s=matrix_s_b, rho=rho,error=error)

    NULLIFY(matrix_s)!sm->dbcsr
    CALL allocate_matrix_set( matrix_s, SIZE(matrix_s_b), error )!sm->dbcsr
    DO ispin=1,SIZE(matrix_s)!sm->dbcsr
       CALL sm_from_dbcsr(matrix_s(ispin)%matrix, matrix_s_b(ispin)%matrix, &
            matrix_s_b(1)%matrix%m%dist2d,error)!sm->dbcsr
    ENDDO!sm->dbcsr


    nspins = dft_control%nspins

    CALL allocate_matrix_set(matrix_w,nspins,error=error)

    DO ispin=1,nspins
      CALL replicate_matrix_structure(matrix_s(1)%matrix,&
                         matrix_w(ispin)%matrix,"W MATRIX",error=error)
      IF (qs_env%dft_control%restricted .AND. ispin>1) THEN
         ! not very elegant, indeed ...
         CALL set_matrix(matrix_w(ispin)%matrix,0.0_dp)
      ELSE
         IF (scf_control%use_ot) THEN
            CALL stop_program(routineP,"KG_GPW is not implemented with  OT method")
         ELSE
         ! Here I should construct the W matrix starting from the molecukar mos
           CALL calculate_w_matrix_per_molecule(kg_fm%kg_fm_mol_set,&
                                  matrix_w(ispin)%matrix,ispin,error=error)
         END IF
      END IF
      IF (BTEST(cp_print_key_should_output(logger%iter_info,&
           qs_env%input,"DFT%PRINT%AO_MATRICES/W_MATRIX",error=error),cp_p_file)) THEN
         iw = cp_print_key_unit_nr(logger,qs_env%input,"DFT%PRINT%AO_MATRICES/W_MATRIX",&
              extension=".Log",error=error)
         CALL write_sparse_matrix(matrix_w(ispin)%matrix,4,6,qs_env,para_env,output_unit=iw,error=error)
         CALL cp_print_key_finished_output(iw,logger,qs_env%input,&
              "DFT%PRINT%AO_MATRICES/W_MATRIX", error=error)
      END IF
    END DO

    NULLIFY(matrix_w_b)!sm->dbcsr
    CALL dbcsr_allocate_matrix_set( matrix_w_b, SIZE(matrix_w), error )!sm->dbcsr
    DO ispin=1,SIZE(matrix_w_b)!sm->dbcsr
       ALLOCATE(matrix_w_b(ispin)%matrix)!sm->dbcsr
       !CALL dbcsr_init(matrix_w_b(ispin)%matrix)
       CALL dbcsr_from_sm(matrix_w_b(ispin)%matrix, matrix_w(ispin)%matrix, error)!sm->dbcsr
    ENDDO!sm->dbcsr


!   *** from an eventual mulliken restraint
    IF (dft_control%qs_control%mulliken_restraint) THEN
        CALL mulliken_restraint(dft_control%qs_control%mulliken_restraint_control, &
                                para_env,matrix_s_b(1)%matrix, rho%rho_ao,w_matrix=matrix_w_b,error=error)
    ENDIF
    CALL set_qs_env(qs_env=qs_env,matrix_w=matrix_w_b,error=error)

    CALL build_core_hamiltonian_matrix(qs_env=qs_env,&
                                       globenv=globenv,&
                                       calculate_forces=.TRUE.,error=error)

! ** here qs_env%rho%rho_r and qs_env%rho%rho_g should be up to date
! *** compute grid-based forces ***
     CALL qs_ks_update_qs_env(ks_env=ks_env,qs_env=qs_env, kg_env=kg_env,&
                              error=error,&
                              calculate_forces=.TRUE.)

!   *** distribute forces ***
    DO ikind=1,SIZE(force)
      CALL mp_sum(force(ikind)%overlap,para_env%group)
      CALL mp_sum(force(ikind)%kinetic,para_env%group)
      CALL mp_sum(force(ikind)%gth_ppl,para_env%group)
      CALL mp_sum(force(ikind)%gth_ppnl,para_env%group)
      CALL mp_sum(force(ikind)%core_overlap,para_env%group)
      CALL mp_sum(force(ikind)%rho_core,para_env%group)
      CALL mp_sum(force(ikind)%rho_elec,para_env%group)
      CALL mp_sum(force(ikind)%fock_4c,para_env%group)
      CALL mp_sum(force(ikind)%kg_gpw_ekin_mol,para_env%group)
      force(ikind)%total(:,:) = force(ikind)%total(:,:) +&
                                force(ikind)%core_overlap(:,:) +&
                                force(ikind)%gth_ppl(:,:) +&
                                force(ikind)%gth_ppnl(:,:) +&
                                force(ikind)%kinetic(:,:) +&
                                force(ikind)%overlap(:,:) +&
                                force(ikind)%rho_core(:,:) +&
                                force(ikind)%rho_elec(:,:) +&
                                force(ikind)%fock_4c(:,:) +&
                                force(ikind)%kg_gpw_ekin_mol(:,:)
    END DO

    DO iatom=1,natom
      ikind = kind_of(iatom)
      i = atom_of_kind(iatom)
      ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      ! the force is - dE/dR, what is called force is actually the gradient
      ! Things should have the right name
      ! The minus sign below is a hack
      ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      force(ikind)%other(1:3,i)=-particle_set(iatom)%f(1:3)
      force(ikind)%total(1:3,i)=force(ikind)%total(1:3,i)+force(ikind)%other(1:3,i)
      particle_set(iatom)%f = -force(ikind)%total(1:3,i)
    END DO
    output_unit = cp_print_key_unit_nr(logger,kg_env%input,"DFT%PRINT%DERIVATIVES",&
         extension=".Log",error=error)
    print_section => section_vals_get_subs_vals(qs_env%input,"DFT%PRINT%DERIVATIVES",error=error)
    CALL write_forces(force,atomic_kind_set,3,output_unit=output_unit,&
         print_section=print_section,error=error)
    CALL cp_print_key_finished_output(output_unit,logger,kg_env%input,&
         "DFT%PRINT%DERIVATIVES",error=error)

    CALL deallocate_matrix_set(matrix_w,error=error)
    CALL dbcsr_deallocate_matrix_set( matrix_w_b, error=error )!sm->dbcsr
    CALL set_qs_env(qs_env=qs_env,matrix_w=matrix_w_b,error=error)

    DEALLOCATE (atom_of_kind,STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)

    DEALLOCATE (kind_of,STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)

    !DO ispin=1,SIZE(matrix_s)!sm->dbcsr
    !   CALL dbcsr_from_sm(matrix_s_b(ispin)%matrix, matrix_s(ispin)%matrix, error)!sm->dbcsr
    !ENDDO!sm->dbcsr
    CALL deallocate_matrix_set( matrix_s, error )!sm->dbcsr

    CALL timestop(handle)

  END SUBROUTINE kg_qs_calculate_forces

END MODULE kg_force
