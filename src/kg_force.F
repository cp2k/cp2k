!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C ) 2002 - 2003  CP2K developers group                         !
!-----------------------------------------------------------------------------!
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
  USE coefficient_types,               ONLY: coeff_allocate,&
                                             coeff_copy,&
                                             coeff_deallocate,&
                                             coeff_sumup,&
                                             coeff_transform_space,&
                                             coeff_type,&
                                             coeff_zero
  USE dft_types,                       ONLY: dft_control_type
  USE external_potential_types,        ONLY: epc_potential_type,&
                                             gth_potential_type,&
                                             get_potential
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
                                             calculate_vp_rspace_forces,&
                                             calculate_epc_density_single,&
                                             calculate_density_single,&
                                             calculate_vp_rspace_forces_s
  USE kg_environment_types,            ONLY: get_kg_env,&
                                             kg_environment_type
  USE kg_energy_types,                 ONLY: init_kg_energy,&
                                             kg_energy_type
  USE kg_force_types,                  ONLY: kg_force_type, &     
                                             init_kg_force
!  USE kg_intra,                        ONLY: calculate_ebond_corr
  USE kg_kxc,                          ONLY: calculate_kxc_potential
  USE kinds,                           ONLY: dbl,dp_size
  USE kg_neighbor_lists,               ONLY: build_neighbor_lists,&
                                             rebuild_neighbor_lists
  USE kg_ppl,                          ONLY: calculate_rho0_ppl,&
                                             calculate_drho_ppl 
  USE kg_rho_types,                    ONLY: kg_rho_type
  USE mathconstants,                   ONLY: pi,&
                                             zero
  USE message_passing,                 ONLY: mp_sum
  USE orbital_pointers,                ONLY: indco,&
                                             init_orbital_pointers,&
                                             nco,&
                                             ncoset
  USE particle_types,                  ONLY: particle_type
  USE pws,                             ONLY: cleanup_pw_poisson_solver,&
                                             init_pw_poisson_solver
  USE simulation_cell,                 ONLY: cell_type,&
                                             get_cell
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
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

SUBROUTINE kg_force_control ( kg_env, globenv )

  IMPLICIT NONE

! Arguments
  TYPE ( kg_environment_type ), POINTER         :: kg_env
  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv
!  TYPE ( debug_variables_type ), INTENT ( OUT ), OPTIONAL :: debug
! locals
  TYPE ( dft_control_type ), POINTER :: dft_control
  LOGICAL :: kgpol

!------------------------------------------------------------------------------

  kgpol = .FALSE.
  CALL get_kg_env( kg_env=kg_env, dft_control=dft_control)
  kgpol = dft_control % kg_control % polarization
 
!  IF (kgpol) THEN
!  CALL kgpol_calculate_forces (kg_env, globenv, debug)
! ELSE
   CALL kg_calculate_forces (kg_env, globenv )
! END IF
  END SUBROUTINE kg_force_control


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

SUBROUTINE kg_calculate_forces ( kg_env, globenv )

  IMPLICIT NONE

! Arguments
  TYPE ( kg_environment_type ), POINTER :: kg_env
  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv
!  TYPE ( debug_variables_type ), INTENT ( OUT ), OPTIONAL :: debug

! Local parameters
    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE kg_calculate_forces (MODULE kg_force)"
! Locals
  TYPE ( atomic_kind_type ), DIMENSION(:), POINTER:: atomic_kind_set
  TYPE ( atomic_kind_type ), POINTER   :: atomic_kind
  TYPE ( cell_type ), POINTER :: cell
  TYPE ( dft_control_type ), POINTER :: dft_control
  TYPE ( particle_type ), POINTER :: part ( : )
  TYPE ( kg_force_type ), DIMENSION(:), POINTER ::force
  TYPE ( kg_energy_type ), POINTER ::energy
  TYPE ( kg_rho_type ), POINTER :: kg_rho_v
  INTEGER :: id, i, ii, natoms, nnodes, handle, isos, n, nkind
  INTEGER :: iatom, ikind 
  REAL ( dbl ) :: ehartree, exc
  REAL ( dbl ) :: total_rho0_rspace
  LOGICAL :: first_time, calculate_forces
  INTEGER, DIMENSION(:), POINTER   :: atom_list
  INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: natom_of_kind
  CHARACTER(LEN=8) :: force_type

!------------------------------------------------------------------------------

  CALL write_checkpoint_information("entering "//routine,globenv)
  CALL timeset ( 'FORCE','I',' ',handle )
  NULLIFY(force,energy,part,atomic_kind_set,atomic_kind,cell,dft_control,&
          kg_rho_v)

  CALL get_kg_env( kg_env=kg_env, atomic_kind_set=atomic_kind_set,&
                    dft_control=dft_control, energy=energy, cell=cell,&
                    force=force,particle_set=part,rho=kg_rho_v)

  calculate_forces=.TRUE. 
  natoms = SIZE ( part )
  nkind = SIZE( atomic_kind_set )
  isos = 0

  IF ( .NOT. ALLOCATED ( natom_of_kind ) ) THEN
     ALLOCATE (natom_of_kind(nkind),STAT=isos)
     IF (isos /= 0) CALL stop_memory(routine,"natom_of_kind",nkind)
     CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                              natom_of_kind=natom_of_kind)
  END IF

  IF (.not.ASSOCIATED(kg_env%sac_ppl)) THEN
     CALL build_neighbor_lists(kg_env,globenv)
  ELSE
     CALL rebuild_neighbor_lists(kg_env,globenv)
  END IF

! allocate density/pot arrays
 
!  CALL allocate_kg_densities (kg_rho_v, globenv )

! initialize forces and energies                    
 
  CALL init_kg_force(force)
  CALL init_kg_energy(energy)

  DO i = 1, natoms
     part ( i ) % f ( 1 ) = 0.0_dbl
     part ( i ) % f ( 2 ) = 0.0_dbl
     part ( i ) % f ( 3 ) = 0.0_dbl
  END DO

! self energy 
  CALL calculate_eself ( kg_env, energy % core_self )

! Calculate the rspace energy of the core charge distribution
! and its force contributions

  CALL calculate_ecore_rspace(kg_env,globenv,energy%core_overlap,force)

! calculate intramolecular bonded_correction energy  

!  CALL calculate_ebond_corr(kg_env,globenv,energy%e_bond_corr,kg_forces)

! Calculate the local pseudopotential energy         

  CALL calculate_rho0_ppl(kg_env,globenv,energy%pseudo, &
                         calculate_forces,force)

! Calculate the density generated by the core charges

  CALL calculate_epc_density(kg_rho_v%rho_core, kg_rho_v%cube_info(1),&
                         kg_rho_v%l_info,&
                         kg_rho_v%total_rho_core_rspace, kg_env, globenv)

! calculate the frozen electronic pseudo density

  CALL calculate_density ( kg_rho_v%rho_elec_rspace, kg_rho_v%rho_gspace, &
            kg_rho_v%mgrid_rspace, kg_rho_v%mgrid_gspace,&
            kg_rho_v%gridlevel_info, kg_rho_v%cube_info,&
            kg_rho_v%l_info, total_rho0_rspace, kg_env, globenv )


  CALL coeff_zero(kg_rho_v%v_rspace)
! from the frozen rho, calculate the thomas-fermi and xc potentials
  CALL calculate_kxc_potential (kg_rho_v% pw_pool, & 
                                kg_rho_v % rho_elec_rspace, kg_rho_v % rho_gspace, &
                                kg_rho_v%v_xc_rspace, energy % exc, dft_control)

! calculate the total density
  kg_rho_v%total_rho_rspace = total_rho0_rspace+&
             kg_rho_v%total_rho_core_rspace

! add the core density to the frozen electronic density

  CALL coeff_sumup(kg_rho_v%rho_core,kg_rho_v%rho_gspace)

  kg_rho_v%total_rho_gspace = calculate_total_density(kg_rho_v%rho_gspace)

! calculate electrostatic potential
  CALL init_pw_poisson_solver(cell%green)
  CALL calculate_hartree(kg_rho_v%rho_gspace, energy % hartree,&
           kg_rho_v%v_gspace)
  CALL cleanup_pw_poisson_solver(cell%green)

!  calclation of the forces on the ions

  CALL coeff_transform_space(kg_rho_v%v_gspace,kg_rho_v%v_hartree_rspace)

  kg_rho_v%v_hartree_rspace%pw%cr3d(:,:,:) =&
               kg_rho_v%v_hartree_rspace%pw%pw_grid%dvol*&
               kg_rho_v%v_hartree_rspace%pw%cr3d(:,:,:)

  CALL calculate_epc_rspace_forces(kg_rho_v % v_hartree_rspace,&
                   kg_rho_v%cube_info(1),&
                   kg_rho_v%l_info, kg_env, force)

! if we are doing debugging calculate separately the contribution
! to the forces due by tfxc and hartree 
! else, sumup all the realspace potentials and get the total
! contribution to the forces

!  IF ( PRESENT ( debug )  ) THEN
!    force_type = 'hartree'
!    CALL calculate_v_rspace_forces(kg_rho_v%v_hartree_rspace,kg_rho_v%mgrid_rspace,&
!               kg_rho_v%mgrid_gspace, kg_rho_v%auxbas_grid,&
!               kg_rho_v%gridlevel_info,&
!               kg_rho_v%cube_info, kg_rho_v%l_info,&
!               kg_env, kg_forces, globenv, force_type)
!    force_type = 'xc'
!    CALL calculate_v_rspace_forces(kg_rho_v%v_xc_rspace,kg_rho_v%mgrid_rspace,&
!               kg_rho_v%mgrid_gspace, kg_rho_v%auxbas_grid,&
!               kg_rho_v%gridlevel_info,&
!               kg_rho_v%cube_info, kg_rho_v%l_info,&
!               kg_env, kg_forces, globenv, force_type)
!  ELSE 
    kg_rho_v%v_rspace%pw%cr3d(:,:,:) =&
               kg_rho_v%v_xc_rspace%pw%cr3d(:,:,:) +&
               kg_rho_v%v_hartree_rspace%pw%cr3d(:,:,:)

    CALL calculate_v_rspace_forces(kg_rho_v%v_rspace,kg_rho_v%mgrid_rspace,&
               kg_rho_v%mgrid_gspace, kg_rho_v%auxbas_grid,&
               kg_rho_v%gridlevel_info,&
               kg_rho_v%cube_info, kg_rho_v%l_info,&
               kg_env, force, globenv)
!  END IF

! add up all the potential energies

  energy % total = energy % core_self + energy % core_overlap + energy % exc + &
                 energy % hartree + energy % pseudo 
     
! add up all the forces

  DO ikind = 1, nkind         
    CALL mp_sum(force(ikind)%f_rho,globenv%group)
    CALL mp_sum(force(ikind)%f_rspace_core,globenv%group)
    CALL mp_sum(force(ikind)%f_hartree_core,globenv%group)
    CALL mp_sum(force(ikind)%f_gth_ppl,globenv%group)
    force(ikind)%f_total(:,:) = force(ikind)%f_rho(:,:) +&
                                force(ikind)%f_rspace_core(:,:) +&
                                force(ikind)%f_hartree_core(:,:) +&
                                force(ikind)%f_gth_ppl(:,:) 
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
! if we are doing debugging, check if variables are present and assign
!      IF ( PRESENT ( debug )  ) THEN
!        debug % f_hartree (:,i) = force(ikind)%f_hartree(:,iatom) + &
!                              force(ikind)%f_hartree_core(:,iatom)
!       debug % f_rspace_core (:,i) = force(ikind)%f_rspace_core(:,iatom)
!       debug % f_bc (:,i) = force(ikind)%f_bc(:,iatom)
!       debug % f_xc (:,i) = force(ikind)%f_xc(:,iatom)
!       debug % f_ppl (:,i) = force(ikind)%f_gth_ppl(:,iatom)
!     END IF
    END DO
  END DO

!  IF ( PRESENT ( debug )  ) THEN
!    debug % pot_hartree = thermo % e_hartree
!    debug % pot_rspace_core = thermo % e_rspace_core
!    debug % pot_bond_corr = thermo % e_bond_corr
!    debug % pot_xc = thermo % e_xc
!    debug % pot_ppl = thermo % e_ppl
!  END IF

! deallocating all local variables
  IF ( ALLOCATED ( natom_of_kind ) ) THEN
     DEALLOCATE (natom_of_kind ,STAT=isos)
     IF (isos /= 0) CALL stop_memory(routine,"natom_of_kind",nkind)
  END IF

! CALL deallocate_kg_densities ( kg_rho_v, globenv )

 CALL timestop ( zero, handle )
 CALL write_checkpoint_information("leaving "//routine,globenv)

END SUBROUTINE kg_calculate_forces

!******************************************************************************

END MODULE kg_force  

!******************************************************************************
