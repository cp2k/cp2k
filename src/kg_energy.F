!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C ) 2002 - 2003  CP2K developers group                         !
!-----------------------------------------------------------------------------!
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
                                             coeff_zero
  USE cp_control_types,                ONLY: dft_control_type
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
!  USE kg_intra,                        ONLY: calculate_ebond_corr

  USE kg_environment_types,            ONLY: get_kg_env,&
                                             kg_environment_type
  USE kg_kxc,                          ONLY: calculate_kxc_derivatives
  USE kg_neighbor_lists,               ONLY: build_kg_neighbor_lists
  USE kg_ppl,                          ONLY: calculate_drho_ppl,&
                                             calculate_rho0_ppl
  USE kg_rho_types,                    ONLY: kg_rho_type
  USE kinds,                           ONLY: dbl
  USE machine,                         ONLY: m_flush,&
                                             m_walltime
  USE mathconstants,                   ONLY: zero
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_max
  USE pws,                             ONLY: cleanup_pw_poisson_solver,&
                                             init_pw_poisson_solver
  USE scf_control_types,               ONLY: scf_control_type
  USE simulation_cell,                 ONLY: cell_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xc,                              ONLY: xc_calculate_pw1

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: kg_energies

!******************************************************************************

CONTAINS

!******************************************************************************
!!****** kg_energy/kg_energies [1.1] *
!!
!!   NAME
!!     kg_energies
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
SUBROUTINE kg_energies ( kg_env, globenv )


    TYPE(kg_environment_type), POINTER       :: kg_env
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "SUBROUTINE kg_energies (MODULE kg_energy)"

    INTEGER                                  :: handle, handle2, i, isos, &
                                                iter, output_unit
    LOGICAL                                  :: ionode, should_stop
    REAL(dbl)                                :: maxerr, maxgrad, t1, t2
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(diis_info_type), POINTER            :: diis_info
    TYPE(dyn_coeff_set_type), POINTER        :: dyn_coeff_set
    TYPE(scf_control_type), POINTER          :: scf_control

!------------------------------------------------------------------------------

  CALL write_checkpoint_information("entering "//routine,globenv)
  CALL timeset ( 'kg_energies','I',' ',handle )

  NULLIFY(atomic_kind_set,dyn_coeff_set,diis_info,scf_control)
  ionode=globenv%ionode
  output_unit=globenv%scr
  isos = 0
  iter = 0

  t1 = m_walltime()

  CALL get_kg_env(kg_env=kg_env,&
                  atomic_kind_set=atomic_kind_set,&
                  dyn_coeff_set=dyn_coeff_set,&
                  scf_control=scf_control)           
   
  CALL build_kg_neighbor_lists(kg_env,globenv)

  CALL diis_info_create(diis_info,scf_control%max_diis,dyn_coeff_set)

  CALL init_energy_calculation(kg_env,globenv)

  CALL calc_c_gradients(kg_env,maxgrad,globenv)

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
    CALL calc_c_gradients(kg_env,maxgrad,globenv)

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
!   *** Write restart file ***
  IF (scf_control%write_restart_each>=0) THEN
    CALL write_coeff_restart(dyn_coeff_set,atomic_kind_set,globenv)  
  END IF
  
  CALL diis_info_destroy(diis_info)

  CALL timestop(0.0_dbl,handle)

  CALL write_checkpoint_information("leaving "//routine,globenv)

END SUBROUTINE kg_energies

!------------------------------------------------------------------------------

SUBROUTINE init_energy_calculation ( kg_env, globenv )


    TYPE(kg_environment_type), POINTER       :: kg_env
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "SUBROUTINE init_energy_calculation (MODULE kg_energy)"

    INTEGER                                  :: auxbas_grid, handle, i, isos, &
                                                output_unit
    LOGICAL                                  :: calculate_forces, ionode
    REAL(dbl)                                :: dvol, total_rho0_rspace
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(dyn_coeff_set_type), POINTER        :: dyn_coeff_set
    TYPE(kg_energy_type), POINTER            :: energy
    TYPE(kg_rho_type), POINTER               :: kg_rho_v
    TYPE(scf_control_type), POINTER          :: scf_control

!------------------------------------------------------------------------------

  CALL write_checkpoint_information("entering "//routine,globenv)
  CALL timeset ( 'init_energy_calculation','I',' ',handle )
  ionode=globenv%ionode
  output_unit=globenv%scr
  NULLIFY(energy,atomic_kind_set,dft_control,&
          kg_rho_v,scf_control,dyn_coeff_set)

  CALL get_kg_env( kg_env=kg_env, atomic_kind_set=atomic_kind_set,&
                    dft_control=dft_control, energy=energy,&
                    rho=kg_rho_v,dyn_coeff_set=dyn_coeff_set,&
                    scf_control=scf_control)

  calculate_forces=.FALSE. 
  isos = 0
  auxbas_grid=kg_rho_v%auxbas_grid
  dvol=kg_rho_v%pw_grid(auxbas_grid)%dvol

  CALL build_kg_neighbor_lists(kg_env,globenv)

! initialize energies                    
 
  CALL init_kg_energy(energy)

! calculate the constant parts of the energy 
! self energy 
  CALL calculate_eself ( kg_env, energy % core_self )

! Calculate the rspace energy of the core charge distribution
! and its force contributions

  CALL calculate_ecore_rspace(kg_env,globenv,energy%core_overlap)

! calculate intramolecular bonded_correction energy  

!  CALL calculate_ebond_corr(kg_env,globenv,energy%e_bond_corr,kg_forces)

! Calculate the local pseudopotential energy         

  CALL calculate_rho0_ppl(kg_env,globenv,calculate_forces)

! Calculate the density generated by the core charges

  CALL calculate_epc_density(kg_rho_v%rho_core, kg_rho_v%cube_info(1),&
                         kg_rho_v%l_info,&
                         kg_rho_v%total_rho_core_rspace, kg_env)

! calculate the frozen electronic pseudo density

  CALL calculate_density ( kg_rho_v%rho_elec_rspace, kg_rho_v%rho_gspace, &
            kg_rho_v%mgrid_rspace, kg_rho_v%mgrid_gspace,&
            kg_rho_v%gridlevel_info, kg_rho_v%cube_info,&
            kg_rho_v%l_info, total_rho0_rspace, kg_env)

! calculate the total density
  kg_rho_v%total_rho_rspace = total_rho0_rspace+&
             kg_rho_v%total_rho_core_rspace

! from the frozen rho, calculate the thomas-fermi and xc potentials
  CALL xc_calculate_pw1 (kg_rho_v%rho_elec_rspace%pw, kg_rho_v%rho_gspace%pw, &
                        kg_rho_v%v_xc_rspace%pw, energy%exc, dft_control,&
                        kg_rho_v%pw_pool)

  kg_rho_v%v_xc_rspace%pw%cr3d(:,:,:) =dvol*kg_rho_v%v_xc_rspace%pw%cr3d(:,:,:)
! add the frozen electronic density to the core density 
  CALL coeff_sumup(kg_rho_v%rho_gspace,kg_rho_v%rho_core)

!  read coef from restart file
  IF (scf_control%density_guess == "RESTART") THEN
    CALL read_coeff_restart(dyn_coeff_set,globenv)
  END IF    
  CALL timestop(0.0_dbl,handle)
  CALL write_checkpoint_information("leaving "//routine,globenv)

END SUBROUTINE init_energy_calculation 
!------------------------------------------------------------------------------

SUBROUTINE calc_c_gradients ( kg_env, maxgrad, globenv )


    TYPE(kg_environment_type), POINTER       :: kg_env
    REAL(dbl), intent(inout)                 :: maxgrad
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "SUBROUTINE calc_c_gradients (MODULE kg_energy)"

    INTEGER                                  :: auxbas_grid, handle, i, &
                                                ikind, nkind, output_unit
    LOGICAL                                  :: calculate_forces, ionode
    REAL(dbl)                                :: dvol, e1_xc, e2_xc, &
                                                total_rhop_rspace
    REAL(dbl), DIMENSION(:, :), POINTER      :: forces
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(dyn_coeff_set_type), POINTER        :: dyn_coeff_set
    TYPE(dyn_coeff_type), POINTER            :: local_coeffs
    TYPE(kg_energy_type), POINTER            :: energy
    TYPE(kg_rho_type), POINTER               :: kg_rho_v

!------------------------------------------------------------------------------

  CALL write_checkpoint_information("entering "//routine,globenv)
  CALL timeset ( 'calc_c_gradients','I',' ',handle )
  ionode=globenv%ionode
  output_unit=globenv%scr
  calculate_forces=.FALSE.

  NULLIFY(energy,atomic_kind_set,cell,dft_control,&
          dyn_coeff_set,kg_rho_v,local_coeffs,forces)

  CALL get_kg_env( kg_env=kg_env, atomic_kind_set=atomic_kind_set,&
                    dft_control=dft_control, energy=energy, cell=cell,&
                    rho=kg_rho_v,dyn_coeff_set=dyn_coeff_set)
        
  auxbas_grid=kg_rho_v%auxbas_grid
  dvol=kg_rho_v%pw_grid(auxbas_grid)%dvol
  nkind=SIZE(dyn_coeff_set%coeffs_of_kind)
! reinitialize forces
  DO ikind= 1, nkind
    local_coeffs => dyn_coeff_set%coeffs_of_kind(ikind)%coeffs
    IF(ASSOCIATED(local_coeffs)) THEN
      CALL get_dyn_coeff(coeffs=local_coeffs,forces=forces)
      forces(:,:)=0.0_dbl
    END IF
  END DO
  CALL coeff_zero(kg_rho_v%rho_gspace)
  CALL coeff_copy(kg_rho_v%rho_core,kg_rho_v%rho_gspace)
  kg_rho_v%total_rho_gspace = calculate_total_density(kg_rho_v%rho_gspace)

  CALL calculate_drho_ppl(kg_env,globenv,calculate_forces)
  CALL calculate_p_density(kg_rho_v%rho_p_rspace, kg_rho_v%rho_p_gspace, &
            kg_rho_v%mgrid_rspace, kg_rho_v%mgrid_gspace,&
            kg_rho_v%gridlevel_info, kg_rho_v%cube_info,&
            kg_rho_v%l_info, total_rhop_rspace, kg_env )

  CALL calculate_kxc_derivatives(kg_rho_v% pw_pool, & 
                                kg_rho_v % rho_elec_rspace, kg_rho_v % rho_gspace, &
                                kg_rho_v % rho_p_rspace, kg_rho_v%v_xc_rspace,&
                                kg_rho_v%v1_xc_rspace, kg_rho_v%v2_xc_rspace,&
                                e1_xc, e2_xc, dft_control,calculate_forces)



! add the polarization density to get the total density
  CALL coeff_sumup(kg_rho_v%rho_p_gspace,kg_rho_v%rho_gspace)

! calculate electrostatic potential
  CALL init_pw_poisson_solver(cell%green)
  CALL calculate_hartree(kg_rho_v%rho_gspace, energy % hartree,&
           kg_rho_v%v_gspace)
  CALL cleanup_pw_poisson_solver(cell%green)

  CALL coeff_transform_space(kg_rho_v%v_gspace,kg_rho_v%v_hartree_rspace)

  kg_rho_v%v_hartree_rspace%pw%cr3d(:,:,:) =&
               dvol*kg_rho_v%v_hartree_rspace%pw%cr3d(:,:,:)

  kg_rho_v%vp_rspace%pw%cr3d(:,:,:) =-1.0_dbl* &
              (kg_rho_v%v_xc_rspace%pw%cr3d(:,:,:) +&
               kg_rho_v%v1_xc_rspace%pw%cr3d(:,:,:)+&
               kg_rho_v%v_hartree_rspace%pw%cr3d(:,:,:))

  CALL calculate_vp_rspace_forces(kg_rho_v%vp_rspace,kg_rho_v%mgrid_rspace,&
               kg_rho_v%mgrid_gspace, kg_rho_v%auxbas_grid,&
               kg_rho_v%gridlevel_info,&
               kg_rho_v%cube_info, kg_rho_v%l_info,&
               kg_env,calculate_forces)

! add up all the potential energies

  energy % total = energy % core_self + energy % core_overlap + energy % exc + &
                 e1_xc + e2_xc+ energy % hartree + energy % pseudo + energy%ppseudo 
     
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
            kg_rho_v%total_rho_core_rspace
       WRITE (UNIT=output_unit,FMT="(T3,A,T60,F20.10)")&
            "Total charge density (r-space):     ",&
            kg_rho_v%total_rho_rspace
       WRITE (UNIT=output_unit,FMT="(T3,A,T60,F20.10)")&
            "Total rho_coefs charge density (r-space):     ",&
            total_rhop_rspace
       WRITE (UNIT=output_unit,FMT="(T3,A,T60,F20.10)")&
            "Maximum coefficent gradient         ",&
            maxgrad                  
       WRITE (UNIT=output_unit,FMT="(/,(T3,A,T55,F25.14))")&
            "Overlap energy of the core charge distribution:",energy%core_overlap,&
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
