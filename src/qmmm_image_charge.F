!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Routines for image charge calculation within QM/MM
!> \par History
!>      12.2011 created
!> \author Dorothea Golze
! *****************************************************************************
MODULE  qmmm_image_charge

  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_generate_filename,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE input_constants,                 ONLY: calc_always,&
                                             calc_once,&
                                             calc_once_done
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_path_length,&
                                             dp,&
                                             int_8
  USE mathconstants,                   ONLY: pi
  USE memory_utilities,                ONLY: reallocate
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_sum
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_methods,                      ONLY: pw_integral_ab,&
                                             pw_scale,&
                                             pw_transfer,&
                                             pw_zero
  USE pw_poisson_methods,              ONLY: pw_poisson_solve
  USE pw_poisson_types,                ONLY: pw_poisson_type
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_p_type,&
                                             pw_release
  USE qmmm_types,                      ONLY: qmmm_env_qm_type
  USE qs_collocate_density,            ONLY: calculate_rho_metal,&
                                             calculate_rho_single_gaussian
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_integrate_potential,          ONLY: integrate_pgf_product_rspace
  USE realspace_grid_types,            ONLY: pw2rs,&
                                             realspace_grid_desc_type,&
                                             realspace_grid_type,&
                                             rs_grid_release,&
                                             rs_grid_retain,&
                                             rs_pw_transfer
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE virial_types,                    ONLY: virial_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qmmm_image_charge'

  PUBLIC :: calculate_image_pot,&
            integrate_potential_devga_rspace,&
            conditional_calc_image_matrix,&
            add_image_pot_to_hartree_pot,&
            print_image_coefficients

!***
CONTAINS
! *****************************************************************************
!> \brief determines coefficients by solving image_matrix*coeff=-pot_const by
!>        Gaussian elimination or in an iterative fashion and calculates
!>        image/metal potential with these coefficients
!> \param v_hartree_rspace Hartree potential in real space
!> \param rho_hartree_gspace Kohn Sham density in reciprocal space
!> \param energy structure where energies are stored
!> \param qmmm_env qmmm environment
!> \param qs_env qs environment
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE calculate_image_pot(v_hartree_rspace,rho_hartree_gspace,energy,&
                                qmmm_env,qs_env,error)

    TYPE(pw_p_type), INTENT(IN)              :: v_hartree_rspace, &
                                                rho_hartree_gspace
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qmmm_env_qm_type), POINTER          :: qmmm_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_image_pot', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure

    CALL timeset(routineN,handle)

    failure=.FALSE.

    IF(qmmm_env%image_charge_pot%coeff_iterative) THEN
       !calculate preconditioner matrix for CG if necessary
       IF(qs_env%calc_image_preconditioner) THEN
        IF(qmmm_env%image_charge_pot%image_restart) THEN
         CALL restart_image_matrix(image_matrix=qs_env%image_matrix,&
                                   qs_env=qs_env,qmmm_env=qmmm_env,&
                                   error=error)
        ELSE
         CALL calculate_image_matrix(image_matrix=qs_env%image_matrix,&
                                     qs_env=qs_env,qmmm_env=qmmm_env,&
                                     error=error)
        ENDIF
       ENDIF
       CALL calc_image_coeff_iterative(v_hartree_rspace=v_hartree_rspace,&
                                    coeff=qs_env%image_coeff,qmmm_env=qmmm_env,&
                                    qs_env=qs_env,error=error)

    ELSE
       CALL calc_image_coeff_gaussalgorithm(v_hartree_rspace=v_hartree_rspace,&
                                    coeff=qs_env%image_coeff,qmmm_env=qmmm_env,&
                                    qs_env=qs_env,error=error)
    ENDIF

    ! calculate the image/metal potential with the optimized coefficients
    CALL calculate_potential_metal(v_metal_rspace=&
         qs_env%ks_qmmm_env%v_metal_rspace,coeff=qs_env%image_coeff,&
         rho_hartree_gspace=rho_hartree_gspace,&
         energy=energy,qs_env=qs_env,error=error)

    CALL timestop(handle)


  END SUBROUTINE calculate_image_pot

! *****************************************************************************
!> \brief determines coefficients by solving the linear set of equations
!>        image_matrix*coeff=-pot_const using a Gaussian elimination scheme
!> \param v_hartree_rspace Hartree potential in real space
!> \param coeff expansion coefficients of the image charge density, i.e.
!>        rho_metal=sum_a c_a*g_a
!> \param qmmm_env qmmm environment
!> \param qs_env qs environment
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE calc_image_coeff_gaussalgorithm(v_hartree_rspace,coeff,qmmm_env,&
                                            qs_env,error)


    TYPE(pw_p_type), INTENT(IN)              :: v_hartree_rspace
    REAL(KIND=dp), DIMENSION(:), POINTER     :: coeff
    TYPE(qmmm_env_qm_type), POINTER          :: qmmm_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'calc_image_coeff_gaussalgorithm', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, info, natom, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: eta, V0
    REAL(KIND=dp), DIMENSION(:), POINTER     :: pot_const

    CALL timeset(routineN,handle)

    NULLIFY(pot_const)
    failure=.FALSE.
    
    !minus sign V0: account for the fact that v_hartree has the opposite sign
    V0=-qmmm_env%image_charge_pot%V0
    eta=qmmm_env%image_charge_pot%eta
    natom=SIZE(qmmm_env%image_charge_pot%image_mm_list)

    ALLOCATE(pot_const(natom),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF(.NOT.ASSOCIATED(coeff)) THEN
       ALLOCATE(coeff(natom),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ENDIF
    coeff=0.0_dp

    CALL integrate_potential_ga_rspace(v_hartree_rspace,qmmm_env,qs_env,&
                                       pot_const,error=error)
    !add integral V0*ga(r)
    pot_const(:)=-pot_const(:)+V0*SQRT((pi/eta)**3)

    !solve linear system of equations T*coeff=-pot_const
    !LU factorization of T  by DGETRF done in calculate_image_matrix
    CALL DGETRS('N',natom,1,qs_env%image_matrix,natom,qs_env%ipiv,&
                pot_const,natom,info)
    CPPrecondition(info==0,cp_failure_level,routineP,error,failure)

    coeff = pot_const

    DEALLOCATE(pot_const,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE calc_image_coeff_gaussalgorithm

! *****************************************************************************
!> \brief determines image coefficients iteratively
!> \param v_hartree_rspace Hartree potential in real space
!> \param coeff expansion coefficients of the image charge density, i.e.
!>        rho_metal=sum_a c_a*g_a
!> \param qmmm_env qmmm environment
!> \param qs_env qs environment
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE calc_image_coeff_iterative(v_hartree_rspace,coeff,qmmm_env,&
                                      qs_env,error)


    TYPE(pw_p_type), INTENT(IN)              :: v_hartree_rspace
    REAL(KIND=dp), DIMENSION(:), POINTER     :: coeff
    TYPE(qmmm_env_qm_type), POINTER          :: qmmm_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_image_coeff_iterative', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iter_steps, natom, &
                                                output_unit, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: alpha, eta, rsnew, rsold, V0
    REAL(KIND=dp), DIMENSION(:), POINTER     :: Ad, d, pot_const, r, &
                                                vmetal_const, z
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(pw_p_type)                          :: auxpot_Ad_rspace, &
                                                v_metal_rspace_guess
    TYPE(section_vals_type), POINTER         :: input

    CALL timeset(routineN,handle)
   
    NULLIFY(pot_const,vmetal_const,logger,input)
    failure=.FALSE.
    logger => cp_error_get_logger(error)
   
    !minus sign V0: account for the fact that v_hartree has the opposite sign
    V0=-qmmm_env%image_charge_pot%V0
    eta=qmmm_env%image_charge_pot%eta
    natom=SIZE(qmmm_env%image_charge_pot%image_mm_list)
   
    ALLOCATE(pot_const(natom),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(vmetal_const(natom),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(r(natom),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(d(natom),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(z(natom),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(Ad(natom),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF(.NOT.ASSOCIATED(coeff)) THEN
       ALLOCATE(coeff(natom),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ENDIF
   
    CALL integrate_potential_ga_rspace(v_hartree_rspace,qmmm_env,qs_env,&
                                       pot_const,error=error)
   
    !add integral V0*ga(r)
    pot_const(:)=-pot_const(:)+V0*SQRT((pi/eta)**3)
   
    !initial guess for coeff
    coeff=1.0_dp
    d=0.0_dp
    z=0.0_dp
    r=0.0_dp
    rsold=0.0_dp
    rsnew=0.0_dp
    iter_steps=0
   
   !calculate first guess of image/metal potential
    CALL calculate_potential_metal(v_metal_rspace=v_metal_rspace_guess,&
         coeff=coeff,qs_env=qs_env,error=error)
    CALL integrate_potential_ga_rspace(potential=v_metal_rspace_guess,&
         qmmm_env=qmmm_env,qs_env=qs_env,int_res=vmetal_const,error=error)
   
    ! modify coefficients iteratively
    r=pot_const-vmetal_const
    z=MATMUL(qs_env%image_matrix,r)
    d=z
    rsold=DOT_PRODUCT(r,z)
   
    DO
      !calculate A*d
      Ad=0.0_dp
      CALL calculate_potential_metal(v_metal_rspace=&
           auxpot_Ad_rspace,coeff=d,qs_env=qs_env,error=error)
      CALL integrate_potential_ga_rspace(potential=&
           auxpot_Ad_rspace,qmmm_env=qmmm_env,&
           qs_env=qs_env,int_res=Ad,error=error)
   
      alpha=rsold/DOT_PRODUCT(d,Ad)
      coeff=coeff+alpha*d
   
      r=r-alpha*Ad
      z=MATMUL(qs_env%image_matrix,r)
      rsnew=DOT_PRODUCT(r,z)
      iter_steps=iter_steps+1
      IF(SQRT(rsnew)<1.0E-08) THEN
         CALL pw_release(auxpot_Ad_rspace%pw,error=error)
         EXIT
      END IF
      d=z+rsnew/rsold*d
      rsold=rsnew
      CALL pw_release(auxpot_Ad_rspace%pw,error=error)
    ENDDO
   
    ! print iteration info
    CALL get_qs_env(qs_env=qs_env,&
                      input=input,&
                      error=error)
    output_unit=cp_print_key_unit_nr(logger,input,&
                   "QMMM%PRINT%PROGRAM_RUN_INFO",&
                    extension=".qmmmLog",error=error)
    IF (output_unit>0)  WRITE (UNIT=output_unit,FMT="(T3,A,T74,I7)")&
      "Number of iteration steps for determination of image coefficients:",iter_steps
    CALL cp_print_key_finished_output(output_unit,logger,input,&
                    "QMMM%PRINT%PROGRAM_RUN_INFO", error=error)
   
   
    IF(iter_steps.lt.25) THEN
     qs_env%calc_image_preconditioner=.FALSE.
    ELSE
     qs_env%calc_image_preconditioner=.TRUE.
    ENDIF
   
    CALL pw_release(v_metal_rspace_guess%pw,error=error)
    DEALLOCATE(pot_const,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(vmetal_const,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(r,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(d,z,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(Ad,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
   
    CALL timestop(handle)

  END SUBROUTINE calc_image_coeff_iterative

! ****************************************************************************
!> \brief calculates the integral V(r)*ga(r)
!> \param potential any potential
!> \param qmmm_env qmmm environment
!> \param qs_env qs environment
!> \param int_res result of the integration
!> \param atom_num atom index, needed when calculating image_matrix 
!> \param atom_num_ref index of reference atom, needed when calculating 
!>        image_matrix 
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE integrate_potential_ga_rspace(potential,qmmm_env,qs_env,int_res,&
                                           atom_num,atom_num_ref,error)


    TYPE(pw_p_type), INTENT(IN)              :: potential
    TYPE(qmmm_env_qm_type), POINTER          :: qmmm_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    REAL(KIND=dp), DIMENSION(:), POINTER     :: int_res
    INTEGER, INTENT(IN), OPTIONAL            :: atom_num, atom_num_ref
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'integrate_potential_ga_rspace', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: atom_a, atom_b, atom_ref, &
                                                handle, iatom, j, k, natom, &
                                                npme, stat
    INTEGER, DIMENSION(:), POINTER           :: cores
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: eps_rho_rspace
    REAL(KIND=dp), DIMENSION(3)              :: ra
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: hab
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(realspace_grid_desc_type), POINTER  :: auxbas_rs_desc
    TYPE(realspace_grid_type), POINTER       :: rs_v

    CALL timeset(routineN,handle)

    NULLIFY(cores,hab,cell,auxbas_rs_desc,pw_env,para_env,&
           dft_control,rs_v)
    ALLOCATE (hab(1,1),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL get_qs_env(qs_env=qs_env,pw_env=pw_env,error=error)
    CALL pw_env_get(pw_env=pw_env,auxbas_rs_desc=auxbas_rs_desc, &
                        auxbas_rs_grid=rs_v,error=error)
    CALL rs_grid_retain(rs_v,error=error)
    CALL rs_pw_transfer(rs_v,potential%pw,pw2rs,error=error)

    CALL get_qs_env(qs_env=qs_env,&
                    cell=cell,&
                    dft_control=dft_control,&
                    para_env=para_env,pw_env=pw_env,&
                    error=error)

    eps_rho_rspace = dft_control%qs_control%eps_rho_rspace


    natom=SIZE(qmmm_env%image_charge_pot%image_mm_list)
    k=1
    IF(PRESENT(atom_num)) k=atom_num

    CALL reallocate ( cores, 1, natom-k+1 )
    int_res = 0.0_dp
    npme = 0
    cores = 0

    DO iatom = k, natom
       IF(rs_v%desc%parallel .AND. .NOT. rs_v%desc%distributed) THEN
          ! replicated realspace grid, split the atoms up between procs
           IF (MODULO(iatom,rs_v%desc%group_size) == rs_v % desc % my_pos ) THEN
              npme = npme + 1
              cores (npme) = iatom
           ENDIF
        ELSE
           npme = npme + 1
           cores (npme) = iatom
        ENDIF
    END DO


    DO j=1,npme

     iatom = cores(j)
     atom_a = qmmm_env%image_charge_pot%image_mm_list(iatom)

     IF(PRESENT(atom_num).AND.PRESENT(atom_num_ref)) THEN
       ! shift the function since potential only calculate for ref atom
       atom_b = qmmm_env%image_charge_pot%image_mm_list(k)
       atom_ref = qmmm_env%image_charge_pot%image_mm_list(atom_num_ref)
       ra(:) = pbc(qmmm_env%image_charge_pot%particles_all(atom_a)%r,cell)&
              -pbc(qmmm_env%image_charge_pot%particles_all(atom_b)%r,cell)&
              +pbc(qmmm_env%image_charge_pot%particles_all(atom_ref)%r,cell)
       
     ELSE
       ra(:) = pbc(qmmm_env%image_charge_pot%particles_all(atom_a)%r,cell)
     ENDIF

     hab(1,1) = 0.0_dp

     CALL integrate_pgf_product_rspace(0,qmmm_env%image_charge_pot%eta,0,&
          0,0.0_dp,0,ra,(/0.0_dp,0.0_dp,0.0_dp/),0.0_dp,&
          rs_v,cell,pw_env%cube_info(1),hab,o1=0,o2=0,&
          eps_gvg_rspace=eps_rho_rspace, calculate_forces=.FALSE.,&
          use_subpatch=.TRUE.,subpatch_pattern=0_int_8,error=error)

     int_res(iatom) = hab(1,1)

    END DO

    CALL mp_sum(int_res,para_env%group)

    DEALLOCATE (hab,cores,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL rs_grid_release(rs_v, error=error)

    CALL timestop(handle)

  END SUBROUTINE integrate_potential_ga_rspace

! *****************************************************************************
!> \brief calculates the image forces on the MM atoms
!> \param potential any potential, in this case: Hartree potential
!> \param coeff expansion coefficients of the image charge density, i.e.
!>        rho_metal=sum_a c_a*g_a
!> \param forces structure storing the force contribution of the image charges 
!>        for the metal (MM) atoms
!> \param qmmm_env qmmm environment
!> \param qs_env qs environment
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE integrate_potential_devga_rspace(potential,coeff,forces,qmmm_env,&
                                           qs_env,error)

    TYPE(pw_p_type), INTENT(IN)              :: potential
    REAL(KIND=dp), DIMENSION(:), POINTER     :: coeff
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: forces
    TYPE(qmmm_env_qm_type), POINTER          :: qmmm_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'integrate_potential_devga_rspace', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: atom_a, handle, iatom, j, &
                                                natom, npme, stat
    INTEGER, DIMENSION(:), POINTER           :: cores
    LOGICAL                                  :: failure, use_virial
    REAL(KIND=dp)                            :: eps_rho_rspace
    REAL(KIND=dp), DIMENSION(3)              :: force_a, force_b, ra
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: hab, pab
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(realspace_grid_desc_type), POINTER  :: auxbas_rs_desc
    TYPE(realspace_grid_type), POINTER       :: rs_v
    TYPE(virial_type), POINTER               :: virial

    CALL timeset(routineN,handle)
   
    NULLIFY(cores,hab,pab,cell,auxbas_rs_desc,pw_env,para_env,&
           dft_control,rs_v,virial)
    failure=.FALSE.
    use_virial=.FALSE.
   
    ALLOCATE (hab(1,1),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (pab(1,1),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
   
    CALL get_qs_env(qs_env=qs_env,pw_env=pw_env,error=error)
    CALL pw_env_get(pw_env=pw_env,auxbas_rs_desc=auxbas_rs_desc, &
                        auxbas_rs_grid=rs_v,error=error)
    CALL rs_grid_retain(rs_v,error=error)
    CALL rs_pw_transfer(rs_v,potential%pw,pw2rs,error=error)
   
    CALL get_qs_env(qs_env=qs_env,&
                    cell=cell,&
                    dft_control=dft_control,&
                    para_env=para_env,pw_env=pw_env,&
                    virial=virial,error=error)
   
    use_virial = virial%pv_availability.AND.(.NOT.virial%pv_numer)
   
    IF (use_virial) THEN
       CALL cp_unimplemented_error(fromWhere=routineP, &
            message="Virial not implemented for image charge method", &
            error=error, error_level=cp_failure_level)
    END IF
   
    eps_rho_rspace = dft_control%qs_control%eps_rho_rspace
   
    natom=SIZE(qmmm_env%image_charge_pot%image_mm_list)
   
    IF(.NOT.ASSOCIATED(forces)) THEN
     ALLOCATE(forces(3,natom), stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ENDIF
   
    forces(:,:)=0.0_dp
   
    CALL reallocate ( cores, 1, natom )
    npme = 0
    cores = 0
   
    DO iatom = 1, natom
       IF(rs_v%desc%parallel .AND. .NOT. rs_v%desc%distributed) THEN
          ! replicated realspace grid, split the atoms up between procs
           IF (MODULO(iatom,rs_v%desc%group_size) == rs_v % desc % my_pos ) THEN
              npme = npme + 1
              cores (npme) = iatom
           ENDIF
        ELSE
           npme = npme + 1
           cores (npme) = iatom
        ENDIF
    END DO
   
   
    DO j=1,npme
   
     iatom = cores(j)
     atom_a = qmmm_env%image_charge_pot%image_mm_list(iatom)
     ra(:) = pbc(qmmm_env%image_charge_pot%particles_all(atom_a)%r,cell)
     hab(1,1) = 0.0_dp
     pab(1,1)=1.0_dp
     force_a(:) = 0.0_dp
     force_b(:) = 0.0_dp
   
     CALL integrate_pgf_product_rspace(0,qmmm_env%image_charge_pot%eta,0,&
          0,0.0_dp,0,ra,(/0.0_dp,0.0_dp,0.0_dp/),0.0_dp,&
          rs_v,cell,pw_env%cube_info(1),hab,pab,o1=0,o2=0,&
          eps_gvg_rspace=eps_rho_rspace, calculate_forces=.TRUE.,&
          force_a=force_a,force_b=force_b,use_subpatch=.TRUE.,&
          subpatch_pattern=0_int_8,error=error)
   
     force_a(:)=coeff(iatom)*force_a(:)
     forces(:,iatom) = force_a(:)
   
    END DO

    CALL mp_sum(forces,para_env%group)

    DEALLOCATE (hab,pab,cores,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL rs_grid_release(rs_v, error=error)

    ! print info on gradients if wanted
    CALL print_gradients_image_atoms(forces, qs_env, error)

    CALL timestop(handle)


  END SUBROUTINE integrate_potential_devga_rspace

!****************************************************************************
!> \brief calculate image matrix T depending on constraints on image atoms
!>        in case coefficients are estimated not iteratively
!> \param qs_env qs environment
!> \param qmmm_env qmmm environment
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE conditional_calc_image_matrix(qs_env,qmmm_env,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qmmm_env_qm_type), POINTER          :: qmmm_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'conditional_calc_image_matrix', &
      routineP = moduleN//':'//routineN

    IF(.NOT.qmmm_env%image_charge_pot%coeff_iterative) THEN
      SELECT CASE(qmmm_env%image_charge_pot%state_image_matrix)
      CASE(calc_always)
          CALL calculate_image_matrix(image_matrix=qs_env%image_matrix,&
          ipiv=qs_env%ipiv,qs_env=qs_env,qmmm_env=qmmm_env,&
          error=error)
      CASE(calc_once)
          !if all image atoms are fully constrained, calculate image matrix
          !only for the first MD or GEO_OPT step
          CALL calculate_image_matrix(image_matrix=qs_env%image_matrix,&
          ipiv=qs_env%ipiv,qs_env=qs_env,qmmm_env=qmmm_env,&
          error=error)
          qmmm_env%image_charge_pot%state_image_matrix = calc_once_done
          CALL cp_assert(.NOT.qmmm_env%center_qm_subsys0,cp_warning_level,&
               cp_assertion_failed,routineP,"The image atoms are fully "//&
               "constrained and the image matrix is only calculated once. "//&
               "To be safe, set CENTER to NEVER " //&
               CPSourceFileRef,only_ionode=.TRUE.)
      CASE(calc_once_done)
          ! do nothing image matrix is stored
      CASE DEFAULT
         CALL cp_assert(.FALSE.,cp_fatal_level,cp_assertion_failed,&
                   routineP,"No initialization for image charges available?",&
                   only_ionode=.TRUE.)
      END SELECT
    ENDIF

  END SUBROUTINE conditional_calc_image_matrix

!****************************************************************************
!> \brief calculate image matrix T
!> \param image_matrix matrix T
!> \param ipiv pivoting prior to DGETRS (for Gaussian elimination)
!> \param qs_env qs environment
!> \param qmmm_env qmmm environment
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE calculate_image_matrix(image_matrix,ipiv,qs_env,qmmm_env,error)

    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: image_matrix
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: ipiv
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qmmm_env_qm_type), POINTER          :: qmmm_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_image_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iatom, iatom_ref, j, &
                                                k, natom, output_unit, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:), POINTER     :: int_res
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type)                          :: rho_gb, vb_gspace, vb_rspace
    TYPE(pw_poisson_type), POINTER           :: poisson_env
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(section_vals_type), POINTER         :: input

    CALL timeset(routineN,handle)
   
    NULLIFY(pw_env, auxbas_pw_pool, poisson_env, para_env, int_res,&
            input, logger)
    failure=.FALSE.
    logger => cp_error_get_logger(error)
   
    natom=SIZE(qmmm_env%image_charge_pot%image_mm_list)
   
    IF (.NOT.ASSOCIATED(image_matrix)) THEN
        ALLOCATE(image_matrix(natom,natom),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ENDIF
    IF (PRESENT(ipiv)) THEN
       IF (.NOT.ASSOCIATED(ipiv)) THEN
           ALLOCATE(ipiv(natom),stat=stat)
           CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ENDIF
    ENDIF
    ALLOCATE(int_res(natom),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
   
    image_matrix=0.0_dp
    IF(PRESENT(ipiv)) ipiv=0
   
   
    IF (.NOT.failure) THEN
        CALL get_qs_env(qs_env, pw_env=pw_env,para_env=para_env,&
                        input=input,error=error)
   
        !print info
        output_unit=cp_print_key_unit_nr(logger,input,&
                    "QMMM%PRINT%PROGRAM_RUN_INFO",&
                    extension=".qmmmLog",error=error)
        IF (qmmm_env%image_charge_pot%coeff_iterative) THEN
           IF (output_unit>0)  WRITE (UNIT=output_unit,FMT="(T3,A)")&
                    "Calculating image matrix"
        ELSE
           IF (output_unit>0)  WRITE (UNIT=output_unit,FMT="(T2,A)")&
                    "Calculating image matrix"
        ENDIF
        CALL cp_print_key_finished_output(output_unit,logger,input,&
                    "QMMM%PRINT%PROGRAM_RUN_INFO", error=error)
   
   
        CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool,&
                        poisson_env=poisson_env, error=error)
        CALL pw_pool_create_pw(auxbas_pw_pool,&
                               rho_gb%pw,&
                               use_data=COMPLEXDATA1D,&
                               in_space=RECIPROCALSPACE,&
                               error=error)
        CALL pw_pool_create_pw(auxbas_pw_pool,&
                               vb_gspace%pw,&
                               use_data=COMPLEXDATA1D,&
                               in_space=RECIPROCALSPACE,&
                               error=error)
        CALL pw_pool_create_pw(auxbas_pw_pool,&
                               vb_rspace%pw,&
                               use_data=REALDATA3D,&
                               in_space=REALSPACE,&
                               error=error)
        ! calculate vb only once for one reference atom 
        iatom_ref = 1 !
        !collocate gaussian of reference MM atom on grid
        CALL pw_zero(rho_gb%pw, error=error)
        CALL calculate_rho_single_gaussian(rho_gb,qs_env,iatom_ref,error)
        !calculate potential vb like hartree potential
        CALL pw_zero(vb_gspace%pw, error=error)
        CALL pw_poisson_solve(poisson_env,rho_gb%pw,vhartree=vb_gspace%pw,&
                              error=error)
        CALL pw_zero(vb_rspace%pw, error=error)
        CALL pw_transfer(vb_gspace%pw,vb_rspace%pw,error=error)
        CALL pw_scale(vb_rspace%pw,vb_rspace%pw%pw_grid%dvol,error=error)
        DO iatom=1,natom
          !calculate integral vb_rspace*ga
          int_res=0.0_dp
          CALL integrate_potential_ga_rspace(vb_rspace,qs_env%qmmm_env_qm,&
                                             qs_env,int_res,atom_num=iatom,&
                                             atom_num_ref=iatom_ref,error=error)
          image_matrix(iatom,iatom:natom)=int_res(iatom:natom)
          image_matrix(iatom+1:natom,iatom)=int_res(iatom+1:natom)
        END DO
   
        IF (qmmm_env%image_charge_pot%coeff_iterative) THEN
           !inversion --> preconditioner matrix for CG
           CALL DPOTRF('L',natom,qs_env%image_matrix,natom,stat)
           CPPostconditionNoFail(stat==0,cp_failure_level,routineP,error)
           CALL DPOTRI('L',natom,qs_env%image_matrix,natom,stat)
           CPPostconditionNoFail(stat==0,cp_failure_level,routineP,error)
           DO j=1,natom
            DO k=j+1,natom
             qs_env%image_matrix(j,k)=qs_env%image_matrix(k,j)
            ENDDO
           ENDDO
           CALL write_image_matrix(qs_env%image_matrix,qs_env,qmmm_env,error)
        ELSE
           !pivoting prior to DGETRS (Gaussian elimination)
           IF(PRESENT(ipiv)) THEN
              CALL DGETRF(natom,natom,image_matrix,natom,ipiv,stat)
              CPPostconditionNoFail(stat==0,cp_failure_level,routineP,error)
           ENDIF
        ENDIF
   
        CALL pw_release(vb_gspace%pw,error=error)
        CALL pw_release(vb_rspace%pw,error=error)
        CALL pw_release(rho_gb%pw,error=error)
   
        DEALLOCATE(int_res,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    ENDIF
   
    CALL timestop(handle)

  END SUBROUTINE calculate_image_matrix

! *****************************************************************************
!> \brief calculates potential of the metal (image potential) given a set of
!>        coefficents coeff
!> \param v_metal_rspace potential generated by rho_metal in real space
!> \param coeff expansion coefficients of the image charge density, i.e.
!>        rho_metal=sum_a c_a*g_a
!> \param rho_hartree_gspace Kohn Sham density in reciprocal space
!> \param energy structure where energies are stored
!> \param qs_env qs environment
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE calculate_potential_metal(v_metal_rspace,coeff,rho_hartree_gspace,energy,&
                                     qs_env,error)

    TYPE(pw_p_type), INTENT(INOUT)           :: v_metal_rspace
    REAL(KIND=dp), DIMENSION(:), POINTER     :: coeff
    TYPE(pw_p_type), INTENT(IN), OPTIONAL    :: rho_hartree_gspace
    TYPE(qs_energy_type), OPTIONAL, POINTER  :: energy
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calculate_potential_metal', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: en_external, &
                                                en_vmetal_rhohartree, &
                                                total_rho_metal
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type)                          :: rho_metal, v_metal_gspace
    TYPE(pw_poisson_type), POINTER           :: poisson_env
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool

    CALL timeset(routineN,handle)
   
    NULLIFY(pw_env,auxbas_pw_pool,poisson_env)
    failure=.FALSE.
    en_vmetal_rhohartree=0.0_dp
    en_external=0.0_dp
   
    CALL get_qs_env(qs_env=qs_env,pw_env=pw_env,error=error)
    CALL pw_env_get(pw_env,auxbas_pw_pool=auxbas_pw_pool,&
                      poisson_env=poisson_env,error=error)
   
    CALL pw_pool_create_pw(auxbas_pw_pool,&
                               rho_metal%pw,&
                               use_data=COMPLEXDATA1D,&
                               in_space=RECIPROCALSPACE,&
                               error=error)
   
    CALL pw_pool_create_pw(auxbas_pw_pool,&
                            v_metal_gspace%pw, &
                            use_data=COMPLEXDATA1D,&
                            in_space=RECIPROCALSPACE,&
                            error=error)
   
    CALL pw_pool_create_pw(auxbas_pw_pool,&
                           v_metal_rspace%pw,&
                           use_data=REALDATA3D,&
                           in_space=REALSPACE,&
                           error=error)
   
    CALL pw_zero(rho_metal%pw,error=error)
    CALL calculate_rho_metal(rho_metal,coeff,total_rho_metal=total_rho_metal,&
                            qs_env=qs_env,error=error)
   
    CALL pw_zero(v_metal_gspace%pw, error=error)
    CALL pw_poisson_solve(poisson_env,rho_metal%pw,&
                         vhartree=v_metal_gspace%pw,error=error)
   
    IF(PRESENT(rho_hartree_gspace)) THEN
      en_vmetal_rhohartree=0.5_dp*pw_integral_ab(v_metal_gspace%pw,&
                                  rho_hartree_gspace%pw,error=error)
      en_external=qs_env%qmmm_env_qm%image_charge_pot%V0*total_rho_metal
      energy%image_charge=en_vmetal_rhohartree-0.5_dp*en_external
      CALL print_image_energy_terms(en_vmetal_rhohartree,en_external,&
                                    total_rho_metal,qs_env,error)
    ENDIF
   
    CALL pw_zero(v_metal_rspace%pw, error=error)
    CALL pw_transfer(v_metal_gspace%pw,v_metal_rspace%pw,error=error)
    CALL pw_scale(v_metal_rspace%pw,v_metal_rspace%pw%pw_grid%dvol,&
                  error=error)
    CALL pw_release(v_metal_gspace%pw,error=error)
    CALL pw_release(rho_metal%pw,error=error)
   
    CALL timestop(handle)

  END SUBROUTINE calculate_potential_metal

! ****************************************************************************
!> \brief Add potential of metal (image charge pot) to Hartree Potential
!> \param v_hartree Hartree potential (in real space)
!> \param v_metal potential generated by rho_metal (in real space)
!> \param qs_env qs environment
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE add_image_pot_to_hartree_pot(v_hartree,v_metal,qs_env,error)

    TYPE(pw_p_type), INTENT(INOUT)           :: v_hartree
    TYPE(pw_p_type), INTENT(IN)              :: v_metal
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'add_image_pot_to_hartree_pot', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, output_unit
    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: input

    CALL timeset(routineN,handle)
   
    NULLIFY(input,logger)
    failure=.FALSE.
    logger => cp_error_get_logger(error)
   
    !add image charge potential
    v_hartree%pw%cr3d = v_hartree%pw%cr3d + v_metal%pw%cr3d
   
    ! print info
    CALL get_qs_env(qs_env=qs_env,&
                      input=input,&
                      error=error)
    output_unit=cp_print_key_unit_nr(logger,input,&
                    "QMMM%PRINT%PROGRAM_RUN_INFO",&
                    extension=".qmmmLog",error=error)
    IF (output_unit>0)  WRITE (UNIT=output_unit,FMT="(T3,A)")&
                    "Adding image charge potential to the Hartree potential."
    CALL cp_print_key_finished_output(output_unit,logger,input,&
                     "QMMM%PRINT%PROGRAM_RUN_INFO", error=error)
   
    CALL timestop(handle)

  END SUBROUTINE add_image_pot_to_hartree_pot

!****************************************************************************
!> \brief writes image matrix T to file when used as preconditioner for
!>        calculating image coefficients iteratively 
!> \param image_matrix matrix T
!> \param qs_env qs environment
!> \param qmmm_env qmmm environment
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE write_image_matrix(image_matrix,qs_env,qmmm_env,error)

    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: image_matrix
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qmmm_env_qm_type), POINTER          :: qmmm_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'write_image_matrix', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: filename
    INTEGER                                  :: handle, rst_unit
    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: print_key, qmmm_section

    CALL timeset(routineN,handle)

    NULLIFY(qmmm_section, print_key, logger, para_env)
    logger => cp_error_get_logger(error)
    failure = .FALSE.
    rst_unit = -1

    CALL get_qs_env(qs_env=qs_env,para_env=para_env,&
                    input=qmmm_section,&
                    error=error)

    print_key => section_vals_get_subs_vals(qmmm_section,&
                                           "QMMM%PRINT%IMAGE_CHARGE_RESTART",&
                                           error=error)

    IF (BTEST(cp_print_key_should_output(logger%iter_info,&
              qmmm_section,"QMMM%PRINT%IMAGE_CHARGE_RESTART",error=error),&
              cp_p_file)) THEN

       rst_unit = cp_print_key_unit_nr(logger,qmmm_section,&
                          "QMMM%PRINT%IMAGE_CHARGE_RESTART",&
                           extension=".Image",&
                           file_status="REPLACE",&
                           file_action="WRITE",&
                           file_form="UNFORMATTED",&
                           error=error)

       IF(rst_unit>0) filename = cp_print_key_generate_filename(logger,&
                                 print_key, extension=".IMAGE", &
                                 my_local=.FALSE.,error=error)

       IF(rst_unit>0) THEN
         WRITE(rst_unit) image_matrix
       ENDIF
       
       CALL cp_print_key_finished_output(rst_unit,logger,qmmm_section,&
                       "QMMM%PRINT%IMAGE_CHARGE_RESTART", error=error)
    ENDIF

    CALL timestop(handle)

  END SUBROUTINE write_image_matrix

!****************************************************************************
!> \brief restarts image matrix T when used as preconditioner for calculating
!>        image coefficients iteratively 
!> \param image_matrix matrix T
!> \param qs_env qs environment
!> \param qmmm_env qmmm environment
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE restart_image_matrix(image_matrix,qs_env,qmmm_env,error)

    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: image_matrix
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qmmm_env_qm_type), POINTER          :: qmmm_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'restart_image_matrix', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: image_filename
    INTEGER                                  :: handle, natom, output_unit, &
                                                rst_unit, stat
    LOGICAL                                  :: exist, failure
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: qmmm_section

    CALL timeset(routineN,handle)

    NULLIFY(qmmm_section, logger, para_env)
    logger => cp_error_get_logger(error)
    failure = .FALSE.
    exist= .FALSE.
    rst_unit=-1

    natom=SIZE(qmmm_env%image_charge_pot%image_mm_list)
 
    IF (.NOT.ASSOCIATED(image_matrix)) THEN
       ALLOCATE(image_matrix(natom,natom),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ENDIF
   
    image_matrix=0.0_dp 
    
    CALL get_qs_env(qs_env=qs_env,para_env=para_env,&
                    input=qmmm_section,&
                    error=error)

    CALL section_vals_val_get(qmmm_section,"QMMM%IMAGE_CHARGE%IMAGE_RESTART_FILE_NAME",&
                               c_val=image_filename, error=error)

    INQUIRE(FILE=image_filename, exist=exist)
 
    IF(exist) THEN
      IF(para_env%ionode) THEN
       CALL open_file(file_name=image_filename,&
                      file_status="OLD",&
                      file_form="UNFORMATTED",& 
                      file_action="READ",&
                      unit_number=rst_unit)

       READ(rst_unit) qs_env%image_matrix
      ENDIF
      
      CALL mp_bcast(qs_env%image_matrix,para_env%source,para_env%group)

      IF (para_env%ionode) CALL close_file(unit_number=rst_unit)

      output_unit=cp_print_key_unit_nr(logger,qmmm_section,&
                             "QMMM%PRINT%PROGRAM_RUN_INFO",&
                            extension=".qmmmLog",error=error)
      IF (output_unit>0)  WRITE (UNIT=output_unit,FMT="(T3,A)")&
                             "Restarted image matrix"
    ELSE
      CALL stop_program(routineN,moduleN,__LINE__,&
           "Restart file for image matrix not found")
    ENDIF

    qmmm_env%image_charge_pot%image_restart=.FALSE.
  
    CALL timestop(handle)

  END SUBROUTINE restart_image_matrix 

! ****************************************************************************
!> \brief Print info on image gradients on image MM atoms
!> \param forces structure storing the force contribution of the image charges 
!>        for the metal (MM) atoms (actually these are only the gradients)
!> \param qs_env qs environment
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE print_gradients_image_atoms(forces,qs_env,error)


    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: forces
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'print_gradients_image_atoms', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: atom_a, iatom, natom, &
                                                output_unit
    REAL(KIND=dp), DIMENSION(3)              :: sum_gradients
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: input

    NULLIFY(input, logger)
    logger => cp_error_get_logger(error)

    sum_gradients=0.0_dp
    natom=SIZE(qs_env%qmmm_env_qm%image_charge_pot%image_mm_list)

    DO iatom=1,natom
      sum_gradients(:)=sum_gradients(:)+forces(:,iatom)
    ENDDO

    CALL get_qs_env(qs_env=qs_env,input=input,error=error)

    output_unit = cp_print_key_unit_nr(logger,input,&
                 "QMMM%PRINT%DERIVATIVES",extension=".Log",error=error)
    IF(output_unit>0) THEN
       WRITE (unit=output_unit,fmt="(/1X,A)")&
        "Image gradients [a.u.] on MM image charge atoms after QMMM calculation: "
       WRITE (unit=output_unit,fmt="(T4,A4,T27,A1,T50,A1,T74,A1)")&
                                    "Atom", "X", "Y", "Z"
           DO iatom=1, natom
              atom_a=qs_env%qmmm_env_qm%image_charge_pot%image_mm_list(iatom)
              WRITE (unit=output_unit,fmt="(T2,I6,T22,ES12.5,T45,ES12.5,T69,ES12.5)")&
               atom_a,forces(:,iatom)
           END DO

       WRITE (unit=output_unit,fmt="(T2,A)") REPEAT("-",79)
       WRITE (unit=output_unit,fmt="(T2,A,T22,ES12.5,T45,ES12.5,T69,ES12.5)")&
        "sum gradients:",sum_gradients
       WRITE (unit=output_unit,fmt="(/)")
    ENDIF

    CALL cp_print_key_finished_output(output_unit,logger,input,&
           "QMMM%PRINT%DERIVATIVES",error=error)

  END SUBROUTINE print_gradients_image_atoms

! ****************************************************************************
!> \brief Print image coefficients
!> \param image_coeff expansion coefficients of the image charge density
!> \param qs_env qs environment
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE print_image_coefficients(image_coeff,qs_env,error)

    REAL(KIND=dp), DIMENSION(:), POINTER     :: image_coeff
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'print_image_coefficients', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: atom_a, iatom, natom, &
                                                output_unit
    REAL(KIND=dp)                            :: normalize_factor, sum_coeff
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: input

    NULLIFY(input, logger)
    logger => cp_error_get_logger(error)
    
    sum_coeff=0.0_dp
    natom=SIZE(qs_env%qmmm_env_qm%image_charge_pot%image_mm_list)
    normalize_factor=SQRT((qs_env%qmmm_env_qm%image_charge_pot%eta/pi)**3)
    
    DO iatom=1,natom
      sum_coeff=sum_coeff+image_coeff(iatom)
    ENDDO
    
    CALL get_qs_env(qs_env=qs_env,input=input,error=error)
    
    output_unit = cp_print_key_unit_nr(logger,input,&
                "QMMM%PRINT%IMAGE_CHARGE_INFO",extension=".Log",error=error)
    IF(output_unit>0) THEN
      WRITE (unit=output_unit,fmt="(/)")
      WRITE (unit=output_unit,fmt="(T2,A)")&
       "Image charges [a.u.] after QMMM calculation: "
      WRITE (unit=output_unit,fmt="(T4,A4,T67,A)") "Atom", "Image charge"
      WRITE(unit=output_unit,fmt="(T4,A,T67,A)") REPEAT("-",4), REPEAT("-",12)
    
         DO iatom=1, natom
            atom_a=qs_env%qmmm_env_qm%image_charge_pot%image_mm_list(iatom)
            !opposite sign for image_coeff; during the calculation they have
            !the 'wrong' sign to ensure consistency with v_hartree which has
            !the opposite sign
            WRITE (unit=output_unit,fmt="(T2,I6,T65,ES16.9)")&
                                atom_a, -image_coeff(iatom)/normalize_factor
         ENDDO
    
      WRITE(unit=output_unit,fmt="(T2,A)") REPEAT("-",79)
      WRITE (unit=output_unit,fmt="(T2,A,T65,ES16.9)") &
                             "sum image charges:",-sum_coeff/normalize_factor
    ENDIF
    
    CALL cp_print_key_finished_output(output_unit,logger,input,&
           "QMMM%PRINT%IMAGE_CHARGE_INFO",error=error)

  END SUBROUTINE print_image_coefficients

! ****************************************************************************
!> \brief Print detailed image charge energies
!> \param en_vmetal_rhohartree energy contribution of the image charges 
!>        without external potential, i.e. 0.5*integral(v_metal*rho_hartree) 
!> \param en_external additional energy contribution of the image charges due 
!>        to an external potential, i.e. V0*total_rho_metal
!> \param total_rho_metal total induced image charge density
!> \param qs_env qs environment
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE print_image_energy_terms(en_vmetal_rhohartree,en_external,&
                                    total_rho_metal,qs_env,error)

    REAL(KIND=dp), INTENT(IN)                :: en_vmetal_rhohartree, &
                                                en_external, total_rho_metal
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'print_image_energy_terms', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: output_unit
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: input

    NULLIFY(input, logger)
    logger => cp_error_get_logger(error)
    
    CALL get_qs_env(qs_env=qs_env,input=input,error=error)
    
    output_unit = cp_print_key_unit_nr(logger,input,&
                "QMMM%PRINT%IMAGE_CHARGE_INFO",extension=".Log",error=error)
    
    IF(output_unit>0) THEN
      WRITE (unit=output_unit,fmt="(T3,A,T56,F25.14)")&
      "Total induced charge density [a.u.]:", total_rho_metal
      WRITE (unit=output_unit,fmt="(T3,A)")"Image energy terms: "
      WRITE (unit=output_unit,fmt="(T3,A,T56,F25.14)")&
      "Coulomb energy of QM and image charge density [a.u.]:", en_vmetal_rhohartree
      WRITE (unit=output_unit,fmt="(T3,A,T56,F25.14)")&
      "External potential energy term [a.u.]:", -0.5_dp*en_external
      WRITE (unit=output_unit,fmt="(T3,A,T56,F25.14)")&
      "Total image charge energy [a.u.]:", en_vmetal_rhohartree-0.5_dp*en_external
    ENDIF
    
    CALL cp_print_key_finished_output(output_unit,logger,input,&
           "QMMM%PRINT%IMAGE_CHARGE_INFO",error=error)

  END SUBROUTINE print_image_energy_terms

END MODULE qmmm_image_charge
