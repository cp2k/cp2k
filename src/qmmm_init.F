!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2004  CP2K developers group                                 !
!-----------------------------------------------------------------------------!


!!****h* cp2k/qmmm_init [1.0] *
!!
!!   NAME
!!     qmmm_init
!!
!!   FUNCTION
!!     Initialize a QM/MM calculation
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     5.2004 created [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE qmmm_init
  USE atomic_kind_list_types,          ONLY: atomic_kind_list_create,&
                                             atomic_kind_list_release,&
                                             atomic_kind_list_type
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             read_atomic_kind_set, &
                                             set_atomic_kind
  USE cp_output_handling, ONLY: cp_print_key_should_output, cp_p_file
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsystem_types,              ONLY: cp_subsys_create,&
                                             cp_subsys_get,&
                                             cp_subsys_release,&
                                             cp_subsys_set,&
                                             cp_subsystem_p_type,&
                                             cp_subsystem_type
  USE force_env_methods,               ONLY: force_env_create
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_p_type,&
                                             force_env_type
  USE force_fields,                    ONLY: force_field_control,&
                                             read_lj_section,&
                                             read_wl_section,&
                                             read_gd_section
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE mol_kind_new_list_types,         ONLY: mol_kind_new_list_create,&
                                             mol_kind_new_list_release,&
                                             mol_kind_new_list_type
  USE mol_new_list_types,              ONLY: mol_new_list_create,&
                                             mol_new_list_release,&
                                             mol_new_list_type
  USE molecule_kind_types,             ONLY: molecule_kind_type
  USE molecule_types_new,              ONLY: molecule_type
  USE particle_list_types,             ONLY: particle_list_create,&
                                             particle_list_release,&
                                             particle_list_type
  USE particle_types,                  ONLY: particle_type
  USE qmmm_types,                      ONLY: qmmm_env_qm_create,&
                                             qmmm_env_qm_release,&
                                             qmmm_env_qm_type,&
                                             qmmm_env_mm_type,&
                                             qmmm_env_mm_create,&
                                             qmmm_env_mm_release,&
                                             qmmm_links_type,&
                                             add_set_type,&
                                             create_add_set_type,&
                                             add_set_release
  USE qmmm_ff_fist,                    ONLY: qmmm_ff_precond_only_qm
  USE qmmm_gaussian_types,             ONLY: qmmm_gaussian_p_type
  USE qmmm_gaussian_init,              ONLY: qmmm_gaussian_initialize
  USE qmmm_pw_grid,                    ONLY: qmmm_pw_grid_init
  USE qmmm_elpot,                      ONLY: qmmm_potential_init
  USE qmmm_per_elpot,                  ONLY: qmmm_per_potential_init
  USE qs_environment,                  ONLY: qs_init_subsys
  USE qs_environment_types,            ONLY: qs_env_create,&
                                             qs_env_release,&
                                             qs_environment_type,&
                                             set_qs_env
  USE simulation_cell,                 ONLY: cell_retain,&
                                             cell_type,&
                                             init_qmmm_cell,&
                                             cell_clone,&
                                             cell_release,&
                                             write_cell,&
                                             get_cell
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE topology,                        ONLY: connectivity_control,&
                                             coordinate_control
  USE topology_pdb,                    ONLY: write_coordinate_pdb
  USE topology_psf,                    ONLY: write_topology_psf
  USE topology_types,                  ONLY: deallocate_topology,&
                                             init_topology,&
                                             topology_parameters_type
  USE topology_util,                   ONLY: topology_connectivity_pack,&
                                             topology_constraint_pack,&
                                             topology_coordinate_pack
  USE external_potential_types,        ONLY: elp_potential_type,&
                                             allocate_potential,&
                                             set_potential,&
                                             get_potential
  USE fist_main,                       ONLY: fist_create_force_env
  USE input_section_types, ONLY: section_type, section_create, &
       section_add_keyword, section_add_subsection, section_release,&
       section_vals_type,&
       section_vals_get_subs_vals, section_vals_get, section_vals_val_get
  USE input_constants, ONLY: RADIUS_QMMM_DEFAULT, do_qmmm_link_imomm, do_qmmm_link_gho,&
                        do_qmmm_link_pseudo,&
                        do_qmmm_gauss, do_center_none, do_center_grid,&
                        do_fist, do_qs, do_qmmm
  USE fist_environment_types,  ONLY: set_fist_env
  USE force_field_types,               ONLY: input_info_type
  USE pair_potential_types,            ONLY: pair_potential_reallocate
  USE cp_subsystem_methods,            ONLY: create_small_subsystem
  USE qs_main,                         ONLY: quickstep_create_force_env
#include "cp_common_uses.h"
  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='qmmm_init'

  PUBLIC :: qmmm_create_force_env
!***
!****************************************************************************
CONTAINS

!!****f* qmmm_init/qmmm_create_force_env [1.0] *
!!
!!   NAME
!!     qmmm_create_force_env
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     05.2004 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE  qmmm_create_force_env(force_env, globenv, error)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_create_force_env', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length), &
      DIMENSION(:), POINTER                  :: qm_atom_type
    INTEGER                                  :: handle, isubf, stat
    INTEGER                                  :: iw
    INTEGER, DIMENSION(:), POINTER           :: qm_atom_index, mm_atom_index, qm_subsys_index
    INTEGER, DIMENSION(:), POINTER           :: mm_link_atoms
    LOGICAL                                  :: failure, qmmm_link
    TYPE(cell_type), POINTER                 :: mm_cell, qm_cell_small, super_cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys_mm
    TYPE(cp_subsystem_type), POINTER         :: subsys_qm
    TYPE(force_env_p_type), DIMENSION(:), &
      POINTER                                :: sub_force_env
    TYPE(qmmm_env_qm_type), POINTER          :: qmmm_env_qm
    TYPE(qmmm_env_mm_type), POINTER          :: qmmm_env_mm
    TYPE(qs_environment_type), POINTER       :: qs_env
    REAL(dp), DIMENSION(:), POINTER          :: mm_atom_chrg, mm_el_pot_radius, mm_el_pot_radius_corr 
    REAL(KIND=dp)                            :: eps_mm_rspace
    INTEGER                                  :: qmmm_coupl_type
    TYPE(section_vals_type), POINTER         :: force_env_section, qmmm_section, poisson_section, print_section,&
                                                qmmm_periodic
    TYPE(qmmm_links_type), pointer           :: qmmm_links
    REAL(KIND=dp),  DIMENSION(:), POINTER    :: mm_link_scale_factor, fist_scale_charge_link
    LOGICAL                                  :: move_mm_charges, add_mm_charges, nocenter
    LOGICAL                                  :: nocompatibility, explicit
    TYPE(add_set_type), POINTER              :: added_charges
    REAL(KIND=dp), DIMENSION(3)              :: abc
    TYPE(section_vals_type), POINTER :: print_gen
    TYPE(cp_logger_type), POINTER :: logger

    CALL timeset(routineN,"I","",handle)
    failure=.FALSE.

    NULLIFY(sub_force_env, qm_atom_index, mm_atom_index, qm_atom_type,   &
         qmmm_env_qm, subsys_mm,subsys_qm, mm_cell, qm_cell_small,       &
         qs_env, para_env, mm_atom_chrg,mm_el_pot_radius, qmmm_env_mm,   &
         mm_link_atoms, mm_link_scale_factor, qmmm_links, added_charges, &
         fist_scale_charge_link, poisson_section, print_section,  &
         print_gen, logger, mm_el_pot_radius_corr, super_cell)
    para_env => globenv%para_env
    iw = globenv%scr
    logger => cp_error_get_logger(error)
    IF (.NOT. failure) THEN
       ! Input section...
       force_env_section => section_vals_get_subs_vals(globenv%input_file,"FORCE_EVAL",error=error)
       qmmm_section => section_vals_get_subs_vals(force_env_section,"QMMM",error=error)
       qmmm_periodic => section_vals_get_subs_vals(qmmm_section,"PERIODIC",error=error)
       print_section => section_vals_get_subs_vals(qmmm_section,"PRINT",error=error)
       poisson_section => section_vals_get_subs_vals(force_env_section,"DFT%POISSON_DFT",error=error)
       print_gen => section_vals_get_subs_vals(print_section,"PROGRAM_RUN_INFO",error=error)
       ALLOCATE(sub_force_env(2),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DO isubf=1,2
          NULLIFY(sub_force_env(isubf)%force_env)
       END DO
       ! Create QM/MM Environments..
       CALL qmmm_env_qm_create(qmmm_env_qm,error=error)
       CALL qmmm_env_mm_create(qmmm_env_mm,error=error)

       ! Set up QM/MM Options
       CALL setup_qmmm_vars_mm(qmmm_section,&
                               qmmm_env_mm,&
                               qm_atom_index,&
                               mm_link_atoms,&
                               mm_link_scale_factor,&
                               fist_scale_charge_link,&
                               qmmm_coupl_type,&
                               qmmm_link,&
                               error=error)

       qmmm_env_mm%qm_atom_index          => qm_atom_index       
       qmmm_env_mm%mm_link_atoms          => mm_link_atoms
       qmmm_env_mm%mm_link_scale_factor   => mm_link_scale_factor
       qmmm_env_mm%fist_scale_charge_link => fist_scale_charge_link
       qmmm_env_mm%qmmm_coupl_type        =  qmmm_coupl_type
       qmmm_env_mm%qmmm_link              =  qmmm_link
       ! Center the qm subsys into the qm box
       CALL section_vals_val_get(qmmm_section,"NOCENTER",l_val=nocenter,error=error)
       qmmm_env_qm%center_qm_subsys  = .NOT.nocenter

       CALL section_vals_val_get(qmmm_section,"NOCENTER0",l_val=nocenter,error=error)
       qmmm_env_qm%center_qm_subsys0 = .NOT.nocenter

       ! Compatibility with the QM/MM in CPMD code
       CALL section_vals_val_get(qmmm_section,"NOCOMPATIBILITY",l_val=nocompatibility,error=error)
       qmmm_env_qm%compatibility = .NOT.nocompatibility       

       ! Periodic boundary condition calculation
       CALL section_vals_get(qmmm_periodic,explicit=explicit,error=error)
       qmmm_env_qm%periodic = explicit

       ! First Initialize Fist...
       globenv%prog_name_id=do_fist
       CALL fist_create_force_env(sub_force_env(1)%force_env, &
            globenv, qmmm=.TRUE., qmmm_env=qmmm_env_mm, error=error)

       CALL force_env_get(sub_force_env(1)%force_env,subsys=subsys_mm,&
            cell=mm_cell,error=error) 
       CPPrecondition(SIZE(subsys_mm)==1,cp_failure_level,routineP,error,failure)

       IF (.NOT.ASSOCIATED(qm_cell_small)) THEN
          ALLOCATE (qm_cell_small,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       CALL cell_clone(mm_cell, qm_cell_small )

       ! Set up QM/MM Options       
       CALL setup_qmmm_vars_qm (qmmm_section,&
                                qmmm_env_qm,&
                                subsys_mm,&
                                qm_atom_type,&
                                qm_atom_index,&
                                mm_atom_index,&
                                qm_subsys_index,&
                                qm_cell_small,&
                                qmmm_coupl_type,&
                                eps_mm_rspace,&
                                qmmm_link,&
                                error)
                               
       qmmm_env_qm%qm_atom_index   => qm_atom_index
       qmmm_env_qm%mm_atom_index   => mm_atom_index
       qmmm_env_qm%qm_subsys_index => qm_subsys_index
       qmmm_env_qm%eps_mm_rspace   = eps_mm_rspace
       qmmm_env_qm%qmmm_coupl_type = qmmm_coupl_type
       qmmm_env_qm%qmmm_link       = qmmm_link
       qmmm_env_qm%num_qm_atoms = size(qm_atom_index)
       qmmm_env_qm%num_mm_atoms = size(mm_atom_index)
       IF (BTEST(cp_print_key_should_output(logger%iter_info,print_gen,error=error),cp_p_file).AND.&
            logger%para_env%mepos==logger%para_env%source) THEN
          WRITE(iw,'(A,I6)')" Number of QM atoms: ",qmmm_env_qm%num_qm_atoms 
          WRITE(iw,'(A,I6)')" Number of MM atoms: ",qmmm_env_qm%num_mm_atoms
          WRITE(iw,'(A)')" QM cell ::"
          CALL write_cell(qm_cell_small, globenv)
       END IF
       CALL get_cell(qm_cell_small, abc=abc)
       !
       ! Assign charges and mm_el_pot_radius from fist_topology
       !
       ALLOCATE(mm_atom_chrg(SIZE(mm_atom_index)),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(mm_el_pot_radius(SIZE(mm_atom_index)),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(mm_el_pot_radius_corr(SIZE(mm_atom_index)),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       mm_atom_chrg     = 0.0_dp
       mm_el_pot_radius = 0.0_dp
       mm_el_pot_radius_corr = 0.0_dp

       CALL assign_mm_charges_and_radius(subsys=subsys_mm,&
                                         mm_atom_chrg=mm_atom_chrg,&
                                         mm_el_pot_radius=mm_el_pot_radius,&
                                         mm_el_pot_radius_corr=mm_el_pot_radius_corr,&
                                         mm_atom_index=mm_atom_index,&
                                         mm_link_atoms=mm_link_atoms,&
                                         mm_link_scale_factor=mm_link_scale_factor,&
                                         error=error)

       qmmm_env_qm%mm_atom_chrg          => mm_atom_chrg
       qmmm_env_qm%mm_el_pot_radius      => mm_el_pot_radius
       qmmm_env_qm%mm_el_pot_radius_corr => mm_el_pot_radius_corr

       IF (qmmm_link) THEN 
          CALL setup_qmmm_links(qmmm_section, qmmm_links, mm_el_pot_radius, &
               mm_el_pot_radius_corr, mm_atom_index, qm_cell_small, iw, error=error)
          qmmm_env_qm%qmmm_links    => qmmm_links

          CALL print_qmmm_links(globenv, qmmm_section, qmmm_links, error)

          CALL add_set_release  (qmmm_env_qm%added_charges, error=error)
          CALL move_or_add_atoms(qmmm_section, move_mm_charges, add_mm_charges, &
                                 mm_atom_chrg, mm_el_pot_radius, mm_el_pot_radius_corr,&
                                 added_charges, mm_atom_index, qm_cell_small, error)
          qmmm_env_qm%move_mm_charges =  move_mm_charges
          qmmm_env_qm%add_mm_charges  =  add_mm_charges
          qmmm_env_qm%added_charges   => added_charges
       END IF
       
       CALL print_qmmm_charges(mm_atom_index, mm_atom_chrg, mm_el_pot_radius,&
                               mm_el_pot_radius_corr, qmmm_env_qm%added_charges,&
                               globenv, qmmm_section, nocompatibility, error)

       globenv%prog_name_id=do_qs
       CALL create_small_subsystem(subsys_qm,&
            big_subsys=subsys_mm,small_para_env=para_env,&
            small_cell=qm_cell_small,sub_atom_index=qm_atom_index,&
            sub_subsys_index=qm_subsys_index,&
            sub_atom_kind_name=qm_atom_type, globenv=globenv,error=error)
       CALL quickstep_create_force_env(sub_force_env(2)%force_env,globenv,subsystem=subsys_qm,&
            cell=qm_cell_small, qmmm=.TRUE., error=error)
       CALL cp_subsys_release(subsys_qm,error=error)
       CALL force_env_get(sub_force_env(2)%force_env,qs_env=qs_env,error=error)
       IF (qmmm_env_qm%periodic) THEN 
          IF (.NOT.ASSOCIATED(super_cell)) THEN
             ALLOCATE (super_cell,STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF
          CALL cell_clone(mm_cell, super_cell)
          CALL set_qs_env(qs_env, super_cell=super_cell, error=error)
          CALL cell_release(super_cell, error=error)
       END IF
       !
       ! Initialize MM Potential fitted with Gaussian
       !
       CALL qmmm_init_gaussian_type (qmmm_env_qm = qmmm_env_qm,&
                                     qm_cell_small=qm_cell_small,&
                                     mm_cell=mm_cell,&
                                     globenv=globenv,&
                                     qs_env=qs_env,&
                                     mm_atom_chrg=mm_atom_chrg,&
                                     added_charges=qmmm_env_qm%added_charges,&
                                     print_section=print_section,&
                                     qmmm_section=qmmm_section,&
                                     error=error)
       !
       ! Initialize the MM potential stored on vector
       !
       CALL qmmm_init_potential(qmmm_env_qm=qmmm_env_qm,&
                                qm_cell_small=qm_cell_small,&
                                mm_cell=mm_cell,&
                                globenv=globenv,&
                                qs_env=qs_env,&
                                added_charges=qmmm_env_qm%added_charges,&
                                print_section=print_section,&
                                error=error)
       !
       ! Initialize the qmmm_pw_grid
       !
       CALL qmmm_pw_grid_init(qmmm_env=qmmm_env_qm,&
                              pw_env=qs_env%pw_env,&
                              cell=qm_cell_small,&
                              para_env=para_env,&
                              error=error)
       !
       ! Initialize the MM periodic potential
       !
       CALL qmmm_init_periodic_potential(qmmm_env_qm=qmmm_env_qm,&
                                qm_cell_small=qm_cell_small,&
                                mm_cell=mm_cell,&
                                globenv=globenv,&
                                qs_env=qs_env,&
                                added_charges=qmmm_env_qm%added_charges,&
                                qmmm_periodic=qmmm_periodic,&
                                print_section=print_section,&
                                mm_atom_chrg=mm_atom_chrg,&
                                error=error)
       !
       ! Preparing for PBC...
       !
       CALL setup_origin_mm_cell(qmmm_section, qmmm_env_qm, mm_cell, qm_cell_small,&
            dr=qs_env%pw_env%pw_pools(qs_env%pw_env%auxbas_grid)%pool%pw_grid%dr,&
            error=error)

       CALL cell_release(qm_cell_small, error)                       

       globenv%prog_name_id=do_qmmm
       CALL force_env_create(force_env, sub_force_env=sub_force_env,&
            qmmm_env=qmmm_env_qm,globenv=globenv,error=error)

       CALL qmmm_env_qm_release(qmmm_env_qm, error)
       CALL qmmm_env_mm_release(qmmm_env_mm, error)
       DEALLOCATE(qm_atom_type,stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    END IF
    CALL timestop(0.0_dp,handle)

  END SUBROUTINE qmmm_create_force_env
!***************************************************************************

!!****f* qmmm_init/assign_mm_charges_and_radius [1.0] *
!!
!!   NAME
!!     assign_mm_charges_and_radius
!!
!!   FUNCTION
!!     Assigns charges and radius to evaluate the MM electrostatic potential
!!
!!   NOTES
!!     
!!
!!   INPUTS
!!     - subsys: the subsystem containing the MM charges
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     06.2004 created [tlaino]
!!
!!*** **********************************************************************
 SUBROUTINE assign_mm_charges_and_radius(subsys, mm_atom_chrg,mm_el_pot_radius,mm_el_pot_radius_corr,&
                                         mm_atom_index, mm_link_atoms, mm_link_scale_factor,&
                                         error)
   IMPLICIT NONE
   ! Arguments
   TYPE(cp_subsystem_p_type), &
        DIMENSION(:), POINTER                  :: subsys
   REAL(dp), DIMENSION(:), POINTER             :: mm_atom_chrg, mm_el_pot_radius,mm_link_scale_factor
   REAL(dp), DIMENSION(:), POINTER             :: mm_el_pot_radius_corr
   INTEGER, DIMENSION(:), POINTER              :: mm_atom_index, mm_link_atoms
   TYPE(cp_error_type), INTENT(inout), &
        OPTIONAL                               :: error
   ! Local Variables
   TYPE(atomic_kind_type), POINTER             :: my_kind
   TYPE(elp_potential_type), POINTER           :: my_potential
   TYPE(particle_list_type), POINTER           :: particles
   TYPE(particle_type), DIMENSION(:), &
        POINTER                                :: particle_set
   INTEGER  :: I, IndMM, ilink, iw
   REAL(dp) :: qi, ri, rc

   NULLIFY(particle_set, my_kind)
   CALL cp_subsys_get(subsys=subsys(1)%subsys,&
                      particles=particles)
   particle_set => particles%els

   DO I = 1, SIZE(mm_atom_index)
      IndMM = mm_atom_index(I)
      my_kind => particle_set(IndMM)%atomic_kind
      CALL get_atomic_kind(atomic_kind=my_kind, elp_potential=my_potential)
      CALL get_potential  (potential=my_potential,& 
                           qeff=qi,&
                           core_charge_radius=ri,&
                           ccore_charge_radius=rc)
      mm_atom_chrg(I)          = qi
      mm_el_pot_radius(I)      = ri
      mm_el_pot_radius_corr(I) = rc
   END DO

   IF (ASSOCIATED(mm_link_atoms)) THEN
      DO ilink = 1, SIZE(mm_link_atoms)
         DO i = 1, SIZE(mm_atom_index)
            IF (mm_atom_index(i) == mm_link_atoms(ilink)) EXIT
         END DO
         IndMM = mm_atom_index(I)
         mm_atom_chrg(i) = mm_atom_chrg(i) * mm_link_scale_factor(ilink)
      END DO
   END IF

 END SUBROUTINE assign_mm_charges_and_radius

!!****f* qmmm_init/print_qmmm_charges [1.0] *
!!
!!   NAME
!!     print_qmmm_charges
!!
!!   FUNCTION
!!     Print info on charges generating the qmmm potential..
!!
!!   NOTES
!!     
!!
!!   INPUTS
!!     - subsys: the subsystem containing the MM charges
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     01.2005 created [tlaino]
!!
!!*** **********************************************************************
 SUBROUTINE print_qmmm_charges(mm_atom_index, mm_atom_chrg, mm_el_pot_radius, mm_el_pot_radius_corr,&
                               added_charges, globenv, qmmm_section, nocompatibility, error)
   IMPLICIT NONE
   ! Arguments
   REAL(dp), DIMENSION(:), POINTER             :: mm_atom_chrg, mm_el_pot_radius
   REAL(dp), DIMENSION(:), POINTER             :: mm_el_pot_radius_corr
   INTEGER, DIMENSION(:), POINTER              :: mm_atom_index
   TYPE(section_vals_type), POINTER            :: qmmm_section
   TYPE(cp_error_type), INTENT(inout), &
        OPTIONAL                               :: error
    TYPE(global_environment_type), POINTER     :: globenv
    TYPE(add_set_type), POINTER                :: added_charges
    LOGICAL, INTENT(IN)                        :: nocompatibility
    ! Local Variables
    INTEGER                                    :: IndMM, I, iw, ind1, ind2
    REAL(KIND=dp)                              :: qi, ri, qtot, rc
    TYPE(cp_logger_type), POINTER :: logger
    TYPE(section_vals_type), POINTER :: print_key

    qtot = 0.0_dp
    iw = globenv%scr
    logger => cp_error_get_logger(error)
    print_key => section_vals_get_subs_vals(qmmm_section,"print%qmmm_charges",error=error)
    IF (BTEST(cp_print_key_should_output(logger%iter_info,print_key,error=error),cp_p_file).AND.&
         logger%para_env%mepos==logger%para_env%source) THEN 
       WRITE(iw,FMT="(/,T2,A)") REPEAT("-",73)
       WRITE(iw,FMT='(/5X,A)')"MM    POINT CHARGES GENERATING THE QM/MM ELECTROSTATIC POTENTIAL"
       WRITE(iw,FMT="(/,T2,A)") REPEAT("-",73)
       DO I = 1, SIZE(mm_atom_index)
          IndMM = mm_atom_index(I)
          qi = mm_atom_chrg(I)    
          qtot = qtot + qi 
          ri = mm_el_pot_radius(I) 
          rc = mm_el_pot_radius_corr(I) 
          IF (nocompatibility) THEN
             WRITE(iw,'(7X,A9,I5,A8,F12.6,A8,F12.6)')' MM ATOM:',IndMM,' RADIUS:',ri,&
                  ' CHARGE:',qi
          ELSE
             WRITE(iw,'(7X,A9,I5,A8,F12.6,A8,F12.6,A,F12.6)')' MM ATOM:',IndMM,' RADIUS:',ri,&
                  ' CHARGE:',qi,'CORR. RADIUS',rc
          END IF
       END DO
       IF (added_charges%num_mm_atoms /= 0) THEN
          WRITE(iw,FMT="(/,T2,A)") REPEAT("-",73)
          WRITE(iw,'(/5X,A)')"ADDED POINT CHARGES GENERATING THE QM/MM ELECTROSTATIC POTENTIAL"
          WRITE(iw,FMT="(/,T2,A)") REPEAT("-",73)
          DO I = 1, SIZE(added_charges%mm_atom_index)
             IndMM = added_charges%mm_atom_index(I)
             qi = added_charges%mm_atom_chrg(I)     
             qtot = qtot + qi
             ri = added_charges%mm_el_pot_radius(I) 
             ind1 = added_charges%add_env(I)%Index1
             ind2 = added_charges%add_env(I)%Index2
             IF (nocompatibility) THEN
                WRITE(iw,'(7X,A9,I5,A8,F12.6,A8,F12.6,I5,I5)')'MM POINT:',IndMM,' RADIUS:',ri,&
                     ' CHARGE:',qi,ind1,ind2
             ELSE
                WRITE(iw,'(7X,A9,I5,A8,F12.6,A8,F12.6,I5,I5,A,F12.6)')'MM POINT:',IndMM,' RADIUS:',ri,&
                     ' CHARGE:',qi,ind1,ind2,' CORR. RADIUS',rc
             END IF
          END DO

       END IF
       WRITE(iw,FMT="(/,T2,A)") REPEAT("-",73)
       WRITE(iw,'(/,T36,A,F12.6)')' TOTAL CHARGE:',qtot
       WRITE(iw,FMT="(/,T2,A,/)") REPEAT("-",73)
    END IF
  END SUBROUTINE print_qmmm_charges

!!****f* qmmm_init/print_qmmm_links [1.0] *
!!
!!   NAME
!!     print_qmmm_links
!!
!!   FUNCTION
!!     Print info on qm/mm links
!!
!!   NOTES
!!     
!!
!!   INPUTS
!!     - subsys: the subsystem containing the MM charges
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     01.2005 created [tlaino]
!!
!!*** **********************************************************************
 SUBROUTINE print_qmmm_links(globenv, qmmm_section, qmmm_links, error)
   IMPLICIT NONE
   ! Arguments
   TYPE(section_vals_type), POINTER            :: qmmm_section
   TYPE(cp_error_type), INTENT(inout), &
        OPTIONAL                               :: error
    TYPE(global_environment_type), POINTER     :: globenv
    TYPE(qmmm_links_type), POINTER             :: qmmm_links
    ! Local Variables
    INTEGER                                    :: mm_index, qm_index, iw, i
    REAL(KIND=dp)                              :: alpha
    TYPE(cp_logger_type), POINTER :: logger
    TYPE(section_vals_type), POINTER :: print_key

    iw = globenv%scr
    logger => cp_error_get_logger(error)
    print_key => section_vals_get_subs_vals(qmmm_section,"print%qmmm_link_info",error=error)
    IF (BTEST(cp_print_key_should_output(logger%iter_info,print_key,error=error),cp_p_file).AND.&
         logger%para_env%mepos==logger%para_env%source) THEN 
       IF (ASSOCIATED(qmmm_links)) THEN
          WRITE(iw,FMT="(/,T2, A)") REPEAT("-",73)
          WRITE(iw,FMT="(/,T31,A)") " QM/MM LINKS "
          WRITE(iw,FMT="(/,T2, A)") REPEAT("-",73)
          IF (ASSOCIATED(qmmm_links%imomm)) THEN
             WRITE(iw,FMT="(/,T31,A)") " IMOMM TYPE LINK "
             DO I = 1, SIZE(qmmm_links%imomm)
                qm_index = qmmm_links%imomm(I)%link%qm_index
                mm_index = qmmm_links%imomm(I)%link%mm_index
                alpha    = qmmm_links%imomm(I)%link%alpha
                WRITE(iw,FMT="(T2,A,T20,A9,I6,1X,A9,I6,T55,A6,F12.6)")"TYPE: IMOMM",&
                     "QM INDEX:",qm_index,"MM INDEX:",mm_index,"ALPHA:",alpha
             END DO
          END IF
          IF (ASSOCIATED(qmmm_links%pseudo)) THEN
             WRITE(iw,FMT="(/,T31,A)") " PSEUDO TYPE LINK "
             DO I = 1, SIZE(qmmm_links%pseudo)
                qm_index = qmmm_links%pseudo(I)%link%qm_index
                mm_index = qmmm_links%pseudo(I)%link%mm_index
                WRITE(iw,FMT="(T2,A,T20,A9,I6,1X,A9,I6,T55)")"TYPE: PSEUDO",&
                     "QM INDEX:",qm_index,"MM INDEX:",mm_index
             END DO
          END IF
          WRITE(iw,FMT="(/,T2,A,/)") REPEAT("-",73)
       ELSE
          WRITE(iw,FMT="(/,T2, A)") REPEAT("-",73)
          WRITE(iw,FMT="(/,T26,A)") " NO QM/MM LINKS DETECTED"
          WRITE(iw,FMT="(/,T2, A)") REPEAT("-",73)          
       END IF
    END IF
  END SUBROUTINE print_qmmm_links

  !!****f* qmmm_init/qmmm_init_gaussian_type [1.0] *
  !!
  !!   NAME
  !!     qmmm_init_gaussian_type
  !!
  !!   FUNCTION
  !!
  !!
  !!   NOTES
  !!
  !!
  !!   INPUTS
  !!     - error: variable to control error logging, stopping,... 
  !!       see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     Teodoro Laino
  !!
  !!   MODIFICATION HISTORY
  !!     1.2005 created [tlaino]
  !!
  !!*************************************************************************
 SUBROUTINE qmmm_init_gaussian_type(qmmm_env_qm, qm_cell_small, mm_cell, globenv,&
                                    mm_atom_chrg, qs_env, added_charges, print_section, &
                                    qmmm_section, error )
   IMPLICIT NONE
   ! Arguments
    TYPE(cell_type), POINTER                 :: mm_cell, qm_cell_small
    TYPE(qmmm_env_qm_type), POINTER          :: qmmm_env_qm
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(qs_environment_type), POINTER       :: qs_env
    REAL(KIND=dp), POINTER, DIMENSION(:)     :: mm_atom_chrg
    TYPE(add_set_type), POINTER              :: added_charges
    TYPE(section_vals_type), POINTER         :: print_section, qmmm_section
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                               :: error
   ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_init_gaussian_type',&
         routineP = moduleN//':'//routineN
    REAL(KIND=dp), DIMENSION(:), POINTER :: maxradius, maxradius2
    REAL(KIND=dp) :: maxchrg
    INTEGER :: i, stat
    LOGICAL :: failure

    failure=.FALSE.
    stat=0
    NULLIFY(maxradius, maxradius2)

    maxchrg = MAXVAL(ABS(mm_atom_chrg(:)))
    IF (qmmm_env_qm%add_mm_charges) maxchrg=MAX(maxchrg,MAXVAL(ABS(added_charges%mm_atom_chrg(:))))
    CALL qmmm_gaussian_initialize(qmmm_gaussian_fns=qmmm_env_qm%pgfs,&
                                  qm_cell_small=qm_cell_small,&
                                  mm_cell=mm_cell,&
                                  globenv=globenv,&
                                  pw_env=qs_env%pw_env,&
                                  mm_el_pot_radius=qmmm_env_qm%mm_el_pot_radius,&
                                  mm_el_pot_radius_corr=qmmm_env_qm%mm_el_pot_radius_corr,&
                                  qmmm_coupl_type=qmmm_env_qm%qmmm_coupl_type,&
                                  eps_mm_rspace=qmmm_env_qm%eps_mm_rspace,&
                                  maxradius=maxradius,&
                                  maxchrg=maxchrg,&
                                  compatibility=qmmm_env_qm%compatibility,&
                                  print_section=print_section,&
                                  qmmm_section=qmmm_section,&
                                  error=error)

    IF (qmmm_env_qm%move_mm_charges.OR.qmmm_env_qm%add_mm_charges) THEN
       CALL qmmm_gaussian_initialize(qmmm_gaussian_fns=added_charges%pgfs,&
                                     qm_cell_small=qm_cell_small,&
                                     mm_cell=mm_cell,&
                                     globenv=globenv,&
                                     pw_env=qs_env%pw_env,&
                                     mm_el_pot_radius=added_charges%mm_el_pot_radius,&
                                     mm_el_pot_radius_corr=added_charges%mm_el_pot_radius_corr,&
                                     qmmm_coupl_type=qmmm_env_qm%qmmm_coupl_type,&
                                     eps_mm_rspace=qmmm_env_qm%eps_mm_rspace,&
                                     maxradius=maxradius2,&
                                     maxchrg=maxchrg,&
                                     compatibility=qmmm_env_qm%compatibility,&
                                     print_section=print_section,&
                                     qmmm_section=qmmm_section,&
                                     error=error)

       IF (qmmm_env_qm%qmmm_coupl_type == do_qmmm_gauss) THEN
          DO i = 1, SIZE(maxradius)
             maxradius(i) = MAX( maxradius(i), maxradius2(i))
          END DO
       END IF
       IF (ASSOCIATED(maxradius2)) DEALLOCATE(maxradius2, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

    qmmm_env_qm%maxradius => maxradius

 END SUBROUTINE qmmm_init_gaussian_type

  !!****f* qmmm_init/qmmm_init_potential [1.0] *
  !!
  !!   NAME
  !!     qmmm_init_potential
  !!
  !!   FUNCTION
  !!
  !!
  !!   NOTES
  !!
  !!
  !!   INPUTS
  !!     - error: variable to control error logging, stopping,... 
  !!       see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     Teodoro Laino
  !!
  !!   MODIFICATION HISTORY
  !!     1.2005 created [tlaino]
  !!
  !!*************************************************************************
 SUBROUTINE qmmm_init_potential(qmmm_env_qm, qm_cell_small, mm_cell, globenv, qs_env,&
                                added_charges, print_section, error)
   IMPLICIT NONE
   ! Arguments
   TYPE(cell_type), POINTER                 :: mm_cell, qm_cell_small
   TYPE(qmmm_env_qm_type), POINTER          :: qmmm_env_qm
   TYPE(global_environment_type), POINTER   :: globenv
   TYPE(qs_environment_type), POINTER       :: qs_env
   TYPE(add_set_type), POINTER              :: added_charges
   TYPE(section_vals_type), POINTER         :: print_section
   TYPE(cp_error_type), INTENT(inout), &
        OPTIONAL                               :: error   
   ! Local Variables
   
   CALL qmmm_potential_init(qmmm_coupl_type=qmmm_env_qm%qmmm_coupl_type,&
                            mm_el_pot_radius=qmmm_env_qm%mm_el_pot_radius,&
                            potentials=qmmm_env_qm%potentials,&
                            pgfs=qmmm_env_qm%pgfs,&
                            qm_cell_small=qm_cell_small,&
                            mm_cell=mm_cell,&
                            globenv=globenv,&
                            pw_env=qs_env%pw_env,&
                            compatibility=qmmm_env_qm%compatibility,&
                            print_section=print_section,&
                            error=error)

   IF (qmmm_env_qm%move_mm_charges.OR.qmmm_env_qm%add_mm_charges) THEN

      CALL qmmm_potential_init(qmmm_coupl_type=qmmm_env_qm%qmmm_coupl_type,&
                               mm_el_pot_radius=added_charges%mm_el_pot_radius,&
                               potentials=added_charges%potentials,&
                               pgfs=added_charges%pgfs,&
                               qm_cell_small=qm_cell_small,&
                               mm_cell=mm_cell,&
                               globenv=globenv,&
                               pw_env=qs_env%pw_env,&
                               compatibility=qmmm_env_qm%compatibility,&
                               print_section=print_section,&
                               error=error)      
   END IF

 END SUBROUTINE qmmm_init_potential

  !!****f* qmmm_init/qmmm_init_periodic_potential [1.0] *
  !!
  !!   NAME
  !!     qmmm_init_periodic_potential
  !!
  !!   FUNCTION
  !!
  !!
  !!   NOTES
  !!
  !!
  !!   INPUTS
  !!     - error: variable to control error logging, stopping,... 
  !!       see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     Teodoro Laino
  !!
  !!   MODIFICATION HISTORY
  !!     7.2005 created [tlaino]
  !!
  !!*************************************************************************
 SUBROUTINE qmmm_init_periodic_potential(qmmm_env_qm, qm_cell_small, mm_cell, globenv, qs_env,&
                                         added_charges, qmmm_periodic, print_section, mm_atom_chrg, error)
   IMPLICIT NONE
   ! Arguments
   TYPE(cell_type), POINTER                 :: mm_cell, qm_cell_small
   TYPE(qmmm_env_qm_type), POINTER          :: qmmm_env_qm
   TYPE(global_environment_type), POINTER   :: globenv
   TYPE(qs_environment_type), POINTER       :: qs_env
   TYPE(add_set_type), POINTER              :: added_charges
   REAL(KIND=dp), POINTER, DIMENSION(:)     :: mm_atom_chrg
   TYPE(section_vals_type), POINTER         :: print_section, qmmm_periodic
   TYPE(cp_error_type), INTENT(inout), &
        OPTIONAL                               :: error   
   ! Local Variables
   REAL(KIND=dp) :: maxchrg

   IF (.NOT.qmmm_env_qm%periodic) RETURN

   maxchrg = MAXVAL(ABS(mm_atom_chrg(:)))
   IF (qmmm_env_qm%add_mm_charges) maxchrg=MAX(maxchrg,MAXVAL(ABS(added_charges%mm_atom_chrg(:))))

   CALL qmmm_per_potential_init(qmmm_coupl_type=qmmm_env_qm%qmmm_coupl_type,&
                                per_potentials=qmmm_env_qm%per_potentials,&
                                potentials=qmmm_env_qm%potentials,&
                                pgfs=qmmm_env_qm%pgfs,&
                                qm_cell_small=qm_cell_small,&
                                mm_cell=mm_cell,&
                                globenv=globenv,&
                                compatibility=qmmm_env_qm%compatibility,&
                                qmmm_periodic=qmmm_periodic,&
                                print_section=print_section,&
                                eps_mm_rspace=qmmm_env_qm%eps_mm_rspace,&
                                maxchrg=maxchrg,&
                                ncp =qmmm_env_qm%aug_pools(SIZE(qmmm_env_qm%aug_pools))%pool%pw_grid%npts,&
                                ncpl=qmmm_env_qm%aug_pools(SIZE(qmmm_env_qm%aug_pools))%pool%pw_grid%npts_local,&
                                error=error)

   IF (qmmm_env_qm%move_mm_charges.OR.qmmm_env_qm%add_mm_charges) THEN

      CALL qmmm_per_potential_init(qmmm_coupl_type=qmmm_env_qm%qmmm_coupl_type,&
                                   per_potentials=added_charges%per_potentials,&
                                   potentials=added_charges%potentials,&
                                   pgfs=added_charges%pgfs,&
                                   qm_cell_small=qm_cell_small,&
                                   mm_cell=mm_cell,&
                                   globenv=globenv,&
                                   compatibility=qmmm_env_qm%compatibility,&
                                   qmmm_periodic=qmmm_periodic,&
                                   print_section=print_section,&
                                   eps_mm_rspace=qmmm_env_qm%eps_mm_rspace,&
                                   maxchrg=maxchrg,&
                                   ncp =qmmm_env_qm%aug_pools(SIZE(qmmm_env_qm%aug_pools))%pool%pw_grid%npts,&
                                   ncpl=qmmm_env_qm%aug_pools(SIZE(qmmm_env_qm%aug_pools))%pool%pw_grid%npts_local,&
                                   error=error)      
   END IF

 END SUBROUTINE qmmm_init_periodic_potential

  !!****f* qmmm_init/setup_qmmm_vars_qm [1.0] *
  !!
  !!   NAME
  !!     setup_qmmm_vars_qm
  !!
  !!   FUNCTION
  !!
  !!
  !!   NOTES
  !!
  !!
  !!   INPUTS
  !!     - error: variable to control error logging, stopping,... 
  !!       see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     Teodoro Laino
  !!
  !!   MODIFICATION HISTORY
  !!     11.2004 created [tlaino]
  !!
  !!*************************************************************************
  SUBROUTINE setup_qmmm_vars_qm( qmmm_section, qmmm_env, subsys_mm, qm_atom_type,&
                              qm_atom_index, mm_atom_index, qm_subsys_index,&
                              qm_cell_small, qmmm_coupl_type, eps_mm_rspace, qmmm_link,&
                              error)
    IMPLICIT NONE
    ! Arguments
    TYPE(section_vals_type), POINTER         :: qmmm_section
    TYPE(qmmm_env_qm_type), POINTER             :: qmmm_env
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys_mm
    CHARACTER(len=default_string_length), &
      DIMENSION(:), POINTER                  :: qm_atom_type
    INTEGER, DIMENSION(:), POINTER           :: qm_atom_index, mm_atom_index, qm_subsys_index
    TYPE(cell_type), POINTER                 :: qm_cell_small
    INTEGER, INTENT(OUT)                     :: qmmm_coupl_type
    REAL(KIND=dp), INTENT(OUT)               :: eps_mm_rspace
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    LOGICAL, INTENT(OUT)                     :: qmmm_link
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN = 'setup_qmmm_vars_qm',&
         routineP = moduleN//':'//routineN
    INTEGER :: i, j, icount, size_mm_system, stat, ikind, ikindr, nkind, n_rep, iatm
    LOGICAL :: is_mm
    INTEGER                                  :: qm_subsys_type, n_rep_val
    LOGICAL                                  :: failure
    CHARACTER(len=default_string_length)     :: qm_atom_kind, mm_atom_kind
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(atomic_kind_type),POINTER           :: atomic_kind
    CHARACTER(LEN=default_string_length)     :: atmname
    REAL(KIND=dp)                            :: tmp_radius, tmp_radius_c
    TYPE(elp_potential_type), POINTER        :: elp_potential
    TYPE(section_vals_type), POINTER         :: mm_kinds
    LOGICAL :: explicit, is_qm
    REAL(KIND=dp), DIMENSION(:), POINTER     :: cell_par
    INTEGER, DIMENSION(:), POINTER           :: mm_link_atoms

    NULLIFY(mm_link_atoms)
    qmmm_link = .FALSE.
    failure   = .FALSE.
    CALL section_vals_get(qmmm_section,explicit=explicit,error=error)
    IF (explicit) THEN
       CALL section_vals_val_get(qmmm_section,"E_COUPL",i_val=qmmm_coupl_type,error=error)
       CALL section_vals_val_get(qmmm_section,"EPS_MM_RSPACE",r_val=eps_mm_rspace,error=error)
       !
       ! Initialization of arrays and core_charge_radius...
       !
       tmp_radius=0.0_dp
       CALL cp_subsys_get(subsys=subsys_mm(1)%subsys,&
                          atomic_kinds=atomic_kinds,&
                          error=error)
       DO Ikind = 1, SIZE(atomic_kinds%els)
          atomic_kind => atomic_kinds%els(Ikind)
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
                               elp_potential=elp_potential)
          CALL set_potential  (potential=elp_potential,&
                               core_charge_radius=tmp_radius,&
                               ccore_charge_radius=tmp_radius)
          CALL set_atomic_kind(atomic_kind=atomic_kind,&
                               elp_potential=elp_potential)
       END DO
       CALL setup_qm_atom_list(qmmm_section=qmmm_section,&
                               qm_atom_index=qm_atom_index,&
                               qm_atom_type=qm_atom_type,&
                               mm_link_atoms=mm_link_atoms,&
                               qm_subsys_index=qm_subsys_index,&
                               qmmm_link=qmmm_link,&
                               error=error)
       !
       ! MM_KINDS
       !
       mm_kinds => section_vals_get_subs_vals(qmmm_section,"MM_KIND",error=error)
       CALL section_vals_get(mm_kinds,explicit=explicit,n_repetition=nkind,error=error)
       !
       ! Default
       !
       tmp_radius = RADIUS_QMMM_DEFAULT * qm_cell_small%unit_of_length
       Set_Radius_Pot_0: DO IkindR = 1, SIZE(atomic_kinds%els)
          atomic_kind => atomic_kinds%els(IkindR)
          CALL get_atomic_kind(atomic_kind=atomic_kind,name=atmname)
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
               elp_potential=elp_potential)
          CALL set_potential  (potential=elp_potential,&
               core_charge_radius=tmp_radius,&
               ccore_charge_radius=tmp_radius)
          CALL set_atomic_kind(atomic_kind=atomic_kind,&
               elp_potential=elp_potential)
       END DO Set_Radius_Pot_0
       !
       ! If present overwrite the kind specified in input file...
       !
       IF (explicit) THEN
          DO ikind=1,nkind
             CALL section_vals_val_get(mm_kinds,"_SECTION_PARAMETERS_",i_rep_section=ikind,&
                  c_val=mm_atom_kind,error=error)
             CALL section_vals_val_get(mm_kinds,"RADIUS",i_rep_section=ikind,r_val=tmp_radius,error=error)
             tmp_radius_c = tmp_radius
             CALL section_vals_val_get(mm_kinds,"CORR_RADIUS",i_rep_section=ikind,n_rep_val=n_rep_val,error=error)
             if (n_rep_val == 1) &
                  CALL section_vals_val_get(mm_kinds,"CORR_RADIUS",i_rep_section=ikind,r_val=tmp_radius_c,error=error)
             tmp_radius   = tmp_radius * qm_cell_small%unit_of_length
             tmp_radius_c = tmp_radius_c * qm_cell_small%unit_of_length
             Set_Radius_Pot_1: DO IkindR = 1, SIZE(atomic_kinds%els)
                atomic_kind => atomic_kinds%els(IkindR)
                CALL get_atomic_kind(atomic_kind=atomic_kind,name=atmname)
                is_qm = qmmm_ff_precond_only_qm(atmname)
                IF (TRIM(mm_atom_kind)==atmname) THEN
                   CALL get_atomic_kind(atomic_kind=atomic_kind,&
                        elp_potential=elp_potential)
                   CALL set_potential  (potential=elp_potential,&
                        core_charge_radius=tmp_radius,&
                        ccore_charge_radius=tmp_radius_c)
                   CALL set_atomic_kind(atomic_kind=atomic_kind,&
                        elp_potential=elp_potential)
                END IF
             END DO Set_Radius_Pot_1
          END DO
       END IF
       !
       ! QM Cell Info
       !
       CALL section_vals_val_get(qmmm_section,"A_QM",n_rep_val=n_rep,error=error)
       IF (n_rep /= 0) THEN
           CALL section_vals_val_get(qmmm_section,"A_QM",r_vals=cell_par,error=error)
          qm_cell_small%hmat(:,1) = cell_par(:)
       END IF
       CALL section_vals_val_get(qmmm_section,"B_QM",n_rep_val=n_rep,error=error)
       IF (n_rep /= 0) THEN
          CALL section_vals_val_get(qmmm_section,"B_QM",r_vals=cell_par,error=error)
          qm_cell_small%hmat(:,2) = cell_par(:)
       END IF
       CALL section_vals_val_get(qmmm_section,"C_QM",n_rep_val=n_rep,error=error)
       IF (n_rep /= 0) THEN
          CALL section_vals_val_get(qmmm_section,"C_QM",r_vals=cell_par,error=error)
          qm_cell_small%hmat(:,3) = cell_par(:)
       END IF
       CALL section_vals_val_get(qmmm_section,"ABC_QM",n_rep_val=n_rep,error=error)
       IF (n_rep /= 0) THEN
          CALL section_vals_val_get(qmmm_section,"ABC_QM",r_vals=cell_par,error=error)
          qm_cell_small%hmat(1,1) = cell_par(1)
          qm_cell_small%hmat(2,2) = cell_par(2)
          qm_cell_small%hmat(3,3) = cell_par(3)
       END IF
       qm_cell_small%hmat=qm_cell_small%hmat(:,:)*qm_cell_small%unit_of_length
    ELSE
     CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
          "QMMM section not present in input file!"//&
          CPSourceFileRef,&
          error,failure)
    ENDIF
    !
    ! Initialize cell parameters..
    !
    CALL init_qmmm_cell(qm_cell_small)
    !
    ! Build MM atoms list
    !
    size_mm_system = SIZE(subsys_mm(1)%subsys%particles%els) - SIZE(qm_atom_index)
    IF (qmmm_link.AND.ASSOCIATED(mm_link_atoms)) size_mm_system = size_mm_system + SIZE(mm_link_atoms)
    ALLOCATE(mm_atom_index(size_mm_system))
    icount = 0
    DO i = 1, SIZE(subsys_mm(1)%subsys%particles%els)
       is_mm = .TRUE.
       IF ( ANY(qm_atom_index == i)) is_mm =.FALSE.
       IF (ASSOCIATED(mm_link_atoms)) THEN 
          IF (ANY( mm_link_atoms == i ).AND.qmmm_link) is_mm =.TRUE.
       END IF
       IF (is_mm) THEN
          icount = icount + 1
          mm_atom_index(icount) = i
       END IF
    END DO
    CPPostcondition(icount==size_mm_system,cp_failure_level,routineP,error,failure)
    IF (ASSOCIATED(mm_link_atoms)) THEN 
       DEALLOCATE(mm_link_atoms, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
  END SUBROUTINE setup_qmmm_vars_qm


  !!****f* qmmm_init/setup_qmmm_vars_mm [1.0] *
  !!
  !!   NAME
  !!     setup_qmmm_vars_mm
  !!
  !!   FUNCTION
  !!
  !!
  !!   NOTES
  !!
  !!
  !!   INPUTS
  !!     - error: variable to control error logging, stopping,... 
  !!       see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     Teodoro Laino
  !!
  !!   MODIFICATION HISTORY
  !!     12.2004 created [tlaino]
  !!
  !!*************************************************************************
  SUBROUTINE setup_qmmm_vars_mm( qmmm_section, qmmm_env, qm_atom_index,  &
                                 mm_link_atoms, mm_link_scale_factor, &
                                 fist_scale_charge_link, qmmm_coupl_type,&
                                 qmmm_link, error)
    IMPLICIT NONE
    ! Arguments
    TYPE(section_vals_type), POINTER         :: qmmm_section
    TYPE(qmmm_env_mm_type), POINTER          :: qmmm_env
    INTEGER, DIMENSION(:), POINTER           :: qm_atom_index, mm_link_atoms
    INTEGER, INTENT(OUT)                     :: qmmm_coupl_type
    LOGICAL, INTENT(OUT)                     :: qmmm_link
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mm_link_scale_factor,fist_scale_charge_link
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                               :: error
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN = 'setup_qmmm_vars_mm',&
         routineP = moduleN//':'//routineN
    LOGICAL                                  :: failure
    LOGICAL :: explicit
    TYPE(section_vals_type), POINTER         :: qmmm_ff_section
    
    NULLIFY(qmmm_ff_section)
    qmmm_link = .FALSE.
    CALL section_vals_get(qmmm_section,explicit=explicit,error=error)
    IF (explicit) THEN
       CALL section_vals_val_get(qmmm_section,"E_COUPL",i_val=qmmm_coupl_type,error=error)
       CALL setup_qm_atom_list(qmmm_section,qm_atom_index=qm_atom_index,qmmm_link=qmmm_link,&
            mm_link_atoms=mm_link_atoms, mm_link_scale_factor = mm_link_scale_factor,&
            fist_scale_charge_link=fist_scale_charge_link, error=error)
       !
       ! Do we want to use a different FF for the non-bonded QM/MM interactions?
       !
       qmmm_ff_section => section_vals_get_subs_vals(qmmm_section,"FORCEFIELD",error=error)
       CALL section_vals_get(qmmm_ff_section,explicit=explicit,error=error)
       IF (explicit) THEN
          CALL read_qmmm_ff_section(qmmm_ff_section, qmmm_env%inp_info, error=error)
          qmmm_env%use_qmmm_ff = .TRUE.
       END IF
    END IF
  END SUBROUTINE setup_qmmm_vars_mm
     
  !!****f* qmmm_init/read_qmmm_ff_section [1.0] *
  !!
  !!   NAME
  !!     read_qmmm_ff_section
  !!
  !!   FUNCTION
  !!     reads information regarding the forcefield specific for the QM/MM
  !!     interactions
  !!
  !!   NOTES
  !!
  !!
  !!   INPUTS
  !!     - error: variable to control error logging, stopping,... 
  !!       see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     Teodoro Laino
  !!
  !!   MODIFICATION HISTORY
  !!     12.2004 created [tlaino]
  !!
  !!*************************************************************************
  SUBROUTINE read_qmmm_ff_section(qmmm_ff_section, inp_info, error)
    IMPLICIT NONE
    ! Arguments
    TYPE(section_vals_type), POINTER                        :: qmmm_ff_section
    TYPE(input_info_type),POINTER                           :: inp_info
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                               :: error    
    ! Local Variables
    INTEGER :: n_lj, n_wl, n_gd, np
    TYPE(section_vals_type), POINTER                        :: lj_section,&
                                                               wl_section,&
                                                               gd_section
    !
    ! NONBONDED
    !
    lj_section => section_vals_get_subs_vals(qmmm_ff_section,"NONBONDED%LENNARD-JONES",error=error)
    wl_section => section_vals_get_subs_vals(qmmm_ff_section,"NONBONDED%WILLIAMS",error=error)
    gd_section => section_vals_get_subs_vals(qmmm_ff_section,"NONBONDED%GOODWIN",error=error)
    CALL section_vals_get(lj_section,n_repetition=n_lj,error=error)
    np = n_lj 
    IF (n_lj /= 0) THEN 
       CALL pair_potential_reallocate(inp_info%nonbonded,1,np,lj=.TRUE.,error=error)
       CALL read_lj_section(inp_info%nonbonded,lj_section, start=0,      error=error)
    END IF
    CALL section_vals_get(wl_section,n_repetition=n_wl,error=error)
    np = n_lj + n_wl 
    IF (n_wl /= 0) THEN 
       CALL pair_potential_reallocate(inp_info%nonbonded,1,np,williams=.TRUE.,error=error)
       CALL read_wl_section(inp_info%nonbonded,wl_section, start=n_lj, error=error)
    END IF
    CALL section_vals_get(gd_section,n_repetition=n_gd,error=error)
    np = n_lj + n_wl + n_gd
    IF (n_gd /= 0) THEN 
       CALL pair_potential_reallocate(inp_info%nonbonded,1,np,goodwin=.TRUE.,error=error)
       CALL read_gd_section(inp_info%nonbonded,gd_section, start=n_lj+n_wl, error=error)
    END IF
    !
    ! NONBONDED14
    !
    lj_section => section_vals_get_subs_vals(qmmm_ff_section,"NONBONDED14%LENNARD-JONES",error=error)
    wl_section => section_vals_get_subs_vals(qmmm_ff_section,"NONBONDED14%WILLIAMS",error=error)
    gd_section => section_vals_get_subs_vals(qmmm_ff_section,"NONBONDED14%GOODWIN",error=error)
    CALL section_vals_get(lj_section,n_repetition=n_lj,error=error)
    np = n_lj
    IF (n_lj /= 0) THEN 
       CALL pair_potential_reallocate(inp_info%nonbonded14,1,np,lj=.TRUE.,error=error)
       CALL read_lj_section(inp_info%nonbonded14,lj_section, start=0,      error=error)
    END IF
    CALL section_vals_get(wl_section,n_repetition=n_wl,error=error)
    np = n_lj + n_wl 
    IF (n_wl /= 0) THEN 
       CALL pair_potential_reallocate(inp_info%nonbonded14,1,np,williams=.TRUE.,error=error)
       CALL read_wl_section(inp_info%nonbonded14,wl_section, start=n_lj, error=error)
    END IF
    CALL section_vals_get(gd_section,n_repetition=n_gd,error=error)
    np = n_lj + n_wl + n_gd
    IF (n_gd /= 0) THEN 
       CALL pair_potential_reallocate(inp_info%nonbonded14,1,np,goodwin=.TRUE.,error=error)
       CALL read_gd_section(inp_info%nonbonded14,gd_section, start=n_lj+n_wl, error=error)
    END IF
  END SUBROUTINE read_qmmm_ff_section

  !!****f* qmmm_init/setup_qm_atom_list [1.0] *
  !!
  !!   NAME
  !!     setup_qm_atom_list
  !!
  !!   FUNCTION
  !!
  !!
  !!   NOTES
  !!
  !!
  !!   INPUTS
  !!     - error: variable to control error logging, stopping,... 
  !!       see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     Teodoro Laino
  !!
  !!   MODIFICATION HISTORY
  !!     12.2004 created [tlaino]
  !!
  !!*************************************************************************
  SUBROUTINE setup_qm_atom_list(qmmm_section,qm_atom_index,qm_atom_type,&
       qm_subsys_index, mm_link_atoms, mm_link_scale_factor, qmmm_link, &
       fist_scale_charge_link, error)
    IMPLICIT NONE
    ! Arguments
    TYPE(section_vals_type), POINTER         :: qmmm_section
    CHARACTER(len=default_string_length), &
      DIMENSION(:), POINTER, OPTIONAL        :: qm_atom_type
    INTEGER, DIMENSION(:), POINTER, OPTIONAL :: qm_atom_index, qm_subsys_index
    INTEGER, DIMENSION(:), POINTER, OPTIONAL :: mm_link_atoms
    REAL(KIND=dp), DIMENSION(:), POINTER, OPTIONAL     :: mm_link_scale_factor
    REAL(KIND=dp), DIMENSION(:), POINTER, OPTIONAL     :: fist_scale_charge_link
    LOGICAL, OPTIONAL, INTENT(OUT) :: qmmm_link
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                               :: error
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN = 'setup_qm_atom_list',&
         routineP = moduleN//':'//routineN
    LOGICAL :: explicit, failure
    INTEGER :: num_qm_atom_tot, nkind, ikind, nlinks, stat, mm_index, n_var, k
    INTEGER :: qm_subsys_type, link_type, link_involv_mm
    INTEGER, DIMENSION(:), POINTER           :: mm_indexes
    TYPE(section_vals_type), POINTER         :: qm_kinds, qmmm_links
    CHARACTER(len=default_string_length)     :: qm_atom_kind, qm_link_element
    REAL(KIND=dp) :: scale_f

    failure         = .FALSE.
    num_qm_atom_tot = 0
    link_involv_mm  = 0
    nlinks          = 0
    stat            = 0
    !
    ! QM_KINDS
    !
    qm_kinds => section_vals_get_subs_vals(qmmm_section,"QM_KIND",error=error)
    CALL section_vals_get(qm_kinds,n_repetition=nkind,error=error)
    DO ikind=1,nkind
       CALL section_vals_val_get(qm_kinds,"MM_INDEX",i_rep_section=ikind,n_rep_val=n_var,&
            error=error)       
       DO k = 1, n_var
          CALL section_vals_val_get(qm_kinds,"MM_INDEX",i_rep_section=ikind,i_rep_val=k,&
               i_vals=mm_indexes,error=error)
          num_qm_atom_tot = num_qm_atom_tot + SIZE(mm_indexes)
       END DO             
    END DO
    !
    ! QM/MM LINKS
    !
    qmmm_links => section_vals_get_subs_vals(qmmm_section,"LINK",error=error)
    CALL section_vals_get(qmmm_links,explicit=explicit,error=error)
    IF (explicit) THEN
       qmmm_link = .TRUE.
       CALL section_vals_get(qmmm_links,n_repetition=nlinks,error=error)
       ! Take care of the various link types
       DO ikind = 1, nlinks
          CALL section_vals_val_get(qmmm_links,"LINK_TYPE",i_rep_section=ikind,&
               i_val=link_type,error=error)
          SELECT CASE(link_type)
          CASE(do_qmmm_link_imomm)
             num_qm_atom_tot = num_qm_atom_tot + 1
             link_involv_mm = link_involv_mm + 1
          CASE(do_qmmm_link_pseudo)
             num_qm_atom_tot = num_qm_atom_tot + 1
          CASE(do_qmmm_link_gho)
             ! do nothing for the moment
          CASE DEFAULT
             CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
          END SELECT
       END DO
    END IF
    IF (PRESENT(mm_link_scale_factor).AND.(link_involv_mm /= 0)) &
         ALLOCATE(mm_link_scale_factor(link_involv_mm), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (PRESENT(fist_scale_charge_link).AND.(link_involv_mm /= 0)) & 
         ALLOCATE(fist_scale_charge_link(link_involv_mm), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (PRESENT(mm_link_atoms).AND.(link_involv_mm /= 0)) &
         ALLOCATE(mm_link_atoms(link_involv_mm),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (PRESENT(qm_atom_index))   ALLOCATE(qm_atom_index(num_qm_atom_tot),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (PRESENT(qm_atom_type))    ALLOCATE(qm_atom_type(num_qm_atom_tot),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (PRESENT(qm_subsys_index)) ALLOCATE(qm_subsys_index(num_qm_atom_tot),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (PRESENT(qm_atom_index))   qm_atom_index    =   0
    IF (PRESENT(qm_atom_type))    qm_atom_type     = " "
    IF (PRESENT(qm_subsys_index)) qm_subsys_index  =   0
    num_qm_atom_tot = 1
    DO ikind=1,nkind
       CALL section_vals_val_get(qm_kinds,"MM_INDEX",i_rep_section=ikind,n_rep_val=n_var,&
            error=error)       
       DO k = 1, n_var
          CALL section_vals_val_get(qm_kinds,"MM_INDEX",i_rep_section=ikind,i_rep_val=k,&
               i_vals=mm_indexes,error=error)
          IF (PRESENT(qm_atom_index)) THEN
             qm_atom_index(num_qm_atom_tot:num_qm_atom_tot + SIZE(mm_indexes)-1)   = mm_indexes(:)
          END IF
          IF (PRESENT(qm_atom_type)) THEN
             CALL section_vals_val_get(qm_kinds,"_SECTION_PARAMETERS_",i_rep_section=ikind,&
                  c_val=qm_atom_kind,error=error)
             qm_atom_type (num_qm_atom_tot:num_qm_atom_tot + SIZE(mm_indexes)-1)   = qm_atom_kind
          END IF
          IF (PRESENT(qm_subsys_index)) THEN
             CALL section_vals_val_get(qm_kinds,"SUBSYS",i_rep_section=ikind,i_val=qm_subsys_type,&
                  error=error)
             qm_subsys_index(num_qm_atom_tot:num_qm_atom_tot + SIZE(mm_indexes)-1) = qm_subsys_type
          END IF
          num_qm_atom_tot = num_qm_atom_tot + SIZE(mm_indexes)
       END DO
    END DO
    IF (PRESENT(mm_link_scale_factor).AND.(link_involv_mm /= 0))     mm_link_scale_factor = 0.0_dp
    IF (PRESENT(fist_scale_charge_link).AND.(link_involv_mm /= 0)) fist_scale_charge_link = 0.0_dp
    IF (PRESENT(mm_link_atoms).AND.(link_involv_mm /= 0))            mm_link_atoms        = 0
    IF (explicit) THEN
       DO ikind = 1, nlinks
          IF (PRESENT(qm_subsys_index)) THEN
             CALL section_vals_val_get(qmmm_links,"SUBSYS",i_rep_section=ikind,&
                  i_val=qm_subsys_type,error=error)
             qm_subsys_index(num_qm_atom_tot:num_qm_atom_tot) = qm_subsys_type
          END IF
          IF (PRESENT(qm_atom_type)) THEN
              CALL section_vals_val_get(qmmm_links,"QM_KIND",i_rep_section=ikind,c_val=qm_link_element,&
                  error=error)
             qm_atom_type (num_qm_atom_tot:num_qm_atom_tot)   = TRIM(qm_link_element)//"_LINK"
          END IF
          IF (PRESENT(qm_atom_index)) THEN
             CALL section_vals_val_get(qmmm_links,"MM_INDEX",i_rep_section=ikind,i_val=mm_index,&
                  error=error)
             qm_atom_index(num_qm_atom_tot:num_qm_atom_tot )  = mm_index
             num_qm_atom_tot = num_qm_atom_tot + 1
          END IF
          IF (PRESENT(mm_link_atoms).AND.(link_involv_mm /= 0))  THEN
             CALL section_vals_val_get(qmmm_links,"MM_INDEX",i_rep_section=ikind,i_val=mm_index,&
                  error=error)
             mm_link_atoms (ikind) = mm_index
          END IF
          IF (PRESENT(mm_link_scale_factor).AND.(link_involv_mm /= 0)) THEN
             CALL section_vals_val_get(qmmm_links,"QMMM_SCALE_FACTOR",i_rep_section=ikind,r_val=scale_f,&
                  error=error)
             mm_link_scale_factor(ikind) = scale_f
          END IF
          IF (PRESENT(fist_scale_charge_link).AND.(link_involv_mm /= 0))  THEN
             CALL section_vals_val_get(qmmm_links,"FIST_SCALE_FACTOR",i_rep_section=ikind,r_val=scale_f,&
                  error=error)
             fist_scale_charge_link(ikind) = scale_f
          END IF
       END DO
    END IF
    CPPostcondition(num_qm_atom_tot-1==SIZE(qm_atom_index),cp_failure_level,routineP,error,failure)

  END SUBROUTINE setup_qm_atom_list


  !!****f* qmmm_init/setup_qmmm_links [1.0] *
  !!
  !!   NAME
  !!     setup_qmmm_links
  !!
  !!   FUNCTION
  !!     this routine sets up all variables to treat qmmm links
  !!
  !!   NOTES
  !!
  !!
  !!   INPUTS
  !!     - error: variable to control error logging, stopping,... 
  !!       see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     Teodoro Laino
  !!
  !!   MODIFICATION HISTORY
  !!     12.2004 created [tlaino]
  !!
  !!*************************************************************************
  SUBROUTINE  setup_qmmm_links(qmmm_section, qmmm_links, mm_el_pot_radius, mm_el_pot_radius_corr,&
       mm_atom_index, qm_cell_small, iw, error)
    IMPLICIT NONE
    ! Arguments
    TYPE(section_vals_type), POINTER         :: qmmm_section, qmmm_link_section
    TYPE(qmmm_links_type), POINTER           :: qmmm_links
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mm_el_pot_radius, mm_el_pot_radius_corr
    INTEGER, DIMENSION(:), POINTER           :: mm_atom_index
    TYPE(cell_type), POINTER                 :: qm_cell_small
    INTEGER, INTENT(IN) :: iw
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                               :: error
    ! Local variables
    CHARACTER(len=*), PARAMETER :: routineN = 'setup_qmmm_links',&
         routineP = moduleN//':'//routineN
    INTEGER  :: ikind, nlinks, n_imomm, n_gho, link_type, qm_index, mm_index, stat
    INTEGER  :: n_rep_val, n_pseudo, n_tot
    LOGICAL  :: failure, charge_scale, explicit
    REAL(KIND=dp) :: alpha, scale_factor, my_radius
    INTEGER, POINTER, DIMENSION(:) :: wrk_tmp

    NULLIFY(wrk_tmp)
    failure = .FALSE.
    n_imomm = 0
    n_gho   = 0
    n_pseudo= 0
    qmmm_link_section => section_vals_get_subs_vals(qmmm_section,"LINK",error=error)
    CALL section_vals_get(qmmm_link_section,n_repetition=nlinks,error=error)    
    CPPostcondition(nlinks /= 0,cp_failure_level,routineP,error,failure)
    DO ikind= 1, nlinks
       CALL section_vals_val_get(qmmm_link_section,"LINK_TYPE",i_rep_section=ikind,i_val=link_type,error=error)
       IF (link_type ==  do_qmmm_link_imomm)  n_imomm  = n_imomm  + 1
       IF (link_type ==  do_qmmm_link_gho  )  n_gho    = n_gho    + 1
       IF (link_type ==  do_qmmm_link_pseudo) n_pseudo = n_pseudo + 1
    END DO
    n_tot = n_imomm + n_gho + n_pseudo
    CPPostcondition(n_tot /= 0,cp_failure_level,routineP,error,failure)
    ALLOCATE(qmmm_links, stat=stat)
    CPPostcondition(stat == 0,cp_failure_level,routineP,error,failure)
    NULLIFY(qmmm_links%imomm,&
            qmmm_links%pseudo)
    ! IMOMM
    IF ( n_imomm /= 0 ) THEN
       ALLOCATE(qmmm_links%imomm(n_imomm), stat=stat)
       CPPostcondition(stat == 0,cp_failure_level,routineP,error,failure)
       ALLOCATE(wrk_tmp(n_imomm), stat=stat)
       CPPostcondition(stat == 0,cp_failure_level,routineP,error,failure)
       DO ikind = 1, n_imomm
          NULLIFY(qmmm_links%imomm(ikind)%link)
          ALLOCATE(qmmm_links%imomm(ikind)%link,stat=stat)
          CPPostcondition(stat == 0,cp_failure_level,routineP,error,failure)
       END DO
       n_imomm = 0
       DO ikind = 1, nlinks
          CALL section_vals_val_get(qmmm_link_section,"LINK_TYPE",i_rep_section=ikind,i_val=link_type,error=error)
          IF (link_type ==  do_qmmm_link_imomm) THEN
             n_imomm = n_imomm + 1
             CALL section_vals_val_get(qmmm_link_section,"QM_INDEX",i_rep_section=ikind,i_val=qm_index,error=error)
             CALL section_vals_val_get(qmmm_link_section,"MM_INDEX",i_rep_section=ikind,i_val=mm_index,error=error)
             CALL section_vals_val_get(qmmm_link_section,"ALPHA_IMOMM",i_rep_section=ikind,r_val=alpha,error=error)
             CALL section_vals_val_get(qmmm_link_section,"RADIUS",i_rep_section=ikind,n_rep_val=n_rep_val,error=error)
             qmmm_links%imomm(n_imomm)%link%qm_index = qm_index
             qmmm_links%imomm(n_imomm)%link%mm_index = mm_index
             qmmm_links%imomm(n_imomm)%link%alpha    = alpha
             wrk_tmp(n_imomm) = mm_index
             IF (n_rep_val == 1) THEN
                CALL section_vals_val_get(qmmm_link_section,"RADIUS",i_rep_section=ikind,r_val=my_radius,error=error)
                WHERE ( mm_atom_index == mm_index ) mm_el_pot_radius = my_radius * qm_cell_small%unit_of_length
                WHERE ( mm_atom_index == mm_index ) mm_el_pot_radius_corr = my_radius * qm_cell_small%unit_of_length
             END IF
             CALL section_vals_val_get(qmmm_link_section,"CORR_RADIUS",i_rep_section=ikind,n_rep_val=n_rep_val,error=error)
             IF (n_rep_val == 1) THEN
                CALL section_vals_val_get(qmmm_link_section,"CORR_RADIUS",i_rep_section=ikind,r_val=my_radius,error=error)
                WHERE ( mm_atom_index == mm_index ) mm_el_pot_radius_corr = my_radius * qm_cell_small%unit_of_length
             END IF
          END IF
       END DO
       !
       ! Checking the link structure
       !
       DO ikind = 1, size(wrk_tmp)
          IF (COUNT(wrk_tmp == wrk_tmp(ikind)) > 1) THEN
             WRITE(iw,'(/A)')"In the IMOMM scheme no more than one QM atom can be bounded to the same MM atom."
             WRITE(iw, '(A)')"Multiple link MM atom not allowed. Check your link sections."
             CPPostcondition(.false.,cp_failure_level,routineP,error,failure)
          END IF
       END DO
       DEALLOCATE(wrk_tmp, stat=stat)
       CPPostcondition(stat == 0,cp_failure_level,routineP,error,failure)
    END IF
    ! PSEUDO
    IF ( n_pseudo /= 0 ) THEN
       ALLOCATE(qmmm_links%pseudo(n_pseudo), stat=stat)
       CPPostcondition(stat == 0,cp_failure_level,routineP,error,failure)
       DO ikind = 1, n_pseudo
          NULLIFY(qmmm_links%pseudo(ikind)%link)
          ALLOCATE(qmmm_links%pseudo(ikind)%link,stat=stat)
          CPPostcondition(stat == 0,cp_failure_level,routineP,error,failure)
       END DO
       n_pseudo = 0
       DO ikind = 1, nlinks
          CALL section_vals_val_get(qmmm_link_section,"LINK_TYPE",i_rep_section=ikind,i_val=link_type,error=error)
          IF (link_type ==  do_qmmm_link_pseudo) THEN
             n_pseudo = n_pseudo + 1
             CALL section_vals_val_get(qmmm_link_section,"QM_INDEX",i_rep_section=ikind,i_val=qm_index,error=error)
             CALL section_vals_val_get(qmmm_link_section,"MM_INDEX",i_rep_section=ikind,i_val=mm_index,error=error)
             qmmm_links%pseudo(n_pseudo)%link%qm_index = qm_index
             qmmm_links%pseudo(n_pseudo)%link%mm_index = mm_index
          END IF
       END DO
    END IF
    ! GHO
    IF ( n_gho /= 0) THEN
       ! not yet implemented
       ! still to define : type, implementation into QS
       CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure) 
    END IF
  END SUBROUTINE setup_qmmm_links


  !!****f* qmmm_init/setup_qmmm_links [1.0] *
  !!
  !!   NAME
  !!     setup_qmmm_links
  !!
  !!   FUNCTION
  !!     this routine sets up all variables to treat qmmm links
  !!
  !!   NOTES
  !!
  !!
  !!   INPUTS
  !!     - error: variable to control error logging, stopping,... 
  !!       see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     Teodoro Laino
  !!
  !!   MODIFICATION HISTORY
  !!     12.2004 created [tlaino]
  !!
  !!*************************************************************************
  SUBROUTINE move_or_add_atoms(qmmm_section, move_mm_charges, add_mm_charges, &
                               mm_atom_chrg, mm_el_pot_radius, mm_el_pot_radius_corr,&
                               added_charges, mm_atom_index, qm_cell_small, error)
    IMPLICIT NONE
    ! Arguments
    TYPE(section_vals_type), POINTER         :: qmmm_section
    LOGICAL, INTENT(OUT)                     :: move_mm_charges, add_mm_charges
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mm_atom_chrg, mm_el_pot_radius
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mm_el_pot_radius_corr
    TYPE(add_set_type), POINTER              :: added_charges
    INTEGER, POINTER, DIMENSION(:)           :: mm_atom_index
    TYPE(cell_type), POINTER                 :: qm_cell_small
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                            :: error
    ! Local variables
    CHARACTER(len=*), PARAMETER :: routineN = 'move_or_add_atoms',&
         routineP = moduleN//':'//routineN
    LOGICAL :: failure, explicit
    INTEGER :: nlinks, n_moves, n_adds, icount, ikind, n_move_tot, n_add_tot
    INTEGER :: Index1, Index2, i_add, ind1, n_rep_val
    TYPE(section_vals_type), POINTER         :: qmmm_link_section, move_section,&
                                                add_section
    REAL(KIND=dp) :: radius, alpha, charge, c_radius

    failure = .FALSE.
    explicit = .FALSE.
    move_mm_charges = .FALSE.
    add_mm_charges  = .FALSE.
    NULLIFY( qmmm_link_section, move_section, add_section)
    IF (.NOT.failure) THEN
       qmmm_link_section => section_vals_get_subs_vals(qmmm_section,"LINK",error=error)
       CALL section_vals_get(qmmm_link_section,n_repetition=nlinks,error=error)
       CPPostcondition(nlinks /= 0,cp_failure_level,routineP,error,failure)
       icount     = 0
       n_move_tot = 0
       n_add_tot  = 0
       DO ikind= 1, nlinks
          move_section => section_vals_get_subs_vals(qmmm_link_section,"MOVE_MM_CHARGE",&
               i_rep_section=ikind, error=error)
          CALL section_vals_get(move_section,n_repetition=n_moves,error=error)
          add_section => section_vals_get_subs_vals(qmmm_link_section,"ADD_MM_CHARGE",&
               i_rep_section=ikind, error=error)
          CALL section_vals_get(add_section,n_repetition=n_adds,error=error)
          n_move_tot = n_move_tot + n_moves
          n_add_tot  = n_add_tot  + n_adds
       END DO
       icount = n_move_tot + n_add_tot
       IF (n_add_tot   /= 0)  add_mm_charges  = .TRUE.
       IF (n_move_tot  /= 0) move_mm_charges  = .TRUE.
       !
       ! create add_set_type
       !
       CALL create_add_set_type(added_charges, ndim=icount, error=error)
       !
       ! Fill in structures
       !
       icount = 0
       DO ikind= 1, nlinks
          move_section => section_vals_get_subs_vals(qmmm_link_section,"MOVE_MM_CHARGE",&
                                                     i_rep_section=ikind, error=error)
          CALL section_vals_get(move_section, explicit=explicit, n_repetition=n_moves, error=error)
          !
          ! Moving charge atoms
          !
          IF (explicit) THEN
             DO i_add = 1, n_moves
                icount = icount + 1
                CALL section_vals_val_get(move_section,"ATOM_INDEX_1",i_val=Index1,i_rep_section=i_add,error=error)
                CALL section_vals_val_get(move_section,"ATOM_INDEX_2",i_val=Index2,i_rep_section=i_add,error=error)
                CALL section_vals_val_get(move_section,"ALPHA",r_val=alpha,i_rep_section=i_add,error=error)
                CALL section_vals_val_get(move_section,"RADIUS",r_val=radius,i_rep_section=i_add,error=error)
                CALL section_vals_val_get(move_section,"CORR_RADIUS",n_rep_val=n_rep_val,i_rep_section=i_add,error=error)
                c_radius = radius
                IF (n_rep_val == 1) &
                     CALL section_vals_val_get(move_section,"CORR_RADIUS",r_val=c_radius,i_rep_section=i_add,error=error)

                radius   = radius * qm_cell_small%unit_of_length
                c_radius = c_radius * qm_cell_small%unit_of_length
                CALL set_add_set_type(added_charges, icount, Index1, Index2, alpha, radius, c_radius,&
                                      mm_atom_chrg=mm_atom_chrg, mm_el_pot_radius=mm_el_pot_radius,&
                                      mm_el_pot_radius_corr=mm_el_pot_radius_corr,&
                                      mm_atom_index=mm_atom_index, move=n_moves, Ind1=ind1,&
                                      error=error)
             END DO
             mm_atom_chrg(ind1) = 0.0_dp
          END IF
          
          add_section => section_vals_get_subs_vals(qmmm_link_section,"ADD_MM_CHARGE",&
                                                    i_rep_section=ikind, error=error)
          CALL section_vals_get(add_section, explicit=explicit, n_repetition=n_adds, error=error)
          !
          ! Adding charge atoms
          !
          IF (explicit) THEN
             DO i_add = 1, n_adds
                icount = icount + 1
                CALL section_vals_val_get(add_section,"ATOM_INDEX_1",i_val=Index1,i_rep_section=i_add,error=error)
                CALL section_vals_val_get(add_section,"ATOM_INDEX_2",i_val=Index2,i_rep_section=i_add,error=error)
                CALL section_vals_val_get(add_section,"ALPHA",r_val=alpha,i_rep_section=i_add,error=error)
                CALL section_vals_val_get(add_section,"RADIUS",r_val=radius,i_rep_section=i_add,error=error)
                CALL section_vals_val_get(add_section,"CHARGE",r_val=charge,i_rep_section=i_add,error=error)
                CALL section_vals_val_get(move_section,"CORR_RADIUS",n_rep_val=n_rep_val,i_rep_section=i_add,error=error)
                c_radius = radius
                IF (n_rep_val == 1) &
                     CALL section_vals_val_get(move_section,"CORR_RADIUS",r_val=c_radius,i_rep_section=i_add,error=error)

                radius   = radius * qm_cell_small%unit_of_length
                c_radius = c_radius * qm_cell_small%unit_of_length
                CALL set_add_set_type(added_charges, icount, Index1, Index2, alpha, radius, c_radius, charge,&
                                      mm_atom_chrg=mm_atom_chrg, mm_el_pot_radius=mm_el_pot_radius,&
                                      mm_el_pot_radius_corr=mm_el_pot_radius_corr,&
                                      mm_atom_index=mm_atom_index, error=error)
             END DO
          END IF
       END DO
    END IF

   END SUBROUTINE move_or_add_atoms

  !!****f* qmmm_init/setup_qmmm_links [1.0] *
  !!
  !!   NAME
  !!     setup_qmmm_links
  !!
  !!   FUNCTION
  !!     this routine sets up all variables of the add_set_type type
  !!
  !!   NOTES
  !!
  !!
  !!   INPUTS
  !!     - error: variable to control error logging, stopping,... 
  !!       see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     Teodoro Laino
  !!
  !!   MODIFICATION HISTORY
  !!     12.2004 created [tlaino]
  !!
  !!*************************************************************************
   SUBROUTINE set_add_set_type(added_charges, icount, Index1, Index2, alpha, radius, c_radius, charge,&
        mm_atom_chrg, mm_el_pot_radius, mm_el_pot_radius_corr, mm_atom_index, move, ind1, error)
     IMPLICIT NONE
     ! Arguments
     REAL(KIND=dp), DIMENSION(:), POINTER     :: mm_atom_chrg, mm_el_pot_radius
     REAL(KIND=dp), DIMENSION(:), POINTER     :: mm_el_pot_radius_corr
     TYPE(add_set_type), POINTER              :: added_charges
     INTEGER, POINTER, DIMENSION(:)           :: mm_atom_index
     INTEGER, INTENT(IN)  :: Index1, Index2, icount
     INTEGER, INTENT(OUT),OPTIONAL  :: ind1
     REAL(KIND=dp), INTENT(IN) :: alpha, radius, c_radius
     REAL(KIND=dp), INTENT(IN), OPTIONAL :: charge
     INTEGER, OPTIONAL, INTENT(in) :: move
     TYPE(cp_error_type), INTENT(inout), &
          OPTIONAL                            :: error
     ! Local Variables
     CHARACTER(len=*), PARAMETER :: routineN = 'set_add_set_type',&
          routineP = moduleN//':'//routineN
     LOGICAL :: failure
     INTEGER :: i, my_move
     REAL(KIND=dp) :: my_charge, my_radius, my_c_radius
     
     failure = .FALSE.
     my_move = 0
     my_radius   = radius
     my_c_radius = c_radius
     IF (PRESENT(charge)) my_charge = charge
     IF (PRESENT(move))   my_move   = move
     IF (.NOT.failure) THEN
        i = 1
        GetId: DO WHILE (i<=SIZE(mm_atom_index))
           IF (Index1 == mm_atom_index(i)) EXIT GetId
           i = i+1
        END DO GetId
        IF (PRESENT(ind1)) ind1 = i
        CPPostcondition(i<=SIZE(mm_atom_index),cp_failure_level,routineP,error,failure)
        IF (.NOT.PRESENT(charge)) my_charge = mm_atom_chrg(i)/REAL(my_move,KIND=dp)
        IF (my_radius   == 0.0_dp) my_radius   = mm_el_pot_radius(i)
        IF (my_c_radius == 0.0_dp) my_c_radius = mm_el_pot_radius_corr(i)

        added_charges%add_env(icount)%Index1        = Index1
        added_charges%add_env(icount)%Index2        = Index2
        added_charges%add_env(icount)%alpha         = alpha 
        added_charges%mm_atom_index(icount)         = icount
        added_charges%mm_atom_chrg(icount)          = my_charge
        added_charges%mm_el_pot_radius(icount)      = my_radius
        added_charges%mm_el_pot_radius_corr(icount) = my_c_radius
     END IF
   END SUBROUTINE set_add_set_type

   !!****f* qmmm_init/setup_origin_mm_cell [1.0] *
   !!
   !!   NAME
   !!     setup_origin_mm_cell
   !!
   !!   FUNCTION
   !!     this routine sets up the origin of the MM cell respect to the 
   !!     origin of the QM cell. The origin of the QM cell is assumed to be 
   !!     in (0.0,0.0,0.0)...
   !!
   !!   NOTES
   !!
   !!
   !!   INPUTS
   !!     - error: variable to control error logging, stopping,... 
   !!       see module cp_error_handling 
   !!
   !!   AUTHOR
   !!     Teodoro Laino
   !!
   !!   MODIFICATION HISTORY
   !!     02.2005 created [tlaino]
   !!
   !!*************************************************************************
   SUBROUTINE setup_origin_mm_cell(qmmm_section, qmmm_env, mm_cell, qm_cell_small, &
        dr, error)
     IMPLICIT NONE
     ! Arguments
     TYPE(section_vals_type), POINTER            :: qmmm_section
     TYPE(cell_type), POINTER :: mm_cell, qm_cell_small
     TYPE(qmmm_env_qm_type), POINTER             :: qmmm_env
     REAL(KIND=dp), DIMENSION(3), INTENT(in) :: dr
     TYPE(cp_error_type), INTENT(inout), &
          OPTIONAL                               :: error     
     ! Local Variabels
     CHARACTER(len=*), PARAMETER :: routineN = 'setup_origin_mm_cell',&
          routineP = moduleN//':'//routineN
     LOGICAL :: failure
     REAL(KIND=dp), DIMENSION(3) :: tmp
     INTEGER :: my_type
     
     failure = .FALSE.
     IF (.NOT.failure) THEN
        ! This is the vector that corrects position to apply properly the PBC
        tmp(1) =  qm_cell_small%hmat(1,1)
        tmp(2) =  qm_cell_small%hmat(2,2)
        tmp(3) =  qm_cell_small%hmat(3,3)
        CPPostcondition(ALL(tmp>0),cp_failure_level,routineP,error,failure)     
        qmmm_env%dOmmOqm = tmp/2.0_dp
        ! This is unit vector to translate the QM system in order to center it
        ! in QM cell
        CALL section_vals_val_get(qmmm_section,"TYP_CENTER",i_val=my_type,error=error)
        SELECT CASE(my_type)
        CASE(do_center_none)
           qmmm_env%utrasl = 1.0_dp
        CASE(do_center_grid)
           qmmm_env%utrasl = dr
        END SELECT
     END IF
   END SUBROUTINE setup_origin_mm_cell

END MODULE qmmm_init
