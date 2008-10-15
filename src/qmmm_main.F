!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Initialize a QM/MM calculation
!> \par History
!>      5.2004 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
MODULE qmmm_main
  USE bibliography,                    ONLY: Laino2005,&
                                             Laino2006,&
                                             cite_reference
  USE cell_types,                      ONLY: cell_clone,&
                                             cell_release,&
                                             cell_type,&
                                             get_cell,&
                                             write_cell
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsys_methods,               ONLY: create_small_subsys
  USE cp_subsys_types,                 ONLY: cp_subsys_release,&
                                             cp_subsys_type
  USE f77_blas
  USE fist_main,                       ONLY: fist_create_force_env
  USE force_env_methods,               ONLY: force_env_create
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_p_type,&
                                             force_env_type
  USE global_types,                    ONLY: global_environment_type
  USE header,                          ONLY: qmmm_header
  USE input_constants,                 ONLY: do_fist,&
                                             do_qmmm,&
                                             do_qmmm_none,&
                                             do_qs
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get,&
                                             section_vals_val_set
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE qmmm_init,                       ONLY: &
       assign_mm_charges_and_radius, move_or_add_atoms, print_qmmm_charges, &
       print_qmmm_links, qmmm_init_gaussian_type, &
       qmmm_init_periodic_potential, qmmm_init_potential, &
       setup_origin_mm_cell, setup_qmmm_links, setup_qmmm_vars_mm, &
       setup_qmmm_vars_qm
  USE qmmm_links_methods,              ONLY: qmmm_link_Imomm_coord
  USE qmmm_pw_grid,                    ONLY: qmmm_pw_grid_init
  USE qmmm_types,                      ONLY: &
       add_set_release, add_set_type, qmmm_env_mm_create, &
       qmmm_env_mm_release, qmmm_env_mm_type, qmmm_env_qm_create, &
       qmmm_env_qm_release, qmmm_env_qm_type, qmmm_links_type
  USE qs_environment_types,            ONLY: qs_environment_type,&
                                             set_qs_env
  USE qs_main,                         ONLY: quickstep_create_force_env
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qmmm_main'

  PUBLIC :: qmmm_create_force_env
!***
CONTAINS

! *****************************************************************************
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      05.2004 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE  qmmm_create_force_env(force_env, root_section, para_env, globenv,&
       force_env_section, subsys_section, use_motion_section, error)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(section_vals_type), POINTER         :: force_env_section, &
                                                subsys_section
    LOGICAL, INTENT(IN)                      :: use_motion_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_create_force_env', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length), &
      DIMENSION(:), POINTER                  :: qm_atom_type
    INTEGER                                  :: handle, isubf, iw, iw2, &
                                                qmmm_coupl_type, stat
    INTEGER, DIMENSION(:), POINTER           :: mm_atom_index, mm_link_atoms, &
                                                qm_atom_index
    LOGICAL :: add_mm_charges, explicit, failure, move_mm_charges, nocenter, &
      nocompatibility, qmmm_link, qmmm_link_Imomm
    REAL(dp), DIMENSION(:), POINTER          :: mm_atom_chrg, &
                                                mm_el_pot_radius, &
                                                mm_el_pot_radius_corr
    REAL(KIND=dp)                            :: eps_mm_rspace
    REAL(KIND=dp), DIMENSION(3)              :: abc
    REAL(KIND=dp), DIMENSION(:), POINTER     :: fist_scale_charge_link, &
                                                mm_link_scale_factor
    TYPE(add_set_type), POINTER              :: added_charges
    TYPE(cell_type), POINTER                 :: mm_cell, qm_cell_small, &
                                                super_cell
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_subsys_type), POINTER            :: subsys_mm, subsys_qm
    TYPE(force_env_p_type), DIMENSION(:), &
      POINTER                                :: sub_force_env
    TYPE(qmmm_env_mm_type), POINTER          :: qmmm_env_mm
    TYPE(qmmm_env_qm_type), POINTER          :: qmmm_env_qm
    TYPE(qmmm_links_type), POINTER           :: qmmm_links
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(section_vals_type), POINTER         :: print_gen, print_section, &
                                                qmmm_periodic, qmmm_section

    CALL timeset(routineN,handle)
    failure=.FALSE.

    NULLIFY(sub_force_env, qm_atom_index, mm_atom_index, qm_atom_type,   &
         qmmm_env_qm, subsys_mm,subsys_qm, mm_cell, qm_cell_small,       &
         qs_env, mm_atom_chrg,mm_el_pot_radius, qmmm_env_mm,   &
         mm_link_atoms, mm_link_scale_factor, qmmm_links, added_charges, &
         fist_scale_charge_link, print_section,  &
         print_gen, logger, mm_el_pot_radius_corr, super_cell)
    logger => cp_error_get_logger(error)
    IF (.NOT. failure) THEN
       ! citations
       CALL cite_reference(Laino2005)
       ! Input section...
       IF (.NOT.ASSOCIATED(subsys_section)) THEN
          subsys_section => section_vals_get_subs_vals(force_env_section,"SUBSYS",&
               error=error)
       END IF
       qmmm_section => section_vals_get_subs_vals(force_env_section,"QMMM",error=error)
       qmmm_periodic => section_vals_get_subs_vals(qmmm_section,"PERIODIC",error=error)
       print_section => section_vals_get_subs_vals(qmmm_section,"PRINT",error=error)
       print_gen => section_vals_get_subs_vals(print_section,"PROGRAM_RUN_INFO",error=error)
       iw = cp_print_key_unit_nr(logger,print_gen,"", extension=".log",error=error)

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

       ! Parallel scheme for the long range
       CALL section_vals_val_get(qmmm_section,"PARALLEL_SCHEME",&
            i_val=qmmm_env_qm%par_scheme,error=error)     

       ! Periodic boundary condition calculation
       CALL section_vals_get(qmmm_periodic,explicit=explicit,error=error)
       qmmm_env_qm%periodic = explicit
       IF (explicit) CALL cite_reference(Laino2006)
       IF (qmmm_coupl_type == do_qmmm_none) THEN
          qmmm_env_qm%periodic = .FALSE.
          CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
               "QMMM PERIODIC CALCULATION WITH NONE COUPLING WAS REQUESTED! "//&
               "SWITCHING OFF THE PERIODIC KEYWORD SINCE PERIODIC AND NON-PERIODIC "//&
               "CALCULATION WITH NONE COUPLING REPRESENT THE SAME METHOD! "//&
CPSourceFileRef,&
               only_ionode=.TRUE.)
       END IF

       ! First Initialize Fist...
       CALL section_vals_val_set(force_env_section,"METHOD",i_val=do_fist,error=error)
       CALL fist_create_force_env(sub_force_env(1)%force_env, &
            root_section, para_env, globenv, qmmm=.TRUE., qmmm_env=qmmm_env_mm,&
            force_env_section=force_env_section, subsys_section=subsys_section,&
            use_motion_section=use_motion_section,error=error)

       CALL force_env_get(sub_force_env(1)%force_env,subsys=subsys_mm,&
            cell=mm_cell,error=error)

       ! Set up QM/MM Options
       CALL setup_qmmm_vars_qm (qmmm_section,&
                                qmmm_env_qm,&
                                subsys_mm,&
                                qm_atom_type,&
                                qm_atom_index,&
                                mm_atom_index,&
                                qm_cell_small,&
                                qmmm_coupl_type,&
                                eps_mm_rspace,&
                                qmmm_link,&
                                para_env,&
                                error)

       qmmm_env_qm%qm_atom_index   => qm_atom_index
       qmmm_env_qm%mm_atom_index   => mm_atom_index
       qmmm_env_qm%eps_mm_rspace   = eps_mm_rspace
       qmmm_env_qm%qmmm_coupl_type = qmmm_coupl_type
       qmmm_env_qm%qmmm_link       = qmmm_link
       qmmm_env_qm%num_qm_atoms = SIZE(qm_atom_index)
       qmmm_env_qm%num_mm_atoms = SIZE(mm_atom_index)
       ! Duplicate structure for link atoms
       IF (qmmm_link) THEN
          IF (ASSOCIATED(mm_link_atoms)) THEN
             ALLOCATE(qmmm_env_qm%mm_link_atoms(SIZE(mm_link_atoms)),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             qmmm_env_qm%mm_link_atoms = mm_link_atoms
          END IF
       END IF
       IF ( iw > 0 ) THEN
          WRITE(iw,'(A,I8)')" Number of QM atoms: ",qmmm_env_qm%num_qm_atoms
          WRITE(iw,'(A,I8)')" Number of MM atoms: ",qmmm_env_qm%num_mm_atoms
          WRITE(iw,'(A)')" QM cell ::"
          CALL write_cell(qm_cell_small, subsys_section, error=error)
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

       qmmm_link_Imomm = .FALSE.
       IF (qmmm_link) THEN
          CALL setup_qmmm_links(qmmm_section, qmmm_links, mm_el_pot_radius, &
               mm_el_pot_radius_corr, mm_atom_index, qm_cell_small, iw, error=error)
          qmmm_env_qm%qmmm_links    => qmmm_links

          CALL print_qmmm_links(qmmm_section, qmmm_links, error)

          CALL add_set_release  (qmmm_env_qm%added_charges, error=error)
          CALL move_or_add_atoms(qmmm_section, move_mm_charges, add_mm_charges, &
                                 mm_atom_chrg, mm_el_pot_radius, mm_el_pot_radius_corr,&
                                 added_charges, mm_atom_index, qm_cell_small, error)
          qmmm_env_qm%move_mm_charges =  move_mm_charges
          qmmm_env_qm%add_mm_charges  =  add_mm_charges
          qmmm_env_qm%added_charges   => added_charges
          IF (ASSOCIATED(qmmm_links%imomm)) qmmm_link_imomm = (SIZE(qmmm_links%imomm) /= 0)
       END IF

       CALL print_qmmm_charges(mm_atom_index, mm_atom_chrg, mm_el_pot_radius,&
                               mm_el_pot_radius_corr, qmmm_env_qm%added_charges,&
                               qmmm_section, nocompatibility, error)

       CALL section_vals_val_set(force_env_section,"METHOD",i_val=do_qs,error=error)
       CALL create_small_subsys(subsys_qm,&
            big_subsys=subsys_mm,small_para_env=para_env,&
            small_cell=qm_cell_small,sub_atom_index=qm_atom_index,&
            sub_atom_kind_name=qm_atom_type, para_env=para_env, &
            force_env_section=force_env_section, subsys_section=subsys_section,error=error)
       IF (qmmm_link_imomm) CALL qmmm_link_Imomm_coord(qmmm_links, subsys_qm%particles%els,&
            qm_atom_index, error)
       CALL quickstep_create_force_env(sub_force_env(2)%force_env,root_section,para_env,&
            globenv, subsys=subsys_qm, cell=qm_cell_small, qmmm=.TRUE.,&
            qmmm_periodic=qmmm_env_qm%periodic, force_env_section=force_env_section,&
            subsys_section=subsys_section,use_motion_section=use_motion_section,&
            error=error)
       CALL cp_subsys_release(subsys_qm,error=error)
       CALL force_env_get(sub_force_env(2)%force_env,qs_env=qs_env,error=error)
       IF (qmmm_env_qm%periodic) THEN
          IF (.NOT.ASSOCIATED(super_cell)) THEN
             ALLOCATE (super_cell,STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF
          CALL cell_clone(mm_cell, super_cell,error=error)
          CALL set_qs_env(qs_env, super_cell=super_cell, qmmm_periodic=qmmm_env_qm%periodic, error=error)
          CALL cell_release(super_cell, error=error)
       END IF
       CALL cp_print_key_finished_output(iw,logger,print_gen,"", error=error)
       iw2 = cp_print_key_unit_nr(logger,qmmm_section,"PRINT%PROGRAM_BANNER",&
            extension=".qmmmLog",error=error)
       CALL qmmm_header(iw2, iw2>0)
       CALL cp_print_key_finished_output(iw2,logger,qmmm_section,&
         "PRINT%PROGRAM_BANNER",error=error)
       !
       ! Initialize MM Potential fitted with Gaussian
       !
       CALL qmmm_init_gaussian_type (qmmm_env_qm = qmmm_env_qm,&
                                     qm_cell_small=qm_cell_small,&
                                     mm_cell=mm_cell,&
                                     para_env=para_env,&
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
                                para_env=para_env,&
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
                              error=error)
       !
       ! Initialize the MM periodic potential
       !
       CALL qmmm_init_periodic_potential(qmmm_env_qm=qmmm_env_qm,&
                                qm_cell_small=qm_cell_small,&
                                mm_cell=mm_cell,&
                                para_env=para_env,&
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

       CALL section_vals_val_set(force_env_section,"METHOD",i_val=do_qmmm,error=error)
       CALL force_env_create(force_env,root_section, sub_force_env=sub_force_env,&
            qmmm_env=qmmm_env_qm,para_env=para_env,globenv=globenv, force_env_section= force_env_section,&
            error=error)
       CALL qmmm_env_qm_release(qmmm_env_qm, error)
       CALL qmmm_env_mm_release(qmmm_env_mm, error)
       DEALLOCATE(qm_atom_type,stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    END IF
    CALL timestop(0.0_dp,handle)

  END SUBROUTINE qmmm_create_force_env

END MODULE qmmm_main
