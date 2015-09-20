!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculates QM/MM energy and forces
!> \par History
!>      2015 Factored out of force_env_methods.F
!> \author Ole Schuett
! *****************************************************************************
MODULE qmmm_force
  USE cell_types,                      ONLY: cell_type
  USE cp_output_handling,              ONLY: cp_iter_string,&
                                             cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_result_methods,               ONLY: cp_results_erase,&
                                             get_results,&
                                             put_results
  USE cp_result_types,                 ONLY: cp_result_type
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE fist_energy_types,               ONLY: fist_energy_type
  USE fist_environment_types,          ONLY: fist_env_get
  USE fist_force,                      ONLY: fist_calc_energy_force
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE particle_types,                  ONLY: particle_type
  USE physcon,                         ONLY: debye
  USE qmmm_gpw_energy,                 ONLY: qmmm_el_coupling
  USE qmmm_gpw_forces,                 ONLY: qmmm_forces
  USE qmmm_links_methods,              ONLY: qmmm_added_chrg_coord,&
                                             qmmm_added_chrg_forces,&
                                             qmmm_link_Imomm_coord,&
                                             qmmm_link_Imomm_forces
  USE qmmm_types,                      ONLY: qmmm_env_type
  USE qmmm_types_low,                  ONLY: qmmm_links_type
  USE qmmm_util,                       ONLY: apply_qmmm_translate,&
                                             apply_qmmm_walls
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env
  USE qs_force,                        ONLY: qs_calc_energy_force
  USE qs_ks_qmmm_methods,              ONLY: ks_qmmm_env_rebuild
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qmmm_force'

  PUBLIC :: qmmm_calc_energy_force

CONTAINS

! *****************************************************************************
!> \brief calculates the qm/mm energy and forces
!> \param qmmm_env ...
!> \param calc_force if also the forces should be calculated
!> \param consistent_energies ...
!> \param linres ...
!> \par History
!>      05.2004 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE qmmm_calc_energy_force(qmmm_env,calc_force,consistent_energies,linres)
    TYPE(qmmm_env_type), POINTER             :: qmmm_env
    LOGICAL, INTENT(IN)                      :: calc_force, &
                                                consistent_energies, linres

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_calc_energy_force', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: description, iter
    INTEGER                                  :: ip, j, nres, output_unit
    INTEGER, DIMENSION(:), POINTER           :: qm_atom_index
    LOGICAL                                  :: check, failure, &
                                                qmmm_added_chrg, qmmm_link, &
                                                qmmm_link_imomm
    REAL(KIND=dp)                            :: energy_mm, energy_qm
    REAL(KIND=dp), DIMENSION(3)              :: dip_mm, dip_qm, dip_qmmm, &
                                                max_coord, min_coord
    TYPE(cell_type), POINTER                 :: mm_cell, qm_cell
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_result_type), POINTER            :: results_mm, results_qm, &
                                                results_qmmm
    TYPE(cp_subsys_type), POINTER            :: subsys, subsys_mm, subsys_qm
    TYPE(fist_energy_type), POINTER          :: fist_energy
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particles_mm, particles_qm
    TYPE(qmmm_links_type), POINTER           :: qmmm_links
    TYPE(qs_energy_type), POINTER            :: qs_energy
    TYPE(section_vals_type), POINTER         :: force_env_section, print_key

    min_coord        =  HUGE(0.0_dp)
    max_coord        = -HUGE(0.0_dp)
    failure          = .FALSE.
    qmmm_link        = .FALSE.
    qmmm_link_imomm  = .FALSE.
    qmmm_added_chrg  = .FALSE.
    logger => cp_get_default_logger()
    output_unit = cp_logger_get_default_io_unit(logger)
    NULLIFY(subsys_mm, subsys_qm, subsys, qm_atom_index,particles_mm,particles_qm, qm_cell, mm_cell)
    NULLIFY(force_env_section, print_key, results_qmmm, results_qm, results_mm)

    CALL get_qs_env(qmmm_env%qs_env, input=force_env_section)
    print_key => section_vals_get_subs_vals(force_env_section,"QMMM%PRINT%DIPOLE")

    CPPrecondition(ASSOCIATED(qmmm_env),cp_failure_level,routineP,failure)
    CPPrecondition(qmmm_env%ref_count>0,cp_failure_level,routineP,failure)

    CALL apply_qmmm_translate(qmmm_env)

    CALL fist_env_get(qmmm_env%fist_env,cell=mm_cell,subsys=subsys_mm)
    CALL get_qs_env(qmmm_env%qs_env,cell=qm_cell,cp_subsys=subsys_qm)
    qm_atom_index   => qmmm_env%qm%qm_atom_index
    qmmm_link       =  qmmm_env%qm%qmmm_link
    qmmm_links      => qmmm_env%qm%qmmm_links
    qmmm_added_chrg =  (qmmm_env%qm%move_mm_charges .OR. qmmm_env%qm%add_mm_charges)
    IF (qmmm_link) THEN
       CPPrecondition(ASSOCIATED(qmmm_links),cp_failure_level,routineP,failure)
       IF (ASSOCIATED(qmmm_links%imomm)) qmmm_link_imomm = (SIZE(qmmm_links%imomm) /= 0)
    END IF
    CPPrecondition(ASSOCIATED(qm_atom_index),cp_failure_level,routineP,failure)

    particles_mm => subsys_mm%particles%els
    particles_qm => subsys_qm%particles%els

    DO j=1,3
       IF (qm_cell%perd(j)==1) CYCLE
       DO ip=1, SIZE(particles_qm)
          check = (DOT_PRODUCT(qm_cell%h_inv(j,:),particles_qm(ip)%r) >= 0.0) .AND. &
                  (DOT_PRODUCT(qm_cell%h_inv(j,:),particles_qm(ip)%r) <= 1.0)
          IF (output_unit>0 .AND. .NOT. check) THEN
            WRITE (unit=output_unit, fmt='("ip, j, pos, lat_pos ",2I6,6F12.5)') ip, j, &
               particles_qm(ip)%r, DOT_PRODUCT(qm_cell%h_inv(j,:),particles_qm(ip)%r)
          ENDIF
          CALL cp_assert(check, cp_failure_level, cp_assertion_failed, routinep,&
                  "QM/MM QM atoms must be fully contained in the same image of the QM box "//&
                  "- No wrapping of coordinates is allowed! "//&
                  CPSourceFileRef)
       END DO
    END DO

    ! If present QM/MM links (just IMOMM) correct the position of the qm-link atom
    IF (qmmm_link_imomm) CALL qmmm_link_Imomm_coord(qmmm_links, particles_qm, qm_atom_index)

    ! If add charges get their position NOW!
    IF (qmmm_added_chrg) CALL qmmm_added_chrg_coord(qmmm_env%qm, particles_mm)

    ! Initialize ks_qmmm_env
    CALL ks_qmmm_env_rebuild(qs_env=qmmm_env%qs_env, qmmm_env=qmmm_env%qm)

    ! Compute the short range QM/MM Electrostatic Potential
    CALL qmmm_el_coupling(qs_env=qmmm_env%qs_env,&
         qmmm_env=qmmm_env%qm,&
         mm_particles=particles_mm,&
         mm_cell=mm_cell)

    ! Fist
    CALL fist_calc_energy_force(qmmm_env%fist_env)

    ! Print Out information on fist energy calculation...
    CALL fist_env_get(qmmm_env%fist_env,thermo=fist_energy)
    energy_mm = fist_energy%pot
    CALL cp_subsys_get(subsys_mm, results=results_mm)

    ! QS
    CALL qs_calc_energy_force(qmmm_env%qs_env,calc_force,consistent_energies,linres)

    ! QM/MM Interaction Potential forces
    CALL qmmm_forces(qmmm_env%qs_env,&
                     qmmm_env%qm,particles_mm,&
                     mm_cell=mm_cell,&
                     calc_force=calc_force)

    ! Forces of quadratic wall on QM atoms
    CALL apply_qmmm_walls(qmmm_env)

    ! Print Out information on QS energy calculation...
    CALL get_qs_env(qmmm_env%qs_env, energy=qs_energy)
    energy_qm = qs_energy%total
    CALL cp_subsys_get(subsys_qm, results=results_qm)

    !TODO: is really results_qm == results_qmmm ???
    CALL cp_subsys_get(subsys_qm, results=results_qmmm)

    IF (calc_force) THEN
       ! If present QM/MM links (just IMOMM) correct the position of the qm-link atom
       IF (qmmm_link_imomm) CALL qmmm_link_Imomm_forces(qmmm_links,particles_qm,qm_atom_index)
       particles_mm => subsys_mm%particles%els
       DO ip=1,SIZE(qm_atom_index)
          particles_mm(qm_atom_index(ip))%f=particles_mm(qm_atom_index(ip))%f+particles_qm(ip)%f
       END DO
       ! If add charges get rid of their derivatives right NOW!
       IF (qmmm_added_chrg) CALL qmmm_added_chrg_forces(qmmm_env%qm, particles_mm)
    END IF

    ! Handle some output
    output_unit = cp_print_key_unit_nr(logger,force_env_section,"QMMM%PRINT%DERIVATIVES",&
             extension=".Log")
    IF (output_unit>0) THEN
       WRITE (unit=output_unit,fmt='(/1X,A,F15.9)')"Energy after QMMM calculation: ",energy_qm
       IF (calc_force) THEN
          WRITE (unit=output_unit,fmt='(/1X,A)')"Derivatives on all atoms after QMMM calculation: "
          DO ip=1,SIZE(particles_mm)
             WRITE (unit=output_unit,fmt='(1X,3F15.9)') particles_mm(ip)%f
          END DO
       END IF
    END IF
    CALL cp_print_key_finished_output(output_unit,logger,force_env_section,&
         "QMMM%PRINT%DERIVATIVES")

    ! Dipole
    print_key => section_vals_get_subs_vals(force_env_section,"QMMM%PRINT%DIPOLE")
    IF (BTEST(cp_print_key_should_output(logger%iter_info,print_key),&
         cp_p_file)) THEN
       description ='[DIPOLE]'
       CALL get_results(results=results_qm,description=description,n_rep=nres)
       CPPrecondition(nres<=1,cp_failure_level,routineP,failure)
       CALL get_results(results=results_mm,description=description,n_rep=nres)
       CPPrecondition(nres<=1,cp_failure_level,routineP,failure)
       CALL get_results(results=results_qm,description=description,values=dip_qm)
       CALL get_results(results=results_mm,description=description,values=dip_mm)
       dip_qmmm = dip_qm + dip_mm
       CALL cp_results_erase(results=results_qmmm,description=description)
       CALL put_results(results=results_qmmm,description=description,values=dip_qmmm)

       output_unit = cp_print_key_unit_nr(logger,force_env_section,"QMMM%PRINT%DIPOLE",&
                extension=".Dipole")
       IF (output_unit>0) THEN
             WRITE(unit=output_unit,fmt="(A)")"QMMM TOTAL DIPOLE"
             WRITE(unit=output_unit,fmt="(A,T31,A,T88,A)")&
                  "# iter_level","dipole(x,y,z)[atomic units]",&
                                 "dipole(x,y,z)[debye]"
          iter=cp_iter_string(logger%iter_info)
          WRITE(unit=output_unit,fmt="(a,6(es18.8))")&
                  iter(1:15), dip_qmmm, dip_qmmm*debye
       END IF
    END IF

  END SUBROUTINE qmmm_calc_energy_force

END MODULE qmmm_force
