!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Initialize a QM/MM calculation
!> \par History
!>      5.2004 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
MODULE qmmm_init
  USE atomic_kind_list_types,          ONLY: atomic_kind_list_type
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             set_atomic_kind
  USE cell_methods,                    ONLY: read_cell
  USE cell_types,                      ONLY: cell_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE cp_units,                        ONLY: cp_unit_from_cp2k,&
                                             cp_unit_to_cp2k
  USE external_potential_types,        ONLY: fist_potential_type,&
                                             get_potential,&
                                             set_potential
  USE force_field_types,               ONLY: input_info_type
  USE force_fields_input,              ONLY: read_gd_section,&
                                             read_gp_section,&
                                             read_lj_section,&
                                             read_wl_section
  USE input_constants,                 ONLY: &
       RADIUS_QMMM_DEFAULT, do_qmmm_gauss, do_qmmm_image_calcmatrix, &
       do_qmmm_image_iter, do_qmmm_link_gho, do_qmmm_link_imomm, &
       do_qmmm_link_pseudo, do_qmmm_pcharge, do_qmmm_swave
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE pair_potential_types,            ONLY: pair_potential_reallocate
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: particle_type
  USE pw_env_types,                    ONLY: pw_env_type
  USE qmmm_elpot,                      ONLY: qmmm_potential_init
  USE qmmm_ff_fist,                    ONLY: qmmm_ff_precond_only_qm
  USE qmmm_gaussian_init,              ONLY: qmmm_gaussian_initialize
  USE qmmm_per_elpot,                  ONLY: qmmm_ewald_potential_init,&
                                             qmmm_per_potential_init
  USE qmmm_types,                      ONLY: add_set_type,&
                                             create_add_set_type,&
                                             qmmm_env_mm_type,&
                                             qmmm_env_qm_type,&
                                             qmmm_links_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE termination,                     ONLY: stop_program
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qmmm_init'

  PUBLIC :: assign_mm_charges_and_radius,&
            print_qmmm_charges,&
            print_qmmm_links,&
            print_image_charge_info,&
            qmmm_init_gaussian_type,&
            qmmm_init_potential,&
            qmmm_init_periodic_potential,&
            setup_qmmm_vars_qm,&
            setup_qmmm_vars_mm,&
            setup_qmmm_links,&
            move_or_add_atoms,&
            setup_origin_mm_cell

!***
CONTAINS

! *****************************************************************************
!> \brief Assigns charges and radius to evaluate the MM electrostatic potential
!> \param subsys the subsys containing the MM charges
!> \param charges ...
!> \param mm_atom_chrg ...
!> \param mm_el_pot_radius ...
!> \param mm_el_pot_radius_corr ...
!> \param mm_atom_index ...
!> \param mm_link_atoms ...
!> \param mm_link_scale_factor ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      06.2004 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
 SUBROUTINE assign_mm_charges_and_radius(subsys, charges, mm_atom_chrg,mm_el_pot_radius,&
                                         mm_el_pot_radius_corr, mm_atom_index, mm_link_atoms,&
                                         mm_link_scale_factor, error)
    TYPE(cp_subsys_type), POINTER            :: subsys
    REAL(KIND=dp), DIMENSION(:), POINTER     :: charges
    REAL(dp), DIMENSION(:), POINTER          :: mm_atom_chrg, &
                                                mm_el_pot_radius, &
                                                mm_el_pot_radius_corr
    INTEGER, DIMENSION(:), POINTER           :: mm_atom_index, mm_link_atoms
    REAL(dp), DIMENSION(:), POINTER          :: mm_link_scale_factor
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'assign_mm_charges_and_radius', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: I, ilink, IndMM
    REAL(dp)                                 :: qi, rc, ri
    TYPE(atomic_kind_type), POINTER          :: my_kind
    TYPE(fist_potential_type), POINTER       :: my_potential
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

    NULLIFY(particle_set, my_kind)
    CALL cp_subsys_get(subsys=subsys,particles=particles,error=error)
    particle_set => particles%els

    DO I = 1, SIZE(mm_atom_index)
       IndMM = mm_atom_index(I)
       my_kind => particle_set(IndMM)%atomic_kind
       CALL get_atomic_kind(atomic_kind=my_kind, fist_potential=my_potential)
       CALL get_potential  (potential=my_potential,&
                            qeff=qi,&
                            qmmm_radius=ri,&
                            qmmm_corr_radius=rc)
       IF (ASSOCIATED(charges)) qi=charges(IndMM)
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

! *****************************************************************************
!> \brief Print info on charges generating the qmmm potential..
!> \param mm_atom_index ...
!> \param mm_atom_chrg ...
!> \param mm_el_pot_radius ...
!> \param mm_el_pot_radius_corr ...
!> \param added_charges ...
!> \param qmmm_section ...
!> \param nocompatibility ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      01.2005 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
 SUBROUTINE print_qmmm_charges(mm_atom_index, mm_atom_chrg, mm_el_pot_radius, mm_el_pot_radius_corr,&
                               added_charges, qmmm_section, nocompatibility, error)
    INTEGER, DIMENSION(:), POINTER           :: mm_atom_index
    REAL(dp), DIMENSION(:), POINTER          :: mm_atom_chrg, &
                                                mm_el_pot_radius, &
                                                mm_el_pot_radius_corr
    TYPE(add_set_type), POINTER              :: added_charges
    TYPE(section_vals_type), POINTER         :: qmmm_section
    LOGICAL, INTENT(IN)                      :: nocompatibility
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'print_qmmm_charges', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: I, ind1, ind2, IndMM, iw
    REAL(KIND=dp)                            :: qi, qtot, rc, ri
    TYPE(cp_logger_type), POINTER            :: logger

    qtot = 0.0_dp
    logger => cp_error_get_logger(error)
    iw=cp_print_key_unit_nr(logger,qmmm_section,"PRINT%QMMM_CHARGES",&
                                         extension=".log",error=error)
    IF (iw>0) THEN
       WRITE(iw,FMT="(/,T2,A)") REPEAT("-",79)
       WRITE(iw,FMT='(/5X,A)')"MM    POINT CHARGES GENERATING THE QM/MM ELECTROSTATIC POTENTIAL"
       WRITE(iw,FMT="(/,T2,A)") REPEAT("-",79)
       DO I = 1, SIZE(mm_atom_index)
          IndMM = mm_atom_index(I)
          qi = mm_atom_chrg(I)
          qtot = qtot + qi
          ri = mm_el_pot_radius(I)
          rc = mm_el_pot_radius_corr(I)
          IF (nocompatibility) THEN
             WRITE(iw,'(5X,A9,T15,I5,T28,A8,T38,F12.6,T60,A8,T69,F12.6)')' MM ATOM:',IndMM,' RADIUS:',ri,&
                  ' CHARGE:',qi
          ELSE
             WRITE(iw,'(5X,A9,T15,I5,T28,A8,T38,F12.6,T60,A8,T69,F12.6,/,T56,A12,T69,F12.6)')&
                  ' MM ATOM:',IndMM,' RADIUS:',ri,' CHARGE:',qi,'CORR. RADIUS',rc
          END IF
       END DO
       IF (added_charges%num_mm_atoms /= 0) THEN
          WRITE(iw,FMT="(/,T2,A)") REPEAT("-",79)
          WRITE(iw,'(/5X,A)')"ADDED POINT CHARGES GENERATING THE QM/MM ELECTROSTATIC POTENTIAL"
          WRITE(iw,FMT="(/,T2,A)") REPEAT("-",79)
          DO I = 1, SIZE(added_charges%mm_atom_index)
             IndMM = added_charges%mm_atom_index(I)
             qi = added_charges%mm_atom_chrg(I)
             qtot = qtot + qi
             ri = added_charges%mm_el_pot_radius(I)
             ind1 = added_charges%add_env(I)%Index1
             ind2 = added_charges%add_env(I)%Index2
             IF (nocompatibility) THEN
                WRITE(iw,'(5X,A9,I5,T25,A8,T35,F12.6,T50,A8,T59,F12.6,I5,I5)')'MM POINT:',IndMM,' RADIUS:',ri,&
                     ' CHARGE:',qi,ind1,ind2
             ELSE
                WRITE(iw,'(5X,A9,I5,T25,A8,T35,F12.6,T50,A8,T59,F12.6,I5,I5,/,T56,A12,T69,F12.6)')&
                     'MM POINT:',IndMM,' RADIUS:',ri,' CHARGE:',qi,ind1,ind2,'CORR. RADIUS',rc
             END IF
          END DO

       END IF
       WRITE(iw,FMT="(/,T2,A)") REPEAT("-",79)
       WRITE(iw,'(/,T50,A,T69,F12.6)')' TOTAL CHARGE:',qtot
       WRITE(iw,FMT="(/,T2,A,/)") REPEAT("-",79)
    END IF
    CALL cp_print_key_finished_output(iw,logger,qmmm_section,&
               "PRINT%QMMM_CHARGES", error=error)
  END SUBROUTINE print_qmmm_charges

! *****************************************************************************
!> \brief Print info on qm/mm links
!> \param qmmm_section ...
!> \param qmmm_links ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      01.2005 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
 SUBROUTINE print_qmmm_links( qmmm_section, qmmm_links, error)
    TYPE(section_vals_type), POINTER         :: qmmm_section
    TYPE(qmmm_links_type), POINTER           :: qmmm_links
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'print_qmmm_links', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iw, mm_index, qm_index
    REAL(KIND=dp)                            :: alpha
    TYPE(cp_logger_type), POINTER            :: logger

    logger => cp_error_get_logger(error)
    iw=cp_print_key_unit_nr(logger,qmmm_section,"PRINT%qmmm_link_info", extension=".log",error=error)
    IF (iw>0) THEN
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
                WRITE(iw,FMT="(T2,A,T20,A9,I8,1X,A9,I8,T62,A6,F12.6)")"TYPE: IMOMM",&
                     "QM INDEX:",qm_index,"MM INDEX:",mm_index,"ALPHA:",alpha
             END DO
          END IF
          IF (ASSOCIATED(qmmm_links%pseudo)) THEN
             WRITE(iw,FMT="(/,T31,A)") " PSEUDO TYPE LINK "
             DO I = 1, SIZE(qmmm_links%pseudo)
                qm_index = qmmm_links%pseudo(I)%link%qm_index
                mm_index = qmmm_links%pseudo(I)%link%mm_index
                WRITE(iw,FMT="(T2,A,T20,A9,I8,1X,A9,I8)")"TYPE: PSEUDO",&
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
    CALL cp_print_key_finished_output(iw,logger,qmmm_section,&
               "PRINT%qmmm_link_info", error=error)
  END SUBROUTINE print_qmmm_links

! *****************************************************************************
!> \brief ...
!> \param qmmm_env_qm ...
!> \param qm_cell_small ...
!> \param mm_cell ...
!> \param para_env ...
!> \param mm_atom_chrg ...
!> \param qs_env ...
!> \param added_charges ...
!> \param print_section ...
!> \param qmmm_section ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      1.2005 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
 SUBROUTINE qmmm_init_gaussian_type(qmmm_env_qm, qm_cell_small, mm_cell, para_env, &
                                    mm_atom_chrg, qs_env, added_charges, print_section, &
                                    qmmm_section, error )
    TYPE(qmmm_env_qm_type), POINTER          :: qmmm_env_qm
    TYPE(cell_type), POINTER                 :: qm_cell_small, mm_cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mm_atom_chrg
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(add_set_type), POINTER              :: added_charges
    TYPE(section_vals_type), POINTER         :: print_section, qmmm_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_init_gaussian_type', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: maxchrg
    REAL(KIND=dp), DIMENSION(:), POINTER     :: maxradius, maxradius2
    TYPE(pw_env_type), POINTER               :: pw_env

    failure=.FALSE.
    stat=0
    NULLIFY(maxradius, maxradius2, pw_env)

    maxchrg = MAXVAL(ABS(mm_atom_chrg(:)))
    CALL get_qs_env(qs_env, pw_env=pw_env, error=error)
    IF (qmmm_env_qm%add_mm_charges) maxchrg=MAX(maxchrg,MAXVAL(ABS(added_charges%mm_atom_chrg(:))))
    CALL qmmm_gaussian_initialize(qmmm_gaussian_fns=qmmm_env_qm%pgfs,&
                                  qm_cell_small=qm_cell_small,&
                                  para_env=para_env,&
                                  pw_env=pw_env,&
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
                                     para_env=para_env,&
                                     pw_env=pw_env,&
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

       SELECT CASE(qmmm_env_qm%qmmm_coupl_type)
       CASE(do_qmmm_gauss,do_qmmm_swave,do_qmmm_pcharge)
          DO i = 1, SIZE(maxradius)
             maxradius(i) = MAX( maxradius(i), maxradius2(i))
          END DO
       END SELECT

       IF (ASSOCIATED(maxradius2)) DEALLOCATE(maxradius2, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

    qmmm_env_qm%maxradius => maxradius

 END SUBROUTINE qmmm_init_gaussian_type

! *****************************************************************************
!> \brief ...
!> \param qmmm_env_qm ...
!> \param qm_cell_small ...
!> \param mm_cell ...
!> \param para_env ...
!> \param added_charges ...
!> \param print_section ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      1.2005 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
 SUBROUTINE qmmm_init_potential(qmmm_env_qm, qm_cell_small, mm_cell, para_env, &
                                added_charges, print_section, error)
    TYPE(qmmm_env_qm_type), POINTER          :: qmmm_env_qm
    TYPE(cell_type), POINTER                 :: qm_cell_small, mm_cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(add_set_type), POINTER              :: added_charges
    TYPE(section_vals_type), POINTER         :: print_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qmmm_init_potential', &
      routineP = moduleN//':'//routineN

    CALL qmmm_potential_init(qmmm_coupl_type=qmmm_env_qm%qmmm_coupl_type,&
                             mm_el_pot_radius=qmmm_env_qm%mm_el_pot_radius,&
                             potentials=qmmm_env_qm%potentials,&
                             pgfs=qmmm_env_qm%pgfs,&
                             mm_cell=mm_cell,&
                             compatibility=qmmm_env_qm%compatibility,&
                             print_section=print_section,&
                             error=error)

    IF (qmmm_env_qm%move_mm_charges.OR.qmmm_env_qm%add_mm_charges) THEN

       CALL qmmm_potential_init(qmmm_coupl_type=qmmm_env_qm%qmmm_coupl_type,&
                                mm_el_pot_radius=added_charges%mm_el_pot_radius,&
                                potentials=added_charges%potentials,&
                                pgfs=added_charges%pgfs,&
                                mm_cell=mm_cell,&
                                compatibility=qmmm_env_qm%compatibility,&
                                print_section=print_section,&
                                error=error)
    END IF

 END SUBROUTINE qmmm_init_potential

! *****************************************************************************
!> \brief ...
!> \param qmmm_env_qm ...
!> \param qm_cell_small ...
!> \param mm_cell ...
!> \param para_env ...
!> \param qs_env ...
!> \param added_charges ...
!> \param qmmm_periodic ...
!> \param print_section ...
!> \param mm_atom_chrg ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      7.2005 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
 SUBROUTINE qmmm_init_periodic_potential(qmmm_env_qm, qm_cell_small, mm_cell, para_env, qs_env,&
                                         added_charges, qmmm_periodic, print_section, mm_atom_chrg, error)
    TYPE(qmmm_env_qm_type), POINTER          :: qmmm_env_qm
    TYPE(cell_type), POINTER                 :: qm_cell_small, mm_cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(add_set_type), POINTER              :: added_charges
    TYPE(section_vals_type), POINTER         :: qmmm_periodic, print_section
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mm_atom_chrg
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qmmm_init_periodic_potential', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp)                            :: maxchrg
    TYPE(dft_control_type), POINTER          :: dft_control

   IF (qmmm_env_qm%periodic) THEN

      NULLIFY(dft_control)
      CALL get_qs_env(qs_env, dft_control=dft_control, error=error)

      IF(dft_control%qs_control%semi_empirical) THEN
         CALL cp_unimplemented_error(fromWhere=routineP, &
            message="QM/MM periodic calculations not implemented for semi empirical methods", &
            error=error, error_level=cp_fatal_level)
      ELSE IF(dft_control%qs_control%dftb) THEN
         CALL qmmm_ewald_potential_init(qmmm_env_qm%ewald_env, qmmm_env_qm%ewald_pw,&
              qmmm_coupl_type=qmmm_env_qm%qmmm_coupl_type, mm_cell=mm_cell,&
              para_env=para_env, qmmm_periodic=qmmm_periodic, print_section=print_section,&
              error=error)
      ELSE IF(dft_control%qs_control%scptb) THEN
         CALL cp_unimplemented_error(fromWhere=routineP, &
            message="QM/MM periodic calculations not implemented for SCPTB method", &
            error=error, error_level=cp_fatal_level)
      ELSE
 
         ! setup for GPW/GPAW
         maxchrg = MAXVAL(ABS(mm_atom_chrg(:)))
         IF (qmmm_env_qm%add_mm_charges) maxchrg=MAX(maxchrg,MAXVAL(ABS(added_charges%mm_atom_chrg(:))))

         CALL qmmm_per_potential_init(qmmm_coupl_type=qmmm_env_qm%qmmm_coupl_type,&
              per_potentials=qmmm_env_qm%per_potentials,&
              potentials=qmmm_env_qm%potentials,&
              pgfs=qmmm_env_qm%pgfs,&
              qm_cell_small=qm_cell_small,&
              mm_cell=mm_cell,&
              para_env=para_env,&
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
                 para_env=para_env,&
                 compatibility=qmmm_env_qm%compatibility,&
                 qmmm_periodic=qmmm_periodic,&
                 print_section=print_section,&
                 eps_mm_rspace=qmmm_env_qm%eps_mm_rspace,&
                 maxchrg=maxchrg,&
                 ncp =qmmm_env_qm%aug_pools(SIZE(qmmm_env_qm%aug_pools))%pool%pw_grid%npts,&
                 ncpl=qmmm_env_qm%aug_pools(SIZE(qmmm_env_qm%aug_pools))%pool%pw_grid%npts_local,&
                 error=error)
         END IF

      END IF

   END IF

 END SUBROUTINE qmmm_init_periodic_potential

! *****************************************************************************
!> \brief ...
!> \param qmmm_section ...
!> \param qmmm_env ...
!> \param subsys_mm ...
!> \param qm_atom_type ...
!> \param qm_atom_index ...
!> \param mm_atom_index ...
!> \param qm_cell_small ...
!> \param qmmm_coupl_type ...
!> \param eps_mm_rspace ...
!> \param qmmm_link ...
!> \param para_env ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      11.2004 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE setup_qmmm_vars_qm( qmmm_section, qmmm_env, subsys_mm, qm_atom_type,&
       qm_atom_index, mm_atom_index, qm_cell_small, qmmm_coupl_type, eps_mm_rspace,&
       qmmm_link, para_env, error)
    TYPE(section_vals_type), POINTER         :: qmmm_section
    TYPE(qmmm_env_qm_type), POINTER          :: qmmm_env
    TYPE(cp_subsys_type), POINTER            :: subsys_mm
    CHARACTER(len=default_string_length), &
      DIMENSION(:), POINTER                  :: qm_atom_type
    INTEGER, DIMENSION(:), POINTER           :: qm_atom_index, mm_atom_index
    TYPE(cell_type), POINTER                 :: qm_cell_small
    INTEGER, INTENT(OUT)                     :: qmmm_coupl_type
    REAL(KIND=dp), INTENT(OUT)               :: eps_mm_rspace
    LOGICAL, INTENT(OUT)                     :: qmmm_link
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'setup_qmmm_vars_qm', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length)     :: atmname, mm_atom_kind
    INTEGER                                  :: i, icount, ikind, ikindr, &
                                                my_type, n_rep_val, nkind, &
                                                size_mm_system, stat
    INTEGER, DIMENSION(:), POINTER           :: mm_link_atoms
    LOGICAL                                  :: explicit, failure, is_mm, &
                                                is_qm
    REAL(KIND=dp)                            :: tmp_radius, tmp_radius_c
    REAL(KIND=dp), DIMENSION(:), POINTER     :: tmp_sph_cut
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(fist_potential_type), POINTER       :: fist_potential
    TYPE(section_vals_type), POINTER         :: cell_section, &
                                                image_charge_section, mm_kinds

    NULLIFY(mm_link_atoms, cell_section, tmp_sph_cut)
    NULLIFY(image_charge_section)
    qmmm_link = .FALSE.
    failure   = .FALSE.

    CALL section_vals_get(qmmm_section,explicit=explicit,error=error)
    IF (explicit) THEN
       !
       ! QM_CELL
       !
       cell_section => section_vals_get_subs_vals(qmmm_section,"CELL",error=error)
       CALL read_cell( qm_cell_small, qm_cell_small, cell_section=cell_section,&
            check_for_ref=.FALSE., para_env=para_env, error=error)
       CALL section_vals_val_get(qmmm_section,"E_COUPL",i_val=qmmm_coupl_type,error=error)
       CALL section_vals_val_get(qmmm_section,"EPS_MM_RSPACE",r_val=eps_mm_rspace,error=error)
       CALL section_vals_val_get(qmmm_section,"SPHERICAL_CUTOFF",r_vals=tmp_sph_cut,error=error)
       CPPostcondition(SIZE(tmp_sph_cut)==2,cp_failure_level,routineP,error,failure)
       qmmm_env%spherical_cutoff = tmp_sph_cut
       IF (qmmm_env%spherical_cutoff(1)<=0.0_dp) THEN
          qmmm_env%spherical_cutoff(2) = 0.0_dp
       ELSE
          IF (qmmm_env%spherical_cutoff(2)<=0.0_dp) qmmm_env%spherical_cutoff(2)=EPSILON(0.0_dp)
          tmp_radius = qmmm_env%spherical_cutoff(1)-20.0_dp*qmmm_env%spherical_cutoff(2)
          CALL cp_assert(tmp_radius>0.0_dp,cp_failure_level,cp_assertion_failed,routineP,&
               "SPHERICAL_CUTOFF(1) > 20*SPHERICAL_CUTOFF(1)! Please correct parameters for "//&
               "the Spherical Cutoff in order to satisfy the previous condition!"//&
CPSourceFileRef,&
               error,failure)
       END IF
       !
       ! Initialization of arrays and core_charge_radius...
       !
       tmp_radius=0.0_dp
       CALL cp_subsys_get(subsys=subsys_mm,atomic_kinds=atomic_kinds,&
                          error=error)
       DO Ikind = 1, SIZE(atomic_kinds%els)
          atomic_kind => atomic_kinds%els(Ikind)
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
                               fist_potential=fist_potential)
          CALL set_potential  (potential=fist_potential,&
                               qmmm_radius=tmp_radius,&
                               qmmm_corr_radius=tmp_radius)
          CALL set_atomic_kind(atomic_kind=atomic_kind,&
                               fist_potential=fist_potential)
       END DO
       CALL setup_qm_atom_list(qmmm_section=qmmm_section,&
                               qm_atom_index=qm_atom_index,&
                               qm_atom_type=qm_atom_type,&
                               mm_link_atoms=mm_link_atoms,&
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
       tmp_radius = cp_unit_to_cp2k(RADIUS_QMMM_DEFAULT,"angstrom",error=error)
       Set_Radius_Pot_0: DO IkindR = 1, SIZE(atomic_kinds%els)
          atomic_kind => atomic_kinds%els(IkindR)
          CALL get_atomic_kind(atomic_kind=atomic_kind,name=atmname)
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
               fist_potential=fist_potential)
          CALL set_potential(potential=fist_potential,qmmm_radius=tmp_radius,&
               qmmm_corr_radius=tmp_radius)
          CALL set_atomic_kind(atomic_kind=atomic_kind,&
               fist_potential=fist_potential)
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
             CALL section_vals_val_get(mm_kinds,"CORR_RADIUS",i_rep_section=ikind,n_rep_val=n_rep_val,&
                  error=error)
             IF (n_rep_val == 1) CALL section_vals_val_get(mm_kinds,"CORR_RADIUS",i_rep_section=ikind,&
                  r_val=tmp_radius_c,error=error)
             Set_Radius_Pot_1: DO IkindR = 1, SIZE(atomic_kinds%els)
                atomic_kind => atomic_kinds%els(IkindR)
                CALL get_atomic_kind(atomic_kind=atomic_kind,name=atmname)
                is_qm = qmmm_ff_precond_only_qm(atmname)
                IF (TRIM(mm_atom_kind)==atmname) THEN
                   CALL get_atomic_kind(atomic_kind=atomic_kind,&
                        fist_potential=fist_potential)
                   CALL set_potential(potential=fist_potential,&
                        qmmm_radius=tmp_radius,&
                        qmmm_corr_radius=tmp_radius_c)
                   CALL set_atomic_kind(atomic_kind=atomic_kind,&
                        fist_potential=fist_potential)
                END IF
             END DO Set_Radius_Pot_1
          END DO
       END IF

       !Image charge section

       image_charge_section => section_vals_get_subs_vals(qmmm_section,"IMAGE_CHARGE",error=error)
       CALL section_vals_get(image_charge_section,explicit=qmmm_env%image_charge,error=error)

    ELSE
     CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
          "QMMM section not present in input file!"//&
CPSourceFileRef,&
          error,failure)
    ENDIF
    !
    ! Build MM atoms list
    !
    size_mm_system = SIZE(subsys_mm%particles%els) - SIZE(qm_atom_index)
    IF (qmmm_link.AND.ASSOCIATED(mm_link_atoms)) size_mm_system = size_mm_system + SIZE(mm_link_atoms)
    ALLOCATE(mm_atom_index(size_mm_system))
    icount = 0

    DO i = 1, SIZE(subsys_mm%particles%els)
       is_mm = .TRUE.
       IF ( ANY(qm_atom_index == i)) THEN
         is_mm =.FALSE.
       END IF
       IF (ASSOCIATED(mm_link_atoms)) THEN
          IF (ANY( mm_link_atoms == i ).AND.qmmm_link) is_mm =.TRUE.
       END IF
       IF (is_mm) THEN
          icount = icount + 1
          IF (icount <= size_mm_system) mm_atom_index(icount) = i
       END IF
    END DO
    CPPostcondition(icount==size_mm_system,cp_failure_level,routineP,error,failure)
    IF (ASSOCIATED(mm_link_atoms)) THEN
       DEALLOCATE(mm_link_atoms, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

    ! Build image charge atom list + set up variables
    !
    IF (qmmm_env%image_charge) THEN
       CALL section_vals_val_get(image_charge_section,"MM_ATOM_LIST",&
                                 explicit=explicit,error=error)
       IF(explicit) qmmm_env%image_charge_pot%all_mm=.FALSE.

       IF(qmmm_env%image_charge_pot%all_mm) THEN
         qmmm_env%image_charge_pot%image_mm_list => mm_atom_index
       ELSE
         CALL setup_image_atom_list(image_charge_section,qmmm_env,&
                                    qm_atom_index,subsys_mm,error)
       END IF

       qmmm_env%image_charge_pot%particles_all => subsys_mm%particles%els

       CALL section_vals_val_get(image_charge_section,"EXT_POTENTIAL",&
                                 r_val=qmmm_env%image_charge_pot%V0,error=error)
       CALL section_vals_val_get(image_charge_section,"WIDTH",&
                                 r_val=qmmm_env%image_charge_pot%eta,error=error)
       CALL section_vals_val_get(image_charge_section,"DETERM_COEFF",&
                                 i_val=my_type,error=error)
       SELECT CASE(my_type)
       CASE(do_qmmm_image_calcmatrix)
         qmmm_env%image_charge_pot%coeff_iterative=.FALSE.
       CASE(do_qmmm_image_iter)
         qmmm_env%image_charge_pot%coeff_iterative=.TRUE.
       END SELECT
     
       CALL section_vals_val_get(image_charge_section,"RESTART_IMAGE_MATRIX",&
                                l_val=qmmm_env%image_charge_pot%image_restart,&
                                error=error)
    END IF

  END SUBROUTINE setup_qmmm_vars_qm

! *****************************************************************************
!> \brief ...
!> \param qmmm_section ...
!> \param qmmm_env ...
!> \param qm_atom_index ...
!> \param mm_link_atoms ...
!> \param mm_link_scale_factor ...
!> \param fist_scale_charge_link ...
!> \param qmmm_coupl_type ...
!> \param qmmm_link ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      12.2004 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE setup_qmmm_vars_mm( qmmm_section, qmmm_env, qm_atom_index,  &
                                 mm_link_atoms, mm_link_scale_factor, &
                                 fist_scale_charge_link, qmmm_coupl_type,&
                                 qmmm_link, error)
    TYPE(section_vals_type), POINTER         :: qmmm_section
    TYPE(qmmm_env_mm_type), POINTER          :: qmmm_env
    INTEGER, DIMENSION(:), POINTER           :: qm_atom_index, mm_link_atoms
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mm_link_scale_factor, &
                                                fist_scale_charge_link
    INTEGER, INTENT(OUT)                     :: qmmm_coupl_type
    LOGICAL, INTENT(OUT)                     :: qmmm_link
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'setup_qmmm_vars_mm', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: explicit
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
       CALL section_vals_get(qmmm_ff_section,explicit=qmmm_env%use_qmmm_ff,error=error)
       IF (qmmm_env%use_qmmm_ff) THEN
          CALL section_vals_val_get(qmmm_ff_section,"MULTIPLE_POTENTIAL",&
               l_val=qmmm_env%multiple_potential,error=error)
          CALL read_qmmm_ff_section(qmmm_ff_section, qmmm_env%inp_info, error=error)
       END IF
    END IF
  END SUBROUTINE setup_qmmm_vars_mm

! *****************************************************************************
!> \brief reads information regarding the forcefield specific for the QM/MM
!>      interactions
!> \param qmmm_ff_section ...
!> \param inp_info ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      12.2004 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE read_qmmm_ff_section(qmmm_ff_section, inp_info, error)
    TYPE(section_vals_type), POINTER         :: qmmm_ff_section
    TYPE(input_info_type), POINTER           :: inp_info
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'read_qmmm_ff_section', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: n_gd, n_gp, n_lj, n_wl, np
    LOGICAL                                  :: failure
    TYPE(section_vals_type), POINTER         :: gd_section, gp_section, &
                                                lj_section, wl_section

    failure = .FALSE.
    !
    ! NONBONDED
    !
    lj_section => section_vals_get_subs_vals(qmmm_ff_section,"NONBONDED%LENNARD-JONES",error=error)
    wl_section => section_vals_get_subs_vals(qmmm_ff_section,"NONBONDED%WILLIAMS",error=error)
    gd_section => section_vals_get_subs_vals(qmmm_ff_section,"NONBONDED%GOODWIN",error=error)
    gp_section => section_vals_get_subs_vals(qmmm_ff_section,"NONBONDED%GENPOT",error=error)
    CALL section_vals_get(lj_section,n_repetition=n_lj,error=error)
    np = n_lj
    IF (n_lj /= 0) THEN
       CALL pair_potential_reallocate(inp_info%nonbonded,1,np,lj_charmm=.TRUE.,error=error)
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
    CALL section_vals_get(gp_section,n_repetition=n_gp,error=error)
    np = n_lj + n_wl + n_gd + n_gp
    IF (n_gp /= 0) THEN
       CALL pair_potential_reallocate(inp_info%nonbonded,1,np,gp=.TRUE.,error=error)
       CALL read_gp_section(inp_info%nonbonded,gp_section, start=n_lj+n_wl+n_gd, error=error)
    END IF
    !
    ! NONBONDED14
    !
    lj_section => section_vals_get_subs_vals(qmmm_ff_section,"NONBONDED14%LENNARD-JONES",error=error)
    wl_section => section_vals_get_subs_vals(qmmm_ff_section,"NONBONDED14%WILLIAMS",error=error)
    gd_section => section_vals_get_subs_vals(qmmm_ff_section,"NONBONDED14%GOODWIN",error=error)
    gp_section => section_vals_get_subs_vals(qmmm_ff_section,"NONBONDED14%GENPOT",error=error)
    CALL section_vals_get(lj_section,n_repetition=n_lj,error=error)
    np = n_lj
    IF (n_lj /= 0) THEN
       CALL pair_potential_reallocate(inp_info%nonbonded14,1,np,lj_charmm=.TRUE.,error=error)
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
    CALL section_vals_get(gp_section,n_repetition=n_gp,error=error)
    np = n_lj + n_wl + n_gd + n_gp
    IF (n_gp /= 0) THEN
       CALL pair_potential_reallocate(inp_info%nonbonded14,1,np,gp=.TRUE.,error=error)
       CALL read_gp_section(inp_info%nonbonded14,gp_section, start=n_lj+n_wl+n_gd, error=error)
    END IF

  END SUBROUTINE read_qmmm_ff_section

! *****************************************************************************
!> \brief ...
!> \param qmmm_section ...
!> \param qm_atom_index ...
!> \param qm_atom_type ...
!> \param mm_link_atoms ...
!> \param mm_link_scale_factor ...
!> \param qmmm_link ...
!> \param fist_scale_charge_link ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      12.2004 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE setup_qm_atom_list(qmmm_section,qm_atom_index,qm_atom_type,&
       mm_link_atoms, mm_link_scale_factor, qmmm_link, fist_scale_charge_link, error)
    TYPE(section_vals_type), POINTER         :: qmmm_section
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: qm_atom_index
    CHARACTER(len=default_string_length), &
      DIMENSION(:), OPTIONAL, POINTER        :: qm_atom_type
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: mm_link_atoms
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: mm_link_scale_factor
    LOGICAL, INTENT(OUT), OPTIONAL           :: qmmm_link
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: fist_scale_charge_link
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'setup_qm_atom_list', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length)     :: qm_atom_kind, qm_link_element
    INTEGER                                  :: ikind, k, link_involv_mm, &
                                                link_type, mm_index, n_var, &
                                                nkind, nlinks, &
                                                num_qm_atom_tot, stat
    INTEGER, DIMENSION(:), POINTER           :: mm_indexes
    LOGICAL                                  :: explicit, failure
    REAL(KIND=dp)                            :: scale_f
    TYPE(section_vals_type), POINTER         :: qm_kinds, qmmm_links

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
    IF (PRESENT(qm_atom_index))   qm_atom_index    =   0
    IF (PRESENT(qm_atom_type))    qm_atom_type     = " "
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
          num_qm_atom_tot = num_qm_atom_tot + SIZE(mm_indexes)
       END DO
    END DO
    IF (PRESENT(mm_link_scale_factor).AND.(link_involv_mm /= 0))     mm_link_scale_factor = 0.0_dp
    IF (PRESENT(fist_scale_charge_link).AND.(link_involv_mm /= 0)) fist_scale_charge_link = 0.0_dp
    IF (PRESENT(mm_link_atoms).AND.(link_involv_mm /= 0))            mm_link_atoms        = 0
    IF (explicit) THEN
       DO ikind = 1, nlinks
          IF (PRESENT(qm_atom_type)) THEN
              CALL section_vals_val_get(qmmm_links,"QM_KIND",i_rep_section=ikind,c_val=qm_link_element,&
                  error=error)
             qm_atom_type (num_qm_atom_tot:num_qm_atom_tot)   = TRIM(qm_link_element)//"_LINK"
          END IF
          IF (PRESENT(qm_atom_index)) THEN
             CALL section_vals_val_get(qmmm_links,"MM_INDEX",i_rep_section=ikind,i_val=mm_index,&
                  error=error)
             CPPostcondition(ALL(qm_atom_index/=mm_index),cp_failure_level,routineP,error,failure)
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

! *****************************************************************************
!> \brief this routine sets up all variables to treat qmmm links
!> \param qmmm_section ...
!> \param qmmm_links ...
!> \param mm_el_pot_radius ...
!> \param mm_el_pot_radius_corr ...
!> \param mm_atom_index ...
!> \param qm_cell_small ...
!> \param iw ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      12.2004 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE  setup_qmmm_links(qmmm_section, qmmm_links, mm_el_pot_radius, mm_el_pot_radius_corr,&
       mm_atom_index, qm_cell_small, iw, error)
    TYPE(section_vals_type), POINTER         :: qmmm_section
    TYPE(qmmm_links_type), POINTER           :: qmmm_links
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mm_el_pot_radius, &
                                                mm_el_pot_radius_corr
    INTEGER, DIMENSION(:), POINTER           :: mm_atom_index
    TYPE(cell_type), POINTER                 :: qm_cell_small
    INTEGER, INTENT(IN)                      :: iw
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'setup_qmmm_links', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikind, link_type, mm_index, &
                                                n_gho, n_imomm, n_pseudo, &
                                                n_rep_val, n_tot, nlinks, &
                                                qm_index, stat
    INTEGER, DIMENSION(:), POINTER           :: wrk_tmp
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: alpha, my_radius
    TYPE(section_vals_type), POINTER         :: qmmm_link_section

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
                WHERE ( mm_atom_index == mm_index ) mm_el_pot_radius = my_radius
                WHERE ( mm_atom_index == mm_index ) mm_el_pot_radius_corr = my_radius
             END IF
             CALL section_vals_val_get(qmmm_link_section,"CORR_RADIUS",i_rep_section=ikind,n_rep_val=n_rep_val,error=error)
             IF (n_rep_val == 1) THEN
                CALL section_vals_val_get(qmmm_link_section,"CORR_RADIUS",i_rep_section=ikind,r_val=my_radius,error=error)
                WHERE ( mm_atom_index == mm_index ) mm_el_pot_radius_corr = my_radius
             END IF
          END IF
       END DO
       !
       ! Checking the link structure
       !
       DO ikind = 1, SIZE(wrk_tmp)
          IF (COUNT(wrk_tmp == wrk_tmp(ikind)) > 1) THEN
             WRITE(iw,'(/A)')"In the IMOMM scheme no more than one QM atom can be bounded to the same MM atom."
             WRITE(iw, '(A)')"Multiple link MM atom not allowed. Check your link sections."
             CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
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

! *****************************************************************************
!> \brief this routine sets up all variables to treat qmmm links
!> \param qmmm_section ...
!> \param move_mm_charges ...
!> \param add_mm_charges ...
!> \param mm_atom_chrg ...
!> \param mm_el_pot_radius ...
!> \param mm_el_pot_radius_corr ...
!> \param added_charges ...
!> \param mm_atom_index ...
!> \param qm_cell_small ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      12.2004 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE move_or_add_atoms(qmmm_section, move_mm_charges, add_mm_charges, &
                               mm_atom_chrg, mm_el_pot_radius, mm_el_pot_radius_corr,&
                               added_charges, mm_atom_index, qm_cell_small, error)
    TYPE(section_vals_type), POINTER         :: qmmm_section
    LOGICAL, INTENT(OUT)                     :: move_mm_charges, &
                                                add_mm_charges
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mm_atom_chrg, &
                                                mm_el_pot_radius, &
                                                mm_el_pot_radius_corr
    TYPE(add_set_type), POINTER              :: added_charges
    INTEGER, DIMENSION(:), POINTER           :: mm_atom_index
    TYPE(cell_type), POINTER                 :: qm_cell_small
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'move_or_add_atoms', &
      routineP = moduleN//':'//routineN

    INTEGER :: i_add, icount, ikind, ind1, Index1, Index2, n_add_tot, n_adds, &
      n_move_tot, n_moves, n_rep_val, nlinks
    LOGICAL                                  :: explicit, failure
    REAL(KIND=dp)                            :: alpha, c_radius, charge, &
                                                radius
    TYPE(section_vals_type), POINTER         :: add_section, move_section, &
                                                qmmm_link_section

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
                CALL section_vals_val_get(add_section,"CORR_RADIUS",n_rep_val=n_rep_val,i_rep_section=i_add,error=error)
                c_radius = radius
                IF (n_rep_val == 1) &
                     CALL section_vals_val_get(add_section,"CORR_RADIUS",r_val=c_radius,i_rep_section=i_add,error=error)

                CALL set_add_set_type(added_charges, icount, Index1, Index2, alpha, radius, c_radius, charge,&
                                      mm_atom_chrg=mm_atom_chrg, mm_el_pot_radius=mm_el_pot_radius,&
                                      mm_el_pot_radius_corr=mm_el_pot_radius_corr,&
                                      mm_atom_index=mm_atom_index, error=error)
             END DO
          END IF
       END DO
    END IF

   END SUBROUTINE move_or_add_atoms

! *****************************************************************************
!> \brief this routine sets up all variables of the add_set_type type
!> \param added_charges ...
!> \param icount ...
!> \param Index1 ...
!> \param Index2 ...
!> \param alpha ...
!> \param radius ...
!> \param c_radius ...
!> \param charge ...
!> \param mm_atom_chrg ...
!> \param mm_el_pot_radius ...
!> \param mm_el_pot_radius_corr ...
!> \param mm_atom_index ...
!> \param move ...
!> \param ind1 ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      12.2004 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
   SUBROUTINE set_add_set_type(added_charges, icount, Index1, Index2, alpha, radius, c_radius, charge,&
        mm_atom_chrg, mm_el_pot_radius, mm_el_pot_radius_corr, mm_atom_index, move, ind1, error)
    TYPE(add_set_type), POINTER              :: added_charges
    INTEGER, INTENT(IN)                      :: icount, Index1, Index2
    REAL(KIND=dp), INTENT(IN)                :: alpha, radius, c_radius
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: charge
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mm_atom_chrg, &
                                                mm_el_pot_radius, &
                                                mm_el_pot_radius_corr
    INTEGER, DIMENSION(:), POINTER           :: mm_atom_index
    INTEGER, INTENT(in), OPTIONAL            :: move
    INTEGER, INTENT(OUT), OPTIONAL           :: ind1
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'set_add_set_type', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, my_move
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: my_c_radius, my_charge, &
                                                my_radius

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

! *****************************************************************************
!> \brief this routine sets up the origin of the MM cell respect to the
!>      origin of the QM cell. The origin of the QM cell is assumed to be
!>      in (0.0,0.0,0.0)...
!> \param qmmm_section ...
!> \param qmmm_env ...
!> \param mm_cell ...
!> \param qm_cell_small ...
!> \param dr ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      02.2005 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
   SUBROUTINE setup_origin_mm_cell(qmmm_section, qmmm_env, mm_cell, qm_cell_small, &
        dr, error)
    TYPE(section_vals_type), POINTER         :: qmmm_section
    TYPE(qmmm_env_qm_type), POINTER          :: qmmm_env
    TYPE(cell_type), POINTER                 :: mm_cell, qm_cell_small
    REAL(KIND=dp), DIMENSION(3), INTENT(in)  :: dr
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'setup_origin_mm_cell', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: center_grid, failure
    REAL(KIND=dp), DIMENSION(3)              :: tmp
    REAL(KINd=dp), DIMENSION(:), POINTER     :: vec

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
        CALL section_vals_val_get(qmmm_section,"CENTER_GRID",l_val=center_grid,error=error)
        IF (center_grid) THEN
           qmmm_env%utrasl = dr
        ELSE
           qmmm_env%utrasl = 1.0_dp
        ENDIF
        CALL section_vals_val_get(qmmm_section,"INITIAL_TRANSLATION_VECTOR",r_vals=vec,error=error)
        qmmm_env%transl_v = vec
     END IF
   END SUBROUTINE setup_origin_mm_cell

! *****************************************************************************
!> \brief this routine sets up list of MM atoms carrying an image charge
!> \param image_charge_section ...
!> \param qmmm_env ...
!> \param qm_atom_index ...
!> \param subsys_mm ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      02.2012 created
!> \author Dorothea Golze
! *****************************************************************************
   SUBROUTINE setup_image_atom_list(image_charge_section, qmmm_env,&
                                    qm_atom_index, subsys_mm, error)

    TYPE(section_vals_type), POINTER         :: image_charge_section
    TYPE(qmmm_env_qm_type), POINTER          :: qmmm_env
    INTEGER, DIMENSION(:), POINTER           :: qm_atom_index
    TYPE(cp_subsys_type), POINTER            :: subsys_mm
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'setup_image_atom_list', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: atom_a, atom_b, i, j, &
                                                max_index, n_var, &
                                                num_image_mm_atom, stat
    INTEGER, DIMENSION(:), POINTER           :: mm_indexes
    LOGICAL                                  :: failure, imageind_in_range

    NULLIFY(mm_indexes)
    failure=.FALSE.
    imageind_in_range=.FALSE.
    num_image_mm_atom=0
    max_index=0

    CALL section_vals_val_get(image_charge_section,"MM_ATOM_LIST",&
                              n_rep_val=n_var,error=error)
    DO i = 1, n_var
       CALL section_vals_val_get(image_charge_section,"MM_ATOM_LIST",&
                                 i_rep_val=i,i_vals=mm_indexes,error=error)
       num_image_mm_atom = num_image_mm_atom + SIZE(mm_indexes)
    END DO

    ALLOCATE(qmmm_env%image_charge_pot%image_mm_list(num_image_mm_atom),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    qmmm_env%image_charge_pot%image_mm_list=0
    num_image_mm_atom=1

    DO i = 1, n_var
       CALL section_vals_val_get(image_charge_section,"MM_ATOM_LIST",&
                                 i_rep_val=i,i_vals=mm_indexes,error=error)
       qmmm_env%image_charge_pot%image_mm_list(num_image_mm_atom:num_image_mm_atom &
                                             + SIZE(mm_indexes)-1) = mm_indexes(:)
       num_image_mm_atom = num_image_mm_atom + SIZE(mm_indexes)
    END DO

   ! checking, if in range, if list contains QM atoms or any atoms doubled
   num_image_mm_atom=num_image_mm_atom-1

   max_index=SIZE(subsys_mm%particles%els)

   CPPrecondition(SIZE(qmmm_env%image_charge_pot%image_mm_list) /= 0,cp_failure_level,routineP,error,failure)
   imageind_in_range=(MAXVAL(qmmm_env%image_charge_pot%image_mm_list)<= max_index)&
                     .AND.(MINVAL(qmmm_env%image_charge_pot%image_mm_list) > 0)
   CPPostcondition(imageind_in_range,cp_failure_level,routineP,error,failure)

   DO i=1,num_image_mm_atom
      atom_a=qmmm_env%image_charge_pot%image_mm_list(i)
      IF(ANY(qm_atom_index==atom_a)) THEN
         CALL stop_program(routineN,moduleN,__LINE__,&
             "Image atom list must only contain MM atoms")
      ENDIF
      DO j=i+1,num_image_mm_atom
         atom_b=qmmm_env%image_charge_pot%image_mm_list(j)
         IF(atom_a==atom_b) &
         CALL stop_program(routineN,moduleN,__LINE__,&
             "There are atoms doubled in image list.")
      ENDDO
   ENDDO


   END SUBROUTINE setup_image_atom_list

! *****************************************************************************
!> \brief Print info on image charges
!> \param qmmm_env ...
!> \param qmmm_section ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      03.2012 created
!> \author Dorothea Golze
! *****************************************************************************
 SUBROUTINE print_image_charge_info(qmmm_env, qmmm_section, error)

    TYPE(qmmm_env_qm_type), POINTER          :: qmmm_env
    TYPE(section_vals_type), POINTER         :: qmmm_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'print_image_charge_info', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iw
    REAL(KIND=dp)                            :: eta, eta_conv, V0, V0_conv
    TYPE(cp_logger_type), POINTER            :: logger

    logger => cp_error_get_logger(error)
    iw=cp_print_key_unit_nr(logger,qmmm_section,"PRINT%PROGRAM_RUN_INFO",&
                                         extension=".log",error=error)
    eta=qmmm_env%image_charge_pot%eta
    eta_conv = cp_unit_from_cp2k(eta,"angstrom",power=-2,error=error)
    V0=qmmm_env%image_charge_pot%V0
    V0_conv = cp_unit_from_cp2k(V0,"volt",error=error)

    IF (iw>0) THEN
       WRITE(iw,FMT="(T25,A)")"IMAGE CHARGE PARAMETERS"
       WRITE(iw,FMT="(T25,A)") REPEAT("-",23)
       WRITE(iw,FMT="(/)")
       WRITE(iw,FMT="(T2,A)")"INDEX OF MM ATOMS CARRYING AN IMAGE CHARGE:"
       WRITE(iw,FMT="(/)")

       WRITE(iw,"(7X,10I6)") qmmm_env%image_charge_pot%image_mm_list
       WRITE(iw,FMT="(/)")
       WRITE(iw,"(T2,A52,T69,F12.8)")&
        "WIDTH OF GAUSSIAN CHARGE DISTRIBUTION [angstrom^-2]:", eta_conv
       WRITE(iw,"(T2,A26,T69,F12.8)") "EXTERNAL POTENTIAL [volt]:", V0_conv
       WRITE(iw,FMT="(/,T2,A,/)") REPEAT("-",79)
    END IF
    CALL cp_print_key_finished_output(iw,logger,qmmm_section,&
               "PRINT%PROGRAM_RUN_INFO", error=error)

  END SUBROUTINE print_image_charge_info

END MODULE qmmm_init
