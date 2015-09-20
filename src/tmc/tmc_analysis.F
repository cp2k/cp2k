!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief module analyses element of the TMC tree element structure
!>        e.g. density, radial distribution function, dipole correlation,... 
!> \par History
!>      02.2013 created [Mandes Schönherr]
!> \author Mandes 
! *****************************************************************************

MODULE tmc_analysis
  USE cell_types,                      ONLY: cell_type,&
                                             get_cell,&
                                             pbc
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE force_fields_input,              ONLY: read_chrg_section
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE mathconstants,                   ONLY: pi
  USE mathlib,                         ONLY: diag
  USE physcon,                         ONLY: a_mass,&
                                             au2a => angstrom,&
                                             boltzmann,&
                                             joule,&
                                             massunit
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE tmc_analysis_types,              ONLY: &
       ana_type_default, ana_type_ice, ana_type_sym_xyz, atom_pairs_type, &
       dipole_moment_type, pair_correl_type, search_pair_in_list, &
       tmc_ana_density_create, tmc_ana_density_file_name, &
       tmc_ana_dipole_analysis_create, tmc_ana_dipole_moment_create, &
       tmc_ana_displacement_create, tmc_ana_env_create, &
       tmc_ana_pair_correl_create, tmc_ana_pair_correl_file_name, &
       tmc_analysis_env
  USE tmc_calculations,                ONLY: get_scaled_cell,&
                                             nearest_distance
  USE tmc_file_io,                     ONLY: analyse_files_close,&
                                             analyse_files_open,&
                                             expand_file_name_char,&
                                             expand_file_name_temp,&
                                             read_element_from_file,&
                                             write_dipoles_in_file
  USE tmc_stati,                       ONLY: TMC_STATUS_OK,&
                                             TMC_STATUS_WAIT_FOR_NEW_TASK,&
                                             tmc_default_restart_in_file_name,&
                                             tmc_default_restart_out_file_name,&
                                             tmc_default_trajectory_file_name,&
                                             tmc_default_unspecified_name
  USE tmc_tree_build,                  ONLY: allocate_new_sub_tree_node,&
                                             deallocate_sub_tree_node
  USE tmc_tree_types,                  ONLY: read_subtree_elem_unformated,&
                                             tree_type,&
                                             write_subtree_elem_unformated
  USE tmc_types,                       ONLY: tmc_atom_type,&
                                             tmc_param_type
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'tmc_analysis'

  PUBLIC :: tmc_read_ana_input
  PUBLIC :: analysis_init, do_tmc_analysis, analyze_file_configurations, finalize_tmc_analysis
  PUBLIC :: analysis_restart_print, analysis_restart_read

CONTAINS

! *****************************************************************************
!> \brief creates a new para environment for tmc analysis
!> \param tmc_ana_section ...
!> \param tmc_ana TMC analysis environment
!> \author Mandes 02.2013
! *****************************************************************************
  SUBROUTINE tmc_read_ana_input(tmc_ana_section, tmc_ana)
    TYPE(section_vals_type), POINTER         :: tmc_ana_section
    TYPE(tmc_analysis_env), POINTER          :: tmc_ana

    CHARACTER(len=*), PARAMETER :: routineN = 'tmc_read_ana_input', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: c_tmp
    CHARACTER(LEN=default_string_length), &
      POINTER                                :: charge_atm(:)
    INTEGER                                  :: i_tmp, ntot
    INTEGER, DIMENSION(3)                    :: nr_bins
    INTEGER, DIMENSION(:), POINTER           :: i_arr_tmp
    LOGICAL                                  :: explicit, explicit_key, flag
    REAL(KIND=dp), POINTER                   :: charge(:)
    TYPE(section_vals_type), POINTER         :: tmp_section

    NULLIFY(tmp_section, charge_atm, i_arr_tmp, charge)

    CPASSERT(ASSOCIATED(tmc_ana_section))
    CPASSERT(.NOT.ASSOCIATED(tmc_ana))


    CALL section_vals_get(tmc_ana_section,explicit=explicit)
    IF(explicit) THEN
      CALL tmc_ana_env_create(tmc_ana=tmc_ana)
      ! restarting
      CALL section_vals_val_get(tmc_ana_section,"RESTART",l_val=tmc_ana%restart)
      ! file name prefix
      CALL section_vals_val_get(tmc_ana_section,"PREFIX_ANA_FILES",&
             c_val=tmc_ana%out_file_prefix)
      IF(tmc_ana%out_file_prefix .NE. "")THEN
        tmc_ana%out_file_prefix = TRIM(tmc_ana%out_file_prefix)//"_"
      END IF

      ! density calculation
      CALL section_vals_val_get(tmc_ana_section,"DENSITY" ,explicit=explicit_key)
      IF(explicit_key)THEN
        CALL section_vals_val_get(tmc_ana_section,"DENSITY",i_vals=i_arr_tmp)

        IF(SIZE(i_arr_tmp(:)).EQ.3)THEN
          CALL cp_assert(ALL(i_arr_tmp(:).GT.0),&
                 cp_failure_level,cp_assertion_failed,routineP,&
                 "The amount of intervals in each direction has to be greater than 0.",&
                 failure=flag, only_ionode=.TRUE.)
          nr_bins(:) = i_arr_tmp(:)
        ELSE IF(SIZE(i_arr_tmp(:)).EQ.1) THEN
          CALL cp_assert(ALL(i_arr_tmp(:).GT.0),&
                 cp_failure_level,cp_assertion_failed,routineP,&
                 "The amount of intervals has to be greater than 0.",&
                 failure=flag, only_ionode=.TRUE.)
          nr_bins(:) = i_arr_tmp(1)
        ELSE IF(SIZE(i_arr_tmp(:)).EQ.0) THEN
          nr_bins(:) = 1
        ELSE
          CALL cp_assert(.FALSE.,&
                 cp_failure_level,cp_assertion_failed,routineP,&
                 "unknown amount of dimentions for the binning.",&
                 failure=flag, only_ionode=.TRUE.)
        END IF
        CALL tmc_ana_density_create(tmc_ana%density_3d, nr_bins)
      END IF

      ! radial distribution function calculation
      CALL section_vals_val_get(tmc_ana_section,"G_R" ,explicit=explicit_key)
      IF(explicit_key)THEN
        CALL section_vals_val_get(tmc_ana_section,"G_R",i_val=i_tmp)
        CALL tmc_ana_pair_correl_create(ana_pair_correl=tmc_ana%pair_correl, &
               nr_bins=i_tmp)
      END IF

      ! radial distribution function calculation
      CALL section_vals_val_get(tmc_ana_section,"CLASSICAL_DIPOLE_MOMENTS" ,explicit=explicit_key)
      IF(explicit_key)THEN
        ! charges for dipoles needed
        tmp_section => section_vals_get_subs_vals(tmc_ana_section,"CHARGE")
        CALL section_vals_get(tmp_section,explicit=explicit,n_repetition=i_tmp)
        IF (explicit) THEN
          ntot=0
          ALLOCATE(charge_atm(i_tmp))
          ALLOCATE(charge(i_tmp))
          CALL read_chrg_section(charge_atm,charge,tmp_section,ntot)
        ELSE 
          CALL cp_assert(.FALSE., cp_failure_level,cp_assertion_failed,routineP,&
                       "to calculate the classical cell dipole moment "// &
                       "the charges has to be specified",&
                       failure=flag)
        END IF
        
        CALL tmc_ana_dipole_moment_create(tmc_ana%dip_mom, charge_atm,charge,&
                                          tmc_ana%dim_per_elem)

        IF(ASSOCIATED(charge_atm)) DEALLOCATE(charge_atm)
        IF(ASSOCIATED(charge))     DEALLOCATE(charge)
      END IF

      ! dipole moment analysis
      CALL section_vals_val_get(tmc_ana_section,"DIPOLE_ANALYSIS" ,explicit=explicit_key)
      IF(explicit_key)THEN
        CALL tmc_ana_dipole_analysis_create(tmc_ana%dip_ana)
        CALL section_vals_val_get(tmc_ana_section,"DIPOLE_ANALYSIS",c_val=c_tmp)
        SELECT CASE (TRIM(c_tmp))
        CASE (TRIM(tmc_default_unspecified_name))
          tmc_ana%dip_ana%ana_type = ana_type_default
        CASE ("ICE")
          tmc_ana%dip_ana%ana_type = ana_type_ice
        CASE ("SYM_XYZ")
          tmc_ana%dip_ana%ana_type = ana_type_sym_xyz
        CASE DEFAULT
          CALL cp_assert(.FALSE.,&
                 cp_warning_level,cp_assertion_failed,routineP,&
                 'unknown analysis type "'//TRIM(c_tmp)//'" specified. '//&
                 " Set to default.",&
                 failure=flag, only_ionode=.TRUE.)
          tmc_ana%dip_ana%ana_type = ana_type_default
        END SELECT
      END IF

    END IF

    ! cell displacement (deviation)
    CALL section_vals_val_get(tmc_ana_section,"DEVIATION",l_val=flag)
    IF(flag) THEN
      CALL tmc_ana_displacement_create(ana_disp=tmc_ana%displace, &
             dim_per_elem=tmc_ana%dim_per_elem)
    END IF
  END SUBROUTINE tmc_read_ana_input

! *****************************************************************************
!> \brief initialize all the necessarry analysis structures
!> \param ana_env ...
!> \param nr_dim dimension of the pos, frc etc. array
!> \author Mandes 02.2013
! *****************************************************************************
  SUBROUTINE analysis_init(ana_env, nr_dim)
    TYPE(tmc_analysis_env), POINTER          :: ana_env
    INTEGER                                  :: nr_dim

    CHARACTER(len=*), PARAMETER :: routineN = 'analysis_init', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: tmp_cell_file, tmp_dip_file, &
                                                tmp_pos_file
    LOGICAL                                  :: failure

    failure = .FALSE.

    CPASSERT(ASSOCIATED(ana_env))
    CPASSERT(nr_dim>0)

    ana_env%nr_dim = nr_dim

    ! save file names
    tmp_pos_file  = ana_env%costum_pos_file_name
    tmp_cell_file = ana_env%costum_cell_file_name
    tmp_dip_file  = ana_env%costum_dip_file_name

    ! unset all filenames
    ana_env%costum_pos_file_name  = tmc_default_unspecified_name
    ana_env%costum_cell_file_name = tmc_default_unspecified_name
    ana_env%costum_dip_file_name  = tmc_default_unspecified_name

    ! set the neccessary files for ...
    ! density
    IF(ASSOCIATED(ana_env%density_3d)) THEN
      ana_env%costum_pos_file_name  = tmp_pos_file
      ana_env%costum_cell_file_name = tmp_cell_file
    END IF
    ! pair correlation
    IF(ASSOCIATED(ana_env%pair_correl)) THEN
      ana_env%costum_pos_file_name  = tmp_pos_file
      ana_env%costum_cell_file_name = tmp_cell_file
    END IF
    ! dipole moment
    IF(ASSOCIATED(ana_env%dip_mom)) THEN
      ana_env%costum_pos_file_name  = tmp_pos_file 
      ana_env%costum_cell_file_name = tmp_cell_file
    END IF
    ! dipole analysis
    IF(ASSOCIATED(ana_env%dip_ana)) THEN
      ana_env%costum_pos_file_name  = tmp_pos_file
      ana_env%costum_cell_file_name = tmp_cell_file
      ana_env%costum_dip_file_name  = tmp_dip_file
    END IF
    ! deviation / displacement
    IF(ASSOCIATED(ana_env%displace)) THEN
      ana_env%costum_pos_file_name  = tmp_pos_file
      ana_env%costum_cell_file_name = tmp_cell_file
    END IF


    ! init radial distribution function
    IF(ASSOCIATED(ana_env%pair_correl)) &
      CALL ana_pair_correl_init(ana_pair_correl=ana_env%pair_correl, &
                                atoms=ana_env%atoms, cell=ana_env%cell)
    ! init classical dipole moment calculations
    IF(ASSOCIATED(ana_env%dip_mom)) &
      CALL ana_dipole_moment_init(ana_dip_mom=ana_env%dip_mom, &
             atoms=ana_env%atoms)
  END SUBROUTINE analysis_init

! *****************************************************************************
!> \brief print analysis restart file
!> \param ana_env ...
!> \param 
!> \author Mandes 02.2013
! *****************************************************************************
  SUBROUTINE analysis_restart_print(ana_env)
    TYPE(tmc_analysis_env), POINTER          :: ana_env

    CHARACTER(len=*), PARAMETER :: routineN = 'analysis_restart_print', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: file_name, file_name_tmp, &
                                                restart_file_name
    INTEGER                                  :: file_ptr
    LOGICAL                                  :: failure, l_tmp

    failure = .FALSE.
    CPASSERT(ASSOCIATED(ana_env))
    CPASSERT(ASSOCIATED(ana_env%last_elem))
    IF(.NOT.ana_env%restart) RETURN

    WRITE(file_name, FMT='(I9.9)') ana_env%last_elem%nr
    file_name_tmp = TRIM(expand_file_name_temp(expand_file_name_char(&
                             TRIM(ana_env%out_file_prefix)//&
                             tmc_default_restart_out_file_name,&
                             "ana"), ana_env%temperature))
    restart_file_name = expand_file_name_char(file_name_tmp, &
                                              file_name)
    CALL open_file(file_name=restart_file_name, file_status="REPLACE", &
                   file_action="WRITE", file_form="UNFORMATTED", &
                   unit_number=file_ptr)
    WRITE(file_ptr) ana_env%temperature
    CALL write_subtree_elem_unformated(ana_env%last_elem, file_ptr)

    ! first mention the different kind of anlysis types initialized
    ! then the variables for each calculation type
    l_tmp = ASSOCIATED(ana_env%density_3d)
    WRITE(file_ptr) l_tmp
    IF(l_tmp)THEN
      WRITE(file_ptr)ana_env%density_3d%conf_counter, &
               ana_env%density_3d%nr_bins, &
               ana_env%density_3d%sum_vol, &
               ana_env%density_3d%sum_vol2, &
               ana_env%density_3d%sum_box_length, &
               ana_env%density_3d%sum_box_length2, &
               ana_env%density_3d%sum_density, &
               ana_env%density_3d%sum_dens2
    END IF

    l_tmp = ASSOCIATED(ana_env%pair_correl)
    WRITE(file_ptr) l_tmp
    IF(l_tmp)THEN
      WRITE(file_ptr)ana_env%pair_correl%conf_counter, &
               ana_env%pair_correl%nr_bins, &
               ana_env%pair_correl%step_lenght, &
               ana_env%pair_correl%pairs, &
               ana_env%pair_correl%g_r
    END IF

    l_tmp = ASSOCIATED(ana_env%dip_mom)
    WRITE(file_ptr) l_tmp
    IF(l_tmp)THEN
      WRITE(file_ptr)ana_env%dip_mom%conf_counter, &
               ana_env%dip_mom%charges, &
               ana_env%dip_mom%last_dip_cl
    END IF

    l_tmp = ASSOCIATED(ana_env%dip_ana)
    WRITE(file_ptr) l_tmp
    IF(l_tmp)THEN
      WRITE(file_ptr)ana_env%dip_ana%conf_counter, &
               ana_env%dip_ana%ana_type,&
               ana_env%dip_ana%mu2_pv_s, &
               ana_env%dip_ana%mu_psv, &
               ana_env%dip_ana%mu_pv, &
               ana_env%dip_ana%mu2_pv_mat, &
               ana_env%dip_ana%mu2_pv_mat
    END IF

    l_tmp = ASSOCIATED(ana_env%displace)
    WRITE(file_ptr) l_tmp
    IF(l_tmp)THEN
      WRITE(file_ptr)ana_env%displace%conf_counter, &
               ana_env%displace%disp
    END IF

    CALL close_file(unit_number=file_ptr)

    file_name_tmp = expand_file_name_char(TRIM(ana_env%out_file_prefix)//&
                      tmc_default_restart_in_file_name, "ana")
    file_name = expand_file_name_temp(file_name_tmp, &
                      ana_env%temperature)
    CALL open_file(file_name=file_name, &
                   file_action="WRITE", file_status="REPLACE", &
                   unit_number=file_ptr)
    WRITE(file_ptr,*)TRIM(restart_file_name)
    CALL close_file(unit_number=file_ptr)
  END SUBROUTINE analysis_restart_print

! *****************************************************************************
!> \brief read analysis restart file
!> \param ana_env ...
!> \param elem ...
!> \param 
!> \author Mandes 02.2013
! *****************************************************************************
  SUBROUTINE analysis_restart_read(ana_env, elem)
    TYPE(tmc_analysis_env), POINTER          :: ana_env
    TYPE(tree_type), POINTER                 :: elem

    CHARACTER(len=*), PARAMETER :: routineN = 'analysis_restart_read', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: file_name, file_name_tmp
    INTEGER                                  :: file_ptr
    LOGICAL                                  :: failure, l_tmp
    REAL(KIND=dp)                            :: temp

    failure = .FALSE.
    CPASSERT(ASSOCIATED(ana_env))
    CPASSERT(ASSOCIATED(elem))
    IF(.NOT.ana_env%restart) RETURN
 
    file_name_tmp = expand_file_name_char(TRIM(ana_env%out_file_prefix)//&
                      tmc_default_restart_in_file_name, "ana")
    file_name = expand_file_name_temp(file_name_tmp, &
                      ana_env%temperature)
    INQUIRE(FILE=file_name, EXIST=l_tmp)
    IF(l_tmp) THEN
      CALL open_file(file_name=file_name, file_status="OLD", &
                     file_action="READ", unit_number=file_ptr)
      READ(file_ptr,*) file_name_tmp
      CALL close_file(unit_number=file_ptr)

      CALL open_file(file_name=file_name_tmp, file_status="OLD", file_form="UNFORMATTED", &
                   file_action="READ", unit_number=file_ptr) 
      READ(file_ptr) temp
      CPASSERT(ana_env%temperature.EQ.temp)
      ana_env%last_elem => elem
      CALL read_subtree_elem_unformated(elem, file_ptr)

      ! first mention the different kind of anlysis types initialized
      ! then the variables for each calculation type
      READ(file_ptr) l_tmp
      CPASSERT(ASSOCIATED(ana_env%density_3d).EQV.l_tmp)
      IF(l_tmp)THEN
        READ(file_ptr) ana_env%density_3d%conf_counter, &
               ana_env%density_3d%nr_bins, &
               ana_env%density_3d%sum_vol, &
               ana_env%density_3d%sum_vol2, &
               ana_env%density_3d%sum_box_length, &
               ana_env%density_3d%sum_box_length2, &
               ana_env%density_3d%sum_density, &
               ana_env%density_3d%sum_dens2
      END IF
  
      READ(file_ptr) l_tmp
      CPASSERT(ASSOCIATED(ana_env%pair_correl).EQV.l_tmp)
      IF(l_tmp)THEN
        READ(file_ptr) ana_env%pair_correl%conf_counter, &
                 ana_env%pair_correl%nr_bins, &
                 ana_env%pair_correl%step_lenght, &
                 ana_env%pair_correl%pairs, &
                 ana_env%pair_correl%g_r
      END IF
  
      READ(file_ptr) l_tmp
      CPASSERT(ASSOCIATED(ana_env%dip_mom).EQV.l_tmp)
      IF(l_tmp)THEN
        READ(file_ptr) ana_env%dip_mom%conf_counter, &
                 ana_env%dip_mom%charges, &
                 ana_env%dip_mom%last_dip_cl
      END IF

      READ(file_ptr) l_tmp
      CPASSERT(ASSOCIATED(ana_env%dip_ana).EQV.l_tmp)
      IF(l_tmp)THEN
        READ(file_ptr) ana_env%dip_ana%conf_counter, &
                 ana_env%dip_ana%ana_type,&
                 ana_env%dip_ana%mu2_pv_s, &
                 ana_env%dip_ana%mu_psv, &
                 ana_env%dip_ana%mu_pv, &
                 ana_env%dip_ana%mu2_pv_mat, &
                 ana_env%dip_ana%mu2_pv_mat
      END IF

      READ(file_ptr) l_tmp
      CPASSERT(ASSOCIATED(ana_env%displace).EQV.l_tmp)
      IF(l_tmp)THEN
        READ(file_ptr)ana_env%displace%conf_counter, &
                ana_env%displace%disp
      END IF

      CALL close_file(unit_number=file_ptr)
      elem => NULL()
    END IF
  END SUBROUTINE analysis_restart_read

! *****************************************************************************
!> \brief call all the necessarry analysis routines
!>         analysis the previous element with the weight of the different 
!>        configuration numbers
!>        and stores the actual in the structur % last_elem
!>        afterwards the previous configuration can be deallocated (outside)
!> \param elem ...
!> \param ana_env ...
!> \param 
!> \author Mandes 02.2013
! *****************************************************************************
  SUBROUTINE do_tmc_analysis(elem, ana_env)
    TYPE(tree_type), POINTER                 :: elem
    TYPE(tmc_analysis_env), POINTER          :: ana_env

    CHARACTER(LEN=*), PARAMETER :: routineN = 'do_tmc_analysis', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, weight_act
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(3)              :: dip_tmp
    TYPE(tree_type), POINTER                 :: elem_tmp

    failure = .FALSE.

    CPASSERT(ASSOCIATED(elem))
    CPASSERT(ASSOCIATED(ana_env))

    ! start the timing
    CALL timeset(routineN,handle)
    IF(ASSOCIATED(ana_env%last_elem) .AND. &
       (ana_env%last_elem%nr .LT. elem%nr)) THEN
        weight_act = elem%nr - ana_env%last_elem%nr
        ! calculates the 3 dimensional distributed density
        IF(ASSOCIATED(ana_env%density_3d)) &
          CALL calc_density_3d(elem=ana_env%last_elem, &
                               weight=weight_act, atoms=ana_env%atoms, &
                               ana_env=ana_env)
        ! calculated the radial distribution function for each atom type
        IF(ASSOCIATED(ana_env%pair_correl)) &
          CALL calc_paircorrelation(elem=ana_env%last_elem, weight=weight_act,&
                                    atoms=ana_env%atoms, ana_env=ana_env)
        ! calculates the classical dipole moments
        IF(ASSOCIATED(ana_env%dip_mom)) &
          CALL calc_dipole_moment(elem=ana_env%last_elem, weight=weight_act,&
                                  ana_env=ana_env)
        ! calculates the dipole moments analysis and dielectric constant
        IF(ASSOCIATED(ana_env%dip_ana)) THEN
          ! in symmetric case use also the dipoles 
          !   (-x,y,z) .. .. (-x,-y,z).... (-x,-y-z) all have the same energy
          IF (ana_env%dip_ana%ana_type .EQ. ana_type_sym_xyz) THEN
            ! (-x,y,z)
            ana_env%last_elem%dipole(1) = -ana_env%last_elem%dipole(1)
            dip_tmp(:) = ana_env%last_elem%dipole(:)
            IF(ASSOCIATED(ana_env%dip_mom)) &
              ana_env%dip_mom%last_dip_cl(1) = -ana_env%dip_mom%last_dip_cl(1)
            CALL calc_dipole_analysis(elem=ana_env%last_elem, weight=weight_act,&
                                      ana_env=ana_env)
            ! (-x,-y,z)
            ana_env%last_elem%dipole(:) = dip_tmp(:)
            ana_env%last_elem%dipole(2) = -ana_env%last_elem%dipole(2)
            dip_tmp(:) = ana_env%last_elem%dipole(:)
            IF(ASSOCIATED(ana_env%dip_mom)) &
              ana_env%dip_mom%last_dip_cl(2) = -ana_env%dip_mom%last_dip_cl(2)
            CALL calc_dipole_analysis(elem=ana_env%last_elem, weight=weight_act,&
                                      ana_env=ana_env)
            ! (-x,-y,-z)
            ana_env%last_elem%dipole(:) = dip_tmp(:)
            ana_env%last_elem%dipole(3) = -ana_env%last_elem%dipole(3)
            dip_tmp(:) = ana_env%last_elem%dipole(:)
            IF(ASSOCIATED(ana_env%dip_mom)) &
              ana_env%dip_mom%last_dip_cl(3) = -ana_env%dip_mom%last_dip_cl(3)
            CALL calc_dipole_analysis(elem=ana_env%last_elem, weight=weight_act,&
                                      ana_env=ana_env)
            ! (x,-y,-z)
            ana_env%last_elem%dipole(:) = dip_tmp(:)
            ana_env%last_elem%dipole(1) = -ana_env%last_elem%dipole(1)
            dip_tmp(:) = ana_env%last_elem%dipole(:)
            IF(ASSOCIATED(ana_env%dip_mom)) &
              ana_env%dip_mom%last_dip_cl(1) = -ana_env%dip_mom%last_dip_cl(1)
            CALL calc_dipole_analysis(elem=ana_env%last_elem, weight=weight_act,&
                                      ana_env=ana_env)
            ! (x,y,-z)
            ana_env%last_elem%dipole(:) = dip_tmp(:)
            ana_env%last_elem%dipole(2) = -ana_env%last_elem%dipole(2)
            dip_tmp(:) = ana_env%last_elem%dipole(:)
            IF(ASSOCIATED(ana_env%dip_mom)) &
              ana_env%dip_mom%last_dip_cl(2) = -ana_env%dip_mom%last_dip_cl(2)
            CALL calc_dipole_analysis(elem=ana_env%last_elem, weight=weight_act,&
                                      ana_env=ana_env)
            ! (-x,y,-z)
            ana_env%last_elem%dipole(:) = dip_tmp(:)
            ana_env%last_elem%dipole(1) = -ana_env%last_elem%dipole(1)
            dip_tmp(:) = ana_env%last_elem%dipole(:)
            IF(ASSOCIATED(ana_env%dip_mom)) &
              ana_env%dip_mom%last_dip_cl(1) = -ana_env%dip_mom%last_dip_cl(1)
            CALL calc_dipole_analysis(elem=ana_env%last_elem, weight=weight_act,&
                                      ana_env=ana_env)
            ! (x,-y,z)
            ana_env%last_elem%dipole(:) = dip_tmp(:)
            ana_env%last_elem%dipole(:) = -ana_env%last_elem%dipole(:)
            dip_tmp(:) = ana_env%last_elem%dipole(:)
            IF(ASSOCIATED(ana_env%dip_mom)) &
              ana_env%dip_mom%last_dip_cl(:) = -ana_env%dip_mom%last_dip_cl(:)
            CALL calc_dipole_analysis(elem=ana_env%last_elem, weight=weight_act,&
                                      ana_env=ana_env)
            ! back to (x,y,z)
            ana_env%last_elem%dipole(:) = dip_tmp(:)
            ana_env%last_elem%dipole(2) = -ana_env%last_elem%dipole(2)
            dip_tmp(:) = ana_env%last_elem%dipole(:)
            IF(ASSOCIATED(ana_env%dip_mom)) &
              ana_env%dip_mom%last_dip_cl(2) = -ana_env%dip_mom%last_dip_cl(2)
          END IF
          CALL calc_dipole_analysis(elem=ana_env%last_elem, weight=weight_act,&
                                    ana_env=ana_env)
          CALL print_act_dipole_analysis(elem=ana_env%last_elem, &
                                         ana_env=ana_env)
        END IF

        ! calculates the cell displacement from last cell
        IF(ASSOCIATED(ana_env%displace).AND.weight_act.GT.0) &
          CALL calc_displacement(elem=elem, ana_env=ana_env)
    END IF
    ! swap elem with last elem, to delete original last element and store the actual one
    elem_tmp => ana_env%last_elem
    ana_env%last_elem => elem
    elem => elem_tmp
    ! end the timing
    CALL timestop(handle)
  END SUBROUTINE do_tmc_analysis

! *****************************************************************************
!> \brief call all the necessarry analysis printing routines
!> \param ana_env ...
!> \param 
!> \author Mandes 02.2013
! *****************************************************************************
  SUBROUTINE finalize_tmc_analysis(ana_env)
    TYPE(tmc_analysis_env), POINTER          :: ana_env

    CHARACTER(LEN=*), PARAMETER :: routineN = 'finalize_tmc_analysis', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CPASSERT(ASSOCIATED(ana_env))

    ! start the timing
    CALL timeset(routineN,handle)
    IF(ASSOCIATED(ana_env%density_3d)) THEN
      IF(ana_env%density_3d%conf_counter.GT.0) &
        CALL print_density_3d(ana_env=ana_env)
    END IF
    IF(ASSOCIATED(ana_env%pair_correl)) THEN
      IF(ana_env%pair_correl%conf_counter.GT.0) &
        CALL print_paircorrelation(ana_env=ana_env)
    END IF
    IF(ASSOCIATED(ana_env%dip_mom)) THEN
      IF(ana_env%dip_mom%conf_counter.GT.0) &
        CALL print_dipole_moment(ana_env)
    END IF
    IF(ASSOCIATED(ana_env%dip_ana)) THEN
      IF(ana_env%dip_ana%conf_counter.GT.0) &
        CALL print_dipole_analysis(ana_env)
    END IF
    IF(ASSOCIATED(ana_env%displace)) THEN
      IF(ana_env%displace%conf_counter.GT.0) &
        CALL print_average_displacement(ana_env)
    END IF

    ! end the timing
    CALL timestop(handle)
  END SUBROUTINE finalize_tmc_analysis

! *****************************************************************************
!> \brief read the files and analyze the configurations
!> \param start_id ...
!> \param end_id ...
!> \param dir_ind ...
!> \param ana_env ...
!> \param tmc_params ...
!> \author Mandes 03.2013
! *****************************************************************************
  SUBROUTINE analyze_file_configurations(start_id, end_id, dir_ind, &
                                         ana_env, tmc_params)
    INTEGER                                  :: start_id, end_id
    INTEGER, OPTIONAL                        :: dir_ind
    TYPE(tmc_analysis_env), POINTER          :: ana_env
    TYPE(tmc_param_type), POINTER            :: tmc_params

    CHARACTER(LEN=*), PARAMETER :: routineN = 'analyze_file_configurations', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: conf_nr, handle, nr_dim, stat
    TYPE(tree_type), POINTER                 :: elem

    NULLIFY(elem)
    conf_nr = -1
    stat = TMC_STATUS_WAIT_FOR_NEW_TASK
    CPASSERT(ASSOCIATED(ana_env))
    CPASSERT(ASSOCIATED(tmc_params))

    ! start the timing
    CALL timeset(routineN,handle)

    ! open the files
    CALL analyse_files_open(tmc_ana=ana_env, stat=stat, dir_ind=dir_ind)
    ! set the existence of exact dipoles (from file)
    IF(ana_env%id_dip .GT. 0) THEN
      tmc_params%print_dipole = .TRUE.
    ELSE
      tmc_params%print_dipole = .FALSE.
    END IF

    ! allocate the actual element structure
    CALL allocate_new_sub_tree_node(tmc_params=tmc_params, next_el=elem,&
           nr_dim=ana_env%nr_dim)

    IF(ASSOCIATED(ana_env%last_elem)) conf_nr = ana_env%last_elem%nr
    nr_dim = SIZE(elem%pos)

    IF(stat .EQ. TMC_STATUS_OK) THEN
      conf_loop: DO
        CALL read_element_from_file(elem=elem, tmc_ana=ana_env, conf_nr=conf_nr, &
                                    stat=stat)
        IF(stat.EQ.TMC_STATUS_WAIT_FOR_NEW_TASK) THEN
          CALL deallocate_sub_tree_node(tree_elem=elem)
          EXIT conf_loop
        END IF
        ! if we want just a certain part of the trajectory
        IF(start_id.LT.0 .OR. conf_nr.GE.start_id) THEN
          IF(end_id.LT.0 .OR. conf_nr.LE.end_id) THEN
            ! do the analysis calculations
            CALL do_tmc_analysis(elem=elem, ana_env=ana_env)
          END IF
        END IF
  
        ! clean temporary element (already analyzed)
        IF(ASSOCIATED(elem)) THEN
          CALL deallocate_sub_tree_node(tree_elem=elem)
        END IF
        ! if there was no previous element, create a new temp element to write in
        IF(.NOT.ASSOCIATED(elem)) &
          CALL allocate_new_sub_tree_node(tmc_params=tmc_params, next_el=elem, &
                 nr_dim=nr_dim)
      END DO conf_loop
    END IF
    ! close the files
    CALL analyse_files_close(tmc_ana=ana_env)

    IF(ASSOCIATED(elem)) &
      CALL deallocate_sub_tree_node(tree_elem=elem)

    ! end the timing
    CALL timestop(handle)
  END SUBROUTINE analyze_file_configurations


  !============================================================================
  ! density calculations
  !============================================================================

! *****************************************************************************
!> \brief calculates the density in rectantangulares 
!>        defined by the number of bins in each direction
!> \param elem ...
!> \param weight ...
!> \param atoms ...
!> \param ana_env ...
!> \param 
!> \author Mandes 02.2013
! *****************************************************************************
  SUBROUTINE calc_density_3d(elem, weight, atoms, ana_env)
    TYPE(tree_type), POINTER                 :: elem
    INTEGER                                  :: weight
    TYPE(tmc_atom_type), DIMENSION(:), &
      POINTER                                :: atoms
    TYPE(tmc_analysis_env), POINTER          :: ana_env

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_density_3d', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: file_name, file_name_tmp
    INTEGER                                  :: atom, bin_x, bin_y, bin_z, &
                                                file_ptr, handle
    LOGICAL                                  :: flag
    REAL(KIND=dp)                            :: mass_total, r_tmp, vol_cell, &
                                                vol_sub_box
    REAL(KIND=dp), DIMENSION(3)              :: atom_pos, cell_size, &
                                                interval_size
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: mass_bin

    NULLIFY(mass_bin)

    CPASSERT(ASSOCIATED(elem))
    CPASSERT(ASSOCIATED(elem%pos))
    CPASSERT(weight.GT.0)
    CPASSERT(ASSOCIATED(atoms))
    CPASSERT(ASSOCIATED(ana_env))
    CPASSERT(ASSOCIATED(ana_env%cell))
    CPASSERT(ASSOCIATED(ana_env%density_3d))
    CPASSERT(ASSOCIATED(ana_env%density_3d%sum_density))
    CPASSERT(ASSOCIATED(ana_env%density_3d%sum_dens2))

    ! start the timing
    CALL timeset(routineN,handle)

    atom_pos(:)      = 0.0_dp
    cell_size(:)     = 0.0_dp
    interval_size(:) = 0.0_dp
    mass_total       = 0.0_dp

    bin_x = SIZE(ana_env%density_3d%sum_density(:,1,1))
    bin_y = SIZE(ana_env%density_3d%sum_density(1,:,1))
    bin_z = SIZE(ana_env%density_3d%sum_density(1,1,:))
    ALLOCATE(mass_bin(bin_x, bin_y, bin_z))
    mass_bin(:,:,:) = 0.0_dp

    ! if NPT -> box_scale.NE.1.0 use the scaled cell 
    ! ATTENTION then the sub box middle points are not correct in the output
    !  espacially if we use multiple sub boxes
    CALL get_scaled_cell(cell=ana_env%cell, box_scale=elem%box_scale, &
                         abc=cell_size, vol=vol_cell)
    ! volume summed over configurations for average volume [A]
    ana_env%density_3d%sum_vol = ana_env%density_3d%sum_vol + &
                                   vol_cell*(au2a)**3 *weight
    ana_env%density_3d%sum_vol2 = ana_env%density_3d%sum_vol2 + &
                                   (vol_cell*(au2a)**3)**2 *weight

    ana_env%density_3d%sum_box_length(:) = ana_env%density_3d%sum_box_length(:)&
                                          + cell_size(:)*(au2a)*weight
    ana_env%density_3d%sum_box_length2(:) = ana_env%density_3d%sum_box_length2(:)&
                                          + (cell_size(:)*(au2a))**2 *weight

    ! sub interval length
    interval_size(1) = cell_size(1)/REAL(bin_x,dp)
    interval_size(2) = cell_size(2)/REAL(bin_y,dp)
    interval_size(3) = cell_size(3)/REAL(bin_z,dp)

    ! volume in [cm^3]
    vol_cell = vol_cell* (au2a*1E-8)**3
    vol_sub_box = interval_size(1)*interval_size(2)*interval_size(3)* &
                  (au2a*1E-8)**3

    ! count every atom
    DO atom=1, SIZE(elem%pos), ana_env%dim_per_elem
      
      atom_pos(:) = elem%pos(atom:atom+2)
      ! fold into box
      CALL get_scaled_cell(cell=ana_env%cell, box_scale=elem%box_scale, &
                           vec=atom_pos)
      ! shifts the box to positive values (before 0,0,0 is the center)
      atom_pos(:) = atom_pos(:)+0.5_dp*cell_size(:)
      ! calculate the index of the sub box
      bin_x = INT(atom_pos(1)/interval_size(1))+1
      bin_y = INT(atom_pos(2)/interval_size(2))+1
      bin_z = INT(atom_pos(3)/interval_size(3))+1
      CPASSERT(bin_x.GT.0.AND.bin_y.GT.0.AND.bin_z.GT.0)
      CPASSERT(bin_x.LE.SIZE(ana_env%density_3d%sum_density(:,1,1)))
      CPASSERT(bin_y.LE.SIZE(ana_env%density_3d%sum_density(1,:,1)))
      CPASSERT(bin_z.LE.SIZE(ana_env%density_3d%sum_density(1,1,:)))

      ! sum mass in [g] (in bins and total)
      mass_bin(bin_x,bin_y,bin_z) = mass_bin(bin_x,bin_y,bin_z) + &
        atoms(INT(atom/REAL(ana_env%dim_per_elem,KIND=dp))+1)%mass/massunit*1000*a_mass
      mass_total = mass_total + &
        atoms(INT(atom/REAL(ana_env%dim_per_elem,KIND=dp))+1)%mass/massunit*1000*a_mass
      !mass_bin(bin_x,bin_y,bin_z) = mass_bin(bin_x,bin_y,bin_z) + &
      !  atoms(INT(atom/REAL(ana_env%dim_per_elem,KIND=dp))+1)%mass/&
      !     massunit/n_avogadro
      !mass_total = mass_total + &
      !  atoms(INT(atom/REAL(ana_env%dim_per_elem,KIND=dp))+1)%mass/&
      !     massunit/n_avogadro
    END DO
    ! check total cell density
    r_tmp = mass_total/vol_cell - SUM(mass_bin(:,:,:))/vol_sub_box/SIZE(mass_bin(:,:,:))
    CPASSERT(ABS(r_tmp).LT.1E-5)

    ! calculate density (mass per volume) and sum up for average value
    ana_env%density_3d%sum_density(:,:,:) = ana_env%density_3d%sum_density(:,:,:) +&
                                             weight* mass_bin(:,:,:)/vol_sub_box

    ! calculate density squared ( (mass per volume)^2 ) for variance and sum up for average value
    ana_env%density_3d%sum_dens2(:,:,:) = ana_env%density_3d%sum_dens2(:,:,:) +&
      weight* (mass_bin(:,:,:)/vol_sub_box)**2

    ana_env%density_3d%conf_counter = ana_env%density_3d%conf_counter + weight

    ! print out the actual and average density in file
    IF(ana_env%density_3d%print_dens) THEN
      file_name_tmp = expand_file_name_temp(TRIM(ana_env%out_file_prefix)//&
                        tmc_default_trajectory_file_name, &
                        ana_env%temperature)
      file_name = TRIM(expand_file_name_char(file_name_tmp,&
                       "dens"))
      INQUIRE(FILE=file_name, EXIST=flag)
      CALL open_file(file_name=file_name, file_status="UNKNOWN", &
                     file_action="WRITE", file_position="APPEND", &
                     unit_number=file_ptr)
      IF(.NOT.flag) &
        WRITE(file_ptr,FMT='(A8,11A20)')"# conf_nr","dens_act[g/cm^3]",&
                   "dens_average[g/cm^3]","density_variance", &
                   "averages:volume", "box_lenth_x", "box_lenth_y","box_lenth_z", &
                   "variances:volume", "box_lenth_x", "box_lenth_y","box_lenth_z"
      WRITE(file_ptr,FMT="(I8,11F20.10)") ana_env%density_3d%conf_counter+1-weight, &
        SUM(mass_bin(:,:,:))/vol_sub_box/SIZE(mass_bin(:,:,:)), &
        SUM(ana_env%density_3d%sum_density(:,:,:))/&
         SIZE(ana_env%density_3d%sum_density(:,:,:))/&
         REAL(ana_env%density_3d%conf_counter, KIND=dp), &
        SUM(ana_env%density_3d%sum_dens2(:,:,:))/&
         SIZE(ana_env%density_3d%sum_dens2(:,:,:))/&
         REAL(ana_env%density_3d%conf_counter, KIND=dp)-&
         (SUM(ana_env%density_3d%sum_density(:,:,:))/&
         SIZE(ana_env%density_3d%sum_density(:,:,:))/&
         REAL(ana_env%density_3d%conf_counter, KIND=dp))**2,&
        ana_env%density_3d%sum_vol/&
         REAL(ana_env%density_3d%conf_counter, KIND=dp),&
        ana_env%density_3d%sum_box_length(:)/&
         REAL(ana_env%density_3d%conf_counter, KIND=dp),&
        ana_env%density_3d%sum_vol2/&
         REAL(ana_env%density_3d%conf_counter, KIND=dp) - &
         (ana_env%density_3d%sum_vol/&
         REAL(ana_env%density_3d%conf_counter, KIND=dp))**2,&
        ana_env%density_3d%sum_box_length2(:)/&
         REAL(ana_env%density_3d%conf_counter, KIND=dp) - &
         (ana_env%density_3d%sum_box_length(:)/&
         REAL(ana_env%density_3d%conf_counter, KIND=dp))**2
      CALL close_file(unit_number=file_ptr)
    END IF
    
    DEALLOCATE(mass_bin)
    ! end the timing
    CALL timestop(handle)
  END SUBROUTINE calc_density_3d

! *****************************************************************************
!> \brief print the density in rectantangulares 
!>        defined by the number of bins in each direction
!> \param ana_env ...
!> \param 
!> \author Mandes 02.2013
! *****************************************************************************
  SUBROUTINE print_density_3d(ana_env)
    TYPE(tmc_analysis_env), POINTER          :: ana_env

    CHARACTER(LEN=*), PARAMETER :: fmt_my = '(T2,A,"| ",A,T41,A40)', &
      plabel = "TMC_ANA", routineN = 'print_density_3d', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: file_name, file_name_vari
    INTEGER                                  :: bin_x, bin_y, bin_z, &
                                                file_ptr_dens, file_ptr_vari, &
                                                handle, i, j, k
    REAL(KIND=dp), DIMENSION(3)              :: cell_size, interval_size

    CPASSERT(ASSOCIATED(ana_env))
    CPASSERT(ASSOCIATED(ana_env%density_3d))
    CPASSERT(ASSOCIATED(ana_env%density_3d%sum_density))
    CPASSERT(ASSOCIATED(ana_env%density_3d%sum_dens2))

    ! start the timing
    CALL timeset(routineN,handle)

    file_name=""
    file_name_vari=""

    bin_x = SIZE(ana_env%density_3d%sum_density(:,1,1))
    bin_y = SIZE(ana_env%density_3d%sum_density(1,:,1))
    bin_z = SIZE(ana_env%density_3d%sum_density(1,1,:))
    CALL get_cell(cell=ana_env%cell, abc=cell_size)
    interval_size(1) = cell_size(1)/REAL(bin_x,KIND=dp)*au2a 
    interval_size(2) = cell_size(2)/REAL(bin_y,KIND=dp)*au2a 
    interval_size(3) = cell_size(3)/REAL(bin_z,KIND=dp)*au2a 
    

    file_name = expand_file_name_temp(TRIM(ana_env%out_file_prefix)//&
                  tmc_ana_density_file_name, &
                  ana_env%temperature)
    CALL open_file(file_name=file_name, file_status="REPLACE", &
                   file_action="WRITE", file_position="APPEND", &
                   unit_number=file_ptr_dens)
    WRITE(file_ptr_dens,FMT='(A,1X,I0,1X,A,3(I0,1X),1X,A,1X,3F10.5)')&
      "# configurations", ana_env%density_3d%conf_counter, "bins", &
      ana_env%density_3d%nr_bins, "interval size", interval_size(:)
    WRITE(file_ptr_dens,FMT='(A,3A10,A20)')"#"," x [A] "," y [A] "," z [A] "," density [g/cm^3] "

    file_name_vari = expand_file_name_temp(expand_file_name_char(&
                       TRIM(ana_env%out_file_prefix)//&
                       tmc_ana_density_file_name,"vari"),&
                       ana_env%temperature)
    CALL open_file(file_name=file_name_vari, file_status="REPLACE", &
                   file_action="WRITE", file_position="APPEND", &
                   unit_number=file_ptr_vari)
    WRITE(file_ptr_vari,FMT='(A,1X,I0,1X,A,3(I0,1X),1X,A,1X,3F10.5)')&
      "# configurations", ana_env%density_3d%conf_counter, "bins", &
      ana_env%density_3d%nr_bins, "interval size", interval_size(:)
    WRITE(file_ptr_vari,FMT='(A,3A10,A20)')"#"," x [A] "," y [A] "," z [A] "," variance"

    DO i=1, SIZE(ana_env%density_3d%sum_density(:,1,1))
       DO j=1, SIZE(ana_env%density_3d%sum_density(1,:,1))
          DO k=1, SIZE(ana_env%density_3d%sum_density(1,1,:))
             WRITE(file_ptr_dens,FMT='(3F10.2,F20.10)') &
               (i-0.5_dp)*interval_size(1), (j-0.5_dp)*interval_size(2), (k-0.5_dp)*interval_size(3), &
               ana_env%density_3d%sum_density(i,j,k)/REAL(ana_env%density_3d%conf_counter, KIND=dp)
             WRITE(file_ptr_vari,FMT='(3F10.2,F20.10)') &
               (i-0.5_dp)*interval_size(1), (j-0.5_dp)*interval_size(2), (k-0.5_dp)*interval_size(3), &
               ana_env%density_3d%sum_dens2(i,j,k)/REAL(ana_env%density_3d%conf_counter, KIND=dp)-&
               (ana_env%density_3d%sum_density(i,j,k)/REAL(ana_env%density_3d%conf_counter, KIND=dp))**2
          END DO
        END DO
    END DO
    CALL close_file(unit_number=file_ptr_dens)
    CALL close_file(unit_number=file_ptr_vari)

    WRITE(ana_env%io_unit,FMT="(/,T2,A)") REPEAT("-",79)
    WRITE(ana_env%io_unit,FMT="(T2,A,T35,A,T80,A)")"-","density calculation","-"
    WRITE(ana_env%io_unit,FMT=fmt_my)plabel,"temperature ", cp_to_string(ana_env%temperature)
    WRITE(ana_env%io_unit,FMT=fmt_my)plabel,"used configurations", &
      cp_to_string(REAL(ana_env%density_3d%conf_counter, KIND=dp))
    WRITE(ana_env%io_unit,FMT=fmt_my)plabel,"average volume", &
      cp_to_string(ana_env%density_3d%sum_vol/&
        REAL(ana_env%density_3d%conf_counter, KIND=dp))
    WRITE(ana_env%io_unit,FMT=fmt_my)plabel,"average density in the cell: ", &
      cp_to_string(SUM(ana_env%density_3d%sum_density(:,:,:))/&
      SIZE(ana_env%density_3d%sum_density(:,:,:))/&
      REAL(ana_env%density_3d%conf_counter, KIND=dp))
    WRITE(ana_env%io_unit,FMT=fmt_my)plabel,"density variance:", &
      cp_to_string(SUM(ana_env%density_3d%sum_dens2(:,:,:))/&
         SIZE(ana_env%density_3d%sum_dens2(:,:,:))/&
         REAL(ana_env%density_3d%conf_counter, KIND=dp)-&
         (SUM(ana_env%density_3d%sum_density(:,:,:))/&
         SIZE(ana_env%density_3d%sum_density(:,:,:))/&
         REAL(ana_env%density_3d%conf_counter, KIND=dp))**2)
    WRITE(ana_env%io_unit,FMT="(/,T2,A)") REPEAT("-",79)
    IF(ana_env%print_test_output) &
      WRITE(ana_env%io_unit,*)"TMC|ANALYSIS_CELL_DENSITY_X= ", &
        SUM(ana_env%density_3d%sum_density(:,:,:))/&
        SIZE(ana_env%density_3d%sum_density(:,:,:))/&
        REAL(ana_env%density_3d%conf_counter, KIND=dp)
    ! end the timing
    CALL timestop(handle)
  END SUBROUTINE print_density_3d

  !============================================================================
  ! radial distribution function
  !============================================================================

! *****************************************************************************
!> \brief init radial distribution function structures
!> \param ana_pair_correl ...
!> \param atoms ...
!> \param cell ...
!> \param
!> \author Mandes 02.2013
! *****************************************************************************
  SUBROUTINE ana_pair_correl_init(ana_pair_correl, atoms, cell)
    TYPE(pair_correl_type), POINTER          :: ana_pair_correl
    TYPE(tmc_atom_type), DIMENSION(:), &
      POINTER                                :: atoms
    TYPE(cell_type), POINTER                 :: cell

    CHARACTER(LEN=*), PARAMETER :: routineN = 'ana_pair_correl_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: counter, f_n, handle, list, &
                                                list_ind, s_n
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(3)              :: cell_size
    TYPE(atom_pairs_type), DIMENSION(:), &
      POINTER                                :: pairs_tmp

    failure = .FALSE.

    CPASSERT(ASSOCIATED(ana_pair_correl))
    CPASSERT(.NOT.ASSOCIATED(ana_pair_correl%g_r))
    CPASSERT(.NOT.ASSOCIATED(ana_pair_correl%pairs))
    CPASSERT(ASSOCIATED(atoms))
    CPASSERT(SIZE(atoms).GT.1)
    CPASSERT(ASSOCIATED(cell))

    ! start the timing
    CALL timeset(routineN,handle)

    CALL get_cell(cell=cell, abc=cell_size)
    IF(ana_pair_correl%nr_bins.LE.0) THEN
      ana_pair_correl%nr_bins= CEILING(MAXVAL(cell_size(:))/2.0_dp/(0.03/au2a))
    END IF
    ana_pair_correl%step_lenght = MAXVAL(cell_size(:))/2.0_dp/&
                                    ana_pair_correl%nr_bins
    ana_pair_correl%conf_counter = 0

    counter = 1
    ! initialise the atom pairs
    ALLOCATE(pairs_tmp(SIZE(atoms)))
    DO f_n=1, SIZE(atoms)
      DO s_n=f_n+1, SIZE(atoms)
        ! search if atom pair is already selected
        list_ind = search_pair_in_list(pair_list=pairs_tmp, n1=atoms(f_n)%name,&
                                       n2=atoms(s_n)%name, list_end=counter-1)
        ! add to list
        IF(list_ind.LT.0) THEN
          pairs_tmp(counter)%f_n = atoms(f_n)%name
          pairs_tmp(counter)%s_n = atoms(s_n)%name
          pairs_tmp(counter)%pair_count = 1
          counter = counter +1
        ELSE
          pairs_tmp(list_ind)%pair_count=pairs_tmp(list_ind)%pair_count+1
        END IF
      END DO
    END DO

    ALLOCATE(ana_pair_correl%pairs(counter-1))
    DO list=1, counter-1
      ana_pair_correl%pairs(list)%f_n = pairs_tmp(list)%f_n
      ana_pair_correl%pairs(list)%s_n = pairs_tmp(list)%s_n
      ana_pair_correl%pairs(list)%pair_count = pairs_tmp(list)%pair_count
    END DO
    DEALLOCATE(pairs_tmp)

    ALLOCATE(ana_pair_correl%g_r(SIZE(ana_pair_correl%pairs(:)),ana_pair_correl%nr_bins))
    ana_pair_correl%g_r = 0.0_dp
    ! end the timing
    CALL timestop(handle)
  END SUBROUTINE ana_pair_correl_init

! *****************************************************************************
!> \brief calculates the radial distribution function  
!> \param elem ...
!> \param weight ...
!> \param atoms ...
!> \param ana_env ...
!> \param 
!> \author Mandes 02.2013
! *****************************************************************************
  SUBROUTINE calc_paircorrelation(elem, weight, atoms, ana_env)
    TYPE(tree_type), POINTER                 :: elem
    INTEGER                                  :: weight
    TYPE(tmc_atom_type), DIMENSION(:), &
      POINTER                                :: atoms
    TYPE(tmc_analysis_env), POINTER          :: ana_env

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_paircorrelation', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ind, j, pair_ind
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: dist
    REAL(KIND=dp), DIMENSION(3)              :: cell_size

    failure = .FALSE.
    CPASSERT(ASSOCIATED(elem))
    CPASSERT(ASSOCIATED(elem%pos))
    CPASSERT(ALL(elem%box_scale(:).GT.0.0_dp))
    CPASSERT(weight.GT.0)
    CPASSERT(ASSOCIATED(atoms))
    CPASSERT(ASSOCIATED(ana_env))
    CPASSERT(ASSOCIATED(ana_env%cell))
    CPASSERT(ASSOCIATED(ana_env%pair_correl))
    CPASSERT(ASSOCIATED(ana_env%pair_correl%g_r))
    CPASSERT(ASSOCIATED(ana_env%pair_correl%pairs))

    ! start the timing
    CALL timeset(routineN,handle)

    dist = -1.0_dp

    first_elem_loop:DO i=1, SIZE(elem%pos), ana_env%dim_per_elem
      second_elem_loop:DO j=i+3, SIZE(elem%pos), ana_env%dim_per_elem
        dist = nearest_distance(x1=elem%pos(i:i+ana_env%dim_per_elem-1), &
                                x2=elem%pos(j:j+ana_env%dim_per_elem-1), &
                                cell=ana_env%cell, box_scale=elem%box_scale)
        ind = CEILING(dist/ana_env%pair_correl%step_lenght)
        IF(ind.LE.ana_env%pair_correl%nr_bins)THEN
          pair_ind = search_pair_in_list(pair_list=ana_env%pair_correl%pairs, &
                       n1=atoms(INT(i/REAL(ana_env%dim_per_elem,KIND=dp))+1)%name, &
                       n2=atoms(INT(j/REAL(ana_env%dim_per_elem,KIND=dp))+1)%name)
          CPASSERT(pair_ind.GT.0)
          ana_env%pair_correl%g_r(pair_ind,ind) = &
            ana_env%pair_correl%g_r(pair_ind,ind) + weight
        END IF
      END DO second_elem_loop
    END DO first_elem_loop
    ana_env%pair_correl%conf_counter = ana_env%pair_correl%conf_counter + weight
    CALL get_cell(cell=ana_env%cell, abc=cell_size)
    ana_env%pair_correl%sum_box_scale = ana_env%pair_correl%sum_box_scale + &
                                (elem%box_scale(:) * weight)
    ! end the timing
    CALL timestop(handle)
  END SUBROUTINE calc_paircorrelation

! *****************************************************************************
!> \brief print the radial distribution function for each pair of atoms 
!> \param ana_env ...
!> \param 
!> \author Mandes 02.2013
! *****************************************************************************
  SUBROUTINE print_paircorrelation(ana_env)
    TYPE(tmc_analysis_env), POINTER          :: ana_env

    CHARACTER(LEN=*), PARAMETER :: routineN = 'print_paircorrelation', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: file_name
    INTEGER                                  :: bin, file_ptr, handle, pair
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: aver_box_scale(3), vol, voldr
    REAL(KIND=dp), DIMENSION(3)              :: cell_size

    failure = .FALSE.

    CPASSERT(ASSOCIATED(ana_env))
    CPASSERT(ASSOCIATED(ana_env%pair_correl))

    ! start the timing
    CALL timeset(routineN,handle)

    CALL get_cell(cell=ana_env%cell, abc=cell_size)
    aver_box_scale(:) = ana_env%pair_correl%sum_box_scale(:) / ana_env%pair_correl%conf_counter
    vol = (cell_size(1)*aver_box_scale(1))*&
          (cell_size(2)*aver_box_scale(2))*&
          (cell_size(3)*aver_box_scale(3))

    DO pair=1, SIZE(ana_env%pair_correl%pairs)
      file_name = expand_file_name_temp(TRIM(ana_env%out_file_prefix)//&
                    tmc_ana_pair_correl_file_name, &
                    ana_env%temperature)
      CALL open_file(file_name=expand_file_name_char(&
                                expand_file_name_char(file_name,&
                                 ana_env%pair_correl%pairs(pair)%f_n),&
                                 ana_env%pair_correl%pairs(pair)%s_n),&
                     file_status="REPLACE", &
                     file_action="WRITE", file_position="APPEND", &
                     unit_number=file_ptr)
      WRITE(file_ptr,*)"# radial distribution function of "//&
                 TRIM(ana_env%pair_correl%pairs(pair)%f_n)//" and "//&
                 TRIM(ana_env%pair_correl%pairs(pair)%s_n)//" of ",&
                 ana_env%pair_correl%conf_counter," configurations"
      WRITE(file_ptr,*)"# using a bin size of ",&
        ana_env%pair_correl%step_lenght*au2a,&
        "[A] (for Vol changes: refering to the reference cell)"
      DO bin=1, ana_env%pair_correl%nr_bins
        voldr = 4.0/3.0* PI* ana_env%pair_correl%step_lenght**3 * &
                  (REAL(bin,KIND=dp)**3-REAL(bin-1,KIND=dp)**3)
      WRITE(file_ptr,*) (bin-0.5)*ana_env%pair_correl%step_lenght*au2a , &
                    (ana_env%pair_correl%g_r(pair,bin)/ana_env%pair_correl%conf_counter)/&
                    (voldr*ana_env%pair_correl%pairs(pair)%pair_count/vol)
      END DO
      CALL close_file(unit_number=file_ptr)

      IF(ana_env%print_test_output) &
        WRITE(*,*)"TMC|ANALYSIS_G_R_"//&
          TRIM(ana_env%pair_correl%pairs(pair)%f_n)//"_"//&
          TRIM(ana_env%pair_correl%pairs(pair)%s_n)//"_X= ",&
          SUM(ana_env%pair_correl%g_r(pair,:)/ana_env%pair_correl%conf_counter/&
             voldr*ana_env%pair_correl%pairs(pair)%pair_count/vol)
    END DO

    ! end the timing
    CALL timestop(handle)
  END SUBROUTINE print_paircorrelation

  !============================================================================
  ! classical cell dipole moment
  !============================================================================

! *****************************************************************************
!> \brief init radial distribution function structures
!> \param ana_dip_mom ...
!> \param atoms ...
!> \param
!> \author Mandes 02.2013
! *****************************************************************************
  SUBROUTINE ana_dipole_moment_init(ana_dip_mom, atoms)
    TYPE(dipole_moment_type), POINTER        :: ana_dip_mom
    TYPE(tmc_atom_type), DIMENSION(:), &
      POINTER                                :: atoms

    CHARACTER(LEN=*), PARAMETER :: routineN = 'ana_dipole_moment_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: atom, charge, handle
    LOGICAL                                  :: failure

    failure = .FALSE.

    CPASSERT(ASSOCIATED(ana_dip_mom))
    CPASSERT(ASSOCIATED(ana_dip_mom%charges_inp))
    CPASSERT(ASSOCIATED(atoms))

    ! start the timing
    CALL timeset(routineN,handle)

    ALLOCATE(ana_dip_mom%charges(SIZE(atoms)))
    ana_dip_mom%charges = 0.0_dp
    ! for every atom searcht the correct charge
    DO atom=1, SIZE(atoms)
      charge_loop:DO charge=1, SIZE(ana_dip_mom%charges_inp)
        IF(atoms(atom)%name .EQ. ana_dip_mom%charges_inp(charge)%name) THEN
          ana_dip_mom%charges(atom) = ana_dip_mom%charges_inp(charge)%mass
          EXIT charge_loop
        END IF
      END DO charge_loop
    END DO

    DEALLOCATE(ana_dip_mom%charges_inp)
    ! end the timing
    CALL timestop(handle)
  END SUBROUTINE ana_dipole_moment_init

! *****************************************************************************
!> \brief calculates the classical cell dipole moment
!> \param elem ...
!> \param weight ...
!> \param ana_env ...
!> \param 
!> \author Mandes 02.2013
! *****************************************************************************
  SUBROUTINE calc_dipole_moment(elem, weight, ana_env)
    TYPE(tree_type), POINTER                 :: elem
    INTEGER                                  :: weight
    TYPE(tmc_analysis_env), POINTER          :: ana_env

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_dipole_moment', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: file_name
    INTEGER                                  :: handle, i
    REAL(KIND=dp), DIMENSION(:), POINTER     :: dip_cl

    CPASSERT(ASSOCIATED(elem))
    CPASSERT(ASSOCIATED(elem%pos))
    CPASSERT(ASSOCIATED(ana_env))
    CPASSERT(ASSOCIATED(ana_env%dip_mom))
    CPASSERT(ASSOCIATED(ana_env%dip_mom%charges))

    ! start the timing
    CALL timeset(routineN,handle)

    ALLOCATE(dip_cl(ana_env%dim_per_elem))
    dip_cl(:) = 0.0_dp

    DO i=1, SIZE(elem%pos,1), ana_env%dim_per_elem
      dip_cl(:) = dip_cl(:)+ elem%pos(i:i+ana_env%dim_per_elem-1)*&
        ana_env%dip_mom%charges(INT(i/REAL(ana_env%dim_per_elem,KIND=dp))+1)
    END DO
   
    ! if there are no exact dipoles save these ones in element structure
    IF(.NOT.ASSOCIATED(elem%dipole))THEN
      ALLOCATE(elem%dipole(ana_env%dim_per_elem))
      elem%dipole(:) = dip_cl(:)
    END IF

    IF(ana_env%dip_mom%print_cl_dip) THEN
      file_name = expand_file_name_temp(tmc_default_trajectory_file_name, &
                    ana_env%temperature)
      CALL write_dipoles_in_file(file_name=file_name, &
             conf_nr=ana_env%dip_mom%conf_counter+1, dip=dip_cl, &
             file_ext="dip_cl")
    END IF
    ana_env%dip_mom%conf_counter = ana_env%dip_mom%conf_counter + weight
    ana_env%dip_mom%last_dip_cl(:)  = dip_cl

    DEALLOCATE(dip_cl)

    ! end the timing
    CALL timestop(handle)
  END SUBROUTINE calc_dipole_moment

! *****************************************************************************
!> \brief prints final values for classical cell dipole moment calculation
!> \param ana_env ...
!> \param 
!> \author Mandes 02.2013
! *****************************************************************************
  SUBROUTINE print_dipole_moment(ana_env)
    TYPE(tmc_analysis_env), POINTER          :: ana_env

    CHARACTER(LEN=*), PARAMETER :: routineN = 'print_dipole_moment', &
      routineP = moduleN//':'//routineN

    IF(ana_env%print_test_output) &
      WRITE(*,*)"TMC|ANALYSIS_FINAL_CLASS_CELL_DIPOLE_MOMENT_X= ", &
        ana_env%dip_mom%last_dip_cl(:)
  END SUBROUTINE print_dipole_moment

! *****************************************************************************
!> \brief calculates the dipole moment analysis
!> \param elem ...
!> \param weight ...
!> \param ana_env ...
!> \param 
!> \author Mandes 03.2013
! *****************************************************************************
  SUBROUTINE calc_dipole_analysis(elem, weight, ana_env)
    TYPE(tree_type), POINTER                 :: elem
    INTEGER                                  :: weight
    TYPE(tmc_analysis_env), POINTER          :: ana_env

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_dipole_analysis', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp)                            :: vol, weight_act
    REAL(KIND=dp), DIMENSION(3, 3)           :: tmp_dip
    TYPE(cell_type), POINTER                 :: scaled_cell

    NULLIFY(scaled_cell)

    CPASSERT(ASSOCIATED(elem))
    CPASSERT(ASSOCIATED(elem%dipole))
    CPASSERT(ASSOCIATED(ana_env))
    CPASSERT(ASSOCIATED(ana_env%dip_ana))

    weight_act = weight
    IF(ana_env%dip_ana%ana_type .EQ. ana_type_sym_xyz) &
      weight_act = weight_act/REAL(8.0,KIND=dp)

    ! get the volume
    ALLOCATE(scaled_cell)
    CALL get_scaled_cell(cell=ana_env%cell, box_scale=elem%box_scale, vol=vol,&
                         scaled_cell=scaled_cell)

    ! fold exact dipole moments using the classical ones
    IF(ASSOCIATED(ana_env%dip_mom)) THEN
      IF(ALL(ana_env%dip_mom%last_dip_cl.NE.elem%dipole)) THEN
        elem%dipole = pbc(r=elem%dipole(:) - ana_env%dip_mom%last_dip_cl, &
                          cell=scaled_cell) + ana_env%dip_mom%last_dip_cl
      END IF
    END IF

    ana_env%dip_ana%conf_counter = ana_env%dip_ana%conf_counter + weight_act

    ! dipole sqared absolut value summed and weight_acted with volume and conf weight_act
    ana_env%dip_ana%mu2_pv_s    = ana_env%dip_ana%mu2_pv_s + &
      DOT_PRODUCT(elem%dipole(:),elem%dipole(:))/vol * weight_act
    
    tmp_dip(:,:) = 0.0_dp
    tmp_dip(:,1) = elem%dipole(:)

    ! dipole sum, weight_acted with volume and conf weight_act
    ana_env%dip_ana%mu_pv(:) = ana_env%dip_ana%mu_pv(:) + &
      tmp_dip(:,1)/vol *weight_act

    ! dipole sum, weight_acted with square root of volume and conf weight_act
    ana_env%dip_ana%mu_psv(:) = ana_env%dip_ana%mu_psv(:) + &
      tmp_dip(:,1)/SQRT(vol) *weight_act

    ! dipole squared sum, weight_acted with volume and conf weight_act
    ana_env%dip_ana%mu2_pv(:) = ana_env%dip_ana%mu2_pv(:) + &
      tmp_dip(:,1)**2 /vol *weight_act

    ! calculate the directional average with componentwise correlation per volume
    tmp_dip(:,:) = MATMUL(tmp_dip(:,:), TRANSPOSE(tmp_dip(:,:)))
    ana_env%dip_ana%mu2_pv_mat(:,:) = ana_env%dip_ana%mu2_pv_mat(:,:) + &
      tmp_dip(:,:)/vol * weight_act

  END SUBROUTINE calc_dipole_analysis

! *****************************************************************************
!> \brief prints the actual dipole moment analysis (trajectories)
!> \param elem ...
!> \param ana_env ...
!> \param 
!> \author Mandes 03.2013
! *****************************************************************************
  SUBROUTINE print_act_dipole_analysis(elem, ana_env)
    TYPE(tree_type), POINTER                 :: elem
    TYPE(tmc_analysis_env), POINTER          :: ana_env

    CHARACTER(LEN=*), PARAMETER :: routineN = 'print_act_dipole_analysis', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: file_name, file_name_tmp
    INTEGER                                  :: counter_tmp, file_ptr
    LOGICAL                                  :: flag
    REAL(KIND=dp)                            :: diel_const, diel_const_norm, &
                                                diel_const_sym, e0, kB
    REAL(KIND=dp), DIMENSION(3, 3)           :: tmp_dip

    kB = boltzmann/joule
    counter_tmp = INT(ana_env%dip_ana%conf_counter)

    ! TODO get correct constant using physcon
    e0 = 0.07957747154594767_dp !e^2*a0*me*hbar^-2
    diel_const_norm = 1/(3.0_dp*e0*kB*ana_env%temperature)

    file_name = expand_file_name_temp(TRIM(ana_env%out_file_prefix)//&
                    tmc_default_trajectory_file_name, &
                    ana_env%temperature)
    CALL write_dipoles_in_file(file_name=file_name, &
           conf_nr=INT(ana_env%dip_ana%conf_counter)+1, dip=elem%dipole, &
           file_ext="dip_folded")

    ! set ouput file name
    file_name_tmp = expand_file_name_temp(TRIM(ana_env%out_file_prefix)//&
                        tmc_default_trajectory_file_name, &
                        ana_env%temperature)

    SELECT CASE (ana_env%dip_ana%ana_type)
    CASE (ana_type_default)
      file_name = TRIM(expand_file_name_char(file_name_tmp,&
                         "diel_const"))
      file_name_tmp = TRIM(expand_file_name_char(file_name_tmp,&
                         "diel_const_tensor"))
    CASE (ana_type_sym_xyz)
      file_name = TRIM(expand_file_name_char(file_name_tmp,&
                         "diel_const_sym"))
      file_name_tmp = TRIM(expand_file_name_char(file_name_tmp,&
                         "diel_const_tensor_sym"))
    CASE DEFAULT
      CALL cp_assert(.FALSE.,&
             cp_warning_level,cp_assertion_failed,routineP,&
             'unknown analysis type "'//cp_to_string(ana_env%dip_ana%ana_type)//&
             '" used. ',&
             failure=flag, only_ionode=.TRUE.)
    END SELECT

    ! calc the dielectric constant
    ! 1+( <M^2> - <M>^2 ) / (3*e_0*V*k*T)
    diel_const = 1.0_dp+(ana_env%dip_ana%mu2_pv_s/(ana_env%dip_ana%conf_counter)-&
      DOT_PRODUCT(ana_env%dip_ana%mu_psv(:)/(ana_env%dip_ana%conf_counter),&
                  ana_env%dip_ana%mu_psv(:)/(ana_env%dip_ana%conf_counter)))*&
      diel_const_norm
    ! symmetrized dielctric constant
    ! 1+( <M^2> ) / (3*e_0*V*k*T)
    diel_const_sym = 1.0_dp+ana_env%dip_ana%mu2_pv_s/(ana_env%dip_ana%conf_counter)*&
      diel_const_norm
    ! print dielectric constant trajectory
    !  if szmetry used print only every 8th configuration, hence every different (not mirrowed)
    INQUIRE(FILE=file_name, EXIST=flag)
    CALL open_file(file_name=file_name, file_status="UNKNOWN", &
                     file_action="WRITE", file_position="APPEND", &
                     unit_number=file_ptr)
    IF(.NOT.flag) THEN
      WRITE(file_ptr,FMT='(A8,5A20)')"# conf", "diel_const", &
                               "diel_const_sym", "diel_const_sym_x", &
                               "diel_const_sym_y", "diel_const_sym_z"
    END IF
    WRITE(file_ptr,FMT="(I8,10F20.10)") counter_tmp, diel_const,&
      diel_const_sym, &
      4.0_dp*PI/(kB*ana_env%temperature) * &
      ana_env%dip_ana%mu2_pv(:)/REAL(ana_env%dip_ana%conf_counter, KIND=dp)
    CALL close_file(unit_number=file_ptr)

    ! print dielectric constant tensor trajectory
    INQUIRE(FILE=file_name_tmp, EXIST=flag)
    CALL open_file(file_name=file_name_tmp, file_status="UNKNOWN", &
                     file_action="WRITE", file_position="APPEND", &
                     unit_number=file_ptr)
    IF(.NOT.flag) THEN
      WRITE(file_ptr,FMT='(A8,9A20)')"# conf", "xx", "xy", "xz", &
                               "yx", "yy", "yz", &
                               "zx", "zy", "zz"
    END IF
    tmp_dip(:,:) = 0.0_dp
    tmp_dip(:,1) = ana_env%dip_ana%mu_psv(:)/REAL(ana_env%dip_ana%conf_counter, KIND=dp)

    WRITE(file_ptr,FMT="(I8,10F20.10)") counter_tmp, &
      4.0_dp*PI/(kB*ana_env%temperature)*&
      (ana_env%dip_ana%mu2_pv_mat(:,:)/REAL(ana_env%dip_ana%conf_counter, KIND=dp) - &
      MATMUL(tmp_dip(:,:), TRANSPOSE(tmp_dip(:,:))))
    CALL close_file(unit_number=file_ptr)
  END SUBROUTINE print_act_dipole_analysis

! *****************************************************************************
!> \brief prints the dipole moment analysis
!> \param ana_env ...
!> \param 
!> \author Mandes 03.2013
! *****************************************************************************
  SUBROUTINE print_dipole_analysis(ana_env)
    TYPE(tmc_analysis_env), POINTER          :: ana_env

    CHARACTER(LEN=*), PARAMETER :: fmt_my = '(T2,A,"| ",A,T41,A40)', &
      plabel = "TMC_ANA", routineN = 'print_dipole_analysis', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    REAL(KIND=dp)                            :: diel_const_scalar, kB
    REAL(KIND=dp), DIMENSION(3)              :: diel_const_sym, dielec_ev
    REAL(KIND=dp), DIMENSION(3, 3)           :: diel_const, tmp_dip, tmp_ev

    kB = boltzmann/joule

    CPASSERT(ASSOCIATED(ana_env))
    CPASSERT(ASSOCIATED(ana_env%dip_ana))

    tmp_dip(:,:)      = 0.0_dp
    diel_const(:,:)   = 0.0_dp
    diel_const_scalar = 0.0_dp
    diel_const_sym    = 0.0_dp

    !dielectric constant 
    tmp_dip(:,1) = ana_env%dip_ana%mu_psv(:)/REAL(ana_env%dip_ana%conf_counter, KIND=dp)
    diel_const(:,:) = 4.0_dp*PI/(kB*ana_env%temperature)*&
      (ana_env%dip_ana%mu2_pv_mat(:,:)/REAL(ana_env%dip_ana%conf_counter, KIND=dp) - &
                        MATMUL(tmp_dip(:,:), TRANSPOSE(tmp_dip(:,:))))

    !dielectric constant for symmetric case
    diel_const_sym(:) = 4.0_dp*PI/(kB*ana_env%temperature) * &
      ana_env%dip_ana%mu2_pv(:)/REAL(ana_env%dip_ana%conf_counter, KIND=dp)

    DO i=1, 3
      diel_const(i,i)   = diel_const(i,i)+1.0_dp ! +1 for unpolarizable models, 1.592 for polarizable
      diel_const_scalar = diel_const_scalar + diel_const(i,i)
    END DO
    diel_const_scalar = diel_const_scalar/REAL(3, KIND=dp)

    tmp_dip(:,:) = diel_const
    CALL diag(3, tmp_dip, dielec_ev, tmp_ev)

    ! print out results
    WRITE(ana_env%io_unit,FMT="(/,T2,A)") REPEAT("-",79)
    WRITE(ana_env%io_unit,FMT="(T2,A,T35,A,T80,A)")"-","average dipoles","-"
    WRITE(ana_env%io_unit,FMT=fmt_my)plabel,"temperature ", cp_to_string(ana_env%temperature)
    WRITE(ana_env%io_unit,FMT=fmt_my)plabel,"used configurations ", &
      cp_to_string(REAL(ana_env%dip_ana%conf_counter, KIND=dp))
    IF(ana_env%dip_ana%ana_type.EQ.ana_type_ice) &
      WRITE(ana_env%io_unit,FMT='(T2,A,"| ",A)')plabel,&
        "ice analysis with directions of hexagonal structure"
    IF(ana_env%dip_ana%ana_type.EQ.ana_type_sym_xyz) &
      WRITE(ana_env%io_unit,FMT='(T2,A,"| ",A)')plabel,&
        "ice analysis with symmetrized dipoles in each direction."

    WRITE(ana_env%io_unit,FMT=fmt_my)plabel,"for product of 2 directions(per vol):"
    DO i=1, 3
      WRITE(ana_env%io_unit,'(A,3F16.8,A)') " |",ana_env%dip_ana%mu2_pv_mat(i, :)/&
        REAL(ana_env%dip_ana%conf_counter, KIND=dp) , " |"
    END DO

    WRITE(ana_env%io_unit,FMT=fmt_my)plabel,"dielectric constant tensor:"
    DO i=1, 3
      WRITE(ana_env%io_unit,'(A,3F16.8,A)') " |", diel_const(i, :), " |"
    END DO

    WRITE(ana_env%io_unit,FMT=fmt_my)plabel,"dielectric tensor eigenvalues", &
      cp_to_string(dielec_ev(1))//" "//&
      cp_to_string(dielec_ev(2))//" "//&
      cp_to_string(dielec_ev(3))
    WRITE(ana_env%io_unit,FMT=fmt_my)plabel,"dielectric constant symm ",&
      cp_to_string(diel_const_sym(1))//" | "//&
      cp_to_string(diel_const_sym(2))//" | "//&
      cp_to_string(diel_const_sym(3))
    WRITE(ana_env%io_unit,FMT=fmt_my)plabel,"dielectric constant ",&
      cp_to_string(diel_const_scalar)
    WRITE(ana_env%io_unit,FMT="(/,T2,A)") REPEAT("-",79)

  END SUBROUTINE print_dipole_analysis

  !============================================================================
  ! particle displacement in cell (from one configuration to the next)
  !============================================================================

! *****************************************************************************
!> \brief calculates the mean square displacement
!> \param elem ...
!> \param ana_env ...
!> \param 
!> \author Mandes 02.2013
! *****************************************************************************
  SUBROUTINE calc_displacement(elem, ana_env)
    TYPE(tree_type), POINTER                 :: elem
    TYPE(tmc_analysis_env), POINTER          :: ana_env

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_displacement', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: file_name, file_name_tmp
    INTEGER                                  :: file_ptr, handle, ind
    LOGICAL                                  :: failure, flag
    REAL(KIND=dp)                            :: disp
    REAL(KIND=dp), DIMENSION(3)              :: atom_disp

    failure = .FALSE.
    disp = 0.0_dp

    CPASSERT(ASSOCIATED(elem))
    CPASSERT(ASSOCIATED(elem%pos))
    CPASSERT(ASSOCIATED(ana_env))
    CPASSERT(ASSOCIATED(ana_env%displace))
    CPASSERT(ASSOCIATED(ana_env%last_elem))

    ! start the timing
    CALL timeset(routineN,handle)

    DO ind=1, SIZE(elem%pos), ana_env%dim_per_elem
      ! fold into box
      atom_disp(:) = elem%pos(ind:ind+2)-ana_env%last_elem%pos(ind:ind+2)
      CALL get_scaled_cell(cell=ana_env%cell, box_scale=elem%box_scale, &
                           vec=atom_disp)
      disp = disp + SUM((atom_disp(:)*au2a)**2)
    END DO
    ana_env%displace%disp = ana_env%displace%disp + disp
    ana_env%displace%conf_counter = ana_env%displace%conf_counter +1

   IF(ana_env%displace%print_disp) THEN
      file_name_tmp = expand_file_name_temp(TRIM(ana_env%out_file_prefix)//&
                        tmc_default_trajectory_file_name, &
                        ana_env%temperature)
      file_name = TRIM(expand_file_name_char(file_name_tmp,&
                          "devi"))
      INQUIRE(FILE=file_name, EXIST=flag)
      CALL open_file(file_name=file_name, file_status="UNKNOWN", &
                     file_action="WRITE", file_position="APPEND", &
                     unit_number=file_ptr)
      IF(.NOT.flag) &
        WRITE(file_ptr,*) "# conf     squared deviation of the cell"
      WRITE(file_ptr,*) elem%nr, disp
      CALL close_file(unit_number=file_ptr)
    END IF

    ! end the timing
    CALL timestop(handle)

  END SUBROUTINE calc_displacement

! *****************************************************************************
!> \brief prints final values for the displacement calculations
!> \param ana_env ...
!> \param 
!> \author Mandes 02.2013
! *****************************************************************************
  SUBROUTINE print_average_displacement(ana_env)
    TYPE(tmc_analysis_env), POINTER          :: ana_env

    CHARACTER(LEN=*), PARAMETER :: fmt_my = '(T2,A,"| ",A,T41,A40)', &
      plabel = "TMC_ANA", routineN = 'print_average_displacement', &
      routineP = moduleN//':'//routineN

    WRITE(ana_env%io_unit,FMT="(/,T2,A)") REPEAT("-",79)
    WRITE(ana_env%io_unit,FMT="(T2,A,T35,A,T80,A)")"-","average displacement","-"
    WRITE(ana_env%io_unit,FMT=fmt_my)plabel,"temperature ", &
      cp_to_string(ana_env%temperature)
    WRITE(ana_env%io_unit,FMT=fmt_my)plabel,"used configurations ", &
      cp_to_string(REAL(ana_env%displace%conf_counter, KIND=dp))
    WRITE(ana_env%io_unit,FMT=fmt_my)plabel,"cell root mean square deviation: ", &
      cp_to_string(SQRT(ana_env%displace%disp/&
      REAL(ana_env%displace%conf_counter, KIND=dp)))
    IF(ana_env%print_test_output) &
      WRITE(*,*)"TMC|ANALYSIS_AVERAGE_CELL_DISPLACEMENT_X= ", &
        SQRT(ana_env%displace%disp/&
        REAL(ana_env%displace%conf_counter, KIND=dp))
  END SUBROUTINE print_average_displacement
END MODULE tmc_analysis
