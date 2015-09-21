!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Chemical shift calculation by dfpt
!>      Initialization of the nmr_env, creation of the special neighbor lists
!>      Perturbation Hamiltonians by application of the p and rxp oprtators to  psi0
!>      Write output
!>      Deallocate everything
!> \note
!>      The psi0 should be localized
!>      the Sebastiani method works within the assumption that the orbitals are
!>      completely contained in the simulation box
!> \par History
!>       created 07-2005 [MI]
!> \author MI
! *****************************************************************************
MODULE qs_linres_nmr_utils

  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE cell_types,                      ONLY: cell_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_parser_methods,               ONLY: parser_get_next_line,&
                                             parser_get_object
  USE cp_parser_types,                 ONLY: cp_parser_type,&
                                             parser_create,&
                                             parser_release
  USE cp_units,                        ONLY: cp_unit_from_cp2k,&
                                             cp_unit_to_cp2k
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_path_length,&
                                             dp
  USE memory_utilities,                ONLY: reallocate
  USE particle_types,                  ONLY: particle_type
  USE pw_env_types,                    ONLY: pw_env_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_linres_types,                 ONLY: linres_control_type,&
                                             nmr_env_create,&
                                             nmr_env_type
  USE qs_matrix_pools,                 ONLY: qs_matrix_pools_type
  USE scf_control_types,               ONLY: scf_control_type
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: nmr_env_cleanup, nmr_env_init

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_linres_nmr_utils'

CONTAINS

! *****************************************************************************
!> \brief Initialize the nmr environment
!> \param nmr_env ...
!> \param qs_env ...
!> \par History
!>      07.2006 created [MI]
!> \author MI
! *****************************************************************************
  SUBROUTINE nmr_env_init(nmr_env,qs_env)
    !
    TYPE(nmr_env_type)                       :: nmr_env
    TYPE(qs_environment_type), POINTER       :: qs_env

    CHARACTER(LEN=*), PARAMETER :: routineN = 'nmr_env_init', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=2)                         :: label
    CHARACTER(LEN=default_path_length)       :: nics_file_name
    INTEGER                                  :: handle, ini, ir, j, n_mo(2), &
                                                n_rep, nao, nat_print, natom, &
                                                nmoloc, nspins, output_unit
    INTEGER, DIMENSION(:), POINTER           :: bounds, list
    LOGICAL                                  :: failure, gapw, my_end
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_parser_type), POINTER            :: parser
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(linres_control_type), POINTER       :: linres_control
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(qs_matrix_pools_type), POINTER      :: mpools
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(section_vals_type), POINTER         :: lr_section, nmr_section

!

    CALL timeset(routineN,handle)

    failure = .FALSE.

    NULLIFY(atomic_kind_set, cell, dft_control, linres_control, scf_control)
    NULLIFY(logger, mpools, nmr_section, particle_set)
    NULLIFY(pw_env, parser)

    n_mo(1:2) = 0
    nao = 0
    nmoloc = 0

    logger => cp_get_default_logger()
    lr_section => section_vals_get_subs_vals(qs_env%input,"PROPERTIES%LINRES")

    output_unit = cp_print_key_unit_nr(logger,lr_section,"PRINT%PROGRAM_RUN_INFO",&
                                       extension=".linresLog")

    IF(nmr_env%ref_count /= 0) THEN
      CALL nmr_env_cleanup(nmr_env)
    END IF

    IF(output_unit>0) THEN
      WRITE(output_unit,"(/,T20,A,/)") "*** Start NMR Chemical Shift Calculation ***"
      WRITE(output_unit,"(T10,A,/)") "Inizialization of the NMR environment"
    ENDIF

    CALL nmr_env_create(nmr_env)
    !
    ! If current_density or full_nmr different allocations are required
    nmr_section => section_vals_get_subs_vals(qs_env%input, &
         &                                    "PROPERTIES%LINRES%NMR")
    CALL section_vals_val_get(nmr_section,"INTERPOLATE_SHIFT",l_val=nmr_env%interpolate_shift)
    CALL section_vals_val_get(nmr_section,"SHIFT_GAPW_RADIUS",r_val=nmr_env%shift_gapw_radius)
    CALL section_vals_val_get(nmr_section,"NICS",l_val=nmr_env%do_nics)
    IF(nmr_env%do_nics) THEN
       CALL section_vals_val_get(nmr_section,"NICS_FILE_NAME",&
                                 c_val=nics_file_name)
       CALL parser_create(parser,nics_file_name)
       CALL parser_get_next_line(parser,1)
       CALL parser_get_object(parser,nmr_env%n_nics)
       ALLOCATE(nmr_env%r_nics(3,nmr_env%n_nics))
       CALL parser_get_next_line(parser,2)
       DO j = 1,nmr_env%n_nics
          CALL parser_get_object(parser,label)
          CALL parser_get_object(parser,nmr_env%r_nics(1,j))
          CALL parser_get_object(parser,nmr_env%r_nics(2,j))
          CALL parser_get_object(parser,nmr_env%r_nics(3,j))
          nmr_env%r_nics(1,j) = cp_unit_to_cp2k(nmr_env%r_nics(1,j),"angstrom")
          nmr_env%r_nics(2,j) = cp_unit_to_cp2k(nmr_env%r_nics(2,j),"angstrom")
          nmr_env%r_nics(3,j) = cp_unit_to_cp2k(nmr_env%r_nics(3,j),"angstrom")
          CALL parser_get_next_line(parser,1,at_end=my_end)
          IF (my_end) EXIT
       ENDDO
       CALL parser_release(parser)
    ENDIF

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    cell=cell,&
                    dft_control=dft_control,&
                    linres_control=linres_control,&
                    mpools=mpools,&
                    particle_set=particle_set,&
                    pw_env=pw_env,&
                    scf_control=scf_control)
    !
    ! Check if restat also psi0 should be restarted
    !IF(nmr_env%restart_nmr .AND. scf_control%density_guess/=restart_guess)THEN
    !   CALL stop_program(routineN,moduleN,__LINE__,"restart_nmr requires density_guess=restart")
    !ENDIF
    !
    ! check that the psi0 are localized and you have all the centers
    CPASSERT(linres_control%localized_psi0)
    IF(failure .AND. (output_unit>0)) THEN
       WRITE(output_unit,'(A)') &
            ' To get NMR parameters within PBC you need localized zero order orbitals '
    ENDIF
    gapw = dft_control%qs_control%gapw
    nspins = dft_control%nspins
    natom = SIZE(particle_set,1)

    !
    ! Conversion factors
    ! factor for the CHEMICAL SHIFTS: alpha^2 *  ppm.
    nmr_env%shift_factor = ( 1.0_dp / 137.03602_dp )**2   * 1.0E+6_dp / cell%deth
    ! factor for the CHEMICAL SHIFTS: alpha^2 *  ppm.
    nmr_env%shift_factor_gapw = ( 1.0_dp / 137.03602_dp )**2   * 1.0E+6_dp
    ! chi_factor =  1/4 * e^2/m * a_0 ^2
    nmr_env%chi_factor = 1.9727566E-29_dp / 1.0E-30_dp ! -> displayed in 10^-30 J/T^2
    ! Factor to convert 10^-30 J/T^2 into ppm cgs = ppm cm^3/mol
    ! = 10^-30 * mu_0/4pi * N_A * 10^6 * 10^6  [one 10^6 for ppm, one for m^3 -> cm^3]
    nmr_env%chi_SI2ppmcgs = 6.022045_dp/1.0E+2_dp
    ! Chi to Shift: 10^-30  *  2/3  mu_0 / Omega  * 1/ppm
    nmr_env%chi_SI2shiftppm = 1.0E-30_dp * 8.37758041E-7_dp/ &
       (cp_unit_from_cp2k(cell%deth,"angstrom^3")*1.0E-30_dp) * 1.0E+6_dp

    IF(output_unit>0) THEN
       WRITE(output_unit,"(T2,A,T65,ES15.6)") "NMR| Shift gapw radius (a.u.) ",nmr_env%shift_gapw_radius
       IF(nmr_env%do_nics) THEN
          WRITE(output_unit,"(T2,A,T50,I5,A)") "NMR| NICS computed in ",nmr_env%n_nics," additional points"
          WRITE(output_unit,"(T2,A,T60,A)") "NMR| NICS coordinates read on file ",TRIM(nics_file_name)
       ENDIF
       WRITE(output_unit,"(T2,A,T65,ES15.6)") "NMR| Shift factor (ppm)",nmr_env%shift_factor
       IF(gapw) THEN
          WRITE(output_unit,"(T2,A,T65,ES15.6)") "NMR| Shift factor gapw (ppm)",nmr_env%shift_factor_gapw
       ENDIF
       WRITE(output_unit,"(T2,A,T65,ES15.6)") "NMR| Chi factor (SI)",nmr_env%chi_factor
       WRITE(output_unit,"(T2,A,T65,ES15.6)") "NMR| Conversion Chi (ppm/cgs)",nmr_env%chi_SI2ppmcgs
       WRITE(output_unit,"(T2,A,T65,ES15.6)") "NMR| Conversion Chi to Shift",nmr_env%chi_SI2shiftppm
    ENDIF

    ALLOCATE (nmr_env%do_calc_cs_atom(natom))
    nmr_env%do_calc_cs_atom = 0

    IF (BTEST(cp_print_key_should_output(logger%iter_info,nmr_section,&
         &    "PRINT%SHIELDING_TENSOR"),cp_p_file)) THEN

       NULLIFY(bounds,list)
       nat_print = 0
       CALL section_vals_val_get(nmr_section,&
            &                    "PRINT%SHIELDING_TENSOR%ATOMS_LU_BOUNDS",&
                                 i_vals=bounds)
       nat_print = bounds(2) - bounds(1)  + 1
       IF(nat_print > 0) THEN
          ALLOCATE(nmr_env%cs_atom_list(nat_print))
          DO ir = 1,nat_print
             nmr_env%cs_atom_list(ir) = bounds(1) + (ir-1)
             nmr_env%do_calc_cs_atom(bounds(1) + (ir-1)) = 1
          ENDDO
       ENDIF

       IF(.NOT. ASSOCIATED(nmr_env%cs_atom_list)) THEN
          CALL section_vals_val_get(nmr_section,"PRINT%SHIELDING_TENSOR%ATOMS_LIST",&
                                    n_rep_val=n_rep)
          nat_print = 0
          DO ir = 1,n_rep
             NULLIFY(list)
             CALL section_vals_val_get(nmr_section,"PRINT%SHIELDING_TENSOR%ATOMS_LIST",&
                                       i_rep_val=ir,i_vals=list)
             IF(ASSOCIATED(list)) THEN
                CALL reallocate(nmr_env%cs_atom_list,1,nat_print + SIZE(list))
                DO ini = 1, SIZE(list)
                   nmr_env%cs_atom_list(ini+nat_print) = list(ini)
                   nmr_env%do_calc_cs_atom(list(ini)) = 1
                ENDDO
                nat_print = nat_print + SIZE(list)
             ENDIF
          ENDDO  ! ir
       ENDIF

       IF(.NOT. ASSOCIATED(nmr_env%cs_atom_list)) THEN
          ALLOCATE(nmr_env%cs_atom_list(natom))
          DO ir = 1,natom
             nmr_env%cs_atom_list(ir) = ir
          ENDDO
          nmr_env%do_calc_cs_atom = 1
       ENDIF
       !
       ! check the list
       CPASSERT(ASSOCIATED(nmr_env%cs_atom_list))
       DO ir = 1,SIZE(nmr_env%cs_atom_list,1)
          IF(nmr_env%cs_atom_list(ir).LT.1.OR.nmr_env%cs_atom_list(ir).GT.natom) THEN
             CALL stop_program(routineN,moduleN,__LINE__,"Unknown atom(s)")
          ENDIF
          DO j = 1,SIZE(nmr_env%cs_atom_list,1)
             IF(j.EQ.ir) CYCLE
             IF(nmr_env%cs_atom_list(ir).EQ.nmr_env%cs_atom_list(j)) THEN
                CALL stop_program(routineN,moduleN,__LINE__,"Duplicate atoms")
             ENDIF
           ENDDO
        ENDDO
    ELSE
       NULLIFY(nmr_env%cs_atom_list)
    ENDIF

    IF(output_unit>0) THEN
       IF(ASSOCIATED(nmr_env%cs_atom_list)) THEN
          WRITE(output_unit,"(T2,A,T69,I5,A)") "NMR| Shielding tensor computed for ",&
                 SIZE(nmr_env%cs_atom_list,1)," atoms"
       ELSE
          WRITE(output_unit,"(T2,A,T50)")&
               & "NMR| Shielding tensor not computed at the atomic positions"
       ENDIF
    ENDIF
    !
    ! Initialize the chemical shift tensor
    ALLOCATE(nmr_env%chemical_shift(3,3,natom),&
             nmr_env%chemical_shift_loc(3,3,natom))
    nmr_env%chemical_shift = 0.0_dp
    nmr_env%chemical_shift_loc = 0.0_dp
    IF(nmr_env%do_nics) THEN
       ALLOCATE(nmr_env%chemical_shift_nics_loc(3,3,nmr_env%n_nics),&
                nmr_env%chemical_shift_nics(3,3,nmr_env%n_nics))
       nmr_env%chemical_shift_nics_loc = 0.0_dp
       nmr_env%chemical_shift_nics = 0.0_dp
    ENDIF

    CALL cp_print_key_finished_output(output_unit,logger,lr_section,&
         &                            "PRINT%PROGRAM_RUN_INFO")

    CALL timestop(handle)

  END SUBROUTINE nmr_env_init

! *****************************************************************************
!> \brief Deallocate the nmr environment
!> \param nmr_env ...
!> \par History
!>      07.2005 created [MI]
!> \author MI
! *****************************************************************************
  SUBROUTINE nmr_env_cleanup(nmr_env)

    TYPE(nmr_env_type)                       :: nmr_env

    CHARACTER(LEN=*), PARAMETER :: routineN = 'nmr_env_cleanup', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.

    nmr_env%ref_count = nmr_env%ref_count - 1
    IF(nmr_env%ref_count == 0 ) THEN
       IF(ASSOCIATED(nmr_env%cs_atom_list)) THEN
          DEALLOCATE(nmr_env%cs_atom_list)
       END IF
       IF(ASSOCIATED(nmr_env%do_calc_cs_atom)) THEN
          DEALLOCATE(nmr_env%do_calc_cs_atom)
       END IF
       !chemical_shift
       IF(ASSOCIATED(nmr_env%chemical_shift)) THEN
          DEALLOCATE(nmr_env%chemical_shift)
       END IF
       IF(ASSOCIATED(nmr_env%chemical_shift_loc)) THEN
          DEALLOCATE(nmr_env%chemical_shift_loc)
       END IF
       ! nics
       IF(ASSOCIATED(nmr_env%r_nics)) THEN
          DEALLOCATE(nmr_env%r_nics)
       END IF
       IF(ASSOCIATED(nmr_env%chemical_shift_nics)) THEN
          DEALLOCATE(nmr_env%chemical_shift_nics)
       END IF
       IF(ASSOCIATED(nmr_env%chemical_shift_nics_loc)) THEN
          DEALLOCATE(nmr_env%chemical_shift_nics_loc)
       END IF
    END IF  ! ref count

  END SUBROUTINE nmr_env_cleanup

END MODULE qs_linres_nmr_utils
