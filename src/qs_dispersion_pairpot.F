!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculation of dispersion using pair potentials
!> \author JGH
! *****************************************************************************
MODULE qs_dispersion_pairpot
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE atprop_types,                    ONLY: atprop_array_init,&
                                             atprop_type
  USE bibliography,                    ONLY: cite_reference,&
                                             grimme2006,&
                                             grimme2010,&
                                             grimme2011
  USE cell_types,                      ONLY: cell_type,&
                                             get_cell,&
                                             pbc,&
                                             plane_distance
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_parser_methods,               ONLY: parser_get_next_line,&
                                             parser_get_object
  USE cp_parser_types,                 ONLY: cp_parser_type,&
                                             parser_create,&
                                             parser_release
  USE input_constants,                 ONLY: vdw_pairpot_dftd2,&
                                             vdw_pairpot_dftd3,&
                                             vdw_pairpot_dftd3bj,&
                                             xc_vdw_fun_none,&
                                             xc_vdw_fun_pairpot
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE mathconstants,                   ONLY: twopi
  USE memory_utilities,                ONLY: reallocate
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_max,&
                                             mp_shift,&
                                             mp_sum
  USE particle_types,                  ONLY: particle_type
  USE physcon,                         ONLY: bohr,&
                                             kcalmol,&
                                             kjmol
  USE qs_dispersion_types,             ONLY: dftd2_pp,&
                                             dftd3_pp,&
                                             qs_atom_dispersion_type,&
                                             qs_dispersion_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_kind_types,                   ONLY: get_qs_kind,&
                                             qs_kind_type,&
                                             set_qs_kind
  USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                             neighbor_list_iterate,&
                                             neighbor_list_iterator_create,&
                                             neighbor_list_iterator_p_type,&
                                             neighbor_list_iterator_release,&
                                             neighbor_list_set_p_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE virial_methods,                  ONLY: virial_pair_force
  USE virial_types,                    ONLY: virial_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_dispersion_pairpot'

  PUBLIC :: qs_dispersion_pairpot_init, calculate_dispersion_pairpot, &
            qs_scaling_init, qs_scaling_dftd3, qs_scaling_dftd3bj

  TYPE dcnum_type
     INTEGER                                :: neighbors
     INTEGER, DIMENSION(:), POINTER         :: nlist
     REAL(KIND=dp), DIMENSION(:), POINTER   :: dvals
     REAL(KIND=dp), DIMENSION(:,:), POINTER :: rik
  END TYPE dcnum_type

! *****************************************************************************

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param atomic_kind_set ...
!> \param qs_kind_set ...
!> \param dispersion_env ...
!> \param pp_section ...
!> \param para_env ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE qs_dispersion_pairpot_init(atomic_kind_set,qs_kind_set,dispersion_env,pp_section,para_env,error)
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(qs_dispersion_type), POINTER        :: dispersion_env
    TYPE(section_vals_type), OPTIONAL, &
      POINTER                                :: pp_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_dispersion_pairpot_init', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=2)                         :: symbol
    CHARACTER(LEN=default_string_length)     :: aname, filename
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: tmpstringlist
    INTEGER                                  :: elem, handle, i, ikind, j, &
                                                max_elem, maxc, n_rep, nkind, &
                                                nl, stat, vdw_pp_type, &
                                                vdw_type
    LOGICAL                                  :: at_end, failure, found, &
                                                is_available
    REAL(KIND=dp)                            :: dum
    TYPE(cp_parser_type), POINTER            :: parser
    TYPE(qs_atom_dispersion_type), POINTER   :: disp

    CALL timeset(routineN,handle)

    failure = .FALSE.

    vdw_type = dispersion_env%type
    SELECT CASE (vdw_type)
       CASE DEFAULT
          ! do nothing
       CASE (xc_vdw_fun_pairpot)
          ! setup information on pair potentials
          vdw_pp_type = dispersion_env%type
          SELECT CASE (dispersion_env%pp_type)
             CASE DEFAULT
               ! do nothing
             CASE (vdw_pairpot_dftd3,vdw_pairpot_dftd3bj)
                !DFT-D3 Method initial setup
                CALL cite_reference(Grimme2010)
                CALL cite_reference(Grimme2011)
                max_elem = 94
                maxc     = 5
                dispersion_env%max_elem = max_elem
                dispersion_env%maxc     = maxc
                ALLOCATE(dispersion_env%maxci(max_elem),STAT=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                ALLOCATE(dispersion_env%c6ab(max_elem,max_elem,maxc,maxc,3),STAT=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                ALLOCATE(dispersion_env%r0ab(max_elem,max_elem),STAT=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                ALLOCATE(dispersion_env%rcov(max_elem),STAT=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                ALLOCATE(dispersion_env%r2r4(max_elem),STAT=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                ALLOCATE(dispersion_env%cn(max_elem),STAT=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

                ! get filename of parameter file
                filename = dispersion_env%parameter_file_name
                CALL dftd3_c6_param(dispersion_env%c6ab,dispersion_env%maxci,filename,para_env,error)
                CALL setr0ab(dispersion_env%r0ab,dispersion_env%rcov,dispersion_env%r2r4)
                ! the default coordination numbers
                CALL setcn(dispersion_env%cn)
                ! scale r4/r2 values of the atoms by sqrt(Z)
                ! sqrt is also globally close to optimum
                ! together with the factor 1/2 this yield reasonable
                ! c8 for he, ne and ar. for larger Z, C8 becomes too large
                ! which effectively mimics higher R^n terms neglected due
                ! to stability reasons
                DO i=1,max_elem
                   dum    = 0.5_dp*dispersion_env%r2r4(i)*REAL(i,dp)**0.5_dp
                   ! store it as sqrt because the geom. av. is taken
                   dispersion_env%r2r4(i)=SQRT(dum)
                ENDDO
                ! parameters
                dispersion_env%k1=16.0_dp
                dispersion_env%k2=4._dp/3._dp
                ! reasonable choices are between 3 and 5
                ! this gives smoth curves with maxima around the integer values
                ! k3=3 give for CN=0 a slightly smaller value than computed
                ! for the free atom. This also yields to larger CN for atoms
                ! in larger molecules but with the same chem. environment
                ! which is physically not right
                ! values >5 might lead to bumps in the potential
                dispersion_env%k3=-4._dp
                dispersion_env%rcov=dispersion_env%k2*dispersion_env%rcov*bohr
                ! alpha default parameter
                dispersion_env%alp=14._dp
          END SELECT
    END SELECT

    nkind  = SIZE(atomic_kind_set)
    DO ikind = 1, nkind
       CALL get_atomic_kind(atomic_kind_set(ikind),element_symbol=symbol,z=elem)
       ALLOCATE(disp,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       SELECT CASE (vdw_type)
       CASE DEFAULT
          CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
       CASE (xc_vdw_fun_none)
          disp%defined = .FALSE.
       CASE (xc_vdw_fun_pairpot)
          ! setup information on pair potentials
          SELECT CASE (dispersion_env%pp_type)
          CASE DEFAULT
             CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
          CASE (vdw_pairpot_dftd2)
             CALL cite_reference(Grimme2006)
             disp%type = dftd2_pp
             ! get filename of parameter file
             filename = dispersion_env%parameter_file_name
             ! check for local parameters
             found = .FALSE.
             IF (PRESENT(pp_section)) THEN
                CALL section_vals_val_get(pp_section,"ATOMPARM",n_rep_val=n_rep,error=error)
                DO i=1,n_rep
                   CALL section_vals_val_get(pp_section,"ATOMPARM", i_rep_val=i,&
                        c_vals=tmpstringlist, error=error)
                   IF ( TRIM(tmpstringlist(1)) == TRIM(symbol) ) THEN
                      ! we assume the parameters are in atomic units!
                      READ(tmpstringlist(2),*) disp%c6
                      READ(tmpstringlist(3),*) disp%vdw_radii
                      found = .TRUE.
                      EXIT
                   END IF
                END DO
             END IF
             IF ( .NOT. found ) THEN
                ! check for internal parameters
                CALL dftd2_param(elem,disp%c6,disp%vdw_radii,found,error)
             END IF
             IF ( .NOT. found ) THEN
                ! check on file
                INQUIRE (FILE=filename,EXIST=is_available)
                IF (is_available) THEN
                   NULLIFY(parser)
                   CALL parser_create(parser,filename,para_env=para_env,error=error)
                   DO
                      at_end = .FALSE.
                      CALL parser_get_next_line(parser,1,at_end,error=error)
                      IF ( at_end ) EXIT
                      CALL parser_get_object(parser,aname,error=error)
                      IF ( TRIM(aname) == TRIM(symbol) ) THEN
                         CALL parser_get_object(parser,disp%c6,error=error)
                         ! we have to change the units J*nm^6*mol^-1 -> Hartree*Bohr^6
                         disp%c6 = disp%c6 * 1000._dp*bohr**6/kjmol
                         CALL parser_get_object(parser,disp%vdw_radii,error=error)
                         disp%vdw_radii = disp%vdw_radii*bohr
                         found = .TRUE.
                         EXIT
                      END IF
                   END DO
                   CALL parser_release(parser,error=error)
                END IF
             END IF
             IF ( found ) THEN
                disp%defined = .TRUE.
             ELSE
                disp%defined = .FALSE.
             END IF
          CASE (vdw_pairpot_dftd3,vdw_pairpot_dftd3bj)
             disp%type = dftd3_pp
             IF ( elem <= 94 ) THEN
                disp%defined = .TRUE.
             ELSE
                disp%defined = .FALSE.
             END IF

             ! Check for coordination numbers
             IF (PRESENT(pp_section)) THEN
                CALL section_vals_val_get(pp_section,"KIND_COORDINATION_NUMBERS", n_rep_val=n_rep, error=error)
                IF (n_rep > 0) THEN
                   ALLOCATE(dispersion_env%cnkind(n_rep),STAT=stat)
                   CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                   DO i=1,n_rep
                      CALL section_vals_val_get(pp_section,"KIND_COORDINATION_NUMBERS", i_rep_val=i,&
                        c_vals=tmpstringlist, error=error)
                      READ(tmpstringlist(1),*) dispersion_env%cnkind(i)%cnum
                      READ(tmpstringlist(2),*) dispersion_env%cnkind(i)%kind
                   END DO
                END IF
                CALL section_vals_val_get(pp_section,"ATOM_COORDINATION_NUMBERS", n_rep_val=n_rep, error=error)
                IF (n_rep > 0) THEN
                   ALLOCATE(dispersion_env%cnlist(n_rep),STAT=stat)
                   CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                   DO i=1,n_rep
                      CALL section_vals_val_get(pp_section,"ATOM_COORDINATION_NUMBERS", i_rep_val=i,&
                        c_vals=tmpstringlist, error=error)
                      nl = SIZE(tmpstringlist)
                      ALLOCATE(dispersion_env%cnlist(i)%atom(nl-1),STAT=stat)
                      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                      dispersion_env%cnlist(i)%natom = nl-1
                      READ(tmpstringlist(1),*) dispersion_env%cnlist(i)%cnum
                      DO j=1,nl-1
                         READ(tmpstringlist(j+1),*) dispersion_env%cnlist(i)%atom(j)
                      END DO
                   END DO
                END IF
             END IF

          END SELECT
       END SELECT
       ! Check if the parameter is defined
       CALL cp_assert(disp%defined,cp_failure_level,cp_assertion_failed,routineP,&
            "Dispersion parameters for element ("//TRIM(symbol)//") are not defined! "//&
            "Please provide a valid set of parameters through the input section or "//&
            "through an external file! "//&
 CPSourceFileRef,&
            error=error,failure=failure)
       CALL set_qs_kind(qs_kind_set(ikind),dispersion=disp,error=error)
    END DO

    CALL timestop(handle)

  END SUBROUTINE qs_dispersion_pairpot_init

! *****************************************************************************

! *****************************************************************************
!> \brief ...
!> \param scaling ...
!> \param vdw_section ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE qs_scaling_init(scaling,vdw_section,error)
    REAL(KIND=dp), INTENT(inout)             :: scaling
    TYPE(section_vals_type), POINTER         :: vdw_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_scaling_init', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: functional
    LOGICAL                                  :: failure

    failure = .FALSE.

    CALL section_vals_val_get(vdw_section,"PAIR_POTENTIAL%REFERENCE_FUNCTIONAL",c_val=functional,error=error)

    SELECT CASE (TRIM(functional))
       CASE DEFAULT
          ! unknown functional
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="No DFT-D2 s6 value available for this functional:"//TRIM(functional), &
               error=error, error_level=cp_failure_level)
       CASE ("BLYP")
          scaling = 1.20_dp
       CASE ("B3LYP")
          scaling = 1.05_dp
       CASE("TPSS")
          scaling = 1.00_dp
       CASE("PBE")
          scaling = 0.75_dp
       CASE("PBE0")
          scaling = 0.6_dp
       CASE("B2PLYP")
          scaling = 0.55_dp
       CASE("BP86")
          scaling = 1.05_dp
       CASE("B97")
          scaling = 1.25_dp
    END SELECT

  END SUBROUTINE qs_scaling_init

! *****************************************************************************

! *****************************************************************************
!> \brief ...
!> \param s6 ...
!> \param sr6 ...
!> \param s8 ...
!> \param vdw_section ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE qs_scaling_dftd3(s6,sr6,s8,vdw_section,error)
    REAL(KIND=dp), INTENT(inout)             :: s6, sr6, s8
    TYPE(section_vals_type), POINTER         :: vdw_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_scaling_dftd3', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: functional
    LOGICAL                                  :: failure

    failure = .FALSE.

    CALL section_vals_val_get(vdw_section,"PAIR_POTENTIAL%REFERENCE_FUNCTIONAL",c_val=functional,error=error)

    ! values for different functionals from:
    ! http://www.thch.uni-bonn.de/tc/downloads/DFT-D3/functionals.html
    SELECT CASE (TRIM(functional))
    CASE DEFAULT
       ! unknown functional
       CALL cp_unimplemented_error(fromWhere=routineP, &
            message="No DFT-D3 values available for this functional:"//TRIM(functional), &
            error=error, error_level=cp_failure_level)
    CASE ("B1B95")
       s6  = 1.000_dp
       sr6 = 1.613_dp
       s8  = 1.868_dp
    CASE ("B2GPPLYP")
       ! L. Goerigk and S. Grimme
       ! J. Chem. Theory Comput. 2011, 7, 291–309; doi:10.1021/ct100466k
       s6  = 0.56_dp
       sr6 = 1.586_dp
       s8  = 0.760_dp
    CASE ("B2PLYP")
       ! L. Goerigk and S. Grimme
       ! J. Chem. Theory Comput. 2011, 7, 291–309; doi:10.1021/ct100466k
       s6  = 0.64_dp
       sr6 = 1.427_dp
       s8  = 1.022_dp
    CASE ("DSD-BLYP")
       ! L. Goerigk and S. Grimme
       ! J. Chem. Theory Comput. 2011, 7, 291–309; doi:10.1021/ct100466k
       s6  = 0.50_dp
       sr6 = 1.569_dp
       s8  = 0.705_dp
    CASE ("B3LYP")
       ! S. Grimme, J. Antony, S. Ehrlich, and H. Krieg
       ! J. Chem. Phys. 132,154104 (2010); doi:10.1063/1.3382344
       s6  = 1.000_dp
       sr6 = 1.261_dp
       s8  = 1.703_dp
    CASE ("B97-D")
       ! S. Grimme, J. Antony, S. Ehrlich, and H. Krieg
       ! J. Chem. Phys. 132,154104 (2010); doi:10.1063/1.3382344
       s6  = 1.000_dp
       sr6 = 0.892_dp
       s8  = 0.909_dp
    CASE ("BHLYP")
       s6  = 1.000_dp
       sr6 = 1.370_dp
       s8  = 1.442_dp
    CASE ("BLYP")
       ! S. Grimme, J. Antony, S. Ehrlich, and H. Krieg
       ! J. Chem. Phys. 132,154104 (2010); doi:10.1063/1.3382344
       s6  = 1.000_dp
       sr6 = 1.094_dp
       s8  = 1.682_dp
    CASE ("BP86")
       ! S. Grimme, J. Antony, S. Ehrlich, and H. Krieg
       ! J. Chem. Phys. 132,154104 (2010); doi:10.1063/1.3382344
       s6  = 1.000_dp
       sr6 = 1.139_dp
       s8  = 1.683_dp
    CASE ("BPBE")
       s6  = 1.000_dp
       sr6 = 1.087_dp
       s8  = 2.033_dp
    CASE ("MPWLYP")
       s6  = 1.000_dp
       sr6 = 1.239_dp
       s8  = 1.098_dp
    CASE ("PBE")
       ! S. Grimme, J. Antony, S. Ehrlich, and H. Krieg
       ! J. Chem. Phys. 132,154104 (2010); doi:10.1063/1.3382344
       s6  = 1.000_dp
       sr6 = 1.217_dp
       s8  = 0.722_dp
    CASE ("PBE0")
       ! S. Grimme, J. Antony, S. Ehrlich, and H. Krieg
       ! J. Chem. Phys. 132,154104 (2010); doi:10.1063/1.3382344
       s6  = 1.000_dp
       sr6 = 1.287_dp
       s8  = 0.928_dp
    CASE ("PW6B95")
       ! S. Grimme, J. Antony, S. Ehrlich, and H. Krieg
       ! J. Chem. Phys. 132,154104 (2010); doi:10.1063/1.3382344
       s6  = 1.000_dp
       sr6 = 1.532_dp
       s8  = 0.862_dp
    CASE ("PWB6K")
       s6  = 1.000_dp
       sr6 = 1.660_dp
       s8  = 0.550_dp
    CASE ("revPBE")
       ! S. Grimme, J. Antony, S. Ehrlich, and H. Krieg
       ! J. Chem. Phys. 132,154104 (2010); doi:10.1063/1.3382344
       s6  = 1.000_dp
       sr6 = 0.923_dp
       s8  = 1.010_dp
    CASE ("TPSS")
       ! S. Grimme, J. Antony, S. Ehrlich, and H. Krieg
       ! J. Chem. Phys. 132,154104 (2010); doi:10.1063/1.3382344
       s6  = 1.000_dp
       sr6 = 1.166_dp
       s8  = 1.105_dp
    CASE ("TPSS0")
       ! S. Grimme, J. Antony, S. Ehrlich, and H. Krieg
       ! J. Chem. Phys. 132,154104 (2010); doi:10.1063/1.3382344
       s6  = 1.000_dp
       sr6 = 1.252_dp
       s8  = 1.242_dp
    CASE ("TPSSh")
       s6  = 1.000_dp
       sr6 = 1.223_dp
       s8  = 1.219_dp
    CASE ("B3PW91")
       s6  = 1.000_dp
       sr6 = 1.176_dp
       s8  = 1.775_dp
    CASE ("BMK")
       s6  = 1.000_dp
       sr6 = 1.931_dp
       s8  = 2.168_dp
    CASE ("CAMB3LYP")
       s6  = 1.000_dp
       sr6 = 1.378_dp
       s8  = 1.217_dp
    CASE ("LCWPBE")
       s6  = 1.000_dp
       sr6 = 1.355_dp
       s8  = 1.279_dp
    CASE ("M052X")
       s6  = 1.000_dp
       sr6 = 1.417_dp
       s8  = 0.000_dp
    CASE ("M05")
       s6  = 1.000_dp
       sr6 = 1.373_dp
       s8  = 0.595_dp
    CASE ("M062X")
       s6  = 1.000_dp
       sr6 = 1.619_dp
       s8  = 0.000_dp
    CASE ("M06HF")
       s6  = 1.000_dp
       sr6 = 1.446_dp
       s8  = 0.000_dp
    CASE ("M06L")
       s6  = 1.000_dp
       sr6 = 1.581_dp
       s8  = 0.000_dp
    CASE ("M06N")
       s6  = 1.000_dp
       sr6 = 1.325_dp
       s8  = 0.000_dp
    CASE ("HCTH120")
       s6  = 1.000_dp
       sr6 = 1.221_dp
       s8  = 1.206_dp
    CASE ("PTPSS")
       s6  = 0.750_dp
       sr6 = 1.541_dp
       s8  = 0.879_dp
    CASE ("PWPB95")
       s6  = 0.820_dp
       sr6 = 1.557_dp
       s8  = 0.705_dp
    CASE ("OLYP")
       s6  = 1.000_dp
       sr6 = 0.806_dp
       s8  = 1.764_dp
    CASE ("OPBE")
       s6  = 1.000_dp
       sr6 = 0.837_dp
       s8  = 2.055_dp
    CASE ("OTPSS")
       s6  = 1.000_dp
       sr6 = 1.128_dp
       s8  = 1.494_dp
    CASE ("PBE38")
       s6  = 1.000_dp
       sr6 = 1.333_dp
       s8  = 0.998_dp
    CASE ("PBEsol")
       s6  = 1.000_dp
       sr6 = 1.345_dp
       s8  = 0.612_dp
    CASE ("revSSB")
       s6  = 1.000_dp
       sr6 = 1.221_dp
       s8  = 0.560_dp
    CASE ("SSB")
       s6  = 1.000_dp
       sr6 = 1.215_dp
       s8  = 0.663_dp
    CASE ("BOP")
       s6  = 1.000_dp
       sr6 = 0.929_dp
       s8  = 1.975_dp
    CASE ("MPW1B95")
       s6  = 1.000_dp
       sr6 = 1.605_dp
       s8  = 1.118_dp
    CASE ("MPWB1K")
       s6  = 1.000_dp
       sr6 = 1.671_dp
       s8  = 1.061_dp
    CASE ("revPBE0")
       s6  = 1.000_dp
       sr6 = 0.949_dp
       s8  = 0.792_dp
    CASE ("revPBE38")
       s6  = 1.000_dp
       sr6 = 1.021_dp
       s8  = 0.862_dp
    CASE ("rPW86PBE")
       s6  = 1.000_dp
       sr6 = 1.224_dp
       s8  = 0.901_dp
    END SELECT

  END SUBROUTINE qs_scaling_dftd3

! *****************************************************************************
!> \brief ...
!> \param s6 ...
!> \param a1 ...
!> \param s8 ...
!> \param a2 ...
!> \param vdw_section ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE qs_scaling_dftd3bj(s6,a1,s8,a2,vdw_section,error)
    REAL(KIND=dp), INTENT(inout)             :: s6, a1, s8, a2
    TYPE(section_vals_type), POINTER         :: vdw_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_scaling_dftd3bj', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: functional
    LOGICAL                                  :: failure

    failure = .FALSE.

    CALL section_vals_val_get(vdw_section,"PAIR_POTENTIAL%REFERENCE_FUNCTIONAL",c_val=functional,error=error)

    ! values for different functionals from:
    ! http://www.thch.uni-bonn.de/tc/downloads/DFT-D3/functionalsbj.html
    SELECT CASE (TRIM(functional))
    CASE DEFAULT
       ! unknown functional
       CALL cp_unimplemented_error(fromWhere=routineP, &
            message="No DFT-D3(BJ) values available for this functional:"//TRIM(functional), &
            error=error, error_level=cp_failure_level)
    CASE ("B1B95")
       s6  = 1.0000_dp
       a1  = 0.2092_dp
       s8  = 1.4507_dp
       a2  = 5.5545_dp
    CASE ("B2GPPLYP")
       s6  = 0.5600_dp
       a1  = 0.0000_dp
       s8  = 0.2597_dp
       a2  = 6.3332_dp
    CASE ("B3PW91")
       s6  = 1.0000_dp
       a1  = 0.4312_dp
       s8  = 2.8524_dp
       a2  = 4.4693_dp
    CASE ("BHLYP")
       s6  = 1.0000_dp
       a1  = 0.2793_dp
       s8  = 1.0354_dp
       a2  = 4.9615_dp
    CASE ("BMK")
       s6  = 1.0000_dp
       a1  = 0.1940_dp
       s8  = 2.0860_dp
       a2  = 5.9197_dp
    CASE ("BOP")
       s6  = 1.0000_dp
       a1  = 0.4870_dp
       s8  = 3.2950_dp
       a2  = 3.5043_dp
    CASE ("BPBE")
       s6  = 1.0000_dp
       a1  = 0.4567_dp
       s8  = 4.0728_dp
       a2  = 4.3908_dp
    CASE ("CAMB3LYP")
       s6  = 1.0000_dp
       a1  = 0.3708_dp
       s8  = 2.0674_dp
       a2  = 5.4743_dp
    CASE ("LCWPBE")
       s6  = 1.0000_dp
       a1  = 0.3919_dp
       s8  = 1.8541_dp
       a2  = 5.0897_dp
    CASE ("MPW1B95")
       s6  = 1.0000_dp
       a1  = 0.1955_dp
       s8  = 1.0508_dp
       a2  = 6.4177_dp
    CASE ("MPWB1K")
       s6  = 1.0000_dp
       a1  = 0.1474_dp
       s8  = 0.9499_dp
       a2  = 6.6223_dp
    CASE ("MPWLYP")
       s6  = 1.0000_dp
       a1  = 0.4831_dp
       s8  = 2.0077_dp
       a2  = 4.5323_dp
    CASE ("OLYP")
       s6  = 1.0000_dp
       a1  = 0.5299_dp
       s8  = 2.6205_dp
       a2  = 2.8065_dp
    CASE ("OPBE")
       s6  = 1.0000_dp
       a1  = 0.5512_dp
       s8  = 3.3816_dp
       a2  = 2.9444_dp
    CASE ("OTPSS")
       s6  = 1.0000_dp
       a1  = 0.4634_dp
       s8  = 2.7495_dp
       a2  = 4.3153_dp
    CASE ("PBE38")
       s6  = 1.0000_dp
       a1  = 0.3995_dp
       s8  = 1.4623_dp
       a2  = 5.1405_dp
    CASE ("PBEsol")
       s6  = 1.0000_dp
       a1  = 0.4466_dp
       s8  = 2.9491_dp
       a2  = 6.1742_dp
    CASE ("PTPSS")
       s6  = 0.7500_dp
       a1  = 0.0000_dp
       s8  = 0.2804_dp
       a2  = 6.5745_dp
    CASE ("PWB6K")
       s6  = 1.0000_dp
       a1  = 0.1805_dp
       s8  = 0.9383_dp
       a2  = 7.7627_dp
    CASE ("revSSB")
       s6  = 1.0000_dp
       a1  = 0.4720_dp
       s8  = 0.4389_dp
       a2  = 4.0986_dp
    CASE ("SSB")
       s6  = 1.0000_dp
       a1  =-0.0952_dp
       s8  =-0.1744_dp
       a2  = 5.2170_dp
    CASE ("TPSSh")
       s6  = 1.0000_dp
       a1  = 0.4529_dp
       s8  = 2.2382_dp
       a2  = 4.6550_dp
    CASE ("HCTH120")
       s6  = 1.0000_dp
       a1  = 0.3563_dp
       s8  = 1.0821_dp
       a2  = 4.3359_dp
    CASE ("B2PLYP")
       s6  = 0.6400_dp
       a1  = 0.3065_dp
       s8  = 0.9147_dp
       a2  = 5.0570_dp
    CASE ("B3LYP")
       s6  = 1.0000_dp
       a1  = 0.3981_dp
       s8  = 1.9889_dp
       a2  = 4.4211_dp
    CASE ("B97-D")
       s6  = 1.0000_dp
       a1  = 0.5545_dp
       s8  = 2.2609_dp
       a2  = 3.2297_dp
    CASE ("BLYP")
       s6  = 1.0000_dp
       a1  = 0.4298_dp
       s8  = 2.6996_dp
       a2  = 4.2359_dp
    CASE ("BP86")
       s6  = 1.0000_dp
       a1  = 0.3946_dp
       s8  = 3.2822_dp
       a2  = 4.8516_dp
    CASE ("DSD-BLYP")
       s6  = 0.5000_dp
       a1  = 0.0000_dp
       s8  = 0.2130_dp
       a2  = 6.0519_dp
    CASE ("PBE0")
       s6  = 1.0000_dp
       a1  = 0.4145_dp
       s8  = 1.2177_dp
       a2  = 4.8593_dp
    CASE ("PBE")
       s6  = 1.0000_dp
       a1  = 0.4289_dp
       s8  = 0.7875_dp
       a2  = 4.4407_dp
    CASE ("PW6B95")
       s6  = 1.0000_dp
       a1  = 0.2076_dp
       s8  = 0.7257_dp
       a2  = 6.3750_dp
    CASE ("PWPB95")
       s6  = 0.8200_dp
       a1  = 0.0000_dp
       s8  = 0.2904_dp
       a2  = 7.3141_dp
    CASE ("revPBE0")
       s6  = 1.0000_dp
       a1  = 0.4679_dp
       s8  = 1.7588_dp
       a2  = 3.7619_dp
    CASE ("revPBE38")
       s6  = 1.0000_dp
       a1  = 0.4309_dp
       s8  = 1.4760_dp
       a2  = 3.9446_dp
    CASE ("revPBE")
       s6  = 1.0000_dp
       a1  = 0.5238_dp
       s8  = 2.3550_dp
       a2  = 3.5016_dp
    CASE ("RPW86PBE")
       s6  = 1.0000_dp
       a1  = 0.4613_dp
       s8  = 1.3845_dp
       a2  = 4.5062_dp
    CASE ("TPSS0")
       s6  = 1.0000_dp
       a1  = 0.3768_dp
       s8  = 1.2576_dp
       a2  = 4.5865_dp
    CASE ("TPSS")
       s6  = 1.0000_dp
       a1  = 0.4535_dp
       s8  = 1.9435_dp
       a2  = 4.4752_dp
    END SELECT

  END SUBROUTINE qs_scaling_dftd3bj

! *****************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param dispersion_env ...
!> \param energy ...
!> \param calculate_forces ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE calculate_dispersion_pairpot(qs_env,dispersion_env,energy,calculate_forces,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qs_dispersion_type), POINTER        :: dispersion_env
    REAL(KIND=dp), INTENT(OUT)               :: energy
    LOGICAL, INTENT(IN)                      :: calculate_forces
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_dispersion_pairpot', &
      routineP = moduleN//':'//routineN

    INTEGER :: atom_a, atom_b, atom_c, atom_d, handle, hashb, hashc, i, ia, &
      iatom, icx, icy, icz, idmp, ikind, ilist, istat, jatom, jkind, katom, &
      kkind, kstart, latom, lkind, max_elem, maxc, na, natom, nb, nc, ni, nj, &
      nkind, unit_nr, za, zb, zc
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind, atomnumber, &
                                                kind_of
    INTEGER, DIMENSION(3)                    :: cell_b, cell_c, ncell, &
                                                periodic
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL :: atenergy, atstress, debug, debugall, debugx, failure, ghost_a, &
      ghost_b, ghost_c, is000, use_virial
    LOGICAL, ALLOCATABLE, DIMENSION(:)       :: dodisp, ghost
    REAL(KIND=dp) :: a1, a2, alp6, alp8, ang, c6, c8, c9, cc6ab, cc6bc, &
      cc6ca, cnab, cnum, dc6a, dc6b, dc8a, dc8b, dcc6aba, dcc6abb, dcc6bcb, &
      dcc6bcc, dcc6caa, dcc6cac, dcnab, dd, de6, de8, de91, de921, de922, &
      dea, devdw, dfdab6, dfdab8, dfdabc, dfdmp, dr, drk, e6, e6tot, e8, &
      e8tot, e9, e9tot, elrc6, elrc8, elrc9, eps_cn, er, evdw, f0ab, fac, &
      fac0, fdab6, fdab8, fdabc, fdmp, gnorm, nab, nabc, r0, r0ab, r2ab, &
      r2abc, r2bc, r2ca, r6, r8, rabc, rc2, rcc, rcova, rcovb, rcut, rs6, &
      rs8, s1, s2, s3, s6, s8, s9, xp
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: c6d2, cnkind, cnumbers, &
                                                cnumfix, radd2
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: rcpbc
    REAL(KIND=dp), DIMENSION(3)              :: fdij, fdik, ra, rab, rb, rb0, &
                                                rbc, rc, rc0, rca, rij, rik, &
                                                sab_max
    REAL(KIND=dp), DIMENSION(3, 3)           :: dvirial
    REAL(KIND=dp), DIMENSION(:), POINTER     :: atener
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: atstr
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atprop_type), POINTER               :: atprop
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dcnum_type), ALLOCATABLE, &
      DIMENSION(:)                           :: dcnum
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: nl_iterator
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_cn, sab_vdw
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_atom_dispersion_type), POINTER   :: disp_a, disp_b, disp_c
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(virial_type), POINTER               :: virial

    CALL timeset(routineN,handle)

    failure = .FALSE.

    energy = 0._dp

    IF ( dispersion_env%type == xc_vdw_fun_pairpot ) THEN

       NULLIFY (atomic_kind_set,qs_kind_set,sab_cn,sab_vdw)

       CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set,qs_kind_set=qs_kind_set,&
            cell=cell, virial=virial, para_env=para_env, atprop=atprop, error=error)

       debugx = dispersion_env%verbose
       debugall = debugx
       debug = debugx .AND. para_env%ionode
       nkind = SIZE(atomic_kind_set)

       NULLIFY(logger)
       logger => cp_error_get_logger(error)
       IF(ASSOCIATED(dispersion_env%dftd_section)) THEN
          unit_nr = cp_print_key_unit_nr(logger,dispersion_env%dftd_section,"PRINT_DFTD",&
                                         extension=".dftd",error=error)
       ELSE
          unit_nr = -1
       END IF

       ! atomic energy and stress arrays
       atenergy = atprop%energy
       IF (atenergy) THEN
          NULLIFY (particle_set)
          CALL get_qs_env(qs_env=qs_env,particle_set=particle_set,error=error)
          natom = SIZE(particle_set)
          CALL atprop_array_init(atprop%atevdw,natom,error)
          atener => atprop%atevdw
       END IF
       atstress = atprop%stress
       IF (atstress) THEN
          atstr => atprop%atstress
       END IF

       IF(unit_nr>0) THEN
        WRITE(unit_nr,*)
        WRITE(unit_nr,*) " Pair potential vdW calculation"
        IF ( dispersion_env%pp_type == vdw_pairpot_dftd2 ) THEN
           WRITE(unit_nr,*) " Dispersion potential type: DFT-D2"
        ELSE IF ( dispersion_env%pp_type == vdw_pairpot_dftd3 ) THEN
           WRITE(unit_nr,*) " Dispersion potential type: DFT-D3"
        ELSE IF ( dispersion_env%pp_type == vdw_pairpot_dftd3bj ) THEN
           WRITE(unit_nr,*) " Dispersion potential type: DFT-D3(BJ)"
        END IF
       END IF

       NULLIFY (particle_set)
       CALL get_qs_env(qs_env=qs_env,particle_set=particle_set,error=error)
       natom = SIZE(particle_set)
       IF (calculate_forces .OR. debugall) THEN
          NULLIFY (force)
          CALL get_qs_env(qs_env=qs_env,force=force,error=error)
          ALLOCATE (atom_of_kind(natom),kind_of(natom),STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          CALL get_atomic_kind_set(atomic_kind_set, atom_of_kind=atom_of_kind,kind_of=kind_of)
          use_virial = virial%pv_availability.AND.(.NOT.virial%pv_numer)
          IF (use_virial .AND. debugall) THEN
            dvirial = virial%pv_virial
          END IF
       ELSE IF ( ( dispersion_env%pp_type == vdw_pairpot_dftd3 .OR.&
                   dispersion_env%pp_type == vdw_pairpot_dftd3bj) .AND. dispersion_env%doabc ) THEN
          ALLOCATE (atom_of_kind(natom),kind_of(natom),STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          CALL get_atomic_kind_set(atomic_kind_set, atom_of_kind=atom_of_kind,kind_of=kind_of)
       END IF

       ALLOCATE (dodisp(nkind),ghost(nkind),atomnumber(nkind),c6d2(nkind),radd2(nkind),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       DO ikind=1,nkind
          CALL get_atomic_kind(atomic_kind_set(ikind), z=za)
          CALL get_qs_kind(qs_kind_set(ikind), dispersion=disp_a, ghost=ghost_a)
          dodisp(ikind)=disp_a%defined
          ghost(ikind)=ghost_a
          atomnumber(ikind)=za
          c6d2(ikind)=disp_a%c6
          radd2(ikind)=disp_a%vdw_radii
       END DO
       ALLOCATE (rcpbc(3,natom),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       DO iatom=1,natom
         rcpbc(:,iatom) = pbc(particle_set(iatom)%r(:),cell)
       END DO

       rcut = 2._dp*dispersion_env%rc_disp
       rc2  = rcut*rcut
       IF ( dispersion_env%pp_type == vdw_pairpot_dftd2 ) THEN
          s6 = dispersion_env%scaling
          dd = dispersion_env%exp_pre
          IF(unit_nr>0) THEN
             WRITE(unit_nr,*) " Scaling parameter (s6) ",s6
             WRITE(unit_nr,*) " Exponential prefactor  ",dd
          END IF
       ELSE IF ( dispersion_env%pp_type == vdw_pairpot_dftd3  .OR. &
                 dispersion_env%pp_type == vdw_pairpot_dftd3bj ) THEN
          maxc = SIZE(dispersion_env%c6ab,3)
          max_elem = SIZE(dispersion_env%c6ab,1)
          alp6 = dispersion_env%alp
          alp8 = alp6 + 2._dp
          s6   = dispersion_env%s6
          s8   = dispersion_env%s8
          s9   = dispersion_env%s6
          rs6  = dispersion_env%sr6
          rs8  = 1._dp
          a1   = dispersion_env%a1
          a2   = dispersion_env%a2
          eps_cn = dispersion_env%eps_cn
          e6tot = 0._dp
          e8tot = 0._dp
          e9tot = 0._dp
          ! damping type
          idmp = 0
          IF ( dispersion_env%pp_type == vdw_pairpot_dftd3 ) THEN
             idmp = 1
          ELSE IF ( dispersion_env%pp_type == vdw_pairpot_dftd3bj ) THEN
             idmp = 2
          END IF

          IF(unit_nr>0) THEN
            WRITE(unit_nr,*) " Scaling parameter (s6) ",s6
            WRITE(unit_nr,*) " Scaling parameter (s8) ",s8
            IF ( dispersion_env%pp_type == vdw_pairpot_dftd3 ) THEN
               WRITE(unit_nr,*) " Zero Damping parameter (sr6)",rs6
               WRITE(unit_nr,*) " Zero Damping parameter (sr8)",rs8
            ELSE IF ( dispersion_env%pp_type == vdw_pairpot_dftd3bj ) THEN
               WRITE(unit_nr,*) " BJ Damping parameter (a1) ",a1
               WRITE(unit_nr,*) " BJ Damping parameter (a2) ",a2
            END IF
            WRITE(unit_nr,*) " Cutoff coordination numbers",eps_cn
            IF (dispersion_env%lrc) THEN
              WRITE(unit_nr,*) " Apply a long range correction"
            END IF
          END IF
          ! Calculate coordination numbers
          NULLIFY (particle_set)
          CALL get_qs_env(qs_env=qs_env,particle_set=particle_set,error=error)
          natom = SIZE (particle_set)
          ALLOCATE (cnumbers(natom),STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          cnumbers = 0._dp
          IF(calculate_forces .OR. debugall) THEN
            ALLOCATE (dcnum(natom),STAT=istat)
            CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
            dcnum(:)%neighbors=0
            DO iatom=1,natom
              ALLOCATE (dcnum(iatom)%nlist(10),dcnum(iatom)%dvals(10),dcnum(iatom)%rik(3,10),STAT=istat)
              CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
            END DO
          END IF
          sab_cn => dispersion_env%sab_cn
          CALL neighbor_list_iterator_create(nl_iterator,sab_cn)
          DO WHILE (neighbor_list_iterate(nl_iterator)==0)
             CALL get_iterator_info(nl_iterator,ikind=ikind,jkind=jkind,iatom=iatom,jatom=jatom,r=rij)
             IF (ghost(ikind).OR.ghost(jkind)) CYCLE
             rcc=SQRT(rij(1)*rij(1)+rij(2)*rij(2)+rij(3)*rij(3))
             IF (rcc > 1.e-6_dp) THEN
                za = atomnumber(ikind)
                zb = atomnumber(jkind)
                rcova = dispersion_env%rcov(za)
                rcovb = dispersion_env%rcov(zb)
                CALL cnparam_d3(rcc,rcova,rcovb,dispersion_env%k1,cnab,dcnab)
                IF (cnab > eps_cn) THEN
                  cnumbers(iatom) = cnumbers(iatom) + cnab
                  IF (iatom /= jatom) cnumbers(jatom) = cnumbers(jatom) + cnab
                END IF
                IF(calculate_forces .OR. debugall .AND. cnab > eps_cn) THEN
                  dcnum(iatom)%neighbors=dcnum(iatom)%neighbors+1
                  ni = dcnum(iatom)%neighbors
                  IF (ni > SIZE(dcnum(iatom)%nlist)) THEN
                    CALL reallocate(dcnum(iatom)%nlist,1,2*ni)
                    CALL reallocate(dcnum(iatom)%dvals,1,2*ni)
                    CALL reallocate(dcnum(iatom)%rik,1,3,1,2*ni)
                  END IF
                  dcnum(iatom)%nlist(ni)=jatom
                  dcnum(iatom)%dvals(ni)=dcnab
                  dcnum(iatom)%rik(1:3,ni)=rij(1:3)
                  IF (iatom/=jatom) THEN
                    dcnum(jatom)%neighbors=dcnum(jatom)%neighbors+1
                    nj = dcnum(jatom)%neighbors
                    IF (nj > SIZE(dcnum(jatom)%dvals)) THEN
                      CALL reallocate(dcnum(jatom)%nlist,1,2*nj)
                      CALL reallocate(dcnum(jatom)%dvals,1,2*nj)
                      CALL reallocate(dcnum(jatom)%rik,1,3,1,2*nj)
                    END IF
                    dcnum(jatom)%nlist(nj)=iatom
                    dcnum(jatom)%dvals(nj)=dcnab
                    dcnum(jatom)%rik(1:3,nj)=-rij(1:3)
                  END IF
                END IF
             END IF
          END DO
          CALL neighbor_list_iterator_release(nl_iterator)
          CALL mp_sum(cnumbers,para_env%group)
          ! for parallel runs we have to update dcnum on all processors
          IF(calculate_forces .OR. debugall) THEN
            CALL dcnum_distribute(dcnum,para_env,error)
            IF (unit_nr>0 .AND. SIZE(dcnum)>0) THEN
              WRITE(unit_nr,*)
              WRITE(unit_nr,*) "  ATOM       Coordination   Neighbors"
              DO i=1,natom
                 WRITE(unit_nr,"(I6,F20.10,I12)") i,cnumbers(i),dcnum(i)%neighbors
              END DO
              WRITE(unit_nr,*)
            END IF
          END IF
       END IF

       nab  = 0._dp
       nabc = 0._dp
       IF (dispersion_env%doabc) THEN
          rcc = 2._dp*dispersion_env%rc_disp
          CALL get_cell(cell=cell,periodic=periodic)
          sab_max(1) = rcc/plane_distance(1,0,0,cell)
          sab_max(2) = rcc/plane_distance(0,1,0,cell)
          sab_max(3) = rcc/plane_distance(0,0,1,cell)
          ncell(:) = (INT(sab_max(:)) + 1)*periodic(:)
          IF (unit_nr>0) THEN
            WRITE(unit_nr,*) " Calculate C9 Terms"
            WRITE(unit_nr,"(A,T20,I3,A,I3)") "  Search in cells ",-ncell(1),":",ncell(1)
            WRITE(unit_nr,"(T20,I3,A,I3)") -ncell(2),":",ncell(2)
            WRITE(unit_nr,"(T20,I3,A,I3)") -ncell(3),":",ncell(3)
            WRITE(unit_nr,*)
          END IF
          IF (dispersion_env%c9cnst) THEN
            IF (unit_nr>0) WRITE(unit_nr,*) " Use reference coordination numbers for C9 term"
            ALLOCATE (cnumfix(natom),STAT=istat)
            CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
            cnumfix = 0._dp
            ! first use the default values
            DO iatom=1,natom
              ikind=kind_of(iatom)
              CALL get_atomic_kind(atomic_kind_set(ikind),z=za)
              cnumfix(iatom) = dispersion_env%cn(za)
            END DO
            ! now check for changes from default
            IF (ASSOCIATED(dispersion_env%cnkind)) THEN
              DO i=1,SIZE(dispersion_env%cnkind)
                ikind=dispersion_env%cnkind(i)%kind
                cnum=dispersion_env%cnkind(i)%cnum
                CPPrecondition(ikind <= nkind,cp_failure_level,routineP,error,failure)
                CPPrecondition(ikind > 0,cp_failure_level,routineP,error,failure)
                CALL get_atomic_kind(atomic_kind_set(ikind), natom=na,atom_list=atom_list)
                DO ia=1,na
                  iatom=atom_list(ia)
                  cnumfix(iatom) = cnum
                END DO
              END DO
            END IF
            IF (ASSOCIATED(dispersion_env%cnlist)) THEN
              DO i=1,SIZE(dispersion_env%cnlist)
                DO ilist=1,dispersion_env%cnlist(i)%natom
                  iatom = dispersion_env%cnlist(i)%atom(ilist)
                  cnumfix(iatom) = dispersion_env%cnlist(i)%cnum
                END DO
              END DO
            END IF
            IF (unit_nr>0) THEN
              DO i=1,natom
                IF (ABS(cnumbers(i)-cnumfix(i)) > 0.5_dp) THEN
                  WRITE(unit_nr,"(A,T20,A,I6,T41,2F10.3)") "  Difference in CN ","Atom:",&
                        i,cnumbers(i),cnumfix(i)
                END IF
              END DO
              WRITE(unit_nr,*)
            END IF
          END IF
       END IF

       evdw = 0._dp

       sab_vdw => dispersion_env%sab_vdw

       nkind = SIZE(atomic_kind_set)

       CALL neighbor_list_iterator_create(nl_iterator,sab_vdw)
       DO WHILE (neighbor_list_iterate(nl_iterator)==0)
          CALL get_iterator_info(nl_iterator,ikind=ikind,jkind=jkind,iatom=iatom,jatom=jatom,r=rij)
          IF (ghost(ikind).OR.ghost(jkind)) CYCLE
          IF (.NOT.(dodisp(ikind).AND.dodisp(jkind))) CYCLE
          za = atomnumber(ikind)
          zb = atomnumber(jkind)
          ! vdW potential
          dr = SQRT(SUM(rij(:)**2))
          IF (dr <= rcut) THEN
             nab = nab+1._dp
             fac = 1._dp
             IF(iatom==jatom) fac = 0.5_dp
             IF ( disp_a%type == dftd2_pp .AND. dr > 0.001_dp ) THEN
                c6 = SQRT(c6d2(ikind)*c6d2(jkind))
                rcc = radd2(ikind) + radd2(jkind)
                er = EXP(-dd*(dr/rcc-1._dp))
                fdmp = 1._dp/(1._dp+er)
                xp = s6*c6/dr**6
                evdw = evdw - xp*fdmp*fac
                IF(calculate_forces) THEN
                   dfdmp = dd/rcc * er * fdmp*fdmp
                   devdw = -xp*(-6._dp*fdmp/dr + dfdmp)
                   fdij(:) = devdw * rij(:)/dr * fac
                   atom_a = atom_of_kind(iatom)
                   atom_b = atom_of_kind(jatom)
                   force(ikind)%dispersion(:,atom_a) = force(ikind)%dispersion(:,atom_a) - fdij(:)
                   force(jkind)%dispersion(:,atom_b) = force(jkind)%dispersion(:,atom_b) + fdij(:)
                   IF ( use_virial ) THEN
                      CALL virial_pair_force ( virial%pv_virial, -1._dp, fdij, rij, error)
                   END IF
                   IF(atstress) THEN
                      CALL virial_pair_force ( atstr(:,:,iatom), -0.5_dp, fdij, rij, error)
                      CALL virial_pair_force ( atstr(:,:,jatom), -0.5_dp, fdij, rij, error)
                   END IF
                END IF
                IF(atenergy) THEN
                  atener(iatom) = atener(iatom) - 0.5_dp*xp*fdmp*fac
                  atener(jatom) = atener(jatom) - 0.5_dp*xp*fdmp*fac
                END IF
             ELSE IF ( disp_a%type == dftd3_pp .AND. dr > 0.001_dp ) THEN
                ! C6 coefficient
                CALL getc6(maxc,max_elem,dispersion_env%c6ab,dispersion_env%maxci,za,zb,&
                           cnumbers(iatom),cnumbers(jatom),dispersion_env%k3,c6,dc6a,dc6b)
                c8 = 3.0d0*c6*dispersion_env%r2r4(za)*dispersion_env%r2r4(zb)
                dc8a = 3.0d0*dc6a*dispersion_env%r2r4(za)*dispersion_env%r2r4(zb)
                dc8b = 3.0d0*dc6b*dispersion_env%r2r4(za)*dispersion_env%r2r4(zb)
                r6 = dr**6
                r8 = r6*dr*dr
                ! damping
                IF(idmp==1) THEN
                   ! zero
                   CALL damping_d3(dr,dispersion_env%r0ab(za,zb),rs6,alp6,rcut,fdab6,dfdab6)
                   CALL damping_d3(dr,dispersion_env%r0ab(za,zb),rs8,alp8,rcut,fdab8,dfdab8)
                   e6 = s6*fac*c6*fdab6/r6
                   e8 = s8*fac*c8*fdab8/r8
                ELSE IF(idmp==2) THEN
                   ! BJ
                   r0ab = SQRT(3.0d0*dispersion_env%r2r4(za)*dispersion_env%r2r4(zb))
                   f0ab = a1*r0ab + a2
                   fdab6 = 1.0_dp/(r6+f0ab**6)
                   fdab8 = 1.0_dp/(r8+f0ab**8)
                   e6 = s6*fac*c6*fdab6
                   e8 = s8*fac*c8*fdab8
                ELSE
                   CALL cp_unimplemented_error(fromWhere=routineP, &
                        message="Unknown DFT-D3 damping function:", &
                        error=error, error_level=cp_failure_level)
                END IF
                evdw = evdw - e6 - e8
                e6tot = e6tot - e6
                e8tot = e8tot - e8
                IF(atenergy) THEN
                  atener(iatom) = atener(iatom) - 0.5_dp*(e6+e8)
                  atener(jatom) = atener(jatom) - 0.5_dp*(e6+e8)
                END IF
                IF(calculate_forces) THEN
                   ! damping
                   IF(idmp==1) THEN
                      ! zero
                      de6 = -s6*c6/r6*(dfdab6-6._dp*fdab6/dr)
                      de8 = -s8*c8/r8*(dfdab8-8._dp*fdab8/dr)
                   ELSE IF(idmp==2) THEN
                      ! BJ  
                      de6 = s6*c6*fdab6*fdab6*6.0_dp*dr**5
                      de8 = s8*c8*fdab8*fdab8*8.0_dp*dr**7
                   ELSE
                      CALL cp_unimplemented_error(fromWhere=routineP, &
                           message="Unknown DFT-D3 damping function:", &
                           error=error, error_level=cp_failure_level)
                   END IF
                   fdij(:) = (de6+de8) * rij(:)/dr * fac
                   atom_a = atom_of_kind(iatom)
                   atom_b = atom_of_kind(jatom)
                   force(ikind)%dispersion(:,atom_a) = force(ikind)%dispersion(:,atom_a) - fdij(:)
                   force(jkind)%dispersion(:,atom_b) = force(jkind)%dispersion(:,atom_b) + fdij(:)
                   IF ( use_virial ) THEN
                      CALL virial_pair_force ( virial%pv_virial, -1._dp, fdij, rij, error)
                   END IF
                   IF(atstress) THEN
                      CALL virial_pair_force ( atstr(:,:,iatom), -0.5_dp, fdij, rij, error)
                      CALL virial_pair_force ( atstr(:,:,jatom), -0.5_dp, fdij, rij, error)
                   END IF
                   ! forces from the r-dependence of the coordination numbers
                   IF(idmp==1) THEN
                      ! zero
                      dc6a = -s6*fac*dc6a*fdab6/r6
                      dc6b = -s6*fac*dc6b*fdab6/r6
                      dc8a = -s8*fac*dc8a*fdab8/r8
                      dc8b = -s8*fac*dc8b*fdab8/r8
                   ELSE IF(idmp==2) THEN
                      ! BJ  
                      dc6a = -s6*fac*dc6a*fdab6
                      dc6b = -s6*fac*dc6b*fdab6
                      dc8a = -s8*fac*dc8a*fdab8
                      dc8b = -s8*fac*dc8b*fdab8
                   ELSE
                      CALL cp_unimplemented_error(fromWhere=routineP, &
                           message="Unknown DFT-D3 damping function:", &
                           error=error, error_level=cp_failure_level)
                   END IF
                   DO i=1,dcnum(iatom)%neighbors
                     katom=dcnum(iatom)%nlist(i)
                     kkind=kind_of(katom)
                     rik = dcnum(iatom)%rik(:,i)
                     drk = SQRT(SUM(rik(:)**2))
                     fdik(:) = (dc6a+dc8a)*dcnum(iatom)%dvals(i)*rik(:)/drk
                     atom_c = atom_of_kind(katom)
                     force(ikind)%dispersion(:,atom_a) = force(ikind)%dispersion(:,atom_a) - fdik(:)
                     force(kkind)%dispersion(:,atom_c) = force(kkind)%dispersion(:,atom_c) + fdik(:)
                     IF ( use_virial ) THEN
                       CALL virial_pair_force ( virial%pv_virial, -1._dp, fdik, rik, error)
                     END IF
                     IF(atstress) THEN
                        CALL virial_pair_force ( atstr(:,:,iatom), -0.5_dp, fdik, rik, error)
                        CALL virial_pair_force ( atstr(:,:,katom), -0.5_dp, fdik, rik, error)
                     END IF
                   END DO
                   DO i=1,dcnum(jatom)%neighbors
                     katom=dcnum(jatom)%nlist(i)
                     kkind=kind_of(katom)
                     rik = dcnum(jatom)%rik(:,i)
                     drk = SQRT(SUM(rik(:)**2))
                     fdik(:) = (dc6b+dc8b)*dcnum(jatom)%dvals(i)*rik(:)/drk
                     atom_c = atom_of_kind(katom)
                     force(jkind)%dispersion(:,atom_b) = force(jkind)%dispersion(:,atom_b) - fdik(:)
                     force(kkind)%dispersion(:,atom_c) = force(kkind)%dispersion(:,atom_c) + fdik(:)
                     IF ( use_virial ) THEN
                       CALL virial_pair_force ( virial%pv_virial, -1._dp, fdik, rik, error)
                     END IF
                     IF(atstress) THEN
                        CALL virial_pair_force ( atstr(:,:,jatom), -0.5_dp, fdik, rik, error)
                        CALL virial_pair_force ( atstr(:,:,katom), -0.5_dp, fdik, rik, error)
                     END IF
                   END DO
                END IF
                IF (dispersion_env%doabc) THEN
                   CALL get_iterator_info(nl_iterator,cell=cell_b)
                   hashb = cellhash(cell_b,ncell,error)
                   is000 = (ALL(cell_b == 0))
                   rb0(:) = MATMUL(cell%hmat,cell_b)
                   ra(:) = pbc(particle_set(iatom)%r(:),cell)
                   rb(:) = pbc(particle_set(jatom)%r(:),cell)+rb0
                   r2ab = SUM( (ra-rb)**2 )
                   DO icx=-ncell(1),ncell(1)
                      DO icy=-ncell(2),ncell(2)
                         DO icz=-ncell(3),ncell(3)
                           cell_c(1) = icx
                           cell_c(2) = icy
                           cell_c(3) = icz
                           hashc = cellhash(cell_c,ncell,error)
                           IF(is000 .AND. (ALL(cell_c == 0))) THEN
                              ! CASE 1: all atoms in (000), use only ordered triples
                              kstart=MAX(jatom+1,iatom+1)
                              fac0 = 1.0_dp
                           ELSE IF(is000) THEN
                              ! CASE 2: AB in (000), C in other cell
                              !         This case covers also all instances with BC in same
                              !         cell not (000)
                              kstart=1
                              fac0 = 1.0_dp
                           ELSE
                              ! These are case 2 again, cycle
                              IF ( hashc == hashb ) CYCLE
                              IF ( ALL(cell_c == 0) ) CYCLE
                              ! CASE 3: A in (000) and B and C in different cells
                              kstart=1
                              fac0 = 1.0_dp/3.0_dp
                           END IF
                           rc0=MATMUL(cell%hmat,cell_c)
                           DO katom=kstart,natom
                              kkind=kind_of(katom)
                              IF (ghost(kkind) .OR. .NOT.dodisp(kkind)) CYCLE
                              rc(:) = rcpbc(:,katom)+rc0(:)
                              r2bc = SUM( (rb-rc)**2 )
                              IF (r2bc >= rc2) CYCLE
                              r2ca = SUM( (rc-ra)**2 )
                              IF (r2ca >= rc2) CYCLE
                              r2abc = r2ab*r2bc*r2ca
                              IF (r2abc <= 0.001_dp) CYCLE
                              ! this is an empirical scaling
                              IF (r2abc <= 0.01*rc2*rc2*rc2) THEN
                                 kkind=kind_of(katom)
                                 atom_c = atom_of_kind(katom)
                                 zc=atomnumber(kkind)
                                 ! avoid double counting!
                                 fac = 1._dp
                                 IF(iatom==jatom .OR. iatom==katom .OR. jatom==katom) fac = 0.5_dp
                                 IF(iatom==jatom .AND. iatom==katom) fac = 1._dp/3._dp
                                 fac = fac*fac0
                                 nabc=nabc+1._dp
                                 IF (dispersion_env%c9cnst) THEN
                                   CALL getc6(maxc,max_elem,dispersion_env%c6ab,dispersion_env%maxci,za,zb,&
                                      cnumfix(iatom),cnumfix(jatom),dispersion_env%k3,cc6ab,dcc6aba,dcc6abb)
                                   CALL getc6(maxc,max_elem,dispersion_env%c6ab,dispersion_env%maxci,zb,zc,&
                                      cnumfix(jatom),cnumfix(katom),dispersion_env%k3,cc6bc,dcc6bcb,dcc6bcc)
                                   CALL getc6(maxc,max_elem,dispersion_env%c6ab,dispersion_env%maxci,zc,za,&
                                      cnumfix(katom),cnumfix(iatom),dispersion_env%k3,cc6ca,dcc6cac,dcc6caa)
                                 ELSE
                                   CALL getc6(maxc,max_elem,dispersion_env%c6ab,dispersion_env%maxci,za,zb,&
                                      cnumbers(iatom),cnumbers(jatom),dispersion_env%k3,cc6ab,dcc6aba,dcc6abb)
                                   CALL getc6(maxc,max_elem,dispersion_env%c6ab,dispersion_env%maxci,zb,zc,&
                                      cnumbers(jatom),cnumbers(katom),dispersion_env%k3,cc6bc,dcc6bcb,dcc6bcc)
                                   CALL getc6(maxc,max_elem,dispersion_env%c6ab,dispersion_env%maxci,zc,za,&
                                      cnumbers(katom),cnumbers(iatom),dispersion_env%k3,cc6ca,dcc6cac,dcc6caa)
                                 END IF
                                 c9 = -SQRT(cc6ab*cc6bc*cc6ca)
                                 rabc = r2abc**(1._dp/6._dp)
                                 r0 = (dispersion_env%r0ab(za,zb)*dispersion_env%r0ab(zb,zc)*&
                                       dispersion_env%r0ab(zc,za))**(1._dp/3._dp)
                                 CALL damping_d3(rabc,r0,4._dp/3._dp,alp6,rcut,fdabc,dfdabc)
                                 s1 = r2ab+r2bc-r2ca
                                 s2 = r2ab+r2ca-r2bc
                                 s3 = r2ca+r2bc-r2ab
                                 ang=0.375_dp*s1*s2*s3/r2abc+1.0_dp

                                 e9=s9*fac*fdabc*c9*ang/r2abc**1.50d0
                                 evdw = evdw - e9
                                 e9tot = e9tot - e9
                                 IF(atenergy) THEN
                                   atener(iatom) = atener(iatom) - e9/3._dp
                                   atener(jatom) = atener(jatom) - e9/3._dp
                                   atener(katom) = atener(katom) - e9/3._dp
                                 END IF

                                 IF(calculate_forces) THEN
                                   rab = rb-ra; rbc = rc-rb; rca = ra-rc
                                   de91 = s9*fac*dfdabc*c9*ang/r2abc**1.50d0
                                   de91 = -de91/3._dp*rabc + 3._dp*e9
                                   dea = s9*fac*fdabc*c9/r2abc**2.50d0 * 0.75_dp
                                   fdij(:) = de91*rab(:)/r2ab
                                   fdij(:) = fdij(:) + dea*s1*s2*s3 * rab(:)/r2ab
                                   fdij(:) = fdij(:) - dea*(s2*s3+s1*s3-s1*s2) * rab(:)
                                   force(ikind)%dispersion(:,atom_a) = force(ikind)%dispersion(:,atom_a) - fdij(:)
                                   force(jkind)%dispersion(:,atom_b) = force(jkind)%dispersion(:,atom_b) + fdij(:)
                                   IF ( use_virial ) THEN
                                      CALL virial_pair_force ( virial%pv_virial, -1._dp, fdij, rab, error)
                                   END IF
                                   IF(atstress) THEN
                                      CALL virial_pair_force ( atstr(:,:,iatom), -0.5_dp, fdij, rab, error)
                                      CALL virial_pair_force ( atstr(:,:,jatom), -0.5_dp, fdij, rab, error)
                                   END IF
                                   fdij(:) = de91*rbc(:)/r2bc
                                   fdij(:) = fdij(:) + dea*s1*s2*s3 * rbc(:)/r2bc
                                   fdij(:) = fdij(:) - dea*(s2*s3-s1*s3+s1*s2) * rbc(:)
                                   force(jkind)%dispersion(:,atom_b) = force(jkind)%dispersion(:,atom_b) - fdij(:)
                                   force(kkind)%dispersion(:,atom_c) = force(kkind)%dispersion(:,atom_c) + fdij(:)
                                   IF ( use_virial ) THEN
                                      CALL virial_pair_force ( virial%pv_virial, -1._dp, fdij, rbc, error)
                                   END IF
                                   IF(atstress) THEN
                                      CALL virial_pair_force ( atstr(:,:,jatom), -0.5_dp, fdij, rbc, error)
                                      CALL virial_pair_force ( atstr(:,:,katom), -0.5_dp, fdij, rbc, error)
                                   END IF
                                   fdij(:) = de91*rca(:)/r2ca
                                   fdij(:) = fdij(:) + dea*s1*s2*s3 * rca(:)/r2ca
                                   fdij(:) = fdij(:) - dea*(-s2*s3+s1*s3+s1*s2) * rca(:)
                                   force(kkind)%dispersion(:,atom_c) = force(kkind)%dispersion(:,atom_c) - fdij(:)
                                   force(ikind)%dispersion(:,atom_a) = force(ikind)%dispersion(:,atom_a) + fdij(:)
                                   IF ( use_virial ) THEN
                                      CALL virial_pair_force ( virial%pv_virial, -1._dp, fdij, rca, error)
                                   END IF
                                   IF(atstress) THEN
                                      CALL virial_pair_force ( atstr(:,:,iatom), -0.5_dp, fdij, rca, error)
                                      CALL virial_pair_force ( atstr(:,:,katom), -0.5_dp, fdij, rca, error)
                                   END IF
                                   IF (.NOT. dispersion_env%c9cnst) THEN
                                     ! forces from the r-dependence of the coordination numbers
                                     ! atomic stress not implemented
                                     de91  = 0.5_dp*e9/cc6ab
                                     de921 = de91*dcc6aba
                                     de922 = de91*dcc6abb
                                     DO i=1,dcnum(iatom)%neighbors
                                       latom=dcnum(iatom)%nlist(i)
                                       lkind=kind_of(latom)
                                       rik(1) = dcnum(iatom)%rik(1,i)
                                       rik(2) = dcnum(iatom)%rik(2,i)
                                       rik(3) = dcnum(iatom)%rik(3,i)
                                       drk = SQRT(rik(1)*rik(1)+rik(2)*rik(2)+rik(3)*rik(3))
                                       fdik(:) = -de921*dcnum(iatom)%dvals(i)*rik(:)/drk
                                       atom_d = atom_of_kind(latom)
                                       force(ikind)%dispersion(:,atom_a) = force(ikind)%dispersion(:,atom_a) - fdik(:)
                                       force(lkind)%dispersion(:,atom_d) = force(lkind)%dispersion(:,atom_d) + fdik(:)
                                       IF ( use_virial ) THEN
                                         CALL virial_pair_force ( virial%pv_virial, -1._dp, fdik, rik, error)
                                       END IF
                                     END DO
                                     DO i=1,dcnum(jatom)%neighbors
                                       latom=dcnum(jatom)%nlist(i)
                                       lkind=kind_of(latom)
                                       rik(1) = dcnum(jatom)%rik(1,i)
                                       rik(2) = dcnum(jatom)%rik(2,i)
                                       rik(3) = dcnum(jatom)%rik(3,i)
                                       drk = SQRT(rik(1)*rik(1)+rik(2)*rik(2)+rik(3)*rik(3))
                                       fdik(:) = -de922*dcnum(jatom)%dvals(i)*rik(:)/drk
                                       atom_d = atom_of_kind(latom)
                                       force(jkind)%dispersion(:,atom_b) = force(jkind)%dispersion(:,atom_b) - fdik(:)
                                       force(lkind)%dispersion(:,atom_d) = force(lkind)%dispersion(:,atom_d) + fdik(:)
                                       IF ( use_virial ) THEN
                                         CALL virial_pair_force ( virial%pv_virial, -1._dp, fdik, rik, error)
                                       END IF
                                     END DO

                                     de91  = 0.5_dp*e9/cc6bc
                                     de921 = de91*dcc6bcb
                                     de922 = de91*dcc6bcc
                                     DO i=1,dcnum(jatom)%neighbors
                                       latom=dcnum(jatom)%nlist(i)
                                       lkind=kind_of(latom)
                                       rik(1) = dcnum(jatom)%rik(1,i)
                                       rik(2) = dcnum(jatom)%rik(2,i)
                                       rik(3) = dcnum(jatom)%rik(3,i)
                                       drk = SQRT(rik(1)*rik(1)+rik(2)*rik(2)+rik(3)*rik(3))
                                       fdik(:) = -de921*dcnum(jatom)%dvals(i)*rik(:)/drk
                                       atom_d = atom_of_kind(latom)
                                       force(jkind)%dispersion(:,atom_b) = force(jkind)%dispersion(:,atom_b) - fdik(:)
                                       force(lkind)%dispersion(:,atom_d) = force(lkind)%dispersion(:,atom_d) + fdik(:)
                                       IF ( use_virial ) THEN
                                         CALL virial_pair_force ( virial%pv_virial, -1._dp, fdik, rik, error)
                                       END IF
                                     END DO
                                     DO i=1,dcnum(katom)%neighbors
                                       latom=dcnum(katom)%nlist(i)
                                       lkind=kind_of(latom)
                                       rik(1) = dcnum(katom)%rik(1,i)
                                       rik(2) = dcnum(katom)%rik(2,i)
                                       rik(3) = dcnum(katom)%rik(3,i)
                                       drk = SQRT(rik(1)*rik(1)+rik(2)*rik(2)+rik(3)*rik(3))
                                       fdik(:) = -de922*dcnum(katom)%dvals(i)*rik(:)/drk
                                       atom_d = atom_of_kind(latom)
                                       force(kkind)%dispersion(:,atom_c) = force(kkind)%dispersion(:,atom_c) - fdik(:)
                                       force(lkind)%dispersion(:,atom_d) = force(lkind)%dispersion(:,atom_d) + fdik(:)
                                       IF ( use_virial ) THEN
                                         CALL virial_pair_force ( virial%pv_virial, -1._dp, fdik, rik, error)
                                       END IF
                                     END DO

                                     de91  = 0.5_dp*e9/cc6ca
                                     de921 = de91*dcc6cac
                                     de922 = de91*dcc6caa
                                     DO i=1,dcnum(katom)%neighbors
                                       latom=dcnum(katom)%nlist(i)
                                       lkind=kind_of(latom)
                                       rik(1) = dcnum(katom)%rik(1,i)
                                       rik(2) = dcnum(katom)%rik(2,i)
                                       rik(3) = dcnum(katom)%rik(3,i)
                                       drk = SQRT(rik(1)*rik(1)+rik(2)*rik(2)+rik(3)*rik(3))
                                       fdik(:) = -de921*dcnum(katom)%dvals(i)*rik(:)/drk
                                       atom_d = atom_of_kind(latom)
                                       force(kkind)%dispersion(:,atom_c) = force(kkind)%dispersion(:,atom_c) - fdik(:)
                                       force(lkind)%dispersion(:,atom_d) = force(lkind)%dispersion(:,atom_d) + fdik(:)
                                       IF ( use_virial ) THEN
                                         CALL virial_pair_force ( virial%pv_virial, -1._dp, fdik, rik, error)
                                       END IF
                                     END DO
                                     DO i=1,dcnum(iatom)%neighbors
                                       latom=dcnum(iatom)%nlist(i)
                                       lkind=kind_of(latom)
                                       rik(1) = dcnum(iatom)%rik(1,i)
                                       rik(2) = dcnum(iatom)%rik(2,i)
                                       rik(3) = dcnum(iatom)%rik(3,i)
                                       drk = SQRT(rik(1)*rik(1)+rik(2)*rik(2)+rik(3)*rik(3))
                                       fdik(:) = -de922*dcnum(iatom)%dvals(i)*rik(:)/drk
                                       atom_d = atom_of_kind(latom)
                                       force(ikind)%dispersion(:,atom_a) = force(ikind)%dispersion(:,atom_a) - fdik(:)
                                       force(lkind)%dispersion(:,atom_d) = force(lkind)%dispersion(:,atom_d) + fdik(:)
                                       IF ( use_virial ) THEN
                                         CALL virial_pair_force ( virial%pv_virial, -1._dp, fdik, rik, error)
                                       END IF
                                     END DO
                                   END IF

                                 END IF

                              END IF
                           END DO
                         END DO
                      END DO
                   END DO

                END IF
             END IF
          END IF
       END DO
       CALL neighbor_list_iterator_release(nl_iterator)

       elrc6 = 0._dp
       elrc8 = 0._dp
       elrc9 = 0._dp
       IF ( dispersion_env%pp_type == vdw_pairpot_dftd3 .OR.&
            dispersion_env%pp_type == vdw_pairpot_dftd3bj ) THEN
          ! Long range correction (atomic contributions not implemented)
          IF (dispersion_env%lrc) THEN
             ALLOCATE (cnkind(nkind),STAT=istat)
             CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
             cnkind = 0._dp
             ! first use the default values
             DO ikind=1,nkind
               CALL get_atomic_kind(atomic_kind_set(ikind), z=za)
               cnkind(ikind) = dispersion_env%cn(za)
             END DO
             ! now check for changes from default
             IF (ASSOCIATED(dispersion_env%cnkind)) THEN
               DO i=1,SIZE(dispersion_env%cnkind)
                 ikind=dispersion_env%cnkind(i)%kind
                 cnkind(ikind)=dispersion_env%cnkind(i)%cnum
               END DO
             END IF
             DO ikind=1,nkind
                CALL get_atomic_kind(atomic_kind_set(ikind), natom=na, z=za)
                CALL get_qs_kind(qs_kind_set(ikind), dispersion=disp_a, ghost=ghost_a)
                IF (.NOT.disp_a%defined .OR. ghost_a) CYCLE
                DO jkind=1,nkind
                   CALL get_atomic_kind(atomic_kind_set(jkind), natom=nb, z=zb)
                   CALL get_qs_kind(qs_kind_set(jkind), dispersion=disp_b, ghost=ghost_b)
                   IF (.NOT.disp_b%defined .OR. ghost_b) CYCLE
                   CALL getc6(maxc,max_elem,dispersion_env%c6ab,dispersion_env%maxci,za,zb,&
                      cnkind(ikind),cnkind(jkind),dispersion_env%k3,cc6ab,dcc6aba,dcc6abb)
                   elrc6 = elrc6 - s6*twopi*REAL(na*nb,KIND=dp)*cc6ab/(3._dp*rcut**3*cell%deth)
                   c8 = 3.0d0*cc6ab*dispersion_env%r2r4(za)*dispersion_env%r2r4(zb)
                   elrc8 = elrc8 - s8*twopi*REAL(na*nb,KIND=dp)*c8/(5._dp*rcut**5*cell%deth)
                   IF (dispersion_env%doabc) THEN
                      DO kkind=1,nkind
                         CALL get_atomic_kind(atomic_kind_set(kkind), natom=nc, z=zc)
                         CALL get_qs_kind(qs_kind_set(kkind), dispersion=disp_c, ghost=ghost_c)
                         IF (.NOT.disp_c%defined .OR. ghost_c) CYCLE
                         CALL getc6(maxc,max_elem,dispersion_env%c6ab,dispersion_env%maxci,za,zb,&
                            cnkind(ikind),cnkind(jkind),dispersion_env%k3,cc6ab,dcc6aba,dcc6abb)
                         CALL getc6(maxc,max_elem,dispersion_env%c6ab,dispersion_env%maxci,zc,za,&
                            cnkind(kkind),cnkind(ikind),dispersion_env%k3,cc6ca,dcc6aba,dcc6abb)
                         CALL getc6(maxc,max_elem,dispersion_env%c6ab,dispersion_env%maxci,zb,zc,&
                            cnkind(jkind),cnkind(kkind),dispersion_env%k3,cc6bc,dcc6aba,dcc6abb)
                         c9 = -SQRT(cc6ab*cc6bc*cc6ca)
                         elrc9 = elrc9 - s9*64._dp*twopi*REAL(na*nb*nc,KIND=dp)*c9/(6._dp*rcut**3*cell%deth**2)
                      END DO
                   END IF
                END DO
             END DO
             IF ( use_virial ) THEN
               DO i=1,3
                 virial%pv_virial(i,i) = virial%pv_virial(i,i) + (elrc6+elrc8+2._dp*elrc9)
               END DO
             END IF
             DEALLOCATE (cnkind,STAT=istat)
             CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          END IF
       END IF

       IF ( dispersion_env%pp_type == vdw_pairpot_dftd3 .OR.&
            dispersion_env%pp_type == vdw_pairpot_dftd3bj ) THEN
          DEALLOCATE (cnumbers,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          IF ( dispersion_env%doabc .AND. dispersion_env%c9cnst ) THEN
            DEALLOCATE (cnumfix,STAT=istat)
            CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          END IF
          IF(calculate_forces .OR. debugall) THEN
            DO iatom=1,natom
              DEALLOCATE (dcnum(iatom)%nlist,dcnum(iatom)%dvals,dcnum(iatom)%rik,STAT=istat)
              CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
            END DO
            DEALLOCATE (dcnum,STAT=istat)
            CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          END IF
       END IF

       ! set dispersion energy
       CALL mp_sum(evdw,para_env%group)
       evdw = evdw + (elrc6+elrc8+elrc9)
       energy = evdw

       IF ( ( dispersion_env%pp_type == vdw_pairpot_dftd3 .OR. &
              dispersion_env%pp_type == vdw_pairpot_dftd3bj ) .AND. debugall ) THEN
         CALL mp_sum(e6tot,para_env%group)
         CALL mp_sum(e8tot,para_env%group)
         CALL mp_sum(e9tot,para_env%group)
         CALL mp_sum(nab,para_env%group)
         CALL mp_sum(nabc,para_env%group)
         e6tot = e6tot+elrc6
         e8tot = e8tot+elrc8
         e9tot = e9tot+elrc9
         IF(unit_nr>0) THEN
            WRITE(unit_nr,"(A,F20.0)") "  E6 vdW terms              :",nab
            WRITE(unit_nr,*) " E6 vdW energy [au/kcal]   :",e6tot,e6tot*kcalmol
            WRITE(unit_nr,*) " E8 vdW energy [au/kcal]   :",e8tot,e8tot*kcalmol
            WRITE(unit_nr,*) " %E8 on total vdW energy   :",e8tot/evdw*100._dp
            WRITE(unit_nr,"(A,F20.0)") "  E9 vdW terms              :",nabc
            WRITE(unit_nr,*) " E9 vdW energy [au/kcal]   :",e9tot,e9tot*kcalmol
            WRITE(unit_nr,*) " %E9 on total vdW energy   :",e9tot/evdw*100._dp
            IF (dispersion_env%lrc) THEN
              WRITE(unit_nr,*) " E LRC C6 [au/kcal]        :",elrc6,elrc6*kcalmol
              WRITE(unit_nr,*) " E LRC C8 [au/kcal]        :",elrc8,elrc8*kcalmol
              WRITE(unit_nr,*) " E LRC C9 [au/kcal]        :",elrc9,elrc9*kcalmol
            END IF
         END IF
       END IF
       IF(unit_nr>0) THEN
         WRITE(unit_nr,*) " Total vdW energy [au]     :",evdw
         WRITE(unit_nr,*) " Total vdW energy [kcal]   :",evdw * kcalmol
         WRITE(unit_nr,*)
       END IF
       IF(calculate_forces .AND. debugall) THEN
         IF(unit_nr>0) THEN
            WRITE(unit_nr,*) " Dispersion Forces         "
            WRITE(unit_nr,*) " Atom   Kind                            Forces    "
         END IF
         gnorm=0._dp
         DO iatom=1,natom
           ikind=kind_of(iatom)
           atom_a = atom_of_kind(iatom)
           fdij(1:3) = force(ikind)%dispersion(:,atom_a)
           CALL mp_sum(fdij,para_env%group)
           gnorm = gnorm + SUM(ABS(fdij))
           IF(unit_nr>0) WRITE(unit_nr,"(i5,i7,3F20.14)") iatom,ikind,fdij
         END DO
         IF(unit_nr>0) THEN
           WRITE(unit_nr,*)
           WRITE(unit_nr,*) "|G| = ",gnorm
           WRITE(unit_nr,*)
         END IF
         IF ( use_virial ) THEN
           dvirial = virial%pv_virial - dvirial
           CALL mp_sum(dvirial,para_env%group)
           IF(unit_nr>0) THEN
              WRITE(unit_nr,*) "Stress Tensor (dispersion)"
              WRITE(unit_nr,"(3G20.12)") dvirial
              WRITE(unit_nr,*) "  Tr(P)/3 :  ",(dvirial(1,1)+dvirial(2,2)+dvirial(3,3) )/3._dp
              WRITE(unit_nr,*)
           END IF
         END IF
       END IF

       IF(calculate_forces .OR. debugall .OR. dispersion_env%doabc) THEN
          DEALLOCATE(atom_of_kind,kind_of,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       END IF

       DEALLOCATE (dodisp,ghost,atomnumber,rcpbc,radd2,c6d2,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

       IF(ASSOCIATED(dispersion_env%dftd_section)) THEN
          CALL cp_print_key_finished_output(unit_nr,logger,dispersion_env%dftd_section,&
                                            "PRINT_DFTD", error=error)
       END IF

    END IF

    CALL timestop(handle)

  END SUBROUTINE calculate_dispersion_pairpot

! *****************************************************************************
!> \brief ...
!> \param cell ...
!> \param ncell ...
!> \param error ...
!> \retval hash ...
! *****************************************************************************
  FUNCTION cellhash(cell,ncell,error) RESULT(hash)
    INTEGER, DIMENSION(3), INTENT(IN)        :: cell, ncell
    TYPE(cp_error_type), INTENT(inout)       :: error
    INTEGER                                  :: hash

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cellhash', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ix, iy, iz, nx, ny, nz
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ALL(ABS(cell)<=ncell),cp_failure_level,routineP,error,failure)

    ix = cell(1)
    IF (ix/=0) THEN
      ix = 2*ABS(ix)-(1+SIGN(1,ix))/2
    END IF
    iy = cell(2)
    IF (iy/=0) THEN
      iy = 2*ABS(iy)-(1+SIGN(1,iy))/2
    END IF
    iz = cell(3)
    IF (iz/=0) THEN
      iz = 2*ABS(iz)-(1+SIGN(1,iz))/2
    END IF

    nx=2*ncell(1)+1
    ny=2*ncell(2)+1
    nz=2*ncell(3)+1

    hash = ix*ny*nz + iy*nz + iz + 1

  END FUNCTION cellhash
! *****************************************************************************
!> \brief ...
!> \param z ...
!> \param c6 ...
!> \param r ...
!> \param found ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE dftd2_param(z,c6,r,found,error)

    INTEGER, INTENT(in)                      :: z
    REAL(KIND=dp), INTENT(inout)             :: c6, r
    LOGICAL, INTENT(inout)                   :: found
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'dftd2_param', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), DIMENSION(54), PARAMETER :: c6val = (/0.14_dp, 0.08_dp, &
      1.61_dp, 1.61_dp, 3.13_dp, 1.75_dp, 1.23_dp, 0.70_dp, 0.75_dp, 0.63_dp, &
      5.71_dp, 5.71_dp, 10.79_dp, 9.23_dp, 7.84_dp, 5.57_dp, 5.07_dp, 4.61_dp,&
      10.80_dp, 10.80_dp, 10.80_dp, 10.80_dp, 10.80_dp, 10.80_dp, 10.80_dp, &
      10.80_dp, 10.80_dp, 10.80_dp, 10.80_dp, 10.80_dp, 16.99_dp, 17.10_dp, &
      16.37_dp, 12.64_dp, 12.47_dp, 12.01_dp, 24.67_dp, 24.67_dp, 24.67_dp, &
      24.67_dp, 24.67_dp, 24.67_dp, 24.67_dp, 24.67_dp, 24.67_dp, 24.67_dp, &
      24.67_dp, 24.67_dp, 37.32_dp, 38.71_dp, 38.44_dp, 31.74_dp, 31.50_dp, &
      29.99_dp /)
    REAL(KIND=dp), DIMENSION(54), PARAMETER :: rval = (/1.001_dp, 1.012_dp, &
      0.825_dp, 1.408_dp, 1.485_dp, 1.452_dp, 1.397_dp, 1.342_dp, 1.287_dp, &
      1.243_dp, 1.144_dp, 1.364_dp, 1.639_dp, 1.716_dp, 1.705_dp, 1.683_dp, &
      1.639_dp, 1.595_dp, 1.485_dp, 1.474_dp, 1.562_dp, 1.562_dp, 1.562_dp, &
      1.562_dp, 1.562_dp, 1.562_dp, 1.562_dp, 1.562_dp, 1.562_dp, 1.562_dp, &
      1.650_dp, 1.727_dp, 1.760_dp, 1.771_dp, 1.749_dp, 1.727_dp, 1.628_dp, &
      1.606_dp, 1.639_dp, 1.639_dp, 1.639_dp, 1.639_dp, 1.639_dp, 1.639_dp, &
      1.639_dp, 1.639_dp, 1.639_dp, 1.639_dp, 1.672_dp, 1.804_dp, 1.881_dp, &
      1.892_dp, 1.892_dp, 1.881_dp /)

!
! GRIMME DISPERSION PARAMETERS
! Stefan Grimme, Semiempirical GGA-Type Density Functional Constructed
!                with a Long-Range Dispersion Correction, J. Comp. Chem. 27: 1787-1799 (2006)
! doi:10.1002/jcc.20495
!
! Conversion factor [Jnm^6mol^-1] -> [a.u.] : 17.34527758021901
! Conversion factor [A] -> [a.u.] : 1.889726132885643
!
! C6 values in [Jnm^6/mol]
! vdW radii [A]

    IF ( z > 0 .AND. z <= 54 ) THEN
      found =.TRUE.
      c6 = c6val(z)* 1000._dp*bohr**6/kjmol
      r  = rval(z)*bohr
    ELSE
      found =.FALSE.
    END IF

  END SUBROUTINE dftd2_param
! *****************************************************************************
!> \brief ...
!> \param c6ab ...
!> \param maxci ...
!> \param filename ...
!> \param para_env ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE dftd3_c6_param(c6ab,maxci,filename,para_env,error)

    REAL(KIND=dp), DIMENSION(:, :, :, :, :)  :: c6ab
    INTEGER, DIMENSION(:)                    :: maxci
    CHARACTER(LEN=*)                         :: filename
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'dftd3_c6_param', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: funit, iadr, iat, jadr, jat, &
                                                kk, nl, nlines, nn, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: pars

    failure = .FALSE.
    IF(para_env%source==para_env%mepos) THEN
      ! Read the DFT-D3 C6AB parameters from file "filename"
      CALL open_file(file_name=filename,unit_number=funit,file_form="FORMATTED")
      READ(funit,*) nl, nlines
    END IF
    CALL mp_bcast(nl,para_env%source,para_env%group)
    CALL mp_bcast(nlines,para_env%source,para_env%group)
    ALLOCATE(pars(nl),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF(para_env%source==para_env%mepos) THEN
      READ(funit,*) pars(1:nl)
      CALL close_file(unit_number=funit)
    END IF
    CALL mp_bcast(pars,para_env%source,para_env%group)

    ! Store C6AB coefficients in an array
    c6ab=-1._dp
    maxci=0
    kk=1
    DO nn=1,nlines
       iat=NINT(pars(kk+1))
       jat=NINT(pars(kk+2))
       CALL limit(iat,jat,iadr,jadr)
       maxci(iat)=MAX(maxci(iat),iadr)
       maxci(jat)=MAX(maxci(jat),jadr)
       c6ab(iat,jat,iadr,jadr,1)=pars(kk)
       c6ab(iat,jat,iadr,jadr,2)=pars(kk+3)
       c6ab(iat,jat,iadr,jadr,3)=pars(kk+4)

       c6ab(jat,iat,jadr,iadr,1)=pars(kk)
       c6ab(jat,iat,jadr,iadr,2)=pars(kk+4)
       c6ab(jat,iat,jadr,iadr,3)=pars(kk+3)
       kk=(nn*5)+1
    ENDDO

    DEALLOCATE(pars,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE dftd3_c6_param

! *****************************************************************************
!> \brief ...
!> \param iat ...
!> \param jat ...
!> \param iadr ...
!> \param jadr ...
! *****************************************************************************
  SUBROUTINE limit(iat,jat,iadr,jadr)
    INTEGER                                  :: iat, jat, iadr, jadr

    INTEGER                                  :: i

      iadr=1
      jadr=1
      i=100
 10   IF(iat.gt.100) THEN
         iat=iat-100
         iadr=iadr+1
         GOTO 10
      ENDIF

      i=100
 20   IF(jat.gt.100) THEN
         jat=jat-100
         jadr=jadr+1
         GOTO 20
      ENDIF
   END SUBROUTINE limit

! *****************************************************************************
!> \brief ...
!> \param rout ...
!> \param rcov ...
!> \param r2r4 ...
! *****************************************************************************
   SUBROUTINE setr0ab(rout,rcov,r2r4)
      ! set cut-off radii
    REAL(KIND=dp), DIMENSION(:, :)           :: rout
    REAL(KIND=dp), DIMENSION(:)              :: rcov, r2r4

    INTEGER                                  :: i, j, k
    REAL(KIND=dp), DIMENSION(4465)           :: r0ab(4465)

      r0ab(   1:  70)=(/&
         2.1823,  1.8547,  1.7347,  2.9086,  2.5732,  3.4956,  2.3550&
      ,  2.5095,  2.9802,  3.0982,  2.5141,  2.3917,  2.9977,  2.9484&
      ,  3.2160,  2.4492,  2.2527,  3.1933,  3.0214,  2.9531,  2.9103&
      ,  2.3667,  2.1328,  2.8784,  2.7660,  2.7776,  2.7063,  2.6225&
      ,  2.1768,  2.0625,  2.6395,  2.6648,  2.6482,  2.5697,  2.4846&
      ,  2.4817,  2.0646,  1.9891,  2.5086,  2.6908,  2.6233,  2.4770&
      ,  2.3885,  2.3511,  2.2996,  1.9892,  1.9251,  2.4190,  2.5473&
      ,  2.4994,  2.4091,  2.3176,  2.2571,  2.1946,  2.1374,  2.9898&
      ,  2.6397,  3.6031,  3.1219,  3.7620,  3.2485,  2.9357,  2.7093&
      ,  2.5781,  2.4839,  3.7082,  2.5129,  2.7321,  3.1052,  3.2962&
      /)
      r0ab(  71: 140)=(/&
         3.1331,  3.2000,  2.9586,  3.0822,  2.8582,  2.7120,  3.2570&
      ,  3.4839,  2.8766,  2.7427,  3.2776,  3.2363,  3.5929,  3.2826&
      ,  3.0911,  2.9369,  2.9030,  2.7789,  3.3921,  3.3970,  4.0106&
      ,  2.8884,  2.6605,  3.7513,  3.1613,  3.3605,  3.3325,  3.0991&
      ,  2.9297,  2.8674,  2.7571,  3.8129,  3.3266,  3.7105,  3.7917&
      ,  2.8304,  2.5538,  3.3932,  3.1193,  3.1866,  3.1245,  3.0465&
      ,  2.8727,  2.7664,  2.6926,  3.4608,  3.2984,  3.5142,  3.5418&
      ,  3.5017,  2.6190,  2.4797,  3.1331,  3.0540,  3.0651,  2.9879&
      ,  2.9054,  2.8805,  2.7330,  2.6331,  3.2096,  3.5668,  3.3684&
      ,  3.3686,  3.3180,  3.3107,  2.4757,  2.4019,  2.9789,  3.1468&
      /)
      r0ab( 141: 210)=(/&
         2.9768,  2.8848,  2.7952,  2.7457,  2.6881,  2.5728,  3.0574&
      ,  3.3264,  3.3562,  3.2529,  3.1916,  3.1523,  3.1046,  2.3725&
      ,  2.3289,  2.8760,  2.9804,  2.9093,  2.8040,  2.7071,  2.6386&
      ,  2.5720,  2.5139,  2.9517,  3.1606,  3.2085,  3.1692,  3.0982&
      ,  3.0352,  2.9730,  2.9148,  3.2147,  2.8315,  3.8724,  3.4621&
      ,  3.8823,  3.3760,  3.0746,  2.8817,  2.7552,  2.6605,  3.9740&
      ,  3.6192,  3.6569,  3.9586,  3.6188,  3.3917,  3.2479,  3.1434&
      ,  4.2411,  2.7597,  3.0588,  3.3474,  3.6214,  3.4353,  3.4729&
      ,  3.2487,  3.3200,  3.0914,  2.9403,  3.4972,  3.7993,  3.6773&
      ,  3.8678,  3.5808,  3.8243,  3.5826,  3.4156,  3.8765,  4.1035&
      /)
      r0ab( 211: 280)=(/&
         2.7361,  2.9765,  3.2475,  3.5004,  3.4185,  3.4378,  3.2084&
      ,  3.2787,  3.0604,  2.9187,  3.4037,  3.6759,  3.6586,  3.8327&
      ,  3.5372,  3.7665,  3.5310,  3.3700,  3.7788,  3.9804,  3.8903&
      ,  2.6832,  2.9060,  3.2613,  3.4359,  3.3538,  3.3860,  3.1550&
      ,  3.2300,  3.0133,  2.8736,  3.4024,  3.6142,  3.5979,  3.5295&
      ,  3.4834,  3.7140,  3.4782,  3.3170,  3.7434,  3.9623,  3.8181&
      ,  3.7642,  2.6379,  2.8494,  3.1840,  3.4225,  3.2771,  3.3401&
      ,  3.1072,  3.1885,  2.9714,  2.8319,  3.3315,  3.5979,  3.5256&
      ,  3.4980,  3.4376,  3.6714,  3.4346,  3.2723,  3.6859,  3.8985&
      ,  3.7918,  3.7372,  3.7211,  2.9230,  2.6223,  3.4161,  2.8999&
      /)
      r0ab( 281: 350)=(/&
         3.0557,  3.3308,  3.0555,  2.8508,  2.7385,  2.6640,  3.5263&
      ,  3.0277,  3.2990,  3.7721,  3.5017,  3.2751,  3.1368,  3.0435&
      ,  3.7873,  3.2858,  3.2140,  3.1727,  3.2178,  3.4414,  2.5490&
      ,  2.7623,  3.0991,  3.3252,  3.1836,  3.2428,  3.0259,  3.1225&
      ,  2.9032,  2.7621,  3.2490,  3.5110,  3.4429,  3.3845,  3.3574&
      ,  3.6045,  3.3658,  3.2013,  3.6110,  3.8241,  3.7090,  3.6496&
      ,  3.6333,  3.0896,  3.5462,  2.4926,  2.7136,  3.0693,  3.2699&
      ,  3.1272,  3.1893,  2.9658,  3.0972,  2.8778,  2.7358,  3.2206&
      ,  3.4566,  3.3896,  3.3257,  3.2946,  3.5693,  3.3312,  3.1670&
      ,  3.5805,  3.7711,  3.6536,  3.5927,  3.5775,  3.0411,  3.4885&
      /)
      r0ab( 351: 420)=(/&
         3.4421,  2.4667,  2.6709,  3.0575,  3.2357,  3.0908,  3.1537&
      ,  2.9235,  3.0669,  2.8476,  2.7054,  3.2064,  3.4519,  3.3593&
      ,  3.2921,  3.2577,  3.2161,  3.2982,  3.1339,  3.5606,  3.7582&
      ,  3.6432,  3.5833,  3.5691,  3.0161,  3.4812,  3.4339,  3.4327&
      ,  2.4515,  2.6338,  3.0511,  3.2229,  3.0630,  3.1265,  2.8909&
      ,  3.0253,  2.8184,  2.6764,  3.1968,  3.4114,  3.3492,  3.2691&
      ,  3.2320,  3.1786,  3.2680,  3.1036,  3.5453,  3.7259,  3.6090&
      ,  3.5473,  3.5327,  3.0018,  3.4413,  3.3907,  3.3593,  3.3462&
      ,  2.4413,  2.6006,  3.0540,  3.1987,  3.0490,  3.1058,  2.8643&
      ,  2.9948,  2.7908,  2.6491,  3.1950,  3.3922,  3.3316,  3.2585&
      /)
      r0ab( 421: 490)=(/&
         3.2136,  3.1516,  3.2364,  3.0752,  3.5368,  3.7117,  3.5941&
      ,  3.5313,  3.5164,  2.9962,  3.4225,  3.3699,  3.3370,  3.3234&
      ,  3.3008,  2.4318,  2.5729,  3.0416,  3.1639,  3.0196,  3.0843&
      ,  2.8413,  2.7436,  2.7608,  2.6271,  3.1811,  3.3591,  3.3045&
      ,  3.2349,  3.1942,  3.1291,  3.2111,  3.0534,  3.5189,  3.6809&
      ,  3.5635,  3.5001,  3.4854,  2.9857,  3.3897,  3.3363,  3.3027&
      ,  3.2890,  3.2655,  3.2309,  2.8502,  2.6934,  3.2467,  3.1921&
      ,  3.5663,  3.2541,  3.0571,  2.9048,  2.8657,  2.7438,  3.3547&
      ,  3.3510,  3.9837,  3.6871,  3.4862,  3.3389,  3.2413,  3.1708&
      ,  3.6096,  3.6280,  3.6860,  3.5568,  3.4836,  3.2868,  3.3994&
      /)
      r0ab( 491: 560)=(/&
         3.3476,  3.3170,  3.2950,  3.2874,  3.2606,  3.9579,  2.9226&
      ,  2.6838,  3.7867,  3.1732,  3.3872,  3.3643,  3.1267,  2.9541&
      ,  2.8505,  2.7781,  3.8475,  3.3336,  3.7359,  3.8266,  3.5733&
      ,  3.3959,  3.2775,  3.1915,  3.9878,  3.8816,  3.5810,  3.5364&
      ,  3.5060,  3.8097,  3.3925,  3.3348,  3.3019,  3.2796,  3.2662&
      ,  3.2464,  3.7136,  3.8619,  2.9140,  2.6271,  3.4771,  3.1774&
      ,  3.2560,  3.1970,  3.1207,  2.9406,  2.8322,  2.7571,  3.5455&
      ,  3.3514,  3.5837,  3.6177,  3.5816,  3.3902,  3.2604,  3.1652&
      ,  3.7037,  3.6283,  3.5858,  3.5330,  3.4884,  3.5789,  3.4094&
      ,  3.3473,  3.3118,  3.2876,  3.2707,  3.2521,  3.5570,  3.6496&
      /)
      r0ab( 561: 630)=(/&
         3.6625,  2.7300,  2.5870,  3.2471,  3.1487,  3.1667,  3.0914&
      ,  3.0107,  2.9812,  2.8300,  2.7284,  3.3259,  3.3182,  3.4707&
      ,  3.4748,  3.4279,  3.4182,  3.2547,  3.1353,  3.5116,  3.9432&
      ,  3.8828,  3.8303,  3.7880,  3.3760,  3.7218,  3.3408,  3.3059&
      ,  3.2698,  3.2446,  3.2229,  3.4422,  3.5023,  3.5009,  3.5268&
      ,  2.6026,  2.5355,  3.1129,  3.2863,  3.1029,  3.0108,  2.9227&
      ,  2.8694,  2.8109,  2.6929,  3.1958,  3.4670,  3.4018,  3.3805&
      ,  3.3218,  3.2815,  3.2346,  3.0994,  3.3937,  3.7266,  3.6697&
      ,  3.6164,  3.5730,  3.2522,  3.5051,  3.4686,  3.4355,  3.4084&
      ,  3.3748,  3.3496,  3.3692,  3.4052,  3.3910,  3.3849,  3.3662&
      /)
      r0ab( 631: 700)=(/&
         2.5087,  2.4814,  3.0239,  3.1312,  3.0535,  2.9457,  2.8496&
      ,  2.7780,  2.7828,  2.6532,  3.1063,  3.3143,  3.3549,  3.3120&
      ,  3.2421,  3.1787,  3.1176,  3.0613,  3.3082,  3.5755,  3.5222&
      ,  3.4678,  3.4231,  3.1684,  3.3528,  3.3162,  3.2827,  3.2527&
      ,  3.2308,  3.2029,  3.3173,  3.3343,  3.3092,  3.2795,  3.2452&
      ,  3.2096,  3.2893,  2.8991,  4.0388,  3.6100,  3.9388,  3.4475&
      ,  3.1590,  2.9812,  2.8586,  2.7683,  4.1428,  3.7911,  3.8225&
      ,  4.0372,  3.7059,  3.4935,  3.3529,  3.2492,  4.4352,  4.0826&
      ,  3.9733,  3.9254,  3.8646,  3.9315,  3.7837,  3.7465,  3.7211&
      ,  3.7012,  3.6893,  3.6676,  3.7736,  4.0660,  3.7926,  3.6158&
      /)
      r0ab( 701: 770)=(/&
         3.5017,  3.4166,  4.6176,  2.8786,  3.1658,  3.5823,  3.7689&
      ,  3.5762,  3.5789,  3.3552,  3.4004,  3.1722,  3.0212,  3.7241&
      ,  3.9604,  3.8500,  3.9844,  3.7035,  3.9161,  3.6751,  3.5075&
      ,  4.1151,  4.2877,  4.1579,  4.1247,  4.0617,  3.4874,  3.9848&
      ,  3.9280,  3.9079,  3.8751,  3.8604,  3.8277,  3.8002,  3.9981&
      ,  3.7544,  4.0371,  3.8225,  3.6718,  4.3092,  4.4764,  2.8997&
      ,  3.0953,  3.4524,  3.6107,  3.6062,  3.5783,  3.3463,  3.3855&
      ,  3.1746,  3.0381,  3.6019,  3.7938,  3.8697,  3.9781,  3.6877&
      ,  3.8736,  3.6451,  3.4890,  3.9858,  4.1179,  4.0430,  3.9563&
      ,  3.9182,  3.4002,  3.8310,  3.7716,  3.7543,  3.7203,  3.7053&
      /)
      r0ab( 771: 840)=(/&
         3.6742,  3.8318,  3.7631,  3.7392,  3.9892,  3.7832,  3.6406&
      ,  4.1701,  4.3016,  4.2196,  2.8535,  3.0167,  3.3978,  3.5363&
      ,  3.5393,  3.5301,  3.2960,  3.3352,  3.1287,  2.9967,  3.6659&
      ,  3.7239,  3.8070,  3.7165,  3.6368,  3.8162,  3.5885,  3.4336&
      ,  3.9829,  4.0529,  3.9584,  3.9025,  3.8607,  3.3673,  3.7658&
      ,  3.7035,  3.6866,  3.6504,  3.6339,  3.6024,  3.7708,  3.7283&
      ,  3.6896,  3.9315,  3.7250,  3.5819,  4.1457,  4.2280,  4.1130&
      ,  4.0597,  3.0905,  2.7998,  3.6448,  3.0739,  3.2996,  3.5262&
      ,  3.2559,  3.0518,  2.9394,  2.8658,  3.7514,  3.2295,  3.5643&
      ,  3.7808,  3.6931,  3.4723,  3.3357,  3.2429,  4.0280,  3.5589&
      /)
      r0ab( 841: 910)=(/&
         3.4636,  3.4994,  3.4309,  3.6177,  3.2946,  3.2376,  3.2050&
      ,  3.1847,  3.1715,  3.1599,  3.5555,  3.8111,  3.7693,  3.5718&
      ,  3.4498,  3.3662,  4.1608,  3.7417,  3.6536,  3.6154,  3.8596&
      ,  3.0301,  2.7312,  3.5821,  3.0473,  3.2137,  3.4679,  3.1975&
      ,  2.9969,  2.8847,  2.8110,  3.6931,  3.2076,  3.4943,  3.5956&
      ,  3.6379,  3.4190,  3.2808,  3.1860,  3.9850,  3.5105,  3.4330&
      ,  3.3797,  3.4155,  3.6033,  3.2737,  3.2145,  3.1807,  3.1596&
      ,  3.1461,  3.1337,  3.4812,  3.6251,  3.7152,  3.5201,  3.3966&
      ,  3.3107,  4.1128,  3.6899,  3.6082,  3.5604,  3.7834,  3.7543&
      ,  2.9189,  2.6777,  3.4925,  2.9648,  3.1216,  3.2940,  3.0975&
      /)
      r0ab( 911: 980)=(/&
         2.9757,  2.8493,  2.7638,  3.6085,  3.1214,  3.4006,  3.4793&
      ,  3.5147,  3.3806,  3.2356,  3.1335,  3.9144,  3.4183,  3.3369&
      ,  3.2803,  3.2679,  3.4871,  3.1714,  3.1521,  3.1101,  3.0843&
      ,  3.0670,  3.0539,  3.3890,  3.5086,  3.5895,  3.4783,  3.3484&
      ,  3.2559,  4.0422,  3.5967,  3.5113,  3.4576,  3.6594,  3.6313&
      ,  3.5690,  2.8578,  2.6334,  3.4673,  2.9245,  3.0732,  3.2435&
      ,  3.0338,  2.9462,  2.8143,  2.7240,  3.5832,  3.0789,  3.3617&
      ,  3.4246,  3.4505,  3.3443,  3.1964,  3.0913,  3.8921,  3.3713&
      ,  3.2873,  3.2281,  3.2165,  3.4386,  3.1164,  3.1220,  3.0761&
      ,  3.0480,  3.0295,  3.0155,  3.3495,  3.4543,  3.5260,  3.4413&
      /)
      r0ab( 981:1050)=(/&
         3.3085,  3.2134,  4.0170,  3.5464,  3.4587,  3.4006,  3.6027&
      ,  3.5730,  3.4945,  3.4623,  2.8240,  2.5960,  3.4635,  2.9032&
      ,  3.0431,  3.2115,  2.9892,  2.9148,  2.7801,  2.6873,  3.5776&
      ,  3.0568,  3.3433,  3.3949,  3.4132,  3.3116,  3.1616,  3.0548&
      ,  3.8859,  3.3719,  3.2917,  3.2345,  3.2274,  3.4171,  3.1293&
      ,  3.0567,  3.0565,  3.0274,  3.0087,  2.9939,  3.3293,  3.4249&
      ,  3.4902,  3.4091,  3.2744,  3.1776,  4.0078,  3.5374,  3.4537&
      ,  3.3956,  3.5747,  3.5430,  3.4522,  3.4160,  3.3975,  2.8004&
      ,  2.5621,  3.4617,  2.9154,  3.0203,  3.1875,  2.9548,  2.8038&
      ,  2.7472,  2.6530,  3.5736,  3.0584,  3.3304,  3.3748,  3.3871&
      /)
      r0ab(1051:1120)=(/&
         3.2028,  3.1296,  3.0214,  3.8796,  3.3337,  3.2492,  3.1883&
      ,  3.1802,  3.4050,  3.0756,  3.0478,  3.0322,  3.0323,  3.0163&
      ,  3.0019,  3.3145,  3.4050,  3.4656,  3.3021,  3.2433,  3.1453&
      ,  3.9991,  3.5017,  3.4141,  3.3520,  3.5583,  3.5251,  3.4243&
      ,  3.3851,  3.3662,  3.3525,  2.7846,  2.5324,  3.4652,  2.8759&
      ,  3.0051,  3.1692,  2.9273,  2.7615,  2.7164,  2.6212,  3.5744&
      ,  3.0275,  3.3249,  3.3627,  3.3686,  3.1669,  3.0584,  2.9915&
      ,  3.8773,  3.3099,  3.2231,  3.1600,  3.1520,  3.4023,  3.0426&
      ,  3.0099,  2.9920,  2.9809,  2.9800,  2.9646,  3.3068,  3.3930&
      ,  3.4486,  3.2682,  3.1729,  3.1168,  3.9952,  3.4796,  3.3901&
      /)
      r0ab(1121:1190)=(/&
         3.3255,  3.5530,  3.5183,  3.4097,  3.3683,  3.3492,  3.3360&
      ,  3.3308,  2.5424,  2.6601,  3.2555,  3.2807,  3.1384,  3.1737&
      ,  2.9397,  2.8429,  2.8492,  2.7225,  3.3875,  3.4910,  3.4520&
      ,  3.3608,  3.3036,  3.2345,  3.2999,  3.1487,  3.7409,  3.8392&
      ,  3.7148,  3.6439,  3.6182,  3.1753,  3.5210,  3.4639,  3.4265&
      ,  3.4075,  3.3828,  3.3474,  3.4071,  3.3754,  3.3646,  3.3308&
      ,  3.4393,  3.2993,  3.8768,  3.9891,  3.8310,  3.7483,  3.3417&
      ,  3.3019,  3.2250,  3.1832,  3.1578,  3.1564,  3.1224,  3.4620&
      ,  2.9743,  2.8058,  3.4830,  3.3474,  3.6863,  3.3617,  3.1608&
      ,  3.0069,  2.9640,  2.8427,  3.5885,  3.5219,  4.1314,  3.8120&
      /)
      r0ab(1191:1260)=(/&
         3.6015,  3.4502,  3.3498,  3.2777,  3.8635,  3.8232,  3.8486&
      ,  3.7215,  3.6487,  3.4724,  3.5627,  3.5087,  3.4757,  3.4517&
      ,  3.4423,  3.4139,  4.1028,  3.8388,  3.6745,  3.5562,  3.4806&
      ,  3.4272,  4.0182,  3.9991,  4.0007,  3.9282,  3.7238,  3.6498&
      ,  3.5605,  3.5211,  3.5009,  3.4859,  3.4785,  3.5621,  4.2623&
      ,  3.0775,  2.8275,  4.0181,  3.3385,  3.5379,  3.5036,  3.2589&
      ,  3.0804,  3.0094,  2.9003,  4.0869,  3.5088,  3.9105,  3.9833&
      ,  3.7176,  3.5323,  3.4102,  3.3227,  4.2702,  4.0888,  3.7560&
      ,  3.7687,  3.6681,  3.6405,  3.5569,  3.4990,  3.4659,  3.4433&
      ,  3.4330,  3.4092,  3.8867,  4.0190,  3.7961,  3.6412,  3.5405&
      /)
      r0ab(1261:1330)=(/&
         3.4681,  4.3538,  4.2136,  3.9381,  3.8912,  3.9681,  3.7909&
      ,  3.6774,  3.6262,  3.5999,  3.5823,  3.5727,  3.5419,  4.0245&
      ,  4.1874,  3.0893,  2.7917,  3.7262,  3.3518,  3.4241,  3.5433&
      ,  3.2773,  3.0890,  2.9775,  2.9010,  3.8048,  3.5362,  3.7746&
      ,  3.7911,  3.7511,  3.5495,  3.4149,  3.3177,  4.0129,  3.8370&
      ,  3.7739,  3.7125,  3.7152,  3.7701,  3.5813,  3.5187,  3.4835&
      ,  3.4595,  3.4439,  3.4242,  3.7476,  3.8239,  3.8346,  3.6627&
      ,  3.5479,  3.4639,  4.1026,  3.9733,  3.9292,  3.8667,  3.9513&
      ,  3.8959,  3.7698,  3.7089,  3.6765,  3.6548,  3.6409,  3.5398&
      ,  3.8759,  3.9804,  4.0150,  2.9091,  2.7638,  3.5066,  3.3377&
      /)
      r0ab(1331:1400)=(/&
         3.3481,  3.2633,  3.1810,  3.1428,  2.9872,  2.8837,  3.5929&
      ,  3.5183,  3.6729,  3.6596,  3.6082,  3.5927,  3.4224,  3.2997&
      ,  3.8190,  4.1865,  4.1114,  4.0540,  3.6325,  3.5697,  3.5561&
      ,  3.5259,  3.4901,  3.4552,  3.4315,  3.4091,  3.6438,  3.6879&
      ,  3.6832,  3.7043,  3.5557,  3.4466,  3.9203,  4.2919,  4.2196&
      ,  4.1542,  3.7573,  3.7039,  3.6546,  3.6151,  3.5293,  3.4849&
      ,  3.4552,  3.5192,  3.7673,  3.8359,  3.8525,  3.8901,  2.7806&
      ,  2.7209,  3.3812,  3.4958,  3.2913,  3.1888,  3.0990,  3.0394&
      ,  2.9789,  2.8582,  3.4716,  3.6883,  3.6105,  3.5704,  3.5059&
      ,  3.4619,  3.4138,  3.2742,  3.7080,  3.9773,  3.9010,  3.8409&
      /)
      r0ab(1401:1470)=(/&
         3.7944,  3.4465,  3.7235,  3.6808,  3.6453,  3.6168,  3.5844&
      ,  3.5576,  3.5772,  3.5959,  3.5768,  3.5678,  3.5486,  3.4228&
      ,  3.8107,  4.0866,  4.0169,  3.9476,  3.6358,  3.5800,  3.5260&
      ,  3.4838,  3.4501,  3.4204,  3.3553,  3.6487,  3.6973,  3.7398&
      ,  3.7405,  3.7459,  3.7380,  2.6848,  2.6740,  3.2925,  3.3386&
      ,  3.2473,  3.1284,  3.0301,  2.9531,  2.9602,  2.8272,  3.3830&
      ,  3.5358,  3.5672,  3.5049,  3.4284,  3.3621,  3.3001,  3.2451&
      ,  3.6209,  3.8299,  3.7543,  3.6920,  3.6436,  3.3598,  3.5701&
      ,  3.5266,  3.4904,  3.4590,  3.4364,  3.4077,  3.5287,  3.5280&
      ,  3.4969,  3.4650,  3.4304,  3.3963,  3.7229,  3.9402,  3.8753&
      /)
      r0ab(1471:1540)=(/&
         3.8035,  3.5499,  3.4913,  3.4319,  3.3873,  3.3520,  3.3209&
      ,  3.2948,  3.5052,  3.6465,  3.6696,  3.6577,  3.6388,  3.6142&
      ,  3.5889,  3.3968,  3.0122,  4.2241,  3.7887,  4.0049,  3.5384&
      ,  3.2698,  3.1083,  2.9917,  2.9057,  4.3340,  3.9900,  4.6588&
      ,  4.1278,  3.8125,  3.6189,  3.4851,  3.3859,  4.6531,  4.3134&
      ,  4.2258,  4.1309,  4.0692,  4.0944,  3.9850,  3.9416,  3.9112&
      ,  3.8873,  3.8736,  3.8473,  4.6027,  4.1538,  3.8994,  3.7419&
      ,  3.6356,  3.5548,  4.8353,  4.5413,  4.3891,  4.3416,  4.3243&
      ,  4.2753,  4.2053,  4.1790,  4.1685,  4.1585,  4.1536,  4.0579&
      ,  4.1980,  4.4564,  4.2192,  4.0528,  3.9489,  3.8642,  5.0567&
      /)
      r0ab(1541:1610)=(/&
         3.0630,  3.3271,  4.0432,  4.0046,  4.1555,  3.7426,  3.5130&
      ,  3.5174,  3.2884,  3.1378,  4.1894,  4.2321,  4.1725,  4.1833&
      ,  3.8929,  4.0544,  3.8118,  3.6414,  4.6373,  4.6268,  4.4750&
      ,  4.4134,  4.3458,  3.8582,  4.2583,  4.1898,  4.1562,  4.1191&
      ,  4.1069,  4.0639,  4.1257,  4.1974,  3.9532,  4.1794,  3.9660&
      ,  3.8130,  4.8160,  4.8272,  4.6294,  4.5840,  4.0770,  4.0088&
      ,  3.9103,  3.8536,  3.8324,  3.7995,  3.7826,  4.2294,  4.3380&
      ,  4.4352,  4.1933,  4.4580,  4.2554,  4.1072,  5.0454,  5.1814&
      ,  3.0632,  3.2662,  3.6432,  3.8088,  3.7910,  3.7381,  3.5093&
      ,  3.5155,  3.3047,  3.1681,  3.7871,  3.9924,  4.0637,  4.1382&
      /)
      r0ab(1611:1680)=(/&
         3.8591,  4.0164,  3.7878,  3.6316,  4.1741,  4.3166,  4.2395&
      ,  4.1831,  4.1107,  3.5857,  4.0270,  3.9676,  3.9463,  3.9150&
      ,  3.9021,  3.8708,  4.0240,  4.1551,  3.9108,  4.1337,  3.9289&
      ,  3.7873,  4.3666,  4.5080,  4.4232,  4.3155,  3.8461,  3.8007&
      ,  3.6991,  3.6447,  3.6308,  3.5959,  3.5749,  4.0359,  4.3124&
      ,  4.3539,  4.1122,  4.3772,  4.1785,  4.0386,  4.7004,  4.8604&
      ,  4.6261,  2.9455,  3.2470,  3.6108,  3.8522,  3.6625,  3.6598&
      ,  3.4411,  3.4660,  3.2415,  3.0944,  3.7514,  4.0397,  3.9231&
      ,  4.0561,  3.7860,  3.9845,  3.7454,  3.5802,  4.1366,  4.3581&
      ,  4.2351,  4.2011,  4.1402,  3.5381,  4.0653,  4.0093,  3.9883&
      /)
      r0ab(1681:1750)=(/&
         3.9570,  3.9429,  3.9112,  3.8728,  4.0682,  3.8351,  4.1054&
      ,  3.8928,  3.7445,  4.3415,  4.5497,  4.3833,  4.3122,  3.8051&
      ,  3.7583,  3.6622,  3.6108,  3.5971,  3.5628,  3.5408,  4.0780&
      ,  4.0727,  4.2836,  4.0553,  4.3647,  4.1622,  4.0178,  4.5802&
      ,  4.9125,  4.5861,  4.6201,  2.9244,  3.2241,  3.5848,  3.8293&
      ,  3.6395,  3.6400,  3.4204,  3.4499,  3.2253,  3.0779,  3.7257&
      ,  4.0170,  3.9003,  4.0372,  3.7653,  3.9672,  3.7283,  3.5630&
      ,  4.1092,  4.3347,  4.2117,  4.1793,  4.1179,  3.5139,  4.0426&
      ,  3.9867,  3.9661,  3.9345,  3.9200,  3.8883,  3.8498,  4.0496&
      ,  3.8145,  4.0881,  3.8756,  3.7271,  4.3128,  4.5242,  4.3578&
      /)
      r0ab(1751:1820)=(/&
         4.2870,  3.7796,  3.7318,  3.6364,  3.5854,  3.5726,  3.5378&
      ,  3.5155,  4.0527,  4.0478,  4.2630,  4.0322,  4.3449,  4.1421&
      ,  3.9975,  4.5499,  4.8825,  4.5601,  4.5950,  4.5702,  2.9046&
      ,  3.2044,  3.5621,  3.8078,  3.6185,  3.6220,  3.4019,  3.4359&
      ,  3.2110,  3.0635,  3.7037,  3.9958,  3.8792,  4.0194,  3.7460&
      ,  3.9517,  3.7128,  3.5474,  4.0872,  4.3138,  4.1906,  4.1593&
      ,  4.0973,  3.4919,  4.0216,  3.9657,  3.9454,  3.9134,  3.8986&
      ,  3.8669,  3.8289,  4.0323,  3.7954,  4.0725,  3.8598,  3.7113&
      ,  4.2896,  4.5021,  4.3325,  4.2645,  3.7571,  3.7083,  3.6136&
      ,  3.5628,  3.5507,  3.5155,  3.4929,  4.0297,  4.0234,  4.2442&
      /)
      r0ab(1821:1890)=(/&
         4.0112,  4.3274,  4.1240,  3.9793,  4.5257,  4.8568,  4.5353&
      ,  4.5733,  4.5485,  4.5271,  2.8878,  3.1890,  3.5412,  3.7908&
      ,  3.5974,  3.6078,  3.3871,  3.4243,  3.1992,  3.0513,  3.6831&
      ,  3.9784,  3.8579,  4.0049,  3.7304,  3.9392,  3.7002,  3.5347&
      ,  4.0657,  4.2955,  4.1705,  4.1424,  4.0800,  3.4717,  4.0043&
      ,  3.9485,  3.9286,  3.8965,  3.8815,  3.8500,  3.8073,  4.0180&
      ,  3.7796,  4.0598,  3.8470,  3.6983,  4.2678,  4.4830,  4.3132&
      ,  4.2444,  3.7370,  3.6876,  3.5935,  3.5428,  3.5314,  3.4958&
      ,  3.4730,  4.0117,  4.0043,  4.2287,  3.9939,  4.3134,  4.1096&
      ,  3.9646,  4.5032,  4.8356,  4.5156,  4.5544,  4.5297,  4.5083&
      /)
      r0ab(1891:1960)=(/&
         4.4896,  2.8709,  3.1737,  3.5199,  3.7734,  3.5802,  3.5934&
      ,  3.3724,  3.4128,  3.1877,  3.0396,  3.6624,  3.9608,  3.8397&
      ,  3.9893,  3.7145,  3.9266,  3.6877,  3.5222,  4.0448,  4.2771&
      ,  4.1523,  4.1247,  4.0626,  3.4530,  3.9866,  3.9310,  3.9115&
      ,  3.8792,  3.8641,  3.8326,  3.7892,  4.0025,  3.7636,  4.0471&
      ,  3.8343,  3.6854,  4.2464,  4.4635,  4.2939,  4.2252,  3.7169&
      ,  3.6675,  3.5739,  3.5235,  3.5126,  3.4768,  3.4537,  3.9932&
      ,  3.9854,  4.2123,  3.9765,  4.2992,  4.0951,  3.9500,  4.4811&
      ,  4.8135,  4.4959,  4.5351,  4.5105,  4.4891,  4.4705,  4.4515&
      ,  2.8568,  3.1608,  3.5050,  3.7598,  3.5665,  3.5803,  3.3601&
      /)
      r0ab(1961:2030)=(/&
         3.4031,  3.1779,  3.0296,  3.6479,  3.9471,  3.8262,  3.9773&
      ,  3.7015,  3.9162,  3.6771,  3.5115,  4.0306,  4.2634,  4.1385&
      ,  4.1116,  4.0489,  3.4366,  3.9732,  3.9176,  3.8983,  3.8659&
      ,  3.8507,  3.8191,  3.7757,  3.9907,  3.7506,  4.0365,  3.8235&
      ,  3.6745,  4.2314,  4.4490,  4.2792,  4.2105,  3.7003,  3.6510&
      ,  3.5578,  3.5075,  3.4971,  3.4609,  3.4377,  3.9788,  3.9712&
      ,  4.1997,  3.9624,  4.2877,  4.0831,  3.9378,  4.4655,  4.7974&
      ,  4.4813,  4.5209,  4.4964,  4.4750,  4.4565,  4.4375,  4.4234&
      ,  2.6798,  3.0151,  3.2586,  3.5292,  3.5391,  3.4902,  3.2887&
      ,  3.3322,  3.1228,  2.9888,  3.4012,  3.7145,  3.7830,  3.6665&
      /)
      r0ab(2031:2100)=(/&
         3.5898,  3.8077,  3.5810,  3.4265,  3.7726,  4.0307,  3.9763&
      ,  3.8890,  3.8489,  3.2706,  3.7595,  3.6984,  3.6772,  3.6428&
      ,  3.6243,  3.5951,  3.7497,  3.6775,  3.6364,  3.9203,  3.7157&
      ,  3.5746,  3.9494,  4.2076,  4.1563,  4.0508,  3.5329,  3.4780&
      ,  3.3731,  3.3126,  3.2846,  3.2426,  3.2135,  3.7491,  3.9006&
      ,  3.8332,  3.8029,  4.1436,  3.9407,  3.7998,  4.1663,  4.5309&
      ,  4.3481,  4.2911,  4.2671,  4.2415,  4.2230,  4.2047,  4.1908&
      ,  4.1243,  2.5189,  2.9703,  3.3063,  3.6235,  3.4517,  3.3989&
      ,  3.2107,  3.2434,  3.0094,  2.8580,  3.4253,  3.8157,  3.7258&
      ,  3.6132,  3.5297,  3.7566,  3.5095,  3.3368,  3.7890,  4.1298&
      /)
      r0ab(2101:2170)=(/&
         4.0190,  3.9573,  3.9237,  3.2677,  3.8480,  3.8157,  3.7656&
      ,  3.7317,  3.7126,  3.6814,  3.6793,  3.6218,  3.5788,  3.8763&
      ,  3.6572,  3.5022,  3.9737,  4.3255,  4.1828,  4.1158,  3.5078&
      ,  3.4595,  3.3600,  3.3088,  3.2575,  3.2164,  3.1856,  3.8522&
      ,  3.8665,  3.8075,  3.7772,  4.1391,  3.9296,  3.7772,  4.2134&
      ,  4.7308,  4.3787,  4.3894,  4.3649,  4.3441,  4.3257,  4.3073&
      ,  4.2941,  4.1252,  4.2427,  3.0481,  2.9584,  3.6919,  3.5990&
      ,  3.8881,  3.4209,  3.1606,  3.1938,  2.9975,  2.8646,  3.8138&
      ,  3.7935,  3.7081,  3.9155,  3.5910,  3.4808,  3.4886,  3.3397&
      ,  4.1336,  4.1122,  3.9888,  3.9543,  3.8917,  3.5894,  3.8131&
      /)
      r0ab(2171:2240)=(/&
         3.7635,  3.7419,  3.7071,  3.6880,  3.6574,  3.6546,  3.9375&
      ,  3.6579,  3.5870,  3.6361,  3.5039,  4.3149,  4.2978,  4.1321&
      ,  4.1298,  3.8164,  3.7680,  3.7154,  3.6858,  3.6709,  3.6666&
      ,  3.6517,  3.8174,  3.8608,  4.1805,  3.9102,  3.8394,  3.8968&
      ,  3.7673,  4.5274,  4.6682,  4.3344,  4.3639,  4.3384,  4.3162&
      ,  4.2972,  4.2779,  4.2636,  4.0253,  4.1168,  4.1541,  2.8136&
      ,  3.0951,  3.4635,  3.6875,  3.4987,  3.5183,  3.2937,  3.3580&
      ,  3.1325,  2.9832,  3.6078,  3.8757,  3.7616,  3.9222,  3.6370&
      ,  3.8647,  3.6256,  3.4595,  3.9874,  4.1938,  4.0679,  4.0430&
      ,  3.9781,  3.3886,  3.9008,  3.8463,  3.8288,  3.7950,  3.7790&
      /)
      r0ab(2241:2310)=(/&
         3.7472,  3.7117,  3.9371,  3.6873,  3.9846,  3.7709,  3.6210&
      ,  4.1812,  4.3750,  4.2044,  4.1340,  3.6459,  3.5929,  3.5036&
      ,  3.4577,  3.4528,  3.4146,  3.3904,  3.9014,  3.9031,  4.1443&
      ,  3.8961,  4.2295,  4.0227,  3.8763,  4.4086,  4.7097,  4.4064&
      ,  4.4488,  4.4243,  4.4029,  4.3842,  4.3655,  4.3514,  4.1162&
      ,  4.2205,  4.1953,  4.2794,  2.8032,  3.0805,  3.4519,  3.6700&
      ,  3.4827,  3.5050,  3.2799,  3.3482,  3.1233,  2.9747,  3.5971&
      ,  3.8586,  3.7461,  3.9100,  3.6228,  3.8535,  3.6147,  3.4490&
      ,  3.9764,  4.1773,  4.0511,  4.0270,  3.9614,  3.3754,  3.8836&
      ,  3.8291,  3.8121,  3.7780,  3.7619,  3.7300,  3.6965,  3.9253&
      /)
      r0ab(2311:2380)=(/&
         3.6734,  3.9733,  3.7597,  3.6099,  4.1683,  4.3572,  4.1862&
      ,  4.1153,  3.6312,  3.5772,  3.4881,  3.4429,  3.4395,  3.4009&
      ,  3.3766,  3.8827,  3.8868,  4.1316,  3.8807,  4.2164,  4.0092&
      ,  3.8627,  4.3936,  4.6871,  4.3882,  4.4316,  4.4073,  4.3858&
      ,  4.3672,  4.3485,  4.3344,  4.0984,  4.2036,  4.1791,  4.2622&
      ,  4.2450,  2.7967,  3.0689,  3.4445,  3.6581,  3.4717,  3.4951&
      ,  3.2694,  3.3397,  3.1147,  2.9661,  3.5898,  3.8468,  3.7358&
      ,  3.9014,  3.6129,  3.8443,  3.6054,  3.4396,  3.9683,  4.1656&
      ,  4.0394,  4.0158,  3.9498,  3.3677,  3.8718,  3.8164,  3.8005&
      ,  3.7662,  3.7500,  3.7181,  3.6863,  3.9170,  3.6637,  3.9641&
      /)
      r0ab(2381:2450)=(/&
         3.7503,  3.6004,  4.1590,  4.3448,  4.1739,  4.1029,  3.6224&
      ,  3.5677,  3.4785,  3.4314,  3.4313,  3.3923,  3.3680,  3.8698&
      ,  3.8758,  4.1229,  3.8704,  4.2063,  3.9987,  3.8519,  4.3832&
      ,  4.6728,  4.3759,  4.4195,  4.3952,  4.3737,  4.3551,  4.3364&
      ,  4.3223,  4.0861,  4.1911,  4.1676,  4.2501,  4.2329,  4.2208&
      ,  2.7897,  3.0636,  3.4344,  3.6480,  3.4626,  3.4892,  3.2626&
      ,  3.3344,  3.1088,  2.9597,  3.5804,  3.8359,  3.7251,  3.8940&
      ,  3.6047,  3.8375,  3.5990,  3.4329,  3.9597,  4.1542,  4.0278&
      ,  4.0048,  3.9390,  3.3571,  3.8608,  3.8056,  3.7899,  3.7560&
      ,  3.7400,  3.7081,  3.6758,  3.9095,  3.6552,  3.9572,  3.7436&
      /)
      r0ab(2451:2520)=(/&
         3.5933,  4.1508,  4.3337,  4.1624,  4.0916,  3.6126,  3.5582&
      ,  3.4684,  3.4212,  3.4207,  3.3829,  3.3586,  3.8604,  3.8658&
      ,  4.1156,  3.8620,  4.1994,  3.9917,  3.8446,  4.3750,  4.6617&
      ,  4.3644,  4.4083,  4.3840,  4.3625,  4.3439,  4.3253,  4.3112&
      ,  4.0745,  4.1807,  4.1578,  4.2390,  4.2218,  4.2097,  4.1986&
      ,  2.8395,  3.0081,  3.3171,  3.4878,  3.5360,  3.5145,  3.2809&
      ,  3.3307,  3.1260,  2.9940,  3.4741,  3.6675,  3.7832,  3.6787&
      ,  3.6156,  3.8041,  3.5813,  3.4301,  3.8480,  3.9849,  3.9314&
      ,  3.8405,  3.8029,  3.2962,  3.7104,  3.6515,  3.6378,  3.6020&
      ,  3.5849,  3.5550,  3.7494,  3.6893,  3.6666,  3.9170,  3.7150&
      /)
      r0ab(2521:2590)=(/&
         3.5760,  4.0268,  4.1596,  4.1107,  3.9995,  3.5574,  3.5103&
      ,  3.4163,  3.3655,  3.3677,  3.3243,  3.2975,  3.7071,  3.9047&
      ,  3.8514,  3.8422,  3.8022,  3.9323,  3.7932,  4.2343,  4.4583&
      ,  4.3115,  4.2457,  4.2213,  4.1945,  4.1756,  4.1569,  4.1424&
      ,  4.0620,  4.0494,  3.9953,  4.0694,  4.0516,  4.0396,  4.0280&
      ,  4.0130,  2.9007,  2.9674,  3.8174,  3.5856,  3.6486,  3.5339&
      ,  3.2832,  3.3154,  3.1144,  2.9866,  3.9618,  3.8430,  3.9980&
      ,  3.8134,  3.6652,  3.7985,  3.5756,  3.4207,  4.4061,  4.2817&
      ,  4.1477,  4.0616,  3.9979,  3.6492,  3.8833,  3.8027,  3.7660&
      ,  3.7183,  3.6954,  3.6525,  3.9669,  3.8371,  3.7325,  3.9160&
      /)
      r0ab(2591:2660)=(/&
         3.7156,  3.5714,  4.6036,  4.4620,  4.3092,  4.2122,  3.8478&
      ,  3.7572,  3.6597,  3.5969,  3.5575,  3.5386,  3.5153,  3.7818&
      ,  4.1335,  4.0153,  3.9177,  3.8603,  3.9365,  3.7906,  4.7936&
      ,  4.7410,  4.5461,  4.5662,  4.5340,  4.5059,  4.4832,  4.4604&
      ,  4.4429,  4.2346,  4.4204,  4.3119,  4.3450,  4.3193,  4.3035&
      ,  4.2933,  4.1582,  4.2450,  2.8559,  2.9050,  3.8325,  3.5442&
      ,  3.5077,  3.4905,  3.2396,  3.2720,  3.0726,  2.9467,  3.9644&
      ,  3.8050,  3.8981,  3.7762,  3.6216,  3.7531,  3.5297,  3.3742&
      ,  4.3814,  4.2818,  4.1026,  4.0294,  3.9640,  3.6208,  3.8464&
      ,  3.7648,  3.7281,  3.6790,  3.6542,  3.6117,  3.8650,  3.8010&
      /)
      r0ab(2661:2730)=(/&
         3.6894,  3.8713,  3.6699,  3.5244,  4.5151,  4.4517,  4.2538&
      ,  4.1483,  3.8641,  3.7244,  3.6243,  3.5589,  3.5172,  3.4973&
      ,  3.4715,  3.7340,  4.0316,  3.9958,  3.8687,  3.8115,  3.8862&
      ,  3.7379,  4.7091,  4.7156,  4.5199,  4.5542,  4.5230,  4.4959&
      ,  4.4750,  4.4529,  4.4361,  4.1774,  4.3774,  4.2963,  4.3406&
      ,  4.3159,  4.3006,  4.2910,  4.1008,  4.1568,  4.0980,  2.8110&
      ,  2.8520,  3.7480,  3.5105,  3.4346,  3.3461,  3.1971,  3.2326&
      ,  3.0329,  2.9070,  3.8823,  3.7928,  3.8264,  3.7006,  3.5797&
      ,  3.7141,  3.4894,  3.3326,  4.3048,  4.2217,  4.0786,  3.9900&
      ,  3.9357,  3.6331,  3.8333,  3.7317,  3.6957,  3.6460,  3.6197&
      /)
      r0ab(2731:2800)=(/&
         3.5779,  3.7909,  3.7257,  3.6476,  3.5729,  3.6304,  3.4834&
      ,  4.4368,  4.3921,  4.2207,  4.1133,  3.8067,  3.7421,  3.6140&
      ,  3.5491,  3.5077,  3.4887,  3.4623,  3.6956,  3.9568,  3.8976&
      ,  3.8240,  3.7684,  3.8451,  3.6949,  4.6318,  4.6559,  4.4533&
      ,  4.4956,  4.4641,  4.4366,  4.4155,  4.3936,  4.3764,  4.1302&
      ,  4.3398,  4.2283,  4.2796,  4.2547,  4.2391,  4.2296,  4.0699&
      ,  4.1083,  4.0319,  3.9855,  2.7676,  2.8078,  3.6725,  3.4804&
      ,  3.3775,  3.2411,  3.1581,  3.1983,  2.9973,  2.8705,  3.8070&
      ,  3.7392,  3.7668,  3.6263,  3.5402,  3.6807,  3.4545,  3.2962&
      ,  4.2283,  4.1698,  4.0240,  3.9341,  3.8711,  3.5489,  3.7798&
      /)
      r0ab(2801:2870)=(/&
         3.7000,  3.6654,  3.6154,  3.5882,  3.5472,  3.7289,  3.6510&
      ,  3.6078,  3.5355,  3.5963,  3.4480,  4.3587,  4.3390,  4.1635&
      ,  4.0536,  3.7193,  3.6529,  3.5512,  3.4837,  3.4400,  3.4191&
      ,  3.3891,  3.6622,  3.8934,  3.8235,  3.7823,  3.7292,  3.8106&
      ,  3.6589,  4.5535,  4.6013,  4.3961,  4.4423,  4.4109,  4.3835&
      ,  4.3625,  4.3407,  4.3237,  4.0863,  4.2835,  4.1675,  4.2272&
      ,  4.2025,  4.1869,  4.1774,  4.0126,  4.0460,  3.9815,  3.9340&
      ,  3.8955,  2.6912,  2.7604,  3.6037,  3.4194,  3.3094,  3.1710&
      ,  3.0862,  3.1789,  2.9738,  2.8427,  3.7378,  3.6742,  3.6928&
      ,  3.5512,  3.4614,  3.4087,  3.4201,  3.2607,  4.1527,  4.0977&
      /)
      r0ab(2871:2940)=(/&
         3.9523,  3.8628,  3.8002,  3.4759,  3.7102,  3.6466,  3.6106&
      ,  3.5580,  3.5282,  3.4878,  3.6547,  3.5763,  3.5289,  3.5086&
      ,  3.5593,  3.4099,  4.2788,  4.2624,  4.0873,  3.9770,  3.6407&
      ,  3.5743,  3.5178,  3.4753,  3.3931,  3.3694,  3.3339,  3.6002&
      ,  3.8164,  3.7478,  3.7028,  3.6952,  3.7669,  3.6137,  4.4698&
      ,  4.5488,  4.3168,  4.3646,  4.3338,  4.3067,  4.2860,  4.2645&
      ,  4.2478,  4.0067,  4.2349,  4.0958,  4.1543,  4.1302,  4.1141&
      ,  4.1048,  3.9410,  3.9595,  3.8941,  3.8465,  3.8089,  3.7490&
      ,  2.7895,  2.5849,  3.6484,  3.0162,  3.1267,  3.2125,  3.0043&
      ,  2.9572,  2.8197,  2.7261,  3.7701,  3.2446,  3.5239,  3.4696&
      /)
      r0ab(2941:3010)=(/&
         3.4261,  3.3508,  3.1968,  3.0848,  4.1496,  3.6598,  3.5111&
      ,  3.4199,  3.3809,  3.5382,  3.2572,  3.2100,  3.1917,  3.1519&
      ,  3.1198,  3.1005,  3.5071,  3.5086,  3.5073,  3.4509,  3.3120&
      ,  3.2082,  4.2611,  3.8117,  3.6988,  3.5646,  3.6925,  3.6295&
      ,  3.5383,  3.4910,  3.4625,  3.4233,  3.4007,  3.2329,  3.6723&
      ,  3.6845,  3.6876,  3.6197,  3.4799,  3.3737,  4.4341,  4.0525&
      ,  3.9011,  3.8945,  3.8635,  3.8368,  3.8153,  3.7936,  3.7758&
      ,  3.4944,  3.4873,  3.9040,  3.7110,  3.6922,  3.6799,  3.6724&
      ,  3.5622,  3.6081,  3.5426,  3.4922,  3.4498,  3.3984,  3.4456&
      ,  2.7522,  2.5524,  3.5742,  2.9508,  3.0751,  3.0158,  2.9644&
      /)
      r0ab(3011:3080)=(/&
         2.8338,  2.7891,  2.6933,  3.6926,  3.1814,  3.4528,  3.4186&
      ,  3.3836,  3.2213,  3.1626,  3.0507,  4.0548,  3.5312,  3.4244&
      ,  3.3409,  3.2810,  3.4782,  3.1905,  3.1494,  3.1221,  3.1128&
      ,  3.0853,  3.0384,  3.4366,  3.4562,  3.4638,  3.3211,  3.2762&
      ,  3.1730,  4.1632,  3.6825,  3.5822,  3.4870,  3.6325,  3.5740&
      ,  3.4733,  3.4247,  3.3969,  3.3764,  3.3525,  3.1984,  3.5989&
      ,  3.6299,  3.6433,  3.4937,  3.4417,  3.3365,  4.3304,  3.9242&
      ,  3.7793,  3.7623,  3.7327,  3.7071,  3.6860,  3.6650,  3.6476&
      ,  3.3849,  3.3534,  3.8216,  3.5870,  3.5695,  3.5584,  3.5508&
      ,  3.4856,  3.5523,  3.4934,  3.4464,  3.4055,  3.3551,  3.3888&
      /)
      r0ab(3081:3150)=(/&
         3.3525,  2.7202,  2.5183,  3.4947,  2.8731,  3.0198,  3.1457&
      ,  2.9276,  2.7826,  2.7574,  2.6606,  3.6090,  3.0581,  3.3747&
      ,  3.3677,  3.3450,  3.1651,  3.1259,  3.0147,  3.9498,  3.3857&
      ,  3.2917,  3.2154,  3.1604,  3.4174,  3.0735,  3.0342,  3.0096&
      ,  3.0136,  2.9855,  2.9680,  3.3604,  3.4037,  3.4243,  3.2633&
      ,  3.1810,  3.1351,  4.0557,  3.5368,  3.4526,  3.3699,  3.5707&
      ,  3.5184,  3.4085,  3.3595,  3.3333,  3.3143,  3.3041,  3.1094&
      ,  3.5193,  3.5745,  3.6025,  3.4338,  3.3448,  3.2952,  4.2158&
      ,  3.7802,  3.6431,  3.6129,  3.5853,  3.5610,  3.5406,  3.5204&
      ,  3.5036,  3.2679,  3.2162,  3.7068,  3.4483,  3.4323,  3.4221&
      /)
      r0ab(3151:3220)=(/&
         3.4138,  3.3652,  3.4576,  3.4053,  3.3618,  3.3224,  3.2711&
      ,  3.3326,  3.2950,  3.2564,  2.5315,  2.6104,  3.2734,  3.2299&
      ,  3.1090,  2.9942,  2.9159,  2.8324,  2.8350,  2.7216,  3.3994&
      ,  3.4475,  3.4354,  3.3438,  3.2807,  3.2169,  3.2677,  3.1296&
      ,  3.7493,  3.8075,  3.6846,  3.6104,  3.5577,  3.2052,  3.4803&
      ,  3.4236,  3.3845,  3.3640,  3.3365,  3.3010,  3.3938,  3.3624&
      ,  3.3440,  3.3132,  3.4035,  3.2754,  3.8701,  3.9523,  3.8018&
      ,  3.7149,  3.3673,  3.3199,  3.2483,  3.2069,  3.1793,  3.1558&
      ,  3.1395,  3.4097,  3.5410,  3.5228,  3.5116,  3.4921,  3.4781&
      ,  3.4690,  4.0420,  4.1759,  4.0078,  4.0450,  4.0189,  3.9952&
      /)
      r0ab(3221:3290)=(/&
         3.9770,  3.9583,  3.9434,  3.7217,  3.8228,  3.7826,  3.8640&
      ,  3.8446,  3.8314,  3.8225,  3.6817,  3.7068,  3.6555,  3.6159&
      ,  3.5831,  3.5257,  3.2133,  3.1689,  3.1196,  3.3599,  2.9852&
      ,  2.7881,  3.5284,  3.3493,  3.6958,  3.3642,  3.1568,  3.0055&
      ,  2.9558,  2.8393,  3.6287,  3.5283,  4.1511,  3.8259,  3.6066&
      ,  3.4527,  3.3480,  3.2713,  3.9037,  3.8361,  3.8579,  3.7311&
      ,  3.6575,  3.5176,  3.5693,  3.5157,  3.4814,  3.4559,  3.4445&
      ,  3.4160,  4.1231,  3.8543,  3.6816,  3.5602,  3.4798,  3.4208&
      ,  4.0542,  4.0139,  4.0165,  3.9412,  3.7698,  3.6915,  3.6043&
      ,  3.5639,  3.5416,  3.5247,  3.5153,  3.5654,  4.2862,  4.0437&
      /)
      r0ab(3291:3360)=(/&
         3.8871,  3.7741,  3.6985,  3.6413,  4.2345,  4.3663,  4.3257&
      ,  4.0869,  4.0612,  4.0364,  4.0170,  3.9978,  3.9834,  3.9137&
      ,  3.8825,  3.8758,  3.9143,  3.8976,  3.8864,  3.8768,  3.9190&
      ,  4.1613,  4.0566,  3.9784,  3.9116,  3.8326,  3.7122,  3.6378&
      ,  3.5576,  3.5457,  4.3127,  3.1160,  2.8482,  4.0739,  3.3599&
      ,  3.5698,  3.5366,  3.2854,  3.1039,  2.9953,  2.9192,  4.1432&
      ,  3.5320,  3.9478,  4.0231,  3.7509,  3.5604,  3.4340,  3.3426&
      ,  4.3328,  3.8288,  3.7822,  3.7909,  3.6907,  3.6864,  3.5793&
      ,  3.5221,  3.4883,  3.4649,  3.4514,  3.4301,  3.9256,  4.0596&
      ,  3.8307,  3.6702,  3.5651,  3.4884,  4.4182,  4.2516,  3.9687&
      /)
      r0ab(3361:3430)=(/&
         3.9186,  3.9485,  3.8370,  3.7255,  3.6744,  3.6476,  3.6295&
      ,  3.6193,  3.5659,  4.0663,  4.2309,  4.0183,  3.8680,  3.7672&
      ,  3.6923,  4.5240,  4.4834,  4.1570,  4.3204,  4.2993,  4.2804&
      ,  4.2647,  4.2481,  4.2354,  3.8626,  3.8448,  4.2267,  4.1799&
      ,  4.1670,  3.8738,  3.8643,  3.8796,  4.0575,  4.0354,  3.9365&
      ,  3.8611,  3.7847,  3.7388,  3.6826,  3.6251,  3.5492,  4.0889&
      ,  4.2764,  3.1416,  2.8325,  3.7735,  3.3787,  3.4632,  3.5923&
      ,  3.3214,  3.1285,  3.0147,  2.9366,  3.8527,  3.5602,  3.8131&
      ,  3.8349,  3.7995,  3.5919,  3.4539,  3.3540,  4.0654,  3.8603&
      ,  3.7972,  3.7358,  3.7392,  3.8157,  3.6055,  3.5438,  3.5089&
      /)
      r0ab(3431:3500)=(/&
         3.4853,  3.4698,  3.4508,  3.7882,  3.8682,  3.8837,  3.7055&
      ,  3.5870,  3.5000,  4.1573,  4.0005,  3.9568,  3.8936,  3.9990&
      ,  3.9433,  3.8172,  3.7566,  3.7246,  3.7033,  3.6900,  3.5697&
      ,  3.9183,  4.0262,  4.0659,  3.8969,  3.7809,  3.6949,  4.2765&
      ,  4.2312,  4.1401,  4.0815,  4.0580,  4.0369,  4.0194,  4.0017&
      ,  3.9874,  3.8312,  3.8120,  3.9454,  3.9210,  3.9055,  3.8951&
      ,  3.8866,  3.8689,  3.9603,  3.9109,  3.9122,  3.8233,  3.7438&
      ,  3.7436,  3.6981,  3.6555,  3.5452,  3.9327,  4.0658,  4.1175&
      ,  2.9664,  2.8209,  3.5547,  3.3796,  3.3985,  3.3164,  3.2364&
      ,  3.1956,  3.0370,  2.9313,  3.6425,  3.5565,  3.7209,  3.7108&
      /)
      r0ab(3501:3570)=(/&
         3.6639,  3.6484,  3.4745,  3.3492,  3.8755,  4.2457,  3.7758&
      ,  3.7161,  3.6693,  3.6155,  3.5941,  3.5643,  3.5292,  3.4950&
      ,  3.4720,  3.4503,  3.6936,  3.7392,  3.7388,  3.7602,  3.6078&
      ,  3.4960,  3.9800,  4.3518,  4.2802,  3.8580,  3.8056,  3.7527&
      ,  3.7019,  3.6615,  3.5768,  3.5330,  3.5038,  3.5639,  3.8192&
      ,  3.8883,  3.9092,  3.9478,  3.7995,  3.6896,  4.1165,  4.5232&
      ,  4.4357,  4.4226,  4.4031,  4.3860,  4.3721,  4.3580,  4.3466&
      ,  4.2036,  4.2037,  3.8867,  4.2895,  4.2766,  4.2662,  4.2598&
      ,  3.8408,  3.9169,  3.8681,  3.8250,  3.7855,  3.7501,  3.6753&
      ,  3.5499,  3.4872,  3.5401,  3.8288,  3.9217,  3.9538,  4.0054&
      /)
      r0ab(3571:3640)=(/&
         2.8388,  2.7890,  3.4329,  3.5593,  3.3488,  3.2486,  3.1615&
      ,  3.1000,  3.0394,  2.9165,  3.5267,  3.7479,  3.6650,  3.6263&
      ,  3.5658,  3.5224,  3.4762,  3.3342,  3.7738,  4.0333,  3.9568&
      ,  3.8975,  3.8521,  3.4929,  3.7830,  3.7409,  3.7062,  3.6786&
      ,  3.6471,  3.6208,  3.6337,  3.6519,  3.6363,  3.6278,  3.6110&
      ,  3.4825,  3.8795,  4.1448,  4.0736,  4.0045,  3.6843,  3.6291&
      ,  3.5741,  3.5312,  3.4974,  3.4472,  3.4034,  3.7131,  3.7557&
      ,  3.7966,  3.8005,  3.8068,  3.8015,  3.6747,  4.0222,  4.3207&
      ,  4.2347,  4.2191,  4.1990,  4.1811,  4.1666,  4.1521,  4.1401&
      ,  3.9970,  3.9943,  3.9592,  4.0800,  4.0664,  4.0559,  4.0488&
      /)
      r0ab(3641:3710)=(/&
         3.9882,  4.0035,  3.9539,  3.9138,  3.8798,  3.8355,  3.5359&
      ,  3.4954,  3.3962,  3.5339,  3.7595,  3.8250,  3.8408,  3.8600&
      ,  3.8644,  2.7412,  2.7489,  3.3374,  3.3950,  3.3076,  3.1910&
      ,  3.0961,  3.0175,  3.0280,  2.8929,  3.4328,  3.5883,  3.6227&
      ,  3.5616,  3.4894,  3.4241,  3.3641,  3.3120,  3.6815,  3.8789&
      ,  3.8031,  3.7413,  3.6939,  3.4010,  3.6225,  3.5797,  3.5443&
      ,  3.5139,  3.4923,  3.4642,  3.5860,  3.5849,  3.5570,  3.5257&
      ,  3.4936,  3.4628,  3.7874,  3.9916,  3.9249,  3.8530,  3.5932&
      ,  3.5355,  3.4757,  3.4306,  3.3953,  3.3646,  3.3390,  3.5637&
      ,  3.7053,  3.7266,  3.7177,  3.6996,  3.6775,  3.6558,  3.9331&
      /)
      r0ab(3711:3780)=(/&
         4.1655,  4.0879,  4.0681,  4.0479,  4.0299,  4.0152,  4.0006&
      ,  3.9883,  3.8500,  3.8359,  3.8249,  3.9269,  3.9133,  3.9025&
      ,  3.8948,  3.8422,  3.8509,  3.7990,  3.7570,  3.7219,  3.6762&
      ,  3.4260,  3.3866,  3.3425,  3.5294,  3.7022,  3.7497,  3.7542&
      ,  3.7494,  3.7370,  3.7216,  3.4155,  3.0522,  4.2541,  3.8218&
      ,  4.0438,  3.5875,  3.3286,  3.1682,  3.0566,  2.9746,  4.3627&
      ,  4.0249,  4.6947,  4.1718,  3.8639,  3.6735,  3.5435,  3.4479&
      ,  4.6806,  4.3485,  4.2668,  4.1690,  4.1061,  4.1245,  4.0206&
      ,  3.9765,  3.9458,  3.9217,  3.9075,  3.8813,  3.9947,  4.1989&
      ,  3.9507,  3.7960,  3.6925,  3.6150,  4.8535,  4.5642,  4.4134&
      /)
      r0ab(3781:3850)=(/&
         4.3688,  4.3396,  4.2879,  4.2166,  4.1888,  4.1768,  4.1660&
      ,  4.1608,  4.0745,  4.2289,  4.4863,  4.2513,  4.0897,  3.9876&
      ,  3.9061,  5.0690,  5.0446,  4.6186,  4.6078,  4.5780,  4.5538&
      ,  4.5319,  4.5101,  4.4945,  4.1912,  4.2315,  4.5534,  4.4373&
      ,  4.4224,  4.4120,  4.4040,  4.2634,  4.7770,  4.6890,  4.6107&
      ,  4.5331,  4.4496,  4.4082,  4.3095,  4.2023,  4.0501,  4.2595&
      ,  4.5497,  4.3056,  4.1506,  4.0574,  3.9725,  5.0796,  3.0548&
      ,  3.3206,  3.8132,  3.9720,  3.7675,  3.7351,  3.5167,  3.5274&
      ,  3.3085,  3.1653,  3.9500,  4.1730,  4.0613,  4.1493,  3.8823&
      ,  4.0537,  3.8200,  3.6582,  4.3422,  4.5111,  4.3795,  4.3362&
      /)
      r0ab(3851:3920)=(/&
         4.2751,  3.7103,  4.1973,  4.1385,  4.1129,  4.0800,  4.0647&
      ,  4.0308,  4.0096,  4.1619,  3.9360,  4.1766,  3.9705,  3.8262&
      ,  4.5348,  4.7025,  4.5268,  4.5076,  3.9562,  3.9065,  3.8119&
      ,  3.7605,  3.7447,  3.7119,  3.6916,  4.1950,  4.2110,  4.3843&
      ,  4.1631,  4.4427,  4.2463,  4.1054,  4.7693,  5.0649,  4.7365&
      ,  4.7761,  4.7498,  4.7272,  4.7076,  4.6877,  4.6730,  4.4274&
      ,  4.5473,  4.5169,  4.5975,  4.5793,  4.5667,  4.5559,  4.3804&
      ,  4.6920,  4.6731,  4.6142,  4.5600,  4.4801,  4.0149,  3.8856&
      ,  3.7407,  4.1545,  4.2253,  4.4229,  4.1923,  4.5022,  4.3059&
      ,  4.1591,  4.7883,  4.9294,  3.3850,  3.4208,  3.7004,  3.8800&
      /)
      r0ab(3921:3990)=(/&
         3.9886,  3.9040,  3.6719,  3.6547,  3.4625,  3.3370,  3.8394&
      ,  4.0335,  4.2373,  4.3023,  4.0306,  4.1408,  3.9297,  3.7857&
      ,  4.1907,  4.3230,  4.2664,  4.2173,  4.1482,  3.6823,  4.0711&
      ,  4.0180,  4.0017,  3.9747,  3.9634,  3.9383,  4.1993,  4.3205&
      ,  4.0821,  4.2547,  4.0659,  3.9359,  4.3952,  4.5176,  4.3888&
      ,  4.3607,  3.9583,  3.9280,  3.8390,  3.7971,  3.7955,  3.7674&
      ,  3.7521,  4.1062,  4.3633,  4.2991,  4.2767,  4.4857,  4.3039&
      ,  4.1762,  4.6197,  4.8654,  4.6633,  4.5878,  4.5640,  4.5422&
      ,  4.5231,  4.5042,  4.4901,  4.3282,  4.3978,  4.3483,  4.4202&
      ,  4.4039,  4.3926,  4.3807,  4.2649,  4.6135,  4.5605,  4.5232&
      /)
      r0ab(3991:4060)=(/&
         4.4676,  4.3948,  4.0989,  3.9864,  3.8596,  4.0942,  4.2720&
      ,  4.3270,  4.3022,  4.5410,  4.3576,  4.2235,  4.6545,  4.7447&
      ,  4.7043,  3.0942,  3.2075,  3.5152,  3.6659,  3.8289,  3.7459&
      ,  3.5156,  3.5197,  3.3290,  3.2069,  3.6702,  3.8448,  4.0340&
      ,  3.9509,  3.8585,  3.9894,  3.7787,  3.6365,  4.1425,  4.1618&
      ,  4.0940,  4.0466,  3.9941,  3.5426,  3.8952,  3.8327,  3.8126&
      ,  3.7796,  3.7635,  3.7356,  4.0047,  3.9655,  3.9116,  4.1010&
      ,  3.9102,  3.7800,  4.2964,  4.3330,  4.2622,  4.2254,  3.8195&
      ,  3.7560,  3.6513,  3.5941,  3.5810,  3.5420,  3.5178,  3.8861&
      ,  4.1459,  4.1147,  4.0772,  4.3120,  4.1207,  3.9900,  4.4733&
      /)
      r0ab(4061:4130)=(/&
         4.6157,  4.4580,  4.4194,  4.3954,  4.3739,  4.3531,  4.3343&
      ,  4.3196,  4.2140,  4.2339,  4.1738,  4.2458,  4.2278,  4.2158&
      ,  4.2039,  4.1658,  4.3595,  4.2857,  4.2444,  4.1855,  4.1122&
      ,  3.7839,  3.6879,  3.5816,  3.8633,  4.1585,  4.1402,  4.1036&
      ,  4.3694,  4.1735,  4.0368,  4.5095,  4.5538,  4.5240,  4.4252&
      ,  3.0187,  3.1918,  3.5127,  3.6875,  3.7404,  3.6943,  3.4702&
      ,  3.4888,  3.2914,  3.1643,  3.6669,  3.8724,  3.9940,  4.0816&
      ,  3.8054,  3.9661,  3.7492,  3.6024,  4.0428,  4.1951,  4.1466&
      ,  4.0515,  4.0075,  3.5020,  3.9158,  3.8546,  3.8342,  3.8008&
      ,  3.7845,  3.7549,  3.9602,  3.8872,  3.8564,  4.0793,  3.8835&
      /)
      r0ab(4131:4200)=(/&
         3.7495,  4.2213,  4.3704,  4.3300,  4.2121,  3.7643,  3.7130&
      ,  3.6144,  3.5599,  3.5474,  3.5093,  3.4853,  3.9075,  4.1115&
      ,  4.0473,  4.0318,  4.2999,  4.1050,  3.9710,  4.4320,  4.6706&
      ,  4.5273,  4.4581,  4.4332,  4.4064,  4.3873,  4.3684,  4.3537&
      ,  4.2728,  4.2549,  4.2032,  4.2794,  4.2613,  4.2491,  4.2375&
      ,  4.2322,  4.3665,  4.3061,  4.2714,  4.2155,  4.1416,  3.7660&
      ,  3.6628,  3.5476,  3.8790,  4.1233,  4.0738,  4.0575,  4.3575&
      ,  4.1586,  4.0183,  4.4593,  4.5927,  4.4865,  4.3813,  4.4594&
      ,  2.9875,  3.1674,  3.4971,  3.6715,  3.7114,  3.6692,  3.4446&
      ,  3.4676,  3.2685,  3.1405,  3.6546,  3.8579,  3.9637,  4.0581&
      /)
      r0ab(4201:4270)=(/&
         3.7796,  3.9463,  3.7275,  3.5792,  4.0295,  4.1824,  4.1247&
      ,  4.0357,  3.9926,  3.4827,  3.9007,  3.8392,  3.8191,  3.7851&
      ,  3.7687,  3.7387,  3.9290,  3.8606,  3.8306,  4.0601,  3.8625&
      ,  3.7269,  4.2062,  4.3566,  4.3022,  4.1929,  3.7401,  3.6888&
      ,  3.5900,  3.5350,  3.5226,  3.4838,  3.4594,  3.8888,  4.0813&
      ,  4.0209,  4.0059,  4.2810,  4.0843,  3.9486,  4.4162,  4.6542&
      ,  4.5005,  4.4444,  4.4196,  4.3933,  4.3741,  4.3552,  4.3406&
      ,  4.2484,  4.2413,  4.1907,  4.2656,  4.2474,  4.2352,  4.2236&
      ,  4.2068,  4.3410,  4.2817,  4.2479,  4.1921,  4.1182,  3.7346&
      ,  3.6314,  3.5168,  3.8582,  4.0927,  4.0469,  4.0313,  4.3391&
      /)
      r0ab(4271:4340)=(/&
         4.1381,  3.9962,  4.4429,  4.5787,  4.4731,  4.3588,  4.4270&
      ,  4.3957,  2.9659,  3.1442,  3.4795,  3.6503,  3.6814,  3.6476&
      ,  3.4222,  3.4491,  3.2494,  3.1209,  3.6324,  3.8375,  3.9397&
      ,  3.8311,  3.7581,  3.9274,  3.7085,  3.5598,  4.0080,  4.1641&
      ,  4.1057,  4.0158,  3.9726,  3.4667,  3.8802,  3.8188,  3.7989&
      ,  3.7644,  3.7474,  3.7173,  3.9049,  3.8424,  3.8095,  4.0412&
      ,  3.8436,  3.7077,  4.1837,  4.3366,  4.2816,  4.1686,  3.7293&
      ,  3.6709,  3.5700,  3.5153,  3.5039,  3.4684,  3.4437,  3.8663&
      ,  4.0575,  4.0020,  3.9842,  4.2612,  4.0643,  3.9285,  4.3928&
      ,  4.6308,  4.4799,  4.4244,  4.3996,  4.3737,  4.3547,  4.3358&
      /)
      r0ab(4341:4410)=(/&
         4.3212,  4.2275,  4.2216,  4.1676,  4.2465,  4.2283,  4.2161&
      ,  4.2045,  4.1841,  4.3135,  4.2562,  4.2226,  4.1667,  4.0932&
      ,  3.7134,  3.6109,  3.4962,  3.8352,  4.0688,  4.0281,  4.0099&
      ,  4.3199,  4.1188,  3.9768,  4.4192,  4.5577,  4.4516,  4.3365&
      ,  4.4058,  4.3745,  4.3539,  2.8763,  3.1294,  3.5598,  3.7465&
      ,  3.5659,  3.5816,  3.3599,  3.4024,  3.1877,  3.0484,  3.7009&
      ,  3.9451,  3.8465,  3.9873,  3.7079,  3.9083,  3.6756,  3.5150&
      ,  4.0829,  4.2780,  4.1511,  4.1260,  4.0571,  3.4865,  3.9744&
      ,  3.9150,  3.8930,  3.8578,  3.8402,  3.8073,  3.7977,  4.0036&
      ,  3.7604,  4.0288,  3.8210,  3.6757,  4.2646,  4.4558,  4.2862&
      /)
      r0ab(4411:4465)=(/&
         4.2122,  3.7088,  3.6729,  3.5800,  3.5276,  3.5165,  3.4783&
      ,  3.4539,  3.9553,  3.9818,  4.2040,  3.9604,  4.2718,  4.0689&
      ,  3.9253,  4.4869,  4.7792,  4.4918,  4.5342,  4.5090,  4.4868&
      ,  4.4680,  4.4486,  4.4341,  4.2023,  4.3122,  4.2710,  4.3587&
      ,  4.3407,  4.3281,  4.3174,  4.1499,  4.3940,  4.3895,  4.3260&
      ,  4.2725,  4.1961,  3.7361,  3.6193,  3.4916,  3.9115,  3.9914&
      ,  3.9809,  3.9866,  4.3329,  4.1276,  3.9782,  4.5097,  4.6769&
      ,  4.5158,  4.3291,  4.3609,  4.3462,  4.3265,  4.4341&
      /)

      k=0
      DO i=1,SIZE(rout,1)
         DO j=1,i
            k=k+1
            rout(i,j)=r0ab(k)*bohr
            rout(j,i)=r0ab(k)*bohr
         ENDDO
      ENDDO

      ! covalent radii (taken from Pyykko and Atsumi, Chem. Eur. J. 15, 2009, 188-197)
      ! values for metals decreased by 10 %
      rcov(1:94)=(/&
        0.32, 0.46, 1.20, 0.94, 0.77, 0.75, 0.71, 0.63, 0.64, 0.67 &
      , 1.40, 1.25, 1.13, 1.04, 1.10, 1.02, 0.99, 0.96, 1.76, 1.54 &
      , 1.33, 1.22, 1.21, 1.10, 1.07, 1.04, 1.00, 0.99, 1.01, 1.09 &
      , 1.12, 1.09, 1.15, 1.10, 1.14, 1.17, 1.89, 1.67, 1.47, 1.39 &
      , 1.32, 1.24, 1.15, 1.13, 1.13, 1.08, 1.15, 1.23, 1.28, 1.26 &
      , 1.26, 1.23, 1.32, 1.31, 2.09, 1.76, 1.62, 1.47, 1.58, 1.57 &
      , 1.56, 1.55, 1.51, 1.52, 1.51, 1.50, 1.49, 1.49, 1.48, 1.53 &
      , 1.46, 1.37, 1.31, 1.23, 1.18, 1.16, 1.11, 1.12, 1.13, 1.32 &
      , 1.30, 1.30, 1.36, 1.31, 1.38, 1.42, 2.01, 1.81, 1.67, 1.58 &
      , 1.52, 1.53, 1.54, 1.55  &
      /)

      ! PBE0/def2-QZVP atomic values
      r2r4(1:94)=(/&
        8.0589,  3.4698, 29.0974, 14.8517, 11.8799,  7.8715,  5.5588,&
        4.7566,  3.8025,  3.1036, 26.1552, 17.2304, 17.7210, 12.7442,&
        9.5361,  8.1652,  6.7463,  5.6004, 29.2012, 22.3934, 19.0598,&
       16.8590, 15.4023, 12.5589, 13.4788, 12.2309, 11.2809, 10.5569,&
       10.1428,  9.4907, 13.4606, 10.8544,  8.9386,  8.1350,  7.1251,&
        6.1971, 30.0162, 24.4103, 20.3537, 17.4780, 13.5528, 11.8451,&
       11.0355, 10.1997,  9.5414,  9.0061,  8.6417,  8.9975, 14.0834,&
       11.8333, 10.0179,  9.3844,  8.4110,  7.5152, 32.7622, 27.5708,&
       23.1671, 21.6003, 20.9615, 20.4562, 20.1010, 19.7475, 19.4828,&
       15.6013, 19.2362, 17.4717, 17.8321, 17.4237, 17.1954, 17.1631,&
       14.5716, 15.8758, 13.8989, 12.4834, 11.4421, 10.2671,  8.3549,&
        7.8496,  7.3278,  7.4820, 13.5124, 11.6554, 10.0959,  9.7340,&
        8.8584,  8.0125, 29.8135, 26.3157, 19.1885, 15.8542, 16.1305,&
       15.6161, 15.1226, 16.1576&
      /)

  END SUBROUTINE setr0ab

! *****************************************************************************
!> \brief ...
!> \param cnout ...
! *****************************************************************************
  SUBROUTINE setcn(cnout)
      ! default coordination numbers
    REAL(KIND=dp), DIMENSION(:)              :: cnout

    INTEGER                                  :: n
    REAL(KIND=dp), DIMENSION(104)            :: cn

      cn(   1: 104)=(/&
         1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 3.0, 2.0, 1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0&
      ,  2.0, 1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 3.0, 2.0, 3.0, 2.0, 2.0, 2.0, 2.0&
      ,  3.0, 4.0, 3.0, 4.0, 1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 4.0, 3.0&
      ,  2.0, 1.0, 2.0, 3.0, 2.0, 3.0, 4.0, 1.0, 0.0, 1.0, 2.0, 3.0, 3.0, 3.0, 3.0&
      ,  3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 4.0, 5.0, 6.0, 7.0&
      ,  4.0, 4.0, 4.0, 3.0, 2.0, 1.0, 2.0, 3.0, 4.0, 1.0, 0.0, 1.0, 2.0, 3.0, 4.0&
      ,  5.0, 6.0, 5.0, 4.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 2.0, 3.0, 3.0&
      /)

      cnout = 0._dp
      n = MIN(SIZE(cn),SIZE(cnout))
      cnout(1:n) = cn(1:n)

  END SUBROUTINE setcn

! *****************************************************************************
!> \brief ...
!> \param rab ...
!> \param rcova ...
!> \param rcovb ...
!> \param k1 ...
!> \param cnab ...
!> \param dcnab ...
! *****************************************************************************
  SUBROUTINE cnparam_d3(rab,rcova,rcovb,k1,cnab,dcnab)

    REAL(KIND=dp), INTENT(IN)                :: rab, rcova, rcovb, k1
    REAL(KIND=dp), INTENT(OUT)               :: cnab, dcnab

    REAL(KIND=dp)                            :: dfz, ee, fz, rco, rr

! covalent distance in Bohr

    rco=rcova+rcovb
    rr=rco/rab
    ! counting function exponential has a better long-range behavior
    ! than MHGs inverse damping
    ! factor k2 already included into rcov
    ee = EXP(-k1*(rr-1.0_dp))
    ! force the function to zero using a second step function
    fz = 0.5_dp*(1.0_dp - TANH(rab-2.0_dp*rco))
    dfz = 0.5_dp*((TANH(rab-2.0_dp*rco))**2 - 1.0_dp)
    cnab = 1.0_dp/(1.0_dp+ee) * fz 
    dcnab = -cnab*cnab*k1*rr/rab*ee + 1.0_dp/(1.0_dp+ee)*dfz

  END SUBROUTINE cnparam_d3

! *****************************************************************************
!> \brief ...
!> \param rab ...
!> \param rcutab ...
!> \param srn ...
!> \param alpn ...
!> \param rcut ...
!> \param fdab ...
!> \param dfdab ...
! *****************************************************************************
  SUBROUTINE damping_d3(rab,rcutab,srn,alpn,rcut,fdab,dfdab)

    REAL(KIND=dp), INTENT(IN)                :: rab, rcutab, srn, alpn, rcut
    REAL(KIND=dp), INTENT(OUT)               :: fdab, dfdab

    REAL(KIND=dp)                            :: a, b, c, d, dd, dfab, dfcc, &
                                                dz, fab, fcc, rl, rr, ru, z, &
                                                zz

    rl = rcut - 1._dp
    ru = rcut
    IF (rab >= ru) THEN
      fcc = 0._dp
      dfcc = 0._dp
    ELSEIF (rab <= rl) THEN
      fcc = 1._dp
      dfcc = 0._dp
    ELSE
      z = rab*rab - rl*rl
      dz = 2._dp*rab
      zz = z*z*z
      d = ru*ru-rl*rl
      dd = d*d*d
      a = -10._dp/dd
      b = 15._dp/(dd*d)
      c = -6._dp/(dd*d*d)
      fcc = 1._dp+zz*(a+b*z+c*z*z)
      dfcc = zz*dz/z * (3._dp*a+4._dp*b*z+5._dp*c*z*z)
    END IF

    rr = 6._dp*(rab/(srn*rcutab))**(-alpn)
    fab = 1._dp/(1._dp + rr)
    dfab = fab*fab*rr*alpn/rab
    fdab = fab*fcc
    dfdab = dfab*fcc + fab*dfcc

  END SUBROUTINE damping_d3

! *****************************************************************************
!> \brief ...
!> \param maxc ...
!> \param max_elem ...
!> \param c6ab ...
!> \param mxc ...
!> \param iat ...
!> \param jat ...
!> \param nci ...
!> \param ncj ...
!> \param k3 ...
!> \param c6 ...
!> \param dc6a ...
!> \param dc6b ...
! *****************************************************************************
  SUBROUTINE getc6(maxc,max_elem,c6ab,mxc,iat,jat,nci,ncj,k3,c6,dc6a,dc6b)

    INTEGER, INTENT(IN)                      :: maxc, max_elem
    REAL(KIND=dp), INTENT(IN) :: c6ab(max_elem,max_elem,maxc,maxc,3)
    INTEGER, INTENT(IN)                      :: mxc(max_elem), iat, jat
    REAL(KIND=dp), INTENT(IN)                :: nci, ncj, k3
    REAL(KIND=dp), INTENT(OUT)               :: c6, dc6a, dc6b

    INTEGER                                  :: i, j
    REAL(KIND=dp)                            :: c6mem, cn1, cn2, csum, dtmpa, &
                                                dtmpb, dwa, dwb, dza, dzb, r, &
                                                rsave, rsum, tmp1

! the exponential is sensitive to numerics
! when nci or ncj is much larger than cn1/cn2

    c6mem=-1.0e+99_dp
    rsave=1.0e+99_dp
    rsum=0.0_dp
    csum=0.0_dp
    dza=0.0_dp
    dzb=0.0_dp
    dwa=0.0_dp
    dwb=0.0_dp
    c6=0.0_dp
    DO i=1,mxc(iat)
      DO j=1,mxc(jat)
         c6=c6ab(iat,jat,i,j,1)
         IF(c6 > 0.0_dp)THEN
            cn1=c6ab(iat,jat,i,j,2)
            cn2=c6ab(iat,jat,i,j,3)
            ! distance
            r=(cn1-nci)**2+(cn2-ncj)**2
            IF (r < rsave) THEN
               rsave=r
               c6mem=c6
            ENDIF
            tmp1=EXP(k3*r)
            dtmpa=-2.0_dp*k3*(cn1-nci)*tmp1
            dtmpb=-2.0_dp*k3*(cn2-ncj)*tmp1
            rsum=rsum+tmp1
            csum=csum+tmp1*c6
            dza=dza+dtmpa*c6
            dwa=dwa+dtmpa
            dzb=dzb+dtmpb*c6
            dwb=dwb+dtmpb
         ENDIF
      ENDDO
    ENDDO

    IF(c6 == 0.0_dp) c6mem=0.0_dp

    IF(rsum > 1.0e-66_dp)THEN
       c6=csum/rsum
       dc6a=(dza-c6*dwa)/rsum
       dc6b=(dzb-c6*dwb)/rsum
    ELSE
       c6=c6mem
       dc6a=0._dp
       dc6b=0._dp
    ENDIF

  END SUBROUTINE getc6

! *****************************************************************************
!> \brief ...
!> \param dcnum ...
!> \param para_env ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE dcnum_distribute(dcnum,para_env,error)

    TYPE(dcnum_type), DIMENSION(:)           :: dcnum
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'dcnum_distribute', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: group, i, ia, ipe, istat, jn, &
                                                n, natom, ntmax, ntot
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: list, nloc
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: dvals
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: rik

    failure = .FALSE.
! we only have to do something if this is a parallel run

    IF (para_env%num_pe > 1) THEN
      group = para_env%group
      natom=SIZE(dcnum)
      !pack my dcnum data
      ALLOCATE (nloc(natom),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      DO ia=1,natom
        nloc(ia)=dcnum(ia)%neighbors
      END DO
      ntot = SUM(nloc)
      ntmax = ntot
      CALL mp_max(ntmax,group)
      ALLOCATE (list(ntmax),dvals(ntmax),rik(3,ntmax),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      list=0
      dvals=0._dp
      rik=0._dp
      i=0
      DO ia=1,natom
        DO jn=1,nloc(ia)
          i=i+1
          list(i)=dcnum(ia)%nlist(jn)
          dvals(i)=dcnum(ia)%dvals(jn)
          rik(1:3,i)=dcnum(ia)%rik(1:3,jn)
        END DO
      END DO
      DO ipe=1,para_env%num_pe-1
        !send/receive packed data
        CALL mp_shift(nloc,group)
        !unpack received data
        CALL mp_shift(list,group)
        CALL mp_shift(dvals,group)
        CALL mp_shift(rik,group)
        !add data to local dcnum
        i=0
        DO ia=1,natom
          n=dcnum(ia)%neighbors+nloc(ia)
          IF ( n > SIZE(dcnum(ia)%nlist) ) THEN
             CALL reallocate(dcnum(ia)%nlist,1,2*n)
             CALL reallocate(dcnum(ia)%dvals,1,2*n)
             CALL reallocate(dcnum(ia)%rik,1,3,1,2*n)
          END IF
          DO jn=1,nloc(ia)
            i=i+1
            n=dcnum(ia)%neighbors+jn
            dcnum(ia)%nlist(n)=list(i)
            dcnum(ia)%dvals(n)=dvals(i)
            dcnum(ia)%rik(1:3,n)=rik(1:3,i)
          END DO
          dcnum(ia)%neighbors=dcnum(ia)%neighbors+nloc(ia)
        END DO
      END DO
      DEALLOCATE (nloc,STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      DEALLOCATE (list,dvals,rik,STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    END IF

  END SUBROUTINE dcnum_distribute

! *****************************************************************************

END MODULE qs_dispersion_pairpot
