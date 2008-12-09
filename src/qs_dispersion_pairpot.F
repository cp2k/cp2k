!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculation of dispersion using pair potentials
!> \author JGH
! *****************************************************************************
MODULE qs_dispersion_pairpot

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set,&
                                             set_atomic_kind
  USE bibliography,                    ONLY: cite_reference,&
                                             grimme2006
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_parser_methods,               ONLY: parser_get_next_line,&
                                             parser_get_object
  USE cp_parser_types,                 ONLY: cp_parser_type,&
                                             parser_create,&
                                             parser_release
  USE f77_blas
  USE input_constants,                 ONLY: vdw_pairpot_grimme,&
                                             xc_vdw_fun_none,&
                                             xc_vdw_fun_pairpot
  USE input_cp2k_check,                ONLY: xc_functionals_expand
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_get_subs_vals2,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE message_passing,                 ONLY: mp_sum
  USE particle_types,                  ONLY: particle_type
  USE physcon,                         ONLY: bohr,&
                                             kjmol
  USE qs_dispersion_types,             ONLY: grimme_pp,&
                                             qs_atom_dispersion_type,&
                                             qs_dispersion_type
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_neighbor_list_types,          ONLY: &
       first_list, first_node, get_neighbor_list, get_neighbor_list_set, &
       get_neighbor_node, neighbor_list_set_p_type, neighbor_list_type, &
       neighbor_node_type, next
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE virial_methods,                  ONLY: virial_pair_force
  USE virial_types,                    ONLY: virial_type
#include "cp_common_uses.h"

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_dispersion_pairpot'

  PUBLIC :: qs_dispersion_pairpot_init, calculate_dispersion_pairpot, &
            qs_scaling_init,  write_dispersion

! *****************************************************************************

CONTAINS

! *****************************************************************************
  SUBROUTINE qs_dispersion_pairpot_init(atomic_kind_set,vdw_section,para_env,error)
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(section_vals_type), POINTER         :: vdw_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_dispersion_pairpot_init', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=2)                         :: symbol
    CHARACTER(LEN=default_string_length)     :: aname, filename
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: tmpstringlist
    INTEGER                                  :: elem, i, ikind, n_rep, nkind, &
                                                stat, vdw_pp_type, vdw_type
    LOGICAL                                  :: at_end, failure, found, &
                                                is_available
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_parser_type), POINTER            :: parser
    TYPE(qs_atom_dispersion_type), POINTER   :: disp
    TYPE(section_vals_type), POINTER         :: pp_section

    failure = .FALSE.
    CALL section_vals_val_get(vdw_section,"POTENTIAL_TYPE",i_val=vdw_type,error=error)
    nkind  = SIZE(atomic_kind_set)
    DO ikind = 1, nkind
       atomic_kind => atomic_kind_set(ikind)
       CALL get_atomic_kind(atomic_kind,element_symbol=symbol,z=elem)
       ALLOCATE(disp,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       SELECT CASE (vdw_type)
       CASE DEFAULT
          CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
       CASE (xc_vdw_fun_none)
          disp%defined = .FALSE.
       CASE (xc_vdw_fun_pairpot)
          ! setup information on pair potentials
          pp_section => section_vals_get_subs_vals(vdw_section,"PAIR_POTENTIAL",error=error)
          CALL section_vals_val_get(pp_section,"TYPE",i_val=vdw_pp_type,error=error)
          SELECT CASE (vdw_pp_type)
          CASE DEFAULT
             CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
          CASE (vdw_pairpot_grimme)
             CALL cite_reference(Grimme2006)
             disp%type = grimme_pp
             ! get filename of parameter file
             CALL section_vals_val_get(pp_section,"PARAMETER_FILE_NAME",c_val=filename,error=error)
             ! check for local parameters
             CALL section_vals_val_get(pp_section,"ATOMPARM",n_rep_val=n_rep,error=error)
             found = .FALSE.
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
             IF ( .NOT. found ) THEN
                ! check for internal parameters
                CALL grimme_param(elem,disp%c6,disp%vdw_radii,found,error)
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
          END SELECT
       END SELECT
       ! Check if the parameter is defined
       CALL cp_assert(disp%defined,cp_failure_level,cp_assertion_failed,routineP,&
            "Dispersion parameters for element ("//TRIM(symbol)//") are not defined! "//&
            "Please provide a valid set of parameters through the input section or "//&
            "through an external file! "//&
 CPSourceFileRef,&
            error=error,failure=failure)
       CALL set_atomic_kind(atomic_kind,dispersion=disp,error=error)
    END DO

  END SUBROUTINE qs_dispersion_pairpot_init

! *****************************************************************************

  SUBROUTINE qs_scaling_init(scaling,dft_section,error)
    REAL(KIND=dp), INTENT(inout)             :: scaling
    TYPE(section_vals_type), POINTER         :: dft_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_scaling_init', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: pptype
    INTEGER                                  :: ifun
    LOGICAL                                  :: failure
    TYPE(section_vals_type), POINTER         :: pp_section, vdw_section, &
                                                xc_fun, xc_fun_section, &
                                                xc_section

    failure = .FALSE.

    NULLIFY(xc_section, vdw_section, xc_fun_section)
    xc_section =>  section_vals_get_subs_vals(dft_section,"XC",error=error)
    vdw_section =>  section_vals_get_subs_vals(dft_section,"XC%vdw_potential",error=error)
    xc_fun_section => section_vals_get_subs_vals(dft_section,"XC%XC_FUNCTIONAL",error=error)
    pp_section => section_vals_get_subs_vals(vdw_section,"PAIR_POTENTIAL",error=error)

    CALL section_vals_val_get(pp_section,"TYPE",c_val=pptype,error=error)

    scaling = 0._dp

    SELECT CASE (TRIM(pptype))
    CASE DEFAULT
       CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
    CASE ("GRIMME")
       ! try to find the correct correlation functional scaling parameter
       CALL xc_functionals_expand(xc_fun_section,xc_section,error=error)
       ifun=0
       DO
          ifun=ifun+1
          NULLIFY(xc_fun)
          xc_fun => section_vals_get_subs_vals2(xc_fun_section,i_section=ifun,error=error)
          IF (.NOT.ASSOCIATED(xc_fun)) EXIT
          SELECT CASE(xc_fun%section%name)
          CASE("LYP")
             IF(scaling==1.00_dp) THEN
                ! assume this is B3LYP
                scaling = 1.05_dp
             ELSE
                scaling = 1.20_dp
             END IF
          CASE("TPSS")
             scaling = 1.00_dp
          CASE("PBE")
             scaling = 0.75_dp
          CASE("P86C")
             scaling = 1.05_dp
          CASE("VWN")
             IF(scaling==1.20_dp) THEN
                ! assume this is B3LYP
                scaling = 1.05_dp
             ELSE
                scaling = 1.00_dp
             END IF
          END SELECT
       END DO
    END SELECT
    ! if the scaling parameter is still 0 use a default of 1
    IF ( scaling == 0._dp ) scaling = 1.0_dp

  END SUBROUTINE qs_scaling_init

! *****************************************************************************
  SUBROUTINE write_dispersion(qs_env,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'write_dispersion', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=2)                         :: symbol
    INTEGER                                  :: ikind, nkind, output_unit
    LOGICAL                                  :: failure
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(qs_atom_dispersion_type), POINTER   :: disp
    TYPE(qs_dispersion_type), POINTER        :: dispersion_env
    TYPE(section_vals_type), POINTER         :: dft_section

    failure=.FALSE.

    NULLIFY(logger)
    logger => cp_error_get_logger(error)

    dft_section => section_vals_get_subs_vals(qs_env%input,"DFT",error=error)
    output_unit = cp_print_key_unit_nr(logger,dft_section,&
         "PRINT%DFT_CONTROL_PARAMETERS",extension=".Log",error=error)

    IF (output_unit>0) THEN
       CALL get_qs_env(qs_env=qs_env,dispersion_env=dispersion_env,error=error)
       ! vdW type specific output
       IF ( dispersion_env%type == xc_vdw_fun_pairpot ) THEN
          WRITE (output_unit, fmt="(' vdW POTENTIAL| ',T67,'Pair Potential')") 
          ! Pair potentials
          IF ( TRIM(dispersion_env%pp_type) == "GRIMME" ) THEN
             WRITE (output_unit, fmt="(' vdW POTENTIAL| ',T35,'Potential Form: S. Grimme, JCC 27: 1787 (2006)')") 
             WRITE (output_unit, fmt="(' vdW POTENTIAL| ',T35,'Cutoff Radius [Bohr]:',T73,F8.2)") dispersion_env%rc_disp
             WRITE (output_unit, fmt="(' vdW POTENTIAL| ',T35,'Scaling Factor:',T73,F8.4)") dispersion_env%scaling
             WRITE (output_unit, fmt="(' vdW POTENTIAL| ',T35,'Exp Prefactor for Damping:',T73,F8.1)") dispersion_env%exp_pre
             CALL get_qs_env(qs_env,atomic_kind_set=atomic_kind_set,error=error)
             nkind  = SIZE(atomic_kind_set)
             DO ikind = 1, nkind
                atomic_kind => atomic_kind_set(ikind)
                CALL get_atomic_kind(atomic_kind,element_symbol=symbol,dispersion=disp)
                IF ( disp%defined ) THEN
                   WRITE (output_unit, fmt="(' vdW PARAMETER| ',T18,'Atom=',A2,"//&
                        "T28,'C6[J*nm^6*mol^-1]=',F8.4,T63,'r(vdW)[A]=',F8.4)") &
                        symbol, disp%c6/(1000._dp*bohr**6/kjmol), disp%vdw_radii/bohr
                ELSE
                   WRITE (output_unit, fmt="(' vdW PARAMETER| ',T20,'Atom=',A2,T70,'not defined')")
                END IF
             END DO
          END IF
       END IF
    END IF
    CALL cp_print_key_finished_output(output_unit,logger,dft_section,&
         "PRINT%DFT_CONTROL_PARAMETERS",error=error)

  END SUBROUTINE write_dispersion
! *****************************************************************************

  SUBROUTINE calculate_dispersion_pairpot(qs_env,calculate_forces,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(IN)                      :: calculate_forces
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_dispersion_pairpot', &
      routineP = moduleN//':'//routineN

    INTEGER :: atom_a, atom_b, handle, iab, iatom, ikind, ilist, inode, &
      istat, jatom, jkind, natom, nkind, nlist, nnode
    INTEGER, DIMENSION(:), POINTER           :: atom_of_kind
    LOGICAL                                  :: defined, failure, use_virial
    REAL(KIND=dp)                            :: c6, c6a, c6b, dd, devdw, &
                                                dfdmp, dr, er, evdw, fdmp, &
                                                rc, rc_a, rc_b, s6, xp
    REAL(KIND=dp), DIMENSION(3)              :: fdij, rij
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_vdw
    TYPE(neighbor_list_type), POINTER        :: sab_vdw_neighbor_list, &
                                                sab_vdw_neighbor_list_local
    TYPE(neighbor_node_type), POINTER        :: sab_vdw_neighbor_node
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_atom_dispersion_type), POINTER   :: disp_a, disp_b
    TYPE(qs_dispersion_type), POINTER        :: dispersion_env
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(virial_type), POINTER               :: virial

    CALL timeset(routineN,handle)

    NULLIFY (atomic_kind_set,sab_vdw)

    CALL get_qs_env(qs_env=qs_env, energy=energy, atomic_kind_set=atomic_kind_set,&
         virial=virial, dft_control=dft_control,para_env=para_env,&
         dispersion_env=dispersion_env,error=error)

    energy%dispersion = 0._dp

    IF ( dispersion_env%type == xc_vdw_fun_pairpot ) THEN

       s6 = dispersion_env%scaling
       dd = dispersion_env%exp_pre

       IF(calculate_forces) THEN
          NULLIFY (force,particle_set)
          CALL get_qs_env(qs_env=qs_env,&
               particle_set=particle_set,&
               force=force,error=error)
          natom = SIZE (particle_set)
          ALLOCATE (atom_of_kind(natom),STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
               atom_of_kind=atom_of_kind)
          use_virial = virial%pv_availability.AND.(.NOT.virial%pv_numer)
       END IF

       evdw = 0._dp

       CALL get_qs_env(qs_env=qs_env,sab_vdw=sab_vdw,error=error)

       nkind = SIZE(atomic_kind_set)

       DO ikind=1,nkind
          atomic_kind => atomic_kind_set(ikind)
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
               natom=natom, dispersion=disp_a)

          IF (.NOT.disp_a%defined) CYCLE
          IF ( disp_a%type == grimme_pp ) THEN
             c6a = disp_a%c6
             rc_a = disp_a%vdw_radii
          END IF

          DO jkind=1,nkind
             atomic_kind => atomic_kind_set(jkind)
             CALL get_atomic_kind(atomic_kind=atomic_kind, dispersion=disp_b)

             IF (.NOT.disp_a%defined) CYCLE

             CPPrecondition(disp_a%type==disp_b%type,cp_failure_level,routineP,error,failure)

             IF ( disp_b%type == grimme_pp ) THEN
                c6b = disp_b%c6
                rc_b = disp_b%vdw_radii

                rc = rc_a + rc_b
                c6 = SQRT(c6a*c6b)
             END IF

             iab = ikind + nkind*(jkind - 1)
             IF (.NOT.ASSOCIATED(sab_vdw(iab)%neighbor_list_set)) CYCLE
             CALL get_neighbor_list_set(neighbor_list_set=&
                  sab_vdw(iab)%neighbor_list_set,&
                  nlist=nlist)
             NULLIFY(sab_vdw_neighbor_list)

             DO ilist=1,nlist
                IF ( .NOT. ASSOCIATED(sab_vdw_neighbor_list) ) THEN
                   sab_vdw_neighbor_list => &
                        first_list(sab_vdw(iab)%neighbor_list_set)
                ELSE
                   sab_vdw_neighbor_list => next(sab_vdw_neighbor_list)
                END IF
                sab_vdw_neighbor_list_local => sab_vdw_neighbor_list
                CALL get_neighbor_list(neighbor_list=sab_vdw_neighbor_list_local,&
                     atom=iatom,nnode=nnode)

                sab_vdw_neighbor_node => first_node(sab_vdw_neighbor_list_local)

                DO inode=1,nnode
                   CALL get_neighbor_node(neighbor_node=sab_vdw_neighbor_node,&
                        neighbor=jatom,r=rij)

                   ! vdW potential
                   dr = SQRT(SUM(rij(:)**2))
                   IF (dr <= dispersion_env%rc_disp) THEN
                      IF ( disp_a%type == grimme_pp .AND. dr > 0.001_dp ) THEN
                         er = EXP(-dd*(dr/rc-1._dp))
                         fdmp = 1._dp/(1._dp+er)
                         xp = s6*c6/dr**6
                         evdw = evdw - xp*fdmp
                         IF(calculate_forces .AND. (dr > 0.001_dp)) THEN
                            dfdmp = dd/rc * er * fdmp*fdmp
                            devdw = -xp*(-6._dp*fdmp/dr + dfdmp)
                         END IF
                      END IF
                      IF(calculate_forces .AND. (dr > 0.001_dp)) THEN
                         fdij(:) = devdw * rij(:)/dr
                         atom_a = atom_of_kind(iatom)
                         atom_b = atom_of_kind(jatom)
                         force(ikind)%dispersion(:,atom_a) =&
                              force(ikind)%dispersion(:,atom_a) - fdij(:)
                         force(jkind)%dispersion(:,atom_b) =&
                              force(jkind)%dispersion(:,atom_b) + fdij(:)
                         IF ( use_virial ) THEN
                            CALL virial_pair_force ( virial%pv_virial, -1._dp, fdij, rij, error)
                         END IF
                      END IF
                   END IF

                   sab_vdw_neighbor_node => next(sab_vdw_neighbor_node)

                END DO ! inode => jatom(atom B)

             END DO ! ilist => iatom(atom A)

          END DO ! jkind

       END DO ! ikind

       IF(calculate_forces) THEN
          DEALLOCATE(atom_of_kind,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       END IF

       ! set dispersion energy
       CALL mp_sum(evdw,para_env%group)
       energy%dispersion = evdw

    END IF

    CALL timestop(handle)

  END SUBROUTINE calculate_dispersion_pairpot

! *****************************************************************************
  SUBROUTINE grimme_param(z,c6,r,found,error)

    INTEGER, INTENT(in)                      :: z
    REAL(KIND=dp), INTENT(inout)             :: c6, r
    LOGICAL, INTENT(inout)                   :: found
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'grimme_param', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), DIMENSION(54), PARAMETER :: c6val = (/0.14_dp, 0.08_dp, &
      1.61_dp, 1.61_dp, 3.13_dp, 1.75_dp, 1.23_dp, 0.70_dp, 0.75_dp, 0.63_dp, &
      5.71_dp, 5.71_dp, 10.79_dp, 9.23_dp, 7.84_dp, 5.57_dp, 5.07_dp, 4.61_dp,&
      10.80_dp, 10.80_dp, 10.80_dp, 10.80_dp, 10.80_dp, 10.80_dp, 10.80_dp, &
      10.80_dp, 10.80_dp, 10.80_dp, 10.80_dp, 10.80_dp, 16.99_dp, 17.10_dp, &
      16.37_dp, 12.64_dp, 12.47_dp, 12.01_dp, 24.67_dp, 24.67_dp, 24.67_dp, &
      24.67_dp, 24.67_dp, 24.67_dp, 24.67_dp, 24.67_dp, 24.67_dp, 24.67_dp, &
      24.67_dp, 24.67_dp, 37.32_dp, 38.71_dp, 38.44_dp, 31.74_dp, 31.50_dp, &
      29.99_dp /), rval = (/1.001_dp, 1.012_dp, 0.825_dp, 1.408_dp, 1.485_dp, &
      1.452_dp, 1.397_dp, 1.342_dp, 1.287_dp, 1.243_dp, 1.144_dp, 1.364_dp, &
      1.639_dp, 1.716_dp, 1.705_dp, 1.683_dp, 1.639_dp, 1.595_dp, 1.485_dp, &
      1.474_dp, 1.562_dp, 1.562_dp, 1.562_dp, 1.562_dp, 1.562_dp, 1.562_dp, &
      1.562_dp, 1.562_dp, 1.562_dp, 1.562_dp, 1.650_dp, 1.727_dp, 1.760_dp, &
      1.771_dp, 1.749_dp, 1.727_dp, 1.628_dp, 1.606_dp, 1.639_dp, 1.639_dp, &
      1.639_dp, 1.639_dp, 1.639_dp, 1.639_dp, 1.639_dp, 1.639_dp, 1.639_dp, &
      1.639_dp, 1.672_dp, 1.804_dp, 1.881_dp, 1.892_dp, 1.892_dp, 1.881_dp /)

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

  END SUBROUTINE grimme_param
! *****************************************************************************

END MODULE qs_dispersion_pairpot

