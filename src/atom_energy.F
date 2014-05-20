!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
MODULE atom_energy
  USE atom_electronic_structure,       ONLY: calculate_atom
  USE atom_fit,                        ONLY: atom_fit_density,&
                                             atom_fit_kgpot
  USE atom_operators,                  ONLY: atom_int_release,&
                                             atom_int_setup,&
                                             atom_ppint_release,&
                                             atom_ppint_setup,&
                                             atom_relint_release,&
                                             atom_relint_setup
  USE atom_output,                     ONLY: atom_print_basis,&
                                             atom_print_info,&
                                             atom_print_method,&
                                             atom_print_orbitals,&
                                             atom_print_potential,&
                                             atom_write_pseudo_param
  USE atom_types,                      ONLY: &
       GTH_PSEUDO, atom_basis_type, atom_gthpot_type, atom_integrals, &
       atom_optimization_type, atom_orbitals, atom_p_type, &
       atom_potential_type, atom_state, atom_type, create_atom_orbs, &
       create_atom_type, init_atom_basis, init_atom_potential, &
       read_atom_opt_section, release_atom_basis, release_atom_potential, &
       release_atom_type, set_atom
! ZMP added atom_read_external_density atom_read_external_vxc to atom_utiles
  USE atom_utils,                      ONLY: &
       atom_consistent_method, atom_core_density, atom_density, &
       atom_local_potential, atom_read_external_density, &
       atom_read_external_vxc, atom_set_occupation, atom_write_zmp_restart, &
       get_maxl_occ, get_maxn_occ
! ZMP added atom_xc module
  USE atom_xc,                         ONLY: calculate_atom_ext_vxc,&
                                             calculate_atom_zmp
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE f77_blas
  USE input_constants,                 ONLY: do_analytic
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE lapack
  USE mathconstants,                   ONLY: dfac,&
                                             gamma1,&
                                             pi
  USE periodic_table,                  ONLY: nelem,&
                                             ptable
  USE physcon,                         ONLY: evolt
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE
  PUBLIC  :: atom_energy_opt

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'atom_energy'

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE atom_energy_opt(atom_section,error)
    TYPE(section_vals_type), POINTER         :: atom_section
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'atom_energy_opt', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=2)                         :: elem
    CHARACTER(LEN=default_string_length)     :: filename
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: tmpstringlist
    INTEGER :: do_eric, do_erie, handle, i, ierr, im, in, iw, k, l, maxl, mb, &
      method, mo, n, n_meth, n_rep, nder, num_gau, num_gto, num_pol, reltyp, &
      zcore, zval, zz
    INTEGER, DIMENSION(0:3)                  :: maxn
    INTEGER, DIMENSION(:), POINTER           :: cn
    LOGICAL                                  :: dm, do_zmp, doread, eri_c, &
                                                eri_e, failure, had_ae, &
                                                had_pp, pp_calc, read_vxc
    REAL(KIND=dp)                            :: delta, lambda
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: ext_density, ext_vxc
    REAL(KIND=dp), DIMENSION(0:3, 10)        :: pocc
    TYPE(atom_basis_type), POINTER           :: ae_basis, pp_basis
    TYPE(atom_integrals), POINTER            :: ae_int, pp_int
    TYPE(atom_optimization_type)             :: optimization
    TYPE(atom_orbitals), POINTER             :: orbitals
    TYPE(atom_p_type), DIMENSION(:, :), &
      POINTER                                :: atom_info
    TYPE(atom_potential_type), POINTER       :: ae_pot, p_pot
    TYPE(atom_state), POINTER                :: state
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER :: basis_section, external_vxc_section, &
      method_section, opt_section, potential_section, powell_section, &
      xc_section, zmp_restart_section, zmp_section

    failure = .FALSE.

    CALL timeset(routineN,handle)

    ! What atom do we calculate
    CALL section_vals_val_get(atom_section,"ATOMIC_NUMBER", i_val=zval, error=error)
    CALL section_vals_val_get(atom_section,"ELEMENT", c_val=elem, error=error)
    zz = 0
    DO i=1,nelem
      IF ( ptable(i)%symbol == elem ) THEN
        zz = i
        EXIT
      END IF
    END DO
    IF ( zz /= 1 ) zval = zz

    ! read and set up inofrmation on the basis sets
    ALLOCATE(ae_basis,pp_basis,STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
    basis_section => section_vals_get_subs_vals(atom_section,"AE_BASIS",error=error)
    NULLIFY(ae_basis%grid)
    CALL init_atom_basis(ae_basis,basis_section,zval,"AE",error)
    NULLIFY(pp_basis%grid)
    basis_section => section_vals_get_subs_vals(atom_section,"PP_BASIS",error=error)
    CALL init_atom_basis(pp_basis,basis_section,zval,"PP",error)

    ! print general and basis set information
    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,atom_section,"PRINT%PROGRAM_BANNER",extension=".log",error=error)
    IF(iw > 0) CALL atom_print_info(zval,"Atomic Energy Calculation",iw,error)
    CALL cp_print_key_finished_output(iw,logger,atom_section,"PRINT%PROGRAM_BANNER",error=error)
    iw = cp_print_key_unit_nr(logger,atom_section,"PRINT%BASIS_SET",extension=".log",error=error)
    IF(iw > 0) THEN
      CALL atom_print_basis(ae_basis,iw," All Electron Basis",error)
      CALL atom_print_basis(pp_basis,iw," Pseudopotential Basis",error)
    END IF
    CALL cp_print_key_finished_output(iw,logger,atom_section,"PRINT%BASIS_SET",error=error)

    ! read and setup information on the pseudopotential
    NULLIFY(potential_section)
    potential_section => section_vals_get_subs_vals(atom_section,"POTENTIAL",error=error)
    ALLOCATE(ae_pot,p_pot,STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
    CALL init_atom_potential(p_pot,potential_section,zval,error)
    CALL init_atom_potential(ae_pot,potential_section,-1,error)

    ! if the ERI's are calculated analytically, we have to precalculate them
    eri_c = .FALSE.
    CALL section_vals_val_get(atom_section,"COULOMB_INTEGRALS", i_val=do_eric, error=error)
    IF(do_eric==do_analytic) eri_c = .TRUE.
    eri_e = .FALSE.
    CALL section_vals_val_get(atom_section,"EXCHANGE_INTEGRALS", i_val=do_erie, error=error)
    IF(do_erie==do_analytic) eri_e = .TRUE.

    ! information on the states to be calculated
    CALL section_vals_val_get(atom_section,"MAX_ANGULAR_MOMENTUM", i_val=maxl, error=error)
    maxn=0
    CALL section_vals_val_get(atom_section,"CALCULATE_STATES", i_vals=cn, error=error)
    DO in = 1, MIN(SIZE(cn),4)
      maxn(in-1) = cn(in)
    END DO
    DO in = 0, 3
      maxn(in) = MIN(maxn(in),ae_basis%nbas(in))
      maxn(in) = MIN(maxn(in),pp_basis%nbas(in))
    END DO

    ! read optimization section
    opt_section => section_vals_get_subs_vals(atom_section,"OPTIMIZATION",error=error)
    CALL read_atom_opt_section(optimization,opt_section,error)

    had_ae = .FALSE.
    had_pp = .FALSE.

    ! Check for the total number of electron configurations to be calculated
    CALL section_vals_val_get(atom_section,"ELECTRON_CONFIGURATION", n_rep_val=n_rep, error=error)
    ! Check for the total number of method types to be calculated
    method_section => section_vals_get_subs_vals(atom_section,"METHOD",error=error)
    CALL section_vals_get(method_section,n_repetition=n_meth,error=error)

    ! integrals
    ALLOCATE(ae_int, pp_int,STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

    ALLOCATE(atom_info(n_rep,n_meth),STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

    DO in = 1, n_rep
      DO im = 1, n_meth

        NULLIFY(atom_info(in,im)%atom)
        CALL create_atom_type(atom_info(in,im)%atom,error)

        atom_info(in,im)%atom%optimization = optimization

        atom_info(in,im)%atom%z = zval
        xc_section => section_vals_get_subs_vals(method_section,"XC",i_rep_section=im,error=error)
        atom_info(in,im)%atom%xc_section => xc_section

        ! ZMP Reading input sections if they are found initialize everything
        do_zmp=.FALSE.
        doread=.FALSE.
        read_vxc=.FALSE.

        zmp_section => section_vals_get_subs_vals(method_section,"ZMP",error=error)
        CALL section_vals_get(zmp_section,explicit=do_zmp,error=error)
        atom_info(in,im)%atom%do_zmp=do_zmp
        CALL section_vals_val_get(zmp_section,"FILE_DENSITY",c_val=filename,error=error)
        atom_info(in,im)%atom%ext_file=filename
        CALL section_vals_val_get(zmp_section,"GRID_TOL",&
             r_val= atom_info(in,im)%atom%zmpgrid_tol,error=error)
        CALL section_vals_val_get(zmp_section,"LAMBDA",r_val=lambda,error=error)
        atom_info(in,im)%atom%lambda=lambda
        CALL section_vals_val_get(zmp_section,"DM",l_val=dm,error=error)
        atom_info(in,im)%atom%dm=dm

        zmp_restart_section => section_vals_get_subs_vals(zmp_section,"RESTART",error=error)
        CALL section_vals_get(zmp_restart_section,explicit=doread,error=error)
        atom_info(in,im)%atom%doread=doread
        CALL section_vals_val_get(zmp_restart_section,"FILE_RESTART",c_val=filename,error=error)
        atom_info(in,im)%atom%zmp_restart_file=filename

        ! ZMP Reading external vxc section, if found initialize
        external_vxc_section => section_vals_get_subs_vals(method_section,"EXTERNAL_VXC",error=error)
        CALL section_vals_get(external_vxc_section,explicit=read_vxc,error=error)
        atom_info(in,im)%atom%read_vxc=read_vxc
        CALL section_vals_val_get(external_vxc_section,"FILE_VXC",c_val=filename,error=error)
        atom_info(in,im)%atom%ext_vxc_file=filename
        CALL section_vals_val_get(external_vxc_section,"GRID_TOL",&
             r_val= atom_info(in,im)%atom%zmpvxcgrid_tol,error=error)

        ALLOCATE(state,STAT=ierr)
        CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

        ! get the electronic configuration
        CALL section_vals_val_get(atom_section,"ELECTRON_CONFIGURATION", i_rep_val=in,&
                                  c_vals=tmpstringlist, error=error)

        ! set occupations
        CALL atom_set_occupation(tmpstringlist,state%occ,state%occupation,state%multiplicity,error)
        state%maxl_occ = get_maxl_occ(state%occ)
        state%maxn_occ = get_maxn_occ(state%occ)

        ! set number of states to be calculated
        state%maxl_calc = MAX(maxl,state%maxl_occ)
        state%maxl_calc = MIN(3,state%maxl_calc)
        state%maxn_calc = 0
        DO k=0,state%maxl_calc
          state%maxn_calc(k) = MAX(maxn(k),state%maxn_occ(k))
        END DO

        ! is there a pseudo potential
        pp_calc = ANY(INDEX(tmpstringlist(1:),"CORE") /= 0)
        IF ( pp_calc ) THEN
          ! get and set the core occupations
          CALL section_vals_val_get(atom_section,"CORE", c_vals=tmpstringlist, error=error)
          CALL atom_set_occupation(tmpstringlist,state%core,pocc,error=error)
          zcore = zval - SUM(state%core)
          CALL set_atom(atom_info(in,im)%atom,zcore=zcore,pp_calc=.TRUE.,error=error)
        ELSE
          state%core=0._dp
          CALL set_atom(atom_info(in,im)%atom,zcore=zval,pp_calc=.FALSE.,error=error)
        END IF

        CALL section_vals_val_get(method_section,"METHOD_TYPE",i_val=method,i_rep_section=im,error=error)
        CALL section_vals_val_get(method_section,"RELATIVISTIC",i_val=reltyp,i_rep_section=im,error=error)
        CALL set_atom(atom_info(in,im)%atom,method_type=method,relativistic=reltyp,error=error)

        IF(atom_consistent_method(method,state%multiplicity)) THEN
          iw = cp_print_key_unit_nr(logger,atom_section,"PRINT%METHOD_INFO",extension=".log",error=error)
          CALL atom_print_method(atom_info(in,im)%atom,iw,error)
          CALL cp_print_key_finished_output(iw,logger,atom_section,"PRINT%METHOD_INFO",error=error)

          iw = cp_print_key_unit_nr(logger,atom_section,"PRINT%POTENTIAL",extension=".log",error=error)
          IF ( pp_calc ) THEN
            IF(iw > 0) CALL atom_print_potential(p_pot,iw,error)
          ELSE
            IF(iw > 0) CALL atom_print_potential(ae_pot,iw,error)
          END IF
          CALL cp_print_key_finished_output(iw,logger,atom_section,"PRINT%POTENTIAL",error=error)
        END IF

        ! calculate integrals
        IF ( pp_calc ) THEN
          ! general integrals
          CALL atom_int_setup(pp_int,pp_basis,&
                              potential=p_pot,eri_coulomb=eri_c,eri_exchange=eri_e,error=error)
          ! potential
          CALL atom_ppint_setup(pp_int,pp_basis,potential=p_pot,error=error)
          !
          NULLIFY(pp_int%tzora,pp_int%hdkh)
          !
          CALL set_atom(atom_info(in,im)%atom,basis=pp_basis,integrals=pp_int,potential=p_pot,error=error)
          state%maxn_calc(:) = MIN( state%maxn_calc(:), pp_basis%nbas(:) )
          CPPostcondition(ALL(state%maxn_calc(:) >= state%maxn_occ), cp_failure_level, routineP, error, failure)
          had_pp = .TRUE.
        ELSE
          ! general integrals
          CALL atom_int_setup(ae_int,ae_basis,potential=ae_pot,&
                              eri_coulomb=eri_c,eri_exchange=eri_e,error=error)
          ! potential
          CALL atom_ppint_setup(ae_int,ae_basis,potential=ae_pot,error=error)
          ! relativistic correction terms
          CALL atom_relint_setup(ae_int,ae_basis,reltyp,zcore=REAL(zval,dp),error=error)
          !
          CALL set_atom(atom_info(in,im)%atom,basis=ae_basis,integrals=ae_int,potential=ae_pot,error=error)
          state%maxn_calc(:) = MIN( state%maxn_calc(:), ae_basis%nbas(:) )
          CPPostcondition(ALL(state%maxn_calc(:) >= state%maxn_occ), cp_failure_level, routineP, error, failure)
          had_ae = .TRUE.
        END IF

        CALL set_atom(atom_info(in,im)%atom,state=state,error=error)

        CALL set_atom(atom_info(in,im)%atom,coulomb_integral_type=do_eric,&
                      exchange_integral_type=do_erie,error=error)

        NULLIFY(orbitals)
        mo = MAXVAL(state%maxn_calc)
        mb = MAXVAL(atom_info(in,im)%atom%basis%nbas)
        CALL create_atom_orbs(orbitals,mb,mo,error)
        CALL set_atom(atom_info(in,im)%atom,orbitals=orbitals,error=error)

        IF(atom_consistent_method(method,state%multiplicity)) THEN
          !Calculate the electronic structure
          iw = cp_print_key_unit_nr(logger,atom_section,"PRINT%SCF_INFO",extension=".log",error=error)
          CALL calculate_atom(atom_info(in,im)%atom,iw,error=error)
          CALL cp_print_key_finished_output(iw,logger,atom_section,"PRINT%SCF_INFO",error=error)

          ! ZMP If we have the external density do zmp
          IF (atom_info(in,im)%atom%do_zmp) THEN
            CALL atom_write_zmp_restart(atom_info(in,im)%atom) 

            ALLOCATE(ext_density(atom_info(in,im)%atom%basis%grid%nr))
            ext_density=0._dp
            CALL atom_read_external_density(ext_density, atom_info(in,im)%atom,iw)
            CALL calculate_atom_zmp(ext_density=ext_density, &
                                    atom=atom_info(in,im)%atom, &
                                    lprint=.TRUE.,error=error)
            DEALLOCATE(ext_density)
          ENDIF
          ! ZMP If we have the external v_xc calculate KS quantities
          IF (atom_info(in,im)%atom%read_vxc) THEN
            ALLOCATE(ext_vxc(atom_info(in,im)%atom%basis%grid%nr))
            ext_vxc=0._dp
            CALL atom_read_external_vxc(ext_vxc, atom_info(in,im)%atom,iw)
            CALL calculate_atom_ext_vxc(vxc=ext_vxc, &
                                        atom=atom_info(in,im)%atom, &
                                        lprint=.TRUE.,error=error)
            DEALLOCATE(ext_vxc)
          ENDIF

          ! rk: Print out the energies needed for the scptb parametrization
          ! (Orbital energies minus Kinetic, etc)
          iw = cp_print_key_unit_nr(logger,atom_section,"PRINT%ENERGIES_MINUS_KINETIC",extension=".log",error=error)
          IF (iw > 0) THEN
                WRITE(iw,'(T36,A,T61,F20.12)') "     Etot-Ekin   ::",atom_info(in,im)%atom%energy%etot &
                         - atom_info(in,im)%atom%energy%ekin

                ! Output of Epsilon-kin, used in SCPTB parameter generation.
                IF ( atom_info(in,im)%atom%state%multiplicity == -1 ) THEN
                   WRITE(iw,'(/,A,T20,A,T30,A,T36,A,T49,A,T71,A,/)') " Orbital Epsilon-T",&
                        "State","L","Occupation","Energy[a.u.]","Energy[eV]"
                   DO l=0,atom_info(in,im)%atom%state%maxl_calc
                     n = atom_info(in,im)%atom%state%maxn_calc(l)
                     DO i=1,n 
                       WRITE(iw,'(T23,I2,T30,I1,T36,F10.3,T46,F15.6,T66,F15.6)') &
                          i,l,atom_info(in,im)%atom%state%occupation(l,i), &
                          atom_info(in,im)%atom%orbitals%ener(i,l)-atom_info(in,im)%atom%integrals%kin(i,i,l),&
                          (atom_info(in,im)%atom%orbitals%ener(i,l)-atom_info(in,im)%atom%integrals%kin(i,i,l))*evolt
                     END DO
                     IF (n>0) WRITE(iw,*)
                   END DO
                ELSE
                   WRITE(iw,'(/,A,T20,A,T30,A,T36,A,T42,A,T55,A,T71,A,/)') " Orbital Epsilon-T",&
                        "State","Spin","L","Occupation","Energy[a.u.]","Energy[eV]"
                   DO l=0,atom_info(in,im)%atom%state%maxl_calc
                     n = atom_info(in,im)%atom%state%maxn_calc(l)
                     DO i=1,n 
                       WRITE(iw,'(T23,I2,T29,A,T36,I1,T42,F10.3,T52,F15.6,T68,F13.6)') &
                          i,"alpha",l,atom_info(in,im)%atom%state%occa(l,i),&
                          atom_info(in,im)%atom%orbitals%enera(i,l)-atom_info(in,im)%atom%integrals%kin(i,i,l),&
                          (atom_info(in,im)%atom%orbitals%ener(i,l)-atom_info(in,im)%atom%integrals%kin(i,i,l))*evolt
                     END DO
                     DO i=1,n 
                       WRITE(iw,'(T23,I2,T29,A,T36,I1,T42,F10.3,T52,F15.6,T68,F13.6)') &
                          i," beta",l,atom_info(in,im)%atom%state%occb(l,i),&
                          atom_info(in,im)%atom%orbitals%enerb(i,l)-atom_info(in,im)%atom%integrals%kin(i,i,l),&
                          (atom_info(in,im)%atom%orbitals%ener(i,l)-atom_info(in,im)%atom%integrals%kin(i,i,l))*evolt
                     END DO
                     IF (n>0) WRITE(iw,*)
                   END DO
                END IF
          END IF


          ! Print out the orbitals if requested
          iw = cp_print_key_unit_nr(logger,atom_section,"PRINT%ORBITALS",extension=".log",error=error)
          IF (iw > 0) THEN
            CALL atom_print_orbitals(atom_info(in,im)%atom,iw,error=error)
          END IF
          CALL cp_print_key_finished_output(iw,logger,atom_section,"PRINT%ORBITALS",error=error)

          ! perform a fit of the total electronic density
          iw = cp_print_key_unit_nr(logger,atom_section,"PRINT%FIT_DENSITY",extension=".log",error=error)
          IF (iw>0) THEN
            CALL section_vals_val_get(atom_section,"PRINT%FIT_DENSITY%NUM_GTO", i_val=num_gto, error=error)
            powell_section => section_vals_get_subs_vals(atom_section,"POWELL",error=error)
            CALL atom_fit_density(atom_info(in,im)%atom,num_gto,0,iw,powell_section=powell_section,error=error)
          END IF
          CALL cp_print_key_finished_output(iw,logger,atom_section,"PRINT%FIT_DENSITY",error=error)

          ! Optimize a local potential for the non-additive kinetic energy term in KG
          iw = cp_print_key_unit_nr(logger,atom_section,"PRINT%FIT_KGPOT",extension=".log",error=error)
          IF (iw>0) THEN
            CALL section_vals_val_get(atom_section,"PRINT%FIT_KGPOT%NUM_GAUSSIAN", i_val=num_gau, error=error)
            CALL section_vals_val_get(atom_section,"PRINT%FIT_KGPOT%NUM_POLYNOM", i_val=num_pol, error=error)
            powell_section => section_vals_get_subs_vals(atom_section,"POWELL",error=error)
            CALL atom_fit_kgpot(atom_info(in,im)%atom,num_gau,num_pol,iw,powell_section=powell_section,error=error)
          END IF
          CALL cp_print_key_finished_output(iw,logger,atom_section,"PRINT%FIT_KGPOT",error=error)

          ! generate a response basis
          iw = cp_print_key_unit_nr(logger,atom_section,"PRINT%RESPONSE_BASIS",extension=".log",error=error)
          IF (iw>0) THEN
            CALL section_vals_val_get(atom_section,"PRINT%RESPONSE_BASIS%DELTA_CHARGE", r_val=delta, error=error)
            CALL section_vals_val_get(atom_section,"PRINT%RESPONSE_BASIS%DERIVATIVES", i_val=nder, error=error)
            CALL atom_response_basis(atom_info(in,im)%atom,delta,nder,iw,error)
          END IF
          CALL cp_print_key_finished_output(iw,logger,atom_section,"PRINT%RESPONSE_BASIS",error=error)

          ! generate a UPF file
         iw = cp_print_key_unit_nr(logger,atom_section,"PRINT%UPF_FILE",extension=".upf",&
              file_position="REWIND",error=error)
          IF (iw>0) THEN
            CALL atom_write_upf(atom_info(in,im)%atom,iw,error)
          END IF
          CALL cp_print_key_finished_output(iw,logger,atom_section,"PRINT%UPF_FILE",error=error)
        END IF

      END DO
    END DO

    ! clean up
    IF ( had_ae ) THEN
      CALL atom_int_release(ae_int,error)
      CALL atom_ppint_release(ae_int,error)
      CALL atom_relint_release(ae_int,error)
    END IF
    IF ( had_pp ) THEN
      CALL atom_int_release(pp_int,error)
      CALL atom_ppint_release(pp_int,error)
      CALL atom_relint_release(pp_int,error)
    END IF
    CALL release_atom_basis(ae_basis,error)
    CALL release_atom_basis(pp_basis,error)

    CALL release_atom_potential(p_pot,error)
    CALL release_atom_potential(ae_pot,error)

    DO in = 1, n_rep
      DO im = 1, n_meth
        CALL release_atom_type(atom_info(in,im)%atom,error)
      END DO
    END DO
    DEALLOCATE(atom_info,STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

    DEALLOCATE(ae_pot,p_pot,ae_basis,pp_basis,ae_int,pp_int,STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

    CALL timestop(handle)

  END SUBROUTINE atom_energy_opt

! *****************************************************************************
  SUBROUTINE atom_response_basis(atom,delta,nder,iw,error)

    TYPE(atom_type), POINTER                 :: atom
    REAL(KIND=dp), INTENT(IN)                :: delta
    INTEGER, INTENT(IN)                      :: nder, iw
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'atom_response_basis', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ider, ierr, j, k, l, &
                                                lhomo, m, n, nhomo, s1, s2
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: dene, emax, expzet, fhomo, o, &
                                                prefac, zeta
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: amat
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: rbasis
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :)                  :: wfn
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: ovlp
    TYPE(atom_state), POINTER                :: state

    failure = .FALSE.
    WRITE(iw,'(/," ",79("*"),/,T34,A,/," ",79("*"))') "RESPONSE BASIS"

    state => atom%state
    ovlp => atom%integrals%ovlp

    ! find HOMO
    lhomo = -1
    nhomo = -1
    emax=-HUGE(1._dp)
    DO l=0,state%maxl_occ
      DO i=1,state%maxn_occ(l)
        IF (atom%orbitals%ener(i,l) > emax) THEN
          lhomo=l
          nhomo=i
          emax=atom%orbitals%ener(i,l)
          fhomo=state%occupation(l,i)
        END IF
      END DO
    END DO

    s1=SIZE(atom%orbitals%wfn,1)
    s2=SIZE(atom%orbitals%wfn,2)
    ALLOCATE(wfn(s1,s2,0:3,-nder:nder),STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
    s2=MAXVAL(state%maxn_occ)+nder
    ALLOCATE(rbasis(s1,s2,0:3),STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
    rbasis = 0._dp

    DO ider=-nder,nder
      dene = REAL(ider,KIND=dp)*delta
      CPPostcondition(fhomo>ABS(dene), cp_failure_level, routineP, error, failure)
      state%occupation(lhomo,nhomo)=fhomo+dene
      CALL calculate_atom(atom,iw=0,noguess=.TRUE.,error=error)
      wfn(:,:,:,ider) = atom%orbitals%wfn
      state%occupation(lhomo,nhomo)=fhomo
    END DO

    DO l=0,state%maxl_occ
      ! occupied states
      DO i=1,MAX(state%maxn_occ(l),1)
        rbasis(:,i,l) = wfn(:,i,l,0)
      END DO
      ! differentiation
      DO ider=1,nder
        i=MAX(state%maxn_occ(l),1)
        SELECT CASE (ider)
          CASE (1)
            rbasis(:,i+1,l) = 0.5_dp*(wfn(:,i,l,1) - wfn(:,i,l,-1))/delta
          CASE (2)
            rbasis(:,i+2,l) = 0.25_dp*(wfn(:,i,l,2) - 2._dp*wfn(:,i,l,0) + wfn(:,i,l,-2))/delta**2
          CASE (3)
            rbasis(:,i+3,l) = 0.125_dp*(wfn(:,i,l,3) - 3._dp*wfn(:,i,l,1) &
                                      + 3._dp*wfn(:,i,l,-1) - wfn(:,i,l,-3))/delta**3
          CASE DEFAULT
            CPPostcondition(.FALSE., cp_failure_level, routineP, error, failure)
        END SELECT
      END DO

      ! orthogonalization, use gram-schmidt in order to keep the natural order (semi-core, valence, response) of the wfn.
      n = state%maxn_occ(l)+nder
      m = atom%basis%nbas(l)
      DO i=1,n
         DO j=1,i-1
            o=DOT_PRODUCT(rbasis(1:m,j,l),RESHAPE(MATMUL(ovlp(1:m,1:m,l),rbasis(1:m,i:i,l)),(/m/)))
            rbasis(1:m,i,l)=rbasis(1:m,i,l)-o*rbasis(1:m,j,l)
         ENDDO
         o=DOT_PRODUCT(rbasis(1:m,i,l),RESHAPE(MATMUL(ovlp(1:m,1:m,l),rbasis(1:m,i:i,l)),(/m/)))
         rbasis(1:m,i,l)=rbasis(1:m,i,l)/SQRT(o)
      ENDDO

      ! check
      ALLOCATE(amat(n,n),STAT=ierr)
      CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)
      amat(1:n,1:n) = MATMUL(TRANSPOSE(rbasis(1:m,1:n,l)),MATMUL(ovlp(1:m,1:m,l),rbasis(1:m,1:n,l)))
      DO i=1,n
        amat(i,i)=amat(i,i) - 1._dp
      END DO
      IF (MAXVAL(ABS(amat)) > 1.e-12) THEN
         WRITE(iw,'(A,G20.10)') " Orthogonality error  ", MAXVAL(ABS(amat))
      END IF
      DEALLOCATE(amat,STAT=ierr)
      CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)

      ! Quickstep normalization
      WRITE(iw,'(/,A,T30,I3)') " Angular momentum :",l

      WRITE(iw,'(/,A,I0,A,I0,A)') "             Exponent     Coef.(Quickstep Normalization), first ",&
                                  n-nder," valence ",nder," response"
      expzet = 0.25_dp*REAL(2*l + 3,dp)
      prefac = SQRT( SQRT(pi)/2._dp**(l+2)*dfac(2*l+1) )
      DO i=1,m
        zeta = (2._dp*atom%basis%am(i,l))**expzet
        WRITE(iw,'(4X,F20.10,4X,15ES20.6)') atom%basis%am(i,l),((prefac*rbasis(i,k,l)/zeta),k=1,n)
      END DO

    END DO

    DEALLOCATE(wfn,rbasis,STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

    WRITE(iw,'(" ",79("*"))')

  END SUBROUTINE atom_response_basis
! *****************************************************************************
  SUBROUTINE atom_write_upf(atom,iw,error)

    TYPE(atom_type), POINTER                 :: atom
    INTEGER, INTENT(IN)                      :: iw
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'atom_write_upf', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: string
    INTEGER                                  :: i, ibeta, ierr, j, k, l, &
                                                lmax, nbeta, nr, nwfc, nwfn
    LOGICAL                                  :: failure, up
    REAL(KIND=dp)                            :: pf, rl, rmax
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: beta, corden, dens, ef, &
                                                locpot, rp
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: dij
    TYPE(atom_gthpot_type), POINTER          :: pot

    failure = .FALSE.
    IF(.NOT.atom%pp_calc) RETURN
    IF(atom%potential%ppot_type /= GTH_PSEUDO) RETURN
    pot => atom%potential%gth_pot
    CPPostcondition(.NOT.pot%lsdpot, cp_failure_level, routineP, error, failure)

    WRITE(iw,'(A)') '<UPF version="2.0.1">'

    WRITE(iw,'(T4,A)') '<PP_INFO>'
    WRITE(iw,'(T8,A)') 'Converted from CP2K GTH format'
    WRITE(iw,'(T8,A)') '<PP_INPUTFILE>'
    CALL atom_write_pseudo_param(pot,iw,error)
    WRITE(iw,'(T8,A)') '</PP_INPUTFILE>'
    WRITE(iw,'(T4,A)') '</PP_INFO>'
    WRITE(iw,'(T4,A)') '<PP_HEADER'
    WRITE(iw,'(T8,A)') 'generated="Generated in analytical, separable form"'
    WRITE(iw,'(T8,A)') 'author="Goedecker/Hartwigsen/Hutter/Teter"'
    WRITE(iw,'(T8,A)') 'date="Phys.Rev.B58, 3641 (1998); B54, 1703 (1996)"'
    WRITE(iw,'(T8,A)') 'comment="This file generated by CP2K ATOM code"'
    CALL compose(string,"element",cval=pot%symbol)
    WRITE(iw,'(T8,A)') TRIM(string)
    WRITE(iw,'(T8,A)') 'pseudo_type="NC"'
    WRITE(iw,'(T8,A)') 'relativistic="no"'
    WRITE(iw,'(T8,A)') 'is_ultrasoft="F"'
    WRITE(iw,'(T8,A)') 'is_paw="F"'
    WRITE(iw,'(T8,A)') 'is_coulomb="F"'
    WRITE(iw,'(T8,A)') 'has_so="F"'
    WRITE(iw,'(T8,A)') 'has_wfc="F"'
    WRITE(iw,'(T8,A)') 'has_gipaw="F"'
    WRITE(iw,'(T8,A)') 'paw_as_gipaw="F"'
    IF(pot%nlcc) THEN
       WRITE(iw,'(T8,A)') 'core_correction="T"'
    ELSE
       WRITE(iw,'(T8,A)') 'core_correction="F"'
    END IF
    WRITE(iw,'(T8,A)') 'functional="DFT"'
    CALL compose(string,"z_valence",rval=pot%zion)
    WRITE(iw,'(T8,A)') TRIM(string)
    CALL compose(string,"total_psenergy",rval=2._dp*atom%energy%etot)
    WRITE(iw,'(T8,A)') TRIM(string)
    WRITE(iw,'(T8,A)') 'wfc_cutoff="0.0E+00"'
    WRITE(iw,'(T8,A)') 'rho_cutoff="0.0E+00"'
    lmax = -1
    DO l=0,3
      IF(pot%nl(l) > 0) lmax = l
    END DO
    CALL compose(string,"l_max",ival=lmax)
    WRITE(iw,'(T8,A)') TRIM(string)
    WRITE(iw,'(T8,A)') 'l_max_rho="0"'
    WRITE(iw,'(T8,A)') 'l_local="-3"'
    nr = atom%basis%grid%nr
    CALL compose(string,"mesh_size",ival=nr)
    WRITE(iw,'(T8,A)') TRIM(string)
    nwfc = SUM(atom%state%maxn_occ)
    CALL compose(string,"number_of_wfc",ival=nwfc)
    WRITE(iw,'(T8,A)') TRIM(string)
    nbeta = SUM(pot%nl)
    CALL compose(string,"number_of_proj",ival=nbeta)
    WRITE(iw,'(T8,A)') TRIM(string)//'/>'

    ! Mesh
    up = atom%basis%grid%rad(1) < atom%basis%grid%rad(nr)
    WRITE(iw,'(T4,A)') '<PP_MESH'
    WRITE(iw,'(T8,A)') 'dx="1.E+00"'
    CALL compose(string,"mesh",ival=nr)
    WRITE(iw,'(T8,A)') TRIM(string)
    WRITE(iw,'(T8,A)') 'xmin="1.E+00"'
    rmax = MAXVAL(atom%basis%grid%rad)
    CALL compose(string,"rmax",rval=rmax)
    WRITE(iw,'(T8,A)') TRIM(string)
    WRITE(iw,'(T8,A)') 'zmesh="1.E+00">'
    WRITE(iw,'(T8,A)') '<PP_R type="real"'
    CALL compose(string,"size",ival=nr)
    WRITE(iw,'(T12,A)') TRIM(string)
    WRITE(iw,'(T12,A)') 'columns="4">'
    IF(up) THEN
       WRITE(iw,'(T12,4ES25.12E3)') (atom%basis%grid%rad(i),i=1,nr)
    ELSE
       WRITE(iw,'(T12,4ES25.12E3)') (atom%basis%grid%rad(i),i=nr,1,-1)
    END IF
    WRITE(iw,'(T8,A)') '</PP_R>'
    WRITE(iw,'(T8,A)') '<PP_RAB type="real"'
    CALL compose(string,"size",ival=nr)
    WRITE(iw,'(T12,A)') TRIM(string)
    WRITE(iw,'(T12,A)') 'columns="4">'
    IF(up) THEN
       WRITE(iw,'(T12,4ES25.12E3)') (atom%basis%grid%wr(i)/atom%basis%grid%rad2(i),i=1,nr)
    ELSE
       WRITE(iw,'(T12,4ES25.12E3)') (atom%basis%grid%wr(i)/atom%basis%grid%rad2(i),i=nr,1,-1)
    END IF
    WRITE(iw,'(T8,A)') '</PP_RAB>'
    WRITE(iw,'(T4,A)') '</PP_MESH>'

    ! NLCC    
    IF(pot%nlcc) THEN
       WRITE(iw,'(T4,A)') '<PP_NLCC type="real"'
       CALL compose(string,"size",ival=nr)
       WRITE(iw,'(T8,A)') TRIM(string)
       WRITE(iw,'(T8,A)') 'columns="4">'
       ALLOCATE(corden(nr),STAT=ierr)
       CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)
       CALL atom_core_density(corden,pot,"RHO",atom%basis%grid%rad,error)
       IF(up) THEN
          WRITE(iw,'(T8,4ES25.12E3)') (corden(i),i=1,nr)
       ELSE
          WRITE(iw,'(T8,4ES25.12E3)') (corden(i),i=nr,1,-1)
       END IF
       DEALLOCATE(corden,STAT=ierr)
       CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)
       WRITE(iw,'(T4,A)') '</PP_NLCC>'
    END IF

    ! local PP
    WRITE(iw,'(T4,A)') '<PP_LOCAL type="real"'
    CALL compose(string,"size",ival=nr)
    WRITE(iw,'(T8,A)') TRIM(string)
    WRITE(iw,'(T8,A)') 'columns="4">'
    ALLOCATE(locpot(nr),STAT=ierr)
    CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)
    CALL atom_local_potential(locpot,pot,atom%basis%grid%rad,error)
    IF(up) THEN
       WRITE(iw,'(T8,4ES25.12E3)') (2.0_dp*locpot(i),i=1,nr)
    ELSE
       WRITE(iw,'(T8,4ES25.12E3)') (2.0_dp*locpot(i),i=nr,1,-1)
    END IF
    DEALLOCATE(locpot,STAT=ierr)
    CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)
    WRITE(iw,'(T4,A)') '</PP_LOCAL>'

    ! nonlocal PP
    WRITE(iw,'(T4,A)') '<PP_NONLOCAL>'
    ALLOCATE(rp(nr),ef(nr),beta(nr),STAT=ierr)
    CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)
    ibeta = 0
    DO l=0,3
       IF(pot%nl(l) == 0) CYCLE
       rl = pot%rcnl(l)
       rp(:) = atom%basis%grid%rad
       ef(:) = EXP(-0.5_dp*rp*rp/(rl*rl))
       DO i=1,pot%nl(l)
          pf = rl**(l+0.5_dp*(4._dp*i-1._dp))
          j = l + 2*i - 1
          pf = SQRT(2._dp) / (pf * SQRT(gamma1(j)))
          beta(:) = pf * rp**(l+2*i-2) * ef
          ibeta=ibeta+1
          CALL compose(string,"<PP_BETA",counter=ibeta)
          WRITE(iw,'(T8,A)') TRIM(string)
          CALL compose(string,"angular_momentum",ival=l)
          WRITE(iw,'(T12,A)') TRIM(string)
          WRITE(iw,'(T12,A)') 'type="real"'
          CALL compose(string,"size",ival=nr)
          WRITE(iw,'(T12,A)') TRIM(string)
          WRITE(iw,'(T12,A)') 'columns="4">'
          beta(:) = beta*rp
          IF(up) THEN
             WRITE(iw,'(T12,4ES25.12E3)') (2._dp*beta(j),j=1,nr)
          ELSE
             WRITE(iw,'(T12,4ES25.12E3)') (2._dp*beta(j),j=nr,1,-1)
          END IF
          CALL compose(string,"</PP_BETA",counter=ibeta,isfinal=.TRUE.)
          WRITE(iw,'(T8,A)') TRIM(string)
       END DO
    END DO
    DEALLOCATE(rp,ef,beta,STAT=ierr)
    CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)
    ! nonlocal PP matrix elements
    ALLOCATE(dij(nbeta,nbeta),STAT=ierr)
    CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)
    dij = 0._dp
    DO l=0,3
       IF(pot%nl(l) == 0) CYCLE
       ibeta = SUM(pot%nl(0:l-1)) + 1
       i = ibeta + pot%nl(l) - 1
       dij(ibeta:i,ibeta:i) = pot%hnl(1:pot%nl(l),1:pot%nl(l),l)
    END DO
    WRITE(iw,'(T8,A)') '<PP_DIJ type="real"'
    WRITE(iw,'(T12,A)') 'columns="4">'
    WRITE(iw,'(T12,4ES25.12E3)') ((0.5_dp*dij(i,j),j=1,nbeta),i=1,nbeta)
    WRITE(iw,'(T8,A)') '</PP_DIJ>'
    DEALLOCATE(dij,STAT=ierr)
    CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)
    WRITE(iw,'(T4,A)') '</PP_NONLOCAL>'

    ! atomic wavefunctions
    ALLOCATE(beta(nr),STAT=ierr)
    CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)
    WRITE(iw,'(T4,A)') '<PP_PSWFC>'
    nwfn=0
    DO l=0,3
       DO i=1,10
          IF(ABS(atom%state%occupation(l,i)) == 0._dp) CYCLE
          nwfn = nwfn+1
          CALL compose(string,"<PP_CHI",counter=nwfn)
          WRITE(iw,'(T8,A)') TRIM(string)
          CALL compose(string,"l",ival=l)
          WRITE(iw,'(T12,A)') TRIM(string)
          CALL compose(string,"occupation",rval=atom%state%occupation(l,i))
          WRITE(iw,'(T12,A)') TRIM(string)
          CALL compose(string,"pseudo_energy",rval=2._dp*atom%orbitals%ener(i,l))
          WRITE(iw,'(T12,A)') TRIM(string)
          WRITE(iw,'(T12,A)') 'columns="4">'
          beta = 0._dp
          DO k=1,atom%basis%nbas(l)
             beta(:) = beta(:) + atom%orbitals%wfn(k,i,l)*atom%basis%bf(:,k,l)
          END DO
          beta(:) = beta*atom%basis%grid%rad
          IF(up) THEN
             WRITE(iw,'(T12,4ES25.12E3)') (beta(j)*atom%basis%grid%rad(j),j=1,nr)
          ELSE
             WRITE(iw,'(T12,4ES25.12E3)') (beta(j)*atom%basis%grid%rad(j),j=nr,1,-1)
          END IF
          CALL compose(string,"</PP_CHI",counter=nwfn,isfinal=.TRUE.)
          WRITE(iw,'(T8,A)') TRIM(string)
       END DO
    END DO
    WRITE(iw,'(T4,A)') '</PP_PSWFC>'
    DEALLOCATE(beta,STAT=ierr)
    CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)

    ! atomic charge
    ALLOCATE(dens(nr),STAT=ierr)
    CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)
    WRITE(iw,'(T4,A)') '<PP_RHOATOM type="real"'
    CALL compose(string,"size",ival=nr)
    WRITE(iw,'(T8,A)') TRIM(string)
    WRITE(iw,'(T8,A)') 'columns="4">'
    CALL atom_density(dens,atom%orbitals%pmat,atom%basis,atom%state%maxl_occ,&
         "RHO",atom%basis%grid%rad,error)
    IF(up) THEN
       WRITE(iw,'(T8,4ES25.12E3)') (4._dp*pi*dens(j)*atom%basis%grid%rad2(j),j=1,nr)
    ELSE
       WRITE(iw,'(T8,4ES25.12E3)') (4._dp*pi*dens(j)*atom%basis%grid%rad2(j),j=nr,1,-1)
    END IF
    WRITE(iw,'(T4,A)') '</PP_RHOATOM>'
    DEALLOCATE(dens,STAT=ierr)
    CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)

    WRITE(iw,'(A)') '</UPF>'

  END SUBROUTINE atom_write_upf
  SUBROUTINE compose(string,tag,counter,rval,ival,cval,isfinal)
    CHARACTER(len=*)                         :: string, tag
    INTEGER, OPTIONAL                        :: counter
    REAL(KIND=dp), OPTIONAL                  :: rval
    INTEGER, OPTIONAL                        :: ival
    CHARACTER(len=*), OPTIONAL               :: cval
    LOGICAL, OPTIONAL                        :: isfinal

    CHARACTER(LEN=default_string_length)     :: str
    LOGICAL                                  :: fin

     IF(PRESENT(counter)) THEN
        WRITE(str,"(I12)") counter
     ELSEIF(PRESENT(rval)) THEN
        WRITE(str,"(G18.8)") rval
     ELSEIF(PRESENT(ival)) THEN
        WRITE(str,"(I12)") ival
     ELSEIF(PRESENT(cval)) THEN
        WRITE(str,"(A)") TRIM(ADJUSTL(cval))
     ELSE
        WRITE(str,"(A)") ""
     END IF
     fin = .FALSE.
     IF(PRESENT(isfinal)) fin=isfinal
     IF(PRESENT(counter)) THEN
       IF(fin) THEN
         WRITE(string,"(A,A1,A,A1)") TRIM(ADJUSTL(tag)),'.',TRIM(ADJUSTL(str)),'>'
       ELSE
         WRITE(string,"(A,A1,A)") TRIM(ADJUSTL(tag)),'.',TRIM(ADJUSTL(str))
       END IF
     ELSE
       IF(fin) THEN
         WRITE(string,"(A,A2,A,A2)") TRIM(ADJUSTL(tag)),'="',TRIM(ADJUSTL(str)),'>"'
       ELSE
         WRITE(string,"(A,A2,A,A1)") TRIM(ADJUSTL(tag)),'="',TRIM(ADJUSTL(str)),'"'
       END IF
     END IF

  END SUBROUTINE compose
! *****************************************************************************

END MODULE atom_energy
