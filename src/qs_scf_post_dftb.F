!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Does all kind of post scf calculations for DFTB
!> \par History
!>      Started as a copy from the GPW file
!> \author JHU (03.2013)
! *****************************************************************************
MODULE qs_scf_post_dftb
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_p_type
  USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm
  USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose,&
                                             cp_fm_cholesky_reduce,&
                                             cp_fm_cholesky_restore
  USE cp_fm_diag,                      ONLY: choose_eigv_solver
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_to_fm_submat,&
                                             cp_fm_type
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE input_section_types,             ONLY: section_get_ival,&
                                             section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE machine,                         ONLY: m_flush
  USE mathconstants,                   ONLY: twopi
  USE moments_utils,                   ONLY: get_reference_point
  USE mulliken,                        ONLY: mulliken_charges
  USE particle_types,                  ONLY: particle_type
  USE physcon,                         ONLY: debye
  USE qs_dftb_types,                   ONLY: qs_dftb_atom_type
  USE qs_dftb_utils,                   ONLY: get_dftb_atom_param
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_kind_types,                   ONLY: get_qs_kind,&
                                             qs_kind_type
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
  USE qs_rho_types,                    ONLY: qs_rho_get,&
                                             qs_rho_type
  USE qs_scf_types,                    ONLY: qs_scf_env_type
  USE qs_scf_wfn_mix,                  ONLY: wfn_mix
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  ! Global parameters
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_scf_post_dftb'
  PUBLIC :: scf_post_calculation_dftb

CONTAINS

! *****************************************************************************
!> \brief collects possible post - scf calculations and prints info / computes properties.
!> \param dft_section ...
!> \param scf_env ...
!> \param qs_env ...
!> \par History
!>      03.2013 copy of scf_post_gpw
!> \author JHU
!> \note
! *****************************************************************************
  SUBROUTINE scf_post_calculation_dftb(dft_section, scf_env, qs_env)

    TYPE(section_vals_type), POINTER         :: dft_section
    TYPE(qs_scf_env_type), OPTIONAL, POINTER :: scf_env
    TYPE(qs_environment_type), POINTER       :: qs_env

    CHARACTER(len=*), PARAMETER :: routineN = 'scf_post_calculation_dftb', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=2)                         :: ana
    CHARACTER(LEN=default_string_length)     :: aname
    INTEGER                                  :: handle, iat, iatom, ikind, &
                                                nat, natom, nkind, nspins, &
                                                unit_nr
    LOGICAL                                  :: explicit, failure
    REAL(KIND=dp)                            :: zeff
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: mcharge
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: charges
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_dbcsr_p_type), DIMENSION(:, :), &
      POINTER                                :: matrix_p, matrix_s
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(qs_dftb_atom_type), POINTER         :: dftb_kind
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(section_vals_type), POINTER         :: moments_section, print_key, &
                                                print_section, wfn_mix_section

    CALL timeset(routineN,handle)

    failure=.FALSE.

    logger => cp_get_default_logger()

    CPASSERT(ASSOCIATED(qs_env))
    NULLIFY(dft_control,rho,para_env,matrix_s,matrix_p)
    CALL get_qs_env(qs_env,atomic_kind_set=atomic_kind_set,qs_kind_set=qs_kind_set,&
         dft_control=dft_control,rho=rho,natom=natom,para_env=para_env,&
         matrix_s_kp=matrix_s)
    nspins = dft_control%nspins
    CALL qs_rho_get(rho, rho_ao_kp=matrix_p)
    ! Mulliken charges
    ALLOCATE(charges(natom,nspins),mcharge(natom))
    !
    CALL mulliken_charges(matrix_p,matrix_s,para_env,charges)
    !
    nkind = SIZE(atomic_kind_set)
    DO ikind=1,nkind
       CALL get_atomic_kind(atomic_kind_set(ikind), natom=nat)
       CALL get_qs_kind(qs_kind_set(ikind),dftb_parameter=dftb_kind)
       CALL get_dftb_atom_param(dftb_kind,zeff=zeff)
       DO iatom=1,nat
          iat = atomic_kind_set(ikind)%atom_list(iatom)
          mcharge(iat) = zeff - SUM(charges(iat,1:nspins))
       END DO
    END DO

    print_section => section_vals_get_subs_vals(dft_section,"PRINT")

    ! Mulliken
    print_key => section_vals_get_subs_vals(print_section,"MULLIKEN")
    IF (BTEST(cp_print_key_should_output(logger%iter_info,print_key),cp_p_file)) THEN
       unit_nr=cp_print_key_unit_nr(logger,print_section,"MULLIKEN",extension=".mulliken",&
                                    middle_name="",log_filename=.FALSE.)
       IF (unit_nr > 0) THEN
          WRITE (UNIT=unit_nr,FMT="(/,/,T2,A)") "MULLIKEN POPULATION ANALYSIS"
          IF (nspins == 1) THEN
             WRITE (UNIT=unit_nr,FMT="(/,T2,A,T70,A)")&
                    " # Atom   Element   Kind        Atomic population"," Net charge"
             DO ikind=1,nkind
                CALL get_atomic_kind(atomic_kind_set(ikind),natom=nat)
                CALL get_qs_kind(qs_kind_set(ikind),dftb_parameter=dftb_kind)
                CALL get_dftb_atom_param(dftb_kind,name=aname)
                ana = ADJUSTR(TRIM(ADJUSTL(aname)))
                DO iatom=1,nat
                   iat = atomic_kind_set(ikind)%atom_list(iatom)
                   WRITE (UNIT=unit_nr,&
                          FMT="(T2,I7,6X,A2,3X,I6,T39,F12.6,T69,F12.6)")&
                          iat,ana,ikind,charges(iat,1),mcharge(iat)
                END DO
             END DO
             WRITE (UNIT=unit_nr,&
                    FMT="(T2,A,T39,F12.6,T69,F12.6,/)")&
                        "# Total charge",SUM(charges(:,1)),SUM(mcharge(:))
          ELSE
             WRITE (UNIT=unit_nr,FMT="(/,T2,A)")&
                    "# Atom  Element  Kind  Atomic population (alpha,beta)   Net charge  Spin moment"
             DO ikind=1,nkind
                CALL get_atomic_kind(atomic_kind_set(ikind),natom=nat)
                CALL get_qs_kind(qs_kind_set(ikind),dftb_parameter=dftb_kind)
                CALL get_dftb_atom_param(dftb_kind,name=aname)
                ana = ADJUSTR(TRIM(ADJUSTL(aname)))
                DO iatom=1,nat
                   iat = atomic_kind_set(ikind)%atom_list(iatom)
                   WRITE (UNIT=unit_nr,&
                          FMT="(T2,I6,5X,A2,2X,I6,T29,4(1X,F12.6))")&
                          iat,ana,ikind,charges(iat,1:2),mcharge(iat),charges(iat,1)-charges(iat,2)
                END DO
             END DO
             WRITE (UNIT=unit_nr,&
                    FMT="(T2,A,T29,4(1X,F12.6),/)")&
                        "# Total charge and spin",SUM(charges(:,1)),SUM(charges(:,2)),SUM(mcharge(:))
          END IF
          CALL m_flush(unit_nr)
       END IF
       CALL cp_print_key_finished_output(unit_nr,logger,print_key)
    END IF

    ! Lowdin
    print_key => section_vals_get_subs_vals(print_section,"LOWDIN")
    IF (BTEST(cp_print_key_should_output(logger%iter_info,print_key),cp_p_file)) THEN
      CALL cp_unimplemented_error(fromWhere=routineP, &
         message="Lowdin population analysis not implemented for DFTB method.", &
         error_level=cp_warning_level)
    END IF

    ! Hirshfeld
    print_key => section_vals_get_subs_vals(print_section,"HIRSHFELD")
    IF (BTEST(cp_print_key_should_output(logger%iter_info,print_key),cp_p_file)) THEN
      CALL cp_unimplemented_error(fromWhere=routineP, &
         message="Hirshfeld charges not available for DFTB method.", &
         error_level=cp_warning_level)
    END IF

    ! Dipole Moments
    print_key => section_vals_get_subs_vals(print_section,"MOMENTS")
    IF(BTEST(cp_print_key_should_output(logger%iter_info,print_key),cp_p_file))THEN
       unit_nr=cp_print_key_unit_nr(logger,print_section,"MOMENTS",&
            extension=".data",middle_name="dftb_dipole",log_filename=.FALSE.)
       moments_section => section_vals_get_subs_vals(print_section,"MOMENTS")
       CALL dftb_dipole(qs_env, moments_section, unit_nr, mcharge)
       CALL cp_print_key_finished_output(unit_nr,logger,print_key)
    END IF

    DEALLOCATE(charges,mcharge)

    ! EPR Hyperfine Coupling
    print_key => section_vals_get_subs_vals(print_section,"HYPERFINE_COUPLING_TENSOR")
    IF (BTEST(cp_print_key_should_output(logger%iter_info,print_key),cp_p_file)) THEN
      CALL cp_unimplemented_error(fromWhere=routineP, &
         message="Hyperfine Coupling not implemented for DFTB method.", &
         error_level=cp_warning_level)
    END IF

    ! EFIELD CUBE FILE 
    print_key => section_vals_get_subs_vals(print_section,"EFIELD_CUBE")
    IF (BTEST(cp_print_key_should_output(logger%iter_info,print_key),cp_p_file)) THEN
      CALL cp_unimplemented_error(fromWhere=routineP, &
         message="Efield cube file not implemented for DFTB method.", &
         error_level=cp_warning_level)
    END IF

    ! DENSITY CUBE FILE 
    print_key => section_vals_get_subs_vals(print_section,"E_DENSITY_CUBE")
    IF (BTEST(cp_print_key_should_output(logger%iter_info,print_key),cp_p_file)) THEN
      CALL cp_unimplemented_error(fromWhere=routineP, &
         message="Electronic density cube file not implemented for DFTB method.", &
         error_level=cp_warning_level)
    END IF

    ! TOTAL DENSITY CUBE FILE 
    print_key => section_vals_get_subs_vals(print_section,"TOT_DENSITY_CUBE")
    IF (BTEST(cp_print_key_should_output(logger%iter_info,print_key),cp_p_file)) THEN
      CALL cp_unimplemented_error(fromWhere=routineP, &
         message="Total density cube file not implemented for DFTB method.", &
         error_level=cp_warning_level)
    END IF

    ! V_Hartree CUBE FILE 
    print_key => section_vals_get_subs_vals(print_section,"V_HARTREE_CUBE")
    IF (BTEST(cp_print_key_should_output(logger%iter_info,print_key),cp_p_file)) THEN
      CALL cp_unimplemented_error(fromWhere=routineP, &
         message="Hartree potential cube file not implemented for DFTB method.", &
         error_level=cp_warning_level)
    END IF

    ! V_XC CUBE FILE 
    print_key => section_vals_get_subs_vals(print_section,"V_XC_CUBE")
    IF (BTEST(cp_print_key_should_output(logger%iter_info,print_key),cp_p_file)) THEN
      CALL cp_unimplemented_error(fromWhere=routineP, &
         message="XC potential cube file not available for DFTB method.", &
         error_level=cp_warning_level)
    END IF

    ! ELF
    print_key => section_vals_get_subs_vals(print_section,"ELF_CUBE")
    IF (BTEST(cp_print_key_should_output(logger%iter_info,print_key),cp_p_file)) THEN
      CALL cp_unimplemented_error(fromWhere=routineP, &
         message="ELF not implemented for DFTB method.", &
         error_level=cp_warning_level)
    END IF

    ! KINETIC ENERGY 
    print_key => section_vals_get_subs_vals(print_section,"KINETIC_ENERGY")
    IF (BTEST(cp_print_key_should_output(logger%iter_info,print_key),cp_p_file)) THEN
      CALL cp_unimplemented_error(fromWhere=routineP, &
         message="Kinetic energy not available for DFTB method.", &
         error_level=cp_warning_level)
    END IF

    ! Xray diffraction spectrum
    print_key => section_vals_get_subs_vals(print_section,"XRAY_DIFFRACTION_SPECTRUM")
    IF (BTEST(cp_print_key_should_output(logger%iter_info,print_key),cp_p_file)) THEN
      CALL cp_unimplemented_error(fromWhere=routineP, &
         message="Xray diffraction spectrum not implemented for DFTB method.", &
         error_level=cp_warning_level)
    END IF

    ! Electric field gradients
    print_key => section_vals_get_subs_vals(print_section,"ELECTRIC_FIELD_GRADIENT")
    IF (BTEST(cp_print_key_should_output(logger%iter_info,print_key),cp_p_file)) THEN
      CALL cp_unimplemented_error(fromWhere=routineP, &
         message="Electric field gradient not implemented for DFTB method.", &
         error_level=cp_warning_level)
    END IF

    ! STM
    print_key => section_vals_get_subs_vals(print_section,"STM")
    IF (BTEST(cp_print_key_should_output(logger%iter_info,print_key),cp_p_file)) THEN
      CALL cp_unimplemented_error(fromWhere=routineP, &
         message="STM not implemented for DFTB method.", &
         error_level=cp_warning_level)
    END IF

    ! MO 
    print_key => section_vals_get_subs_vals(print_section,"MO")
    IF (BTEST(cp_print_key_should_output(logger%iter_info,print_key),cp_p_file)) THEN
      CALL cp_unimplemented_error(fromWhere=routineP, &
         message="Printing of MO properties not implemented for DFTB method.", &
         error_level=cp_warning_level)
    END IF

    ! MO CUBES
    print_key => section_vals_get_subs_vals(print_section,"MO_CUBES")
    IF (BTEST(cp_print_key_should_output(logger%iter_info,print_key),cp_p_file)) THEN
      CALL cp_unimplemented_error(fromWhere=routineP, &
         message="Printing of MO cube files not implemented for DFTB method.", &
         error_level=cp_warning_level)
    END IF

    ! Wavefunction mixing
    wfn_mix_section => section_vals_get_subs_vals(dft_section,"PRINT%WFN_MIX")
    CALL section_vals_get(wfn_mix_section,explicit=explicit)
    IF(explicit.AND..NOT.qs_env%run_rtp) CALL wfn_mix_dftb(qs_env,dft_section,scf_env)

    ! PLUS_U
    print_key => section_vals_get_subs_vals(print_section,"PLUS_U")
    IF (BTEST(cp_print_key_should_output(logger%iter_info,print_key),cp_p_file)) THEN
      CALL cp_unimplemented_error(fromWhere=routineP, &
         message="DFT+U method not implemented for DFTB method.", &
         error_level=cp_warning_level)
    END IF

    ! PDOS
    print_key => section_vals_get_subs_vals(print_section,"PDOS")
    IF (BTEST(cp_print_key_should_output(logger%iter_info,print_key),cp_p_file)) THEN
      CALL cp_unimplemented_error(fromWhere=routineP, &
         message="Projected DOS not implemented for DFTB method.", &
         error_level=cp_warning_level)
    END IF

    CALL timestop(handle)

  END SUBROUTINE scf_post_calculation_dftb

! *****************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param input ...
!> \param unit_nr ...
!> \param charges ...
! *****************************************************************************
  SUBROUTINE dftb_dipole(qs_env, input, unit_nr, charges)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(section_vals_type), POINTER         :: input
    INTEGER, INTENT(in)                      :: unit_nr
    REAL(KIND=dp), DIMENSION(:), INTENT(in)  :: charges

    CHARACTER(len=*), PARAMETER :: routineN = 'dftb_dipole', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: description, dipole_type
    COMPLEX(KIND=dp)                         :: dzeta, dzphase(3), zeta, &
                                                zphase(3)
    COMPLEX(KIND=dp), DIMENSION(3)           :: dggamma, ggamma
    INTEGER                                  :: i, iat, ikind, j, nat, &
                                                reference
    LOGICAL                                  :: do_berry
    REAL(KIND=dp) :: charge_tot, ci(3), dci(3), dipole(3), dipole_deriv(3), &
      drcc(3), dria(3), dtheta, gvec(3), q, rcc(3), ria(3), theta, tmp(3), &
      via(3)
    REAL(KIND=dp), DIMENSION(:), POINTER     :: ref_point
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

    NULLIFY(atomic_kind_set,cell)
    CALL get_qs_env(qs_env,atomic_kind_set=atomic_kind_set,&
                    particle_set=particle_set,cell=cell)

    ! Reference point
    reference = section_get_ival(input,keyword_name="REFERENCE")
    NULLIFY(ref_point)
    description='[DIPOLE]'
    CALL section_vals_val_get(input,"REF_POINT",r_vals=ref_point)
    CALL section_vals_val_get(input,"PERIODIC",l_val=do_berry)

    CALL get_reference_point(rcc,drcc,qs_env=qs_env,reference=reference,ref_point=ref_point)

    ! Dipole deriv will be the derivative of the Dipole(dM/dt=\sum e_j v_j)
    dipole_deriv = 0.0_dp
    dipole       = 0.0_dp
    IF(do_berry)THEN
       dipole_type="[BERRY PHASE]"
       rcc = pbc(rcc,cell)
       charge_tot = 0._dp
       charge_tot=SUM(charges)
       ria  = twopi * MATMUL(cell%h_inv, rcc)
       zphase  = CMPLX(COS(ria), SIN(ria), dp)**charge_tot

       dria = twopi * MATMUL(cell%h_inv, drcc)
       dzphase = charge_tot * CMPLX(-SIN(ria), COS(ria), dp)**(charge_tot-1.0_dp)* dria

       ggamma  = CMPLX(1.0_dp,0.0_dp,KIND=dp)
       dggamma = CMPLX(0.0_dp,0.0_dp,KIND=dp)
       DO ikind = 1, SIZE(atomic_kind_set)
          CALL get_atomic_kind(atomic_kind_set(ikind),natom=nat)
          DO i = 1, nat
             iat = atomic_kind_set(ikind)%atom_list(i)
             ria = particle_set(iat)%r(:)
             ria = pbc(ria,cell)
             via = particle_set(iat)%v(:)
             q=charges(iat)
             DO j = 1, 3
                gvec  = twopi*cell%h_inv(j,:)
                theta = SUM(ria(:)*gvec(:))
                dtheta= SUM(via(:)*gvec(:))
                zeta  = CMPLX(COS(theta),SIN(theta),KIND=dp)**(-q)
                dzeta = -q * CMPLX(-SIN(theta),COS(theta),KIND=dp)**(-q-1.0_dp)* dtheta
                dggamma(j)= dggamma(j)* zeta + ggamma(j)* dzeta
                ggamma(j) = ggamma(j)* zeta
             END DO
          ENDDO
       END DO
       dggamma = dggamma * zphase + ggamma * dzphase
       ggamma = ggamma * zphase
       IF(ALL(REAL(ggamma,KIND=dp)/= 0.0_dp))THEN
          tmp = AIMAG(ggamma)/REAL(ggamma,KIND=dp)
          ci  = ATAN(tmp)
          dci =(1.0_dp /(1.0_dp + tmp**2))* &
              (AIMAG(dggamma)*REAL(ggamma,KIND=dp)-AIMAG(ggamma)*REAL(dggamma,KIND=dp))/(REAL(ggamma,KIND=dp))**2
          dipole       = MATMUL(cell%hmat, ci)/twopi
          dipole_deriv = MATMUL(cell%hmat, dci)/twopi
       END IF
    ELSE
       dipole_type="[Non Periodic]"
       DO i = 1, SIZE(particle_set)
          ! no pbc(particle_set(i)%r(:),cell) so that the total dipole is the sum of the molecular dipoles
          ria = particle_set(i)%r(:) 
          q=charges(i)
          dipole = dipole - q *(ria-rcc)
          dipole_deriv(:)= dipole_deriv(:)- q *(particle_set(i)%v(:)- drcc)
       END DO
    END IF
    IF(unit_nr>0)THEN
       WRITE(unit_nr,'(1X,A,T48,3F11.6)')"DIPOLE "//TRIM(dipole_type)//"(A.U.)|",dipole
       WRITE(unit_nr,'(1X,A,T48,3F11.6)')"DIPOLE "//TRIM(dipole_type)//"(Debye)|",dipole*debye
       WRITE(unit_nr,'(1X,A,T48,3F11.6)')"DIPOLE "//TRIM(dipole_type)//" DERIVATIVE(A.U.)|",dipole_deriv
    END IF

  END SUBROUTINE dftb_dipole

! *****************************************************************************
!> \brief computes the MOs and calls the wavefunction mixing routine.
!> \param qs_env ...
!> \param dft_section ...
!> \param scf_env ...
!> \author Florian Schiffmann
!> \note
! *****************************************************************************

  SUBROUTINE wfn_mix_dftb(qs_env,dft_section,scf_env)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(section_vals_type), POINTER         :: dft_section
    TYPE(qs_scf_env_type), POINTER           :: scf_env

    CHARACTER(len=*), PARAMETER :: routineN = 'wfn_mix_dftb', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ispin, nao, nmo, output_unit
    REAL(dp), DIMENSION(:), POINTER          :: mo_eigenvalues
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_s
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: lumos
    TYPE(cp_fm_struct_type), POINTER         :: ao_ao_fmstruct, ao_lumo_struct
    TYPE(cp_fm_type), POINTER                :: KS_tmp, mo_coeff, MO_tmp, &
                                                S_tmp, work
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(section_vals_type), POINTER         :: wfn_mix_section

    logger => cp_get_default_logger()
    CALL get_qs_env(qs_env=qs_env,matrix_s=matrix_s,matrix_ks=matrix_ks,&
                    particle_set=particle_set,atomic_kind_set=atomic_kind_set,&
                    qs_kind_set=qs_kind_set, mos=mos)

    wfn_mix_section => section_vals_get_subs_vals(dft_section,"PRINT%WFN_MIX") 

    CALL get_mo_set(mos(1)%mo_set,mo_coeff=mo_coeff,nao=nao)
    
    CALL cp_fm_struct_create(fmstruct=ao_ao_fmstruct, nrow_global=nao, ncol_global=nao,&
         template_fmstruct=mo_coeff%matrix_struct)
    CALL cp_fm_create(S_tmp,matrix_struct=ao_ao_fmstruct)
    CALL cp_fm_create(KS_tmp,matrix_struct=ao_ao_fmstruct)
    CALL cp_fm_create(MO_tmp,matrix_struct=ao_ao_fmstruct)
    CALL cp_fm_create(work,matrix_struct=ao_ao_fmstruct)
    ALLOCATE(lumos(SIZE(mos)))

    CALL copy_dbcsr_to_fm(matrix_s(1)%matrix,S_tmp)
    CALL cp_fm_cholesky_decompose(S_tmp)

    DO ispin=1,SIZE(mos)
       CALL get_mo_set(mos(ispin)%mo_set,mo_coeff=mo_coeff,eigenvalues=mo_eigenvalues,nmo=nmo)
       CALL cp_fm_struct_create(fmstruct=ao_lumo_struct, nrow_global=nao, ncol_global=nao-nmo,&
            template_fmstruct=mo_coeff%matrix_struct)
       
       CALL cp_fm_create(lumos(ispin)%matrix,matrix_struct=ao_lumo_struct)
       CALL copy_dbcsr_to_fm(matrix_ks(ispin)%matrix,KS_tmp)
       CALL cp_fm_cholesky_reduce(KS_tmp,S_tmp)
       CALL choose_eigv_solver(KS_tmp,work,mo_eigenvalues)
       CALL cp_fm_cholesky_restore(work,nao,S_tmp,MO_tmp,"SOLVE")
       CALL cp_fm_to_fm_submat(MO_tmp, mo_coeff, nao, nmo, 1, 1, 1, 1)
       CALL cp_fm_to_fm_submat(MO_tmp, lumos(ispin)%matrix, nao, nao-nmo, 1, nmo+1, 1, 1)
       
       CALL cp_fm_struct_release(ao_lumo_struct)
    END DO

    output_unit= cp_logger_get_default_io_unit(logger)
    CALL wfn_mix(mos, particle_set, dft_section, atomic_kind_set, qs_kind_set, &
                 lumos, scf_env, matrix_s, output_unit)

    DO ispin=1,SIZE(mos)
       CALL cp_fm_release(lumos(ispin)%matrix)
    END DO
    DEALLOCATE(lumos)   
    CALL cp_fm_release(S_tmp)
    CALL cp_fm_release(MO_tmp)
    CALL cp_fm_release(KS_tmp)
    CALL cp_fm_release(work)
    CALL cp_fm_struct_release(ao_ao_fmstruct)

  END SUBROUTINE wfn_mix_dftb

! *****************************************************************************

END MODULE qs_scf_post_dftb
