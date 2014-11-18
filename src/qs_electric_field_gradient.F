!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculates electric field gradients
!>      H.M. Petrili, P.E. Blochl, P. Blaha, K. Schwarz, PRB 57, 14690 (1998)
!> \par History
!>      12.2007 Added checksum for interpolation regtest [rdeclerck]
!> \author JGH (03-05-2006)
! *****************************************************************************
MODULE qs_electric_field_gradient
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_output_handling,              ONLY: cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE eigenvalueproblems,              ONLY: diagonalise
  USE input_section_types,             ONLY: section_get_lval,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: fac,&
                                             fourpi
  USE message_passing,                 ONLY: mp_sum
  USE orbital_pointers,                ONLY: indso,&
                                             nsoset
  USE particle_types,                  ONLY: particle_type
  USE paw_proj_set_types,              ONLY: get_paw_proj_set,&
                                             paw_proj_set_type
  USE physcon,                         ONLY: a_bohr,&
                                             e_charge,&
                                             joule
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_methods,                      ONLY: pw_dr2,&
                                             pw_integral_ab,&
                                             pw_smoothing,&
                                             pw_structure_factor,&
                                             pw_transfer
  USE pw_poisson_methods,              ONLY: pw_poisson_solve
  USE pw_poisson_types,                ONLY: pw_poisson_type
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_type
  USE pw_spline_utils,                 ONLY: &
       Eval_Interp_Spl3_pbc, Eval_d_Interp_Spl3_pbc, find_coeffs, &
       pw_spline_do_precond, pw_spline_precond_create, &
       pw_spline_precond_release, pw_spline_precond_set_kind, &
       pw_spline_precond_type, spl3_pbc
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_p_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_gapw_densities,               ONLY: prepare_gapw_den
  USE qs_grid_atom,                    ONLY: grid_atom_type
  USE qs_harmonics_atom,               ONLY: get_none0_cg_list,&
                                             harmonics_atom_type
  USE qs_kind_types,                   ONLY: get_qs_kind,&
                                             qs_kind_type
  USE qs_ks_methods,                   ONLY: calc_rho_tot_gspace
  USE qs_rho_atom_types,               ONLY: rho_atom_type
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: get_limit,&
                                             sort
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: qs_efg_calc

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_electric_field_gradient'

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE qs_efg_calc(qs_env,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_efg_calc', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=2)                         :: element_symbol
    INTEGER :: aint_precond, handle, i, iat, iatom, ij, ikind, istat, j, &
      max_iter, natom, natomkind, nkind, nspins, precond_kind, unit_nr
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: efg_debug, efg_interpolation, &
                                                failure, gapw, paw_atom, &
                                                smoothing, success
    REAL(KIND=dp)                            :: chk_spl, ecut, efg_units, &
                                                efg_val, ehartree, eps_r, &
                                                eps_x, f1, f2, sigma
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: efg_diagval, vh0
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: efg_pw, efg_tensor
    REAL(KIND=dp), DIMENSION(3)              :: eigenvalues, ra
    REAL(KIND=dp), DIMENSION(3, 3)           :: eigenvectors
    REAL(KIND=dp), DIMENSION(:), POINTER     :: rvals
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(grid_atom_type), POINTER            :: grid_atom
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type)                          :: dvr2rs, rho_tot_gspace, &
                                                structure_factor, &
                                                v_hartree_gspace
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: dvr2, dvspl
    TYPE(pw_poisson_type), POINTER           :: poisson_env
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(pw_spline_precond_type), POINTER    :: precond
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(rho_atom_type), DIMENSION(:), &
      POINTER                                :: rho_atom_set
    TYPE(section_vals_type), POINTER         :: dft_section, input, &
                                                interp_section

    NULLIFY(atomic_kind_set, qs_kind_set, dft_control, para_env, particle_set, rho, &
         rho_atom_set, input, dft_section, interp_section)

    CALL timeset(routineN,handle)

    logger => cp_error_get_logger(error)

    chk_spl = 0.0_dp
    efg_units = Joule/a_bohr**2/e_charge * 1.e-21_dp
    f1 = SQRT(15._dp/fourpi)
    f2 = SQRT(5._dp/fourpi)

    CALL get_qs_env(qs_env=qs_env,dft_control=dft_control,&
                    rho=rho,qs_kind_set=qs_kind_set,&
                    atomic_kind_set=atomic_kind_set,&
                    rho_atom_set=rho_atom_set,pw_env=pw_env,&
                    particle_set=particle_set,para_env=para_env,&
                    input=input,error=error)

    dft_section => section_vals_get_subs_vals(input,"DFT",error=error)

    efg_interpolation = section_get_lval(section_vals=dft_section,&
                keyword_name="PRINT%ELECTRIC_FIELD_GRADIENT%INTERPOLATION",&
                error=error)
    efg_debug = section_get_lval(section_vals=dft_section,&
                keyword_name="PRINT%ELECTRIC_FIELD_GRADIENT%DEBUG",&
                error=error)
    CALL section_vals_val_get(dft_section,&
                "PRINT%ELECTRIC_FIELD_GRADIENT%GSPACE_SMOOTHING",&
                r_vals=rvals,error=error)
    ecut = rvals(1)
    sigma = rvals(2)
    IF (ecut == 0._dp .AND. sigma <= 0._dp) THEN
      smoothing = .FALSE.
      ecut = 1.e10_dp   ! not used, just to have vars defined
      sigma = 1._dp     ! not used, just to have vars defined
    ELSEIF (ecut == -1._dp .AND. sigma == -1._dp) THEN
      smoothing = .TRUE.
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool,error=error)
      CALL pw_pool_create_pw(auxbas_pw_pool,dvr2rs%pw, &
                              use_data=REALDATA3D,in_space=REALSPACE,&
                              error=error)
      ecut  = 2._dp * dvr2rs%pw%pw_grid%cutoff * 0.875_dp
      sigma = 2._dp * dvr2rs%pw%pw_grid%cutoff * 0.125_dp
      CALL pw_pool_give_back_pw(auxbas_pw_pool,dvr2rs%pw,error=error)
    ELSE
      smoothing = .TRUE.
    END IF
    CPPostcondition(ecut > 0._dp,cp_failure_level,routineP,error,failure)
    CPPostcondition(sigma > 0._dp,cp_failure_level,routineP,error,failure)

    unit_nr=cp_print_key_unit_nr(logger,dft_section,"PRINT%ELECTRIC_FIELD_GRADIENT",&
               extension=".efg",log_filename=.FALSE.,error=error)

    IF (unit_nr > 0) THEN
      WRITE(unit_nr,"(/,A,/)") " ELECTRIC FIELD GRADIENTS [10**21 V/m^2]"
      IF (efg_interpolation) THEN
        WRITE(unit_nr,"(T16,A)") &
            " Smooth potential contribution calculated by spline interpolation"
      ELSE
        WRITE(unit_nr,"(T12,A)") &
            " Smooth potential contribution calculated by plane wave interpolation"
      END IF
      IF (smoothing) THEN
        WRITE(unit_nr,"(T36,A)") &
            " G-Space potential smoothed by Fermi function"
        WRITE(unit_nr,"(T36,A,T71,F10.4)") " Cutoff (eV) ",ecut
        WRITE(unit_nr,"(T36,A,T71,F10.4)") " Width (eV) ",sigma
      END IF
      WRITE(unit_nr,*)
    END IF

    gapw = dft_control%qs_control%gapw
    nspins = dft_control%nspins

    natom = SIZE(particle_set,1)
    ALLOCATE (efg_tensor(3,3,natom),STAT=istat )
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    efg_tensor = 0._dp
    IF ( efg_debug ) THEN
      ALLOCATE (efg_pw(3,3,natom),STAT=istat )
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      efg_pw = 0._dp
    END IF
    ALLOCATE (efg_diagval(3,natom),STAT=istat )
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    efg_diagval = 0._dp

    ALLOCATE (vh0(1:natom,-2:2),STAT=istat )
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    vh0 = 0._dp

    ALLOCATE (dvr2(6),STAT=istat )
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (dvspl(6),STAT=istat )
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    !prepare calculation
    CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool,&
                    poisson_env=poisson_env,error=error)
    IF (gapw) CALL prepare_gapw_den(qs_env,do_rho0=.TRUE.,error=error)

    !calculate electrostatic potential
    CALL pw_pool_create_pw(auxbas_pw_pool,v_hartree_gspace%pw, &
                            use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,&
                            error=error)
    CALL pw_pool_create_pw(auxbas_pw_pool,rho_tot_gspace%pw,&
                            use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,&
                            error=error)
    CALL calc_rho_tot_gspace(rho_tot_gspace, qs_env, rho,error=error)

    CALL pw_poisson_solve(poisson_env,rho_tot_gspace%pw,ehartree,&
                          v_hartree_gspace%pw,error=error)

    ! smoothing of potential
    IF ( smoothing ) CALL pw_smoothing(v_hartree_gspace%pw,ecut,sigma,error=error)

    CALL pw_pool_give_back_pw(auxbas_pw_pool,rho_tot_gspace%pw,error=error)

    DO i=1,3
      DO j=1,i
        ij=(i*(i-1))/2+j
        CALL pw_pool_create_pw(auxbas_pw_pool,dvr2(ij)%pw, &
                                use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,&
                                error=error)
        CALL pw_dr2 ( v_hartree_gspace%pw, dvr2(ij)%pw, i, j, error=error)
      END DO
    END DO

    IF ( .NOT. efg_interpolation ) THEN
      CALL pw_pool_create_pw(auxbas_pw_pool,structure_factor%pw, &
                              use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,&
                              error=error)
    ELSE

      interp_section =>  section_vals_get_subs_vals(dft_section,&
                           "PRINT%ELECTRIC_FIELD_GRADIENT%INTERPOLATOR",error=error)
      CALL section_vals_val_get(interp_section,"aint_precond", &
                                i_val=aint_precond, error=error)
      CALL section_vals_val_get(interp_section,"precond",i_val=precond_kind, error=error)
      CALL section_vals_val_get(interp_section,"max_iter",i_val=max_iter, error=error)
      CALL section_vals_val_get(interp_section,"eps_r",r_val=eps_r, error=error)
      CALL section_vals_val_get(interp_section,"eps_x",r_val=eps_x, error=error)

      CALL pw_pool_create_pw(auxbas_pw_pool,dvr2rs%pw, &
                              use_data=REALDATA3D,in_space=REALSPACE,&
                              error=error)
      DO i=1,6
        CALL pw_pool_create_pw(auxbas_pw_pool,dvspl(i)%pw, &
                            use_data=REALDATA3D,in_space=REALSPACE,&
                            error=error)
        CALL pw_transfer ( dvr2(i)%pw, dvr2rs%pw , error=error)
        ! calculate spline coefficients
        CALL pw_spline_precond_create(precond,precond_kind=aint_precond,&
          pool=auxbas_pw_pool,pbc=.TRUE.,transpose=.FALSE.,error=error)
        CALL pw_spline_do_precond(precond,dvr2rs%pw,dvspl(i)%pw,error=error)
        CALL pw_spline_precond_set_kind(precond,precond_kind,error=error)
        success=find_coeffs(values=dvr2rs%pw,coeffs=dvspl(i)%pw,&
             linOp=spl3_pbc,preconditioner=precond,pool=auxbas_pw_pool, &
             eps_r=eps_r,eps_x=eps_x,max_iter=max_iter, &
             error=error)
        CPPostconditionNoFail(success,cp_warning_level,routineP,error)
        CALL pw_spline_precond_release(precond,error=error)
        CALL pw_pool_give_back_pw(auxbas_pw_pool,dvr2(i)%pw,error=error)
      END DO
      CALL pw_pool_give_back_pw(auxbas_pw_pool,dvr2rs%pw,error=error)
    END IF

    nkind = SIZE(qs_kind_set)

    DO ikind = 1,nkind
      NULLIFY(grid_atom,atom_list)
      CALL get_atomic_kind(atomic_kind_set(ikind), atom_list=atom_list,natom=natomkind)
      CALL get_qs_kind(qs_kind_set(ikind),grid_atom=grid_atom,paw_atom=paw_atom)
      DO iat = 1,natomkind
        iatom = atom_list(iat)
        ra = particle_set(iatom)%r
        IF ( efg_interpolation ) THEN
          DO i=1,3
            DO j=1,i
              ij=(i*(i-1))/2+j
              efg_val = Eval_Interp_Spl3_pbc(ra,dvspl(ij)%pw,error)
              efg_tensor(i,j,iatom)=-efg_val
              efg_tensor(j,i,iatom)=efg_tensor(i,j,iatom)
              IF ( efg_debug ) THEN
                 chk_spl = chk_spl + efg_val + &
                    SUM(Eval_d_Interp_Spl3_pbc(ra,dvspl(ij)%pw,error))
              END IF
            END DO
          END DO
        ELSE
          CALL pw_structure_factor(structure_factor%pw,ra,error=error)
          DO i=1,3
            DO j=1,i
              ij=(i*(i-1))/2+j
              efg_tensor(i,j,iatom)=-pw_integral_ab(dvr2(ij)%pw,structure_factor%pw,error=error)
              efg_tensor(j,i,iatom)=efg_tensor(i,j,iatom)
            END DO
          END DO
          efg_tensor(:,:,iatom)=efg_tensor(:,:,iatom)/structure_factor%pw%pw_grid%vol
        END IF
        IF ( efg_debug ) THEN
          efg_pw(:,:,iatom)=efg_tensor(:,:,iatom)
        END IF
      END DO

      IF( paw_atom) THEN
        CALL vlimit_atom(para_env,vh0,rho_atom_set,qs_kind_set(ikind),&
                         atom_list,grid_atom,natomkind,nspins,error)
        DO iat = 1,natomkind
          iatom = atom_list(iat)
          efg_tensor(1,1,iatom)=efg_tensor(1,1,iatom) &
                                   +f1*(vh0(iatom,2))-f2*(vh0(iatom,0))
          efg_tensor(2,2,iatom)=efg_tensor(2,2,iatom) &
                                   -f1*(vh0(iatom,2))-f2*(vh0(iatom,0))
          efg_tensor(3,3,iatom)=efg_tensor(3,3,iatom)+2._dp*f2*(vh0(iatom,0))
          efg_tensor(1,2,iatom)=efg_tensor(1,2,iatom)+f1*(vh0(iatom,-2))
          efg_tensor(2,1,iatom)=efg_tensor(2,1,iatom)+f1*(vh0(iatom,-2))
          efg_tensor(1,3,iatom)=efg_tensor(1,3,iatom)+f1*(vh0(iatom,1))
          efg_tensor(3,1,iatom)=efg_tensor(3,1,iatom)+f1*(vh0(iatom,1))
          efg_tensor(2,3,iatom)=efg_tensor(2,3,iatom)+f1*(vh0(iatom,-1))
          efg_tensor(3,2,iatom)=efg_tensor(3,2,iatom)+f1*(vh0(iatom,-1))
        END DO
      END IF

      DO iat = 1,natomkind
        iatom = atom_list(iat)
        CALL diagonalise ( efg_tensor(:,:,iatom), 3, "UPPER", &
                           eigenvalues, eigenvectors )
        CALL efgsort(eigenvalues,efg_diagval(:,iatom))
      END DO
    END DO ! ikind

    efg_tensor(:,:,:) = efg_tensor(:,:,:)*efg_units
    efg_diagval(:,:) = efg_diagval(:,:)*efg_units

    IF ( efg_debug ) THEN
      efg_pw(:,:,:) = efg_pw(:,:,:)*efg_units
      DO iatom=1,natom
        IF (unit_nr > 0) THEN
          CALL get_atomic_kind(particle_set(iatom)%atomic_kind,&
                            element_symbol=element_symbol)
          WRITE (UNIT=unit_nr,FMT="(T2,I5,T8,A,T12,A,T15,6F11.5)")&
           iatom, element_symbol,"PW",((efg_pw(i,j,iatom),i=1,j),j=1,3)
          WRITE (UNIT=unit_nr,FMT="(T12,A,T15,6F11.5)")&
           "AT",((efg_tensor(i,j,iatom)-efg_pw(i,j,iatom),i=1,j),j=1,3)
        END IF
      ENDDO
      IF (unit_nr > 0) THEN
          WRITE (UNIT=unit_nr,FMT=*)
      END IF
      IF ( efg_interpolation ) THEN
         IF (unit_nr > 0 ) THEN
             WRITE (UNIT=unit_nr,FMT="(T2,A,E24.16)") "CheckSum splines =",&
              chk_spl
             WRITE (UNIT=unit_nr,FMT=*)
         END IF
      END IF
    END IF

    DO iatom=1,natom
      IF (unit_nr > 0) THEN
        CALL get_atomic_kind(particle_set(iatom)%atomic_kind,&
                            element_symbol=element_symbol)
          WRITE (UNIT=unit_nr,FMT="(T2,I5,T8,A,T12,A,3(T39,3F14.7,/))")&
           iatom, element_symbol,"EFG Tensor",(efg_tensor(i,:,iatom),i=1,3)
          WRITE (UNIT=unit_nr,FMT="(T12,A,T39,3F14.7)")&
           "EFG Tensor eigenvalues",efg_diagval(:,iatom)
          WRITE (UNIT=unit_nr,FMT="(T12,A,T67,F14.7)") "EFG Tensor anisotropy",&
             (efg_diagval(1,iatom)-efg_diagval(2,iatom))/efg_diagval(3,iatom)
          WRITE (UNIT=unit_nr,FMT=*)
      END IF
    ENDDO

    CALL pw_pool_give_back_pw(auxbas_pw_pool,v_hartree_gspace%pw,error=error)
    IF ( .NOT. efg_interpolation ) THEN
      CALL pw_pool_give_back_pw(auxbas_pw_pool,structure_factor%pw,error=error)
      DO i=1,6
        CALL pw_pool_give_back_pw(auxbas_pw_pool,dvr2(i)%pw,error=error)
      END DO
    ELSE
      DO i=1,6
        CALL pw_pool_give_back_pw(auxbas_pw_pool,dvspl(i)%pw,error=error)
      END DO
    END IF

    DEALLOCATE (efg_tensor,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    IF ( efg_debug ) THEN
      DEALLOCATE (efg_pw,STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    END IF

    DEALLOCATE (vh0,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    DEALLOCATE (dvr2,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (dvspl,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE qs_efg_calc

! *****************************************************************************
!> \brief ...
!> \param para_env ...
!> \param vlimit ...
!> \param rho_atom_set ...
!> \param qs_kind ...
!> \param atom_list ...
!> \param grid_atom ...
!> \param natom ...
!> \param nspins ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE vlimit_atom(para_env,vlimit,rho_atom_set,qs_kind,&
                         atom_list,grid_atom,natom,nspins,error)

  ! calculate : Limit(r->0) V_hartree(r)/r^2

    TYPE(cp_para_env_type), POINTER          :: para_env
    REAL(dp), DIMENSION(:, -2:), &
      INTENT(inout)                          :: vlimit
    TYPE(rho_atom_type), DIMENSION(:), &
      POINTER                                :: rho_atom_set
    TYPE(qs_kind_type), INTENT(IN)           :: qs_kind
    INTEGER, DIMENSION(:), INTENT(IN)        :: atom_list
    TYPE(grid_atom_type), INTENT(IN)         :: grid_atom
    INTEGER, INTENT(IN)                      :: natom, nspins
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'vlimit_atom', &
      routineP = moduleN//':'//routineN

    INTEGER :: i, i1, i2, iat, iatom, icg, ipgf1, ipgf2, iset1, iset2, iso, &
      iso1, iso1_first, iso1_last, iso2, iso2_first, iso2_last, istat, l, &
      l_iso, llmax, m1s, m2s, m_iso, max_iso_not0, max_iso_not0_local, &
      max_s_harm, maxl, maxso, mepos, n1s, n2s, nset, num_pe, size1, size2
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: cg_n_list
    INTEGER, ALLOCATABLE, DIMENSION(:, :, :) :: cg_list
    INTEGER, DIMENSION(2)                    :: bo
    INTEGER, DIMENSION(:), POINTER           :: lmax, lmin, npgf, o2nindex
    LOGICAL                                  :: failure
    REAL(dp)                                 :: zet12
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: CPC_sphere
    REAL(dp), DIMENSION(20)                  :: vgg
    REAL(dp), DIMENSION(:, :), POINTER       :: coeff, zet
    REAL(dp), DIMENSION(:, :, :), POINTER    :: my_CG
    TYPE(gto_basis_set_type), POINTER        :: orb_basis
    TYPE(harmonics_atom_type), POINTER       :: harmonics
    TYPE(paw_proj_set_type), POINTER         :: paw_proj

    failure = .FALSE.
    NULLIFY(orb_basis)
    NULLIFY(harmonics)
    NULLIFY(lmin,lmax,npgf,zet,my_CG,coeff)

    CALL get_qs_kind(qs_kind=qs_kind, orb_basis_set=orb_basis,&
                         paw_proj_set=paw_proj,harmonics=harmonics)

    CALL get_gto_basis_set(gto_basis_set=orb_basis,lmax=lmax,lmin=lmin,&
                           maxl=maxl,npgf=npgf,nset=nset,zet=zet,&
                           maxso=maxso)

    CALL get_paw_proj_set(paw_proj_set=paw_proj,o2nindex=o2nindex)

    max_iso_not0 = harmonics%max_iso_not0
    max_s_harm = harmonics%max_s_harm
    llmax = harmonics%llmax

    ! Distribute the atoms of this kind
    num_pe = para_env%num_pe
    mepos  = para_env%mepos
    bo = get_limit( natom, num_pe, mepos )

    my_CG        => harmonics%my_CG

    ALLOCATE(CPC_sphere(nsoset(maxl),nsoset(maxl)),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(cg_list(2,nsoset(maxl)**2,max_s_harm),cg_n_list(max_s_harm),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    m1s = 0
    DO iset1 = 1,nset
       m2s = 0
       DO iset2 = 1,nset
          CALL get_none0_cg_list(my_CG,lmin(iset1),lmax(iset1),lmin(iset2),lmax(iset2),&
                                 max_s_harm,llmax,cg_list,cg_n_list,max_iso_not0_local,error)
          CPPrecondition(max_iso_not0_local.LE.max_iso_not0,cp_failure_level,routineP,error,failure)

          n1s = nsoset(lmax(iset1))
          DO ipgf1  = 1,npgf(iset1)
             iso1_first = nsoset(lmin(iset1)-1)+1+n1s*(ipgf1-1)+m1s
             iso1_last  = nsoset(lmax(iset1))+n1s*(ipgf1-1)+m1s
             size1 = iso1_last - iso1_first + 1
             iso1_first = o2nindex(iso1_first)
             iso1_last  = o2nindex(iso1_last)
             i1 = iso1_last - iso1_first + 1
             CPPrecondition(size1==i1,cp_failure_level,routineP,error,failure)
             i1 = nsoset(lmin(iset1)-1)+1

          n2s=nsoset(lmax(iset2))
          DO ipgf2 = 1,npgf(iset2)
            iso2_first = nsoset(lmin(iset2)-1)+1+n2s*(ipgf2-1)+m2s
            iso2_last  = nsoset(lmax(iset2))+n2s*(ipgf2-1)+m2s
            size2 = iso2_last - iso2_first + 1
            iso2_first = o2nindex(iso2_first)
            iso2_last  = o2nindex(iso2_last)
            i2 = iso2_last - iso2_first + 1
            CPPrecondition(size2==i2,cp_failure_level,routineP,error,failure)
            i2 = nsoset(lmin(iset2)-1)+1

            zet12 = zet(ipgf1,iset1)+zet(ipgf2,iset2)

            vgg = 0.0_dp

            DO iso = 1,max_iso_not0_local
               l_iso = indso(1,iso)
               IF ( l_iso /= 2 ) CYCLE
               DO icg = 1,cg_n_list(iso)
                  iso1 = cg_list(1,icg,iso)
                  iso2 = cg_list(2,icg,iso)
                  l = indso(1,iso1) + indso(1,iso2)
                  IF ( MOD(l,2) == 0 .AND. l > 0 ) THEN
                     vgg(l/2) = fourpi/10._dp * fac(l-2) / zet12**(l/2)
                  END IF
               END DO
            END DO

            DO iat = bo(1),bo(2)
              iatom = atom_list(iat)

              CPC_sphere = 0.0_dp
              DO i=1,nspins
                coeff => rho_atom_set(iatom)%cpc_h(i)%r_coef
                CPC_sphere(i1:i1+size1-1,i2:i2+size2-1) = &
                  CPC_sphere(i1:i1+size1-1,i2:i2+size2-1) + &
                  coeff(iso1_first:iso1_last,iso2_first:iso2_last)
                coeff => rho_atom_set(iatom)%cpc_s(i)%r_coef
                CPC_sphere(i1:i1+size1-1,i2:i2+size2-1) = &
                  CPC_sphere(i1:i1+size1-1,i2:i2+size2-1) - &
                  coeff(iso1_first:iso1_last,iso2_first:iso2_last)
              ENDDO ! i

              DO iso = 1,max_iso_not0_local
                 l_iso = indso(1,iso)
                 m_iso = indso(1,iso)
                 IF ( l_iso /= 2 ) CYCLE
                 DO icg = 1,cg_n_list(iso)
                    iso1 = cg_list(1,icg,iso)
                    iso2 = cg_list(2,icg,iso)
                    l = indso(1,iso1) + indso(1,iso2)
                    IF ( MOD(l,2) == 0 .AND. l > 0 ) THEN
                       vlimit(iatom,m_iso) = vlimit(iatom,m_iso) +&
                            vgg(l/2)*CPC_sphere(iso1,iso2)*my_CG(iso1,iso2,iso)
                    END IF
                 ENDDO  ! icg
              ENDDO  ! iso

            ENDDO  ! iat

          ENDDO ! ipgf2
          ENDDO   ! ipgf1
          m2s = m2s+maxso
       ENDDO  ! iset2
       m1s = m1s+maxso
    ENDDO    ! iset1

    CALL mp_sum(vlimit,para_env%group)

    DEALLOCATE(CPC_sphere,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(cg_list,cg_n_list,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE vlimit_atom

! *****************************************************************************
!> \brief ...
!> \param ein ...
!> \param eout ...
! *****************************************************************************
  SUBROUTINE efgsort(ein,eout)
    REAL(dp), DIMENSION(3), INTENT(in)       :: ein
    REAL(dp), DIMENSION(3), INTENT(inout)    :: eout

    INTEGER                                  :: i
    INTEGER, DIMENSION(3)                    :: ind
    REAL(dp), DIMENSION(3)                   :: eab

    eab = ABS(ein)
    CALL sort(eab,3,ind)
    DO i=1,3
      eout(i) = ein(ind(i))
    END DO

  END SUBROUTINE efgsort

END MODULE qs_electric_field_gradient

