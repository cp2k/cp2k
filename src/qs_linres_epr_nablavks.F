!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculates Nabla V_KS (local part if PSP) on the different grids
!> \par History
!>       created 06-2007 [RD]
!> \author RD
! *****************************************************************************
MODULE qs_linres_epr_nablavks
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cell_types,                      ONLY: cell_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_realspace_grid_cube,          ONLY: cp_pw_to_cube
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE erf_fn,                          ONLY: erf,&
                                             erfc
  USE external_potential_types,        ONLY: all_potential_type,&
                                             get_potential,&
                                             gth_potential_type
  USE f77_blas
  USE hartree_local_methods,           ONLY: calculate_Vh_1center
  USE input_section_types,             ONLY: section_get_ivals,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE mathconstants,                   ONLY: rootpi,&
                                             twopi
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: particle_type
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_methods,                      ONLY: pw_axpy,&
                                             pw_copy,&
                                             pw_derive,&
                                             pw_transfer,&
                                             pw_zero
  USE pw_poisson_methods,              ONLY: pw_poisson_solve
  USE pw_poisson_types,                ONLY: pw_poisson_type
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_p_type,&
                                             pw_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_gapw_densities,               ONLY: prepare_gapw_den
  USE qs_grid_atom,                    ONLY: grid_atom_type
  USE qs_harmonics_atom,               ONLY: harmonics_atom_type
  USE qs_ks_methods,                   ONLY: calc_rho_tot_gspace
  USE qs_linres_types,                 ONLY: epr_env_type,&
                                             get_epr_env,&
                                             nablavks_atom_type
! R0
  USE qs_local_rho_types,              ONLY: rhoz_type
  USE qs_rho0_types,                   ONLY: rho0_atom_type
  USE qs_rho_atom_methods,             ONLY: calculate_rho_atom_coeff
  USE qs_rho_atom_types,               ONLY: get_rho_atom,&
                                             rho_atom_coeff,&
                                             rho_atom_type
! R0
  USE qs_rho_types,                    ONLY: qs_rho_p_type,&
                                             qs_rho_type
  USE qs_vxc,                          ONLY: qs_vxc_create
  USE qs_vxc_atom,                     ONLY: calculate_vxc_atom
  USE termination,                     ONLY: stop_memory
  USE util,                            ONLY: get_limit
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: epr_nablavks

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_linres_epr_nablavks'

CONTAINS

! *****************************************************************************
!> \brief Evaluates Nabla V_KS on the grids
!> \par History
!>      06.2006 created [RD]
!> \author RD
! *****************************************************************************
  SUBROUTINE epr_nablavks(epr_env,qs_env,error)

    TYPE(epr_env_type)                       :: epr_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'epr_nablavks', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: ext, filename
    COMPLEX(KIND=dp)                         :: gtemp
    INTEGER :: bo_atom(2), ia, iat, iatom, idir, iexp, ig, ikind, ir, iso, &
      ispin, istat, ix, iy, iz, natom, nexp_ppl, nkind, nspins, output_unit, &
      unit_nr
    INTEGER, DIMENSION(2, 3)                 :: bo
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: failure, gapw, gapw_xc, &
                                                gth_gspace, ionode, &
                                                make_soft, paw_atom
    REAL(KIND=dp) :: alpha, alpha_core, arg, charge, ehartree, exc, exp_rap, &
      gapw_max_alpha, hard_value, soft_value, sqrt_alpha, sqrt_rap
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: vh1_rad_h, vh1_rad_s
    REAL(KIND=dp), DIMENSION(3)              :: rap, ratom, roffset, rpoint
    REAL(KIND=dp), DIMENSION(:), POINTER     :: cexp_ppl, rho_rad_z
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: rho_rad_0
    TYPE(all_potential_type), POINTER        :: all_potential
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(grid_atom_type), POINTER            :: grid_atom
    TYPE(gth_potential_type), POINTER        :: gth_potential
    TYPE(harmonics_atom_type), POINTER       :: harmonics
    TYPE(nablavks_atom_type), DIMENSION(:), &
      POINTER                                :: nablavks_atom_set
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: v_rspace_new, v_tau_rspace
    TYPE(pw_p_type), POINTER :: rho_tot_gspace, v_coulomb_gspace, &
      v_coulomb_gtemp, v_coulomb_rtemp, v_hartree_gspace, v_hartree_gtemp, &
      v_hartree_rtemp, v_xc_gtemp, v_xc_rtemp, wf_r
    TYPE(pw_poisson_type), POINTER           :: poisson_env
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(pw_type), POINTER                   :: pwx, pwy, pwz
    TYPE(qs_rho_p_type), DIMENSION(:, :), &
      POINTER                                :: nablavks_set
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(rho0_atom_type), DIMENSION(:), &
      POINTER                                :: rho0_atom_set
    TYPE(rho_atom_coeff), DIMENSION(:), &
      POINTER                                :: rho_rad_h, rho_rad_s
    TYPE(rho_atom_coeff), DIMENSION(:, :), &
      POINTER                                :: nablavks_vec_rad_h, &
                                                nablavks_vec_rad_s
    TYPE(rho_atom_type), DIMENSION(:), &
      POINTER                                :: rho_atom_set
    TYPE(rho_atom_type), POINTER             :: rho_atom
    TYPE(rhoz_type), DIMENSION(:), POINTER   :: rhoz_set
    TYPE(section_vals_type), POINTER         :: g_section, input, lr_section, &
                                                xc_section

! R0
! R0
! R0
! R0
! R0
! R0

    NULLIFY(auxbas_pw_pool)
    NULLIFY(cell)
    NULLIFY(dft_control)
    NULLIFY(g_section)
    NULLIFY(logger)
    NULLIFY(lr_section)
    NULLIFY(nablavks_set)
    NULLIFY(nablavks_atom_set) ! R0
    NULLIFY(nablavks_vec_rad_h) ! R0
    NULLIFY(nablavks_vec_rad_s) ! R0
    NULLIFY(para_env)
    NULLIFY(particle_set)
    NULLIFY(particles)
    NULLIFY(poisson_env)
    NULLIFY(pw_env)
    NULLIFY(pwx) ! R0
    NULLIFY(pwy) ! R0
    NULLIFY(pwz) ! R0
    NULLIFY(rho)
    NULLIFY(rho0_atom_set)
    NULLIFY(rho_atom_set)
    NULLIFY(rhoz_set)
    NULLIFY(subsys)
    NULLIFY(v_rspace_new)
    NULLIFY(v_tau_rspace)
    NULLIFY(xc_section)
    NULLIFY(input)

    logger => cp_error_get_logger(error)
    lr_section => section_vals_get_subs_vals(qs_env%input,"PROPERTIES%LINRES",error=error)
    ionode = logger%para_env%mepos==logger%para_env%source

    output_unit = cp_print_key_unit_nr(logger,lr_section,"PRINT%PROGRAM_RUN_INFO",&
         extension=".linresLog",error=error)

    failure = .FALSE.

!   -------------------------------------
!   Read settings
!   -------------------------------------

    g_section => section_vals_get_subs_vals(lr_section,&
                 "EPR%PRINT%G_TENSOR",error=error)

    CALL section_vals_val_get(g_section,"gapw_max_alpha",r_val=gapw_max_alpha, error=error)

    gth_gspace = .TRUE.

!   -------------------------------------
!   Get nablavks arrays
!   -------------------------------------

    CALL get_epr_env(epr_env,nablavks_set=nablavks_set,& ! R0
                     nablavks_atom_set=nablavks_atom_set,& ! R0
                     error=error) ! R0

    DO ispin = 1,SIZE(nablavks_set,2)
       DO idir = 1,SIZE(nablavks_set,1)
          CALL pw_zero(nablavks_set(idir,ispin)%rho%rho_r(1)%pw, error=error)
       END DO
    END DO

    pwx => nablavks_set(1,1)%rho%rho_r(1)%pw
    pwy => nablavks_set(2,1)%rho%rho_r(1)%pw
    pwz => nablavks_set(3,1)%rho%rho_r(1)%pw

    roffset = -REAL(MODULO(pwx%pw_grid%npts,2),dp)*pwx%pw_grid%dr/2.0_dp

!   -------------------------------------
!   Get grids / atom info
!   -------------------------------------

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    input=input,&
                    cell=cell,&
                    dft_control=dft_control,&
                    para_env=para_env,&
                    particle_set=particle_set,&
                    pw_env=pw_env,&
                    rho=rho,&
                    rho_atom_set=rho_atom_set,&
                    rho0_atom_set=rho0_atom_set,&
                    rhoz_set=rhoz_set,&
                    subsys=subsys,&
                    error=error)

    CALL pw_env_get(pw_env,auxbas_pw_pool=auxbas_pw_pool,&
                    poisson_env=poisson_env,&
                    error=error)

    CALL cp_subsys_get(subsys,particles=particles,error=error)

    gapw    = dft_control%qs_control%gapw
    gapw_xc = dft_control%qs_control%gapw_xc
    nkind = SIZE(atomic_kind_set)
    nspins  = dft_control%nspins

!   -------------------------------------
!   Add Hartree potential
!   -------------------------------------

    ALLOCATE(v_hartree_gspace,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "v_hartree_gspace",0)
    ALLOCATE(v_hartree_gtemp,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "v_hartree_gtemp",0)
    ALLOCATE(v_hartree_rtemp,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "v_hartree_rtemp",0)
    ALLOCATE(rho_tot_gspace,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "rho_tot_gspace",0)

    CALL pw_pool_create_pw(auxbas_pw_pool,v_hartree_gspace%pw, &
                           use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,&
                           error=error)
    CALL pw_pool_create_pw(auxbas_pw_pool,v_hartree_gtemp%pw, &
                           use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,&
                           error=error)
    CALL pw_pool_create_pw(auxbas_pw_pool,v_hartree_rtemp%pw,&
                           use_data=REALDATA3D,in_space=REALSPACE,&
                           error=error)
    CALL pw_pool_create_pw(auxbas_pw_pool,rho_tot_gspace%pw,&
                           use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,&
                           error=error)

    IF (gapw) THEN
       ! need to rebuild the coeff !
       CALL calculate_rho_atom_coeff(qs_env,rho%rho_ao,error=error)
       CALL prepare_gapw_den(qs_env,do_rho0=.TRUE.,error=error)
    END IF

    CALL pw_zero(rho_tot_gspace%pw, error=error)

    CALL calc_rho_tot_gspace(rho_tot_gspace,qs_env,rho,&
                             skip_nuclear_density=.NOT. gapw,error=error)

    CALL pw_poisson_solve(poisson_env,rho_tot_gspace%pw,ehartree,&
                          v_hartree_gspace%pw,error=error)

    !   -------------------------------------
    !   Atomic grids part
    !   -------------------------------------

    IF (gapw) THEN

       DO ikind = 1,nkind ! loop over atom types

          NULLIFY(atomic_kind)
          NULLIFY(atom_list)
          NULLIFY(grid_atom)
          NULLIFY(harmonics)
          NULLIFY(rho_rad_z)

          atomic_kind => atomic_kind_set(ikind)
          rho_rad_z  => rhoz_set(ikind)%r_coef

          CALL get_atomic_kind(atomic_kind=atomic_kind,&
                               atom_list=atom_list,&
                               grid_atom=grid_atom,&
                               harmonics=harmonics,&
                               natom=natom,&
                               paw_atom=paw_atom,&
                               zeff=charge,&
                               alpha_core_charge=alpha_core)

          IF (paw_atom) THEN

             ALLOCATE(vh1_rad_h(grid_atom%nr,harmonics%max_iso_not0),STAT=istat)
             ALLOCATE(vh1_rad_s(grid_atom%nr,harmonics%max_iso_not0),STAT=istat)
             CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

             ! DO iat = 1, natom ! natom = # atoms for ikind
             !
             !    iatom = atom_list(iat)
             !    ratom = particle_set(iatom)%r
             !
             !    DO ig = v_hartree_gspace%pw%pw_grid%first_gne0,v_hartree_gspace%pw%pw_grid%ngpts_cut_local
             !
             !       gtemp = fourpi * charge / cell%deth * &
             !               EXP ( - v_hartree_gspace%pw%pw_grid%gsq(ig) / (4.0_dp * alpha_core) ) &
             !               / v_hartree_gspace%pw%pw_grid%gsq(ig)
             !
             !       arg = DOT_PRODUCT(v_hartree_gspace%pw%pw_grid%g(:,ig),ratom)
             !
             !       gtemp = gtemp * CMPLX(COS(arg),-SIN(arg),KIND=dp)
             !
             !       v_hartree_gspace%pw%cc(ig) = v_hartree_gspace%pw%cc(ig) + gtemp
             !    END DO
             !    IF ( v_hartree_gspace%pw%pw_grid%have_g0 ) v_hartree_gspace%pw%cc(1) = 0.0_dp
             !
             ! END DO

             bo_atom = get_limit(natom, para_env%num_pe, para_env%mepos)

             DO iat = bo_atom(1),bo_atom(2) ! natomkind = # atoms for ikind

                iatom = atom_list(iat)
                ratom = particle_set(iatom)%r

                nablavks_vec_rad_h => nablavks_atom_set(iatom)%nablavks_vec_rad_h
                nablavks_vec_rad_s => nablavks_atom_set(iatom)%nablavks_vec_rad_s

                DO ispin = 1,nspins
                   DO idir = 1,3
                      nablavks_vec_rad_h(idir,ispin)%r_coef(:,:) = 0.0_dp
                      nablavks_vec_rad_s(idir,ispin)%r_coef(:,:) = 0.0_dp
                   END DO ! idir
                END DO ! ispin

                rho_atom => rho_atom_set(iatom)
                NULLIFY(rho_rad_h,rho_rad_s,rho_rad_0)
                CALL get_rho_atom(rho_atom=rho_atom,rho_rad_h=rho_rad_h,&
                                  rho_rad_s=rho_rad_s)
                rho_rad_0 => rho0_atom_set(iatom)%rho0_rad_h%r_coef
                vh1_rad_h = 0.0_dp
                vh1_rad_s = 0.0_dp

                CALL calculate_Vh_1center(vh1_rad_h,vh1_rad_s,rho_rad_h,rho_rad_s,rho_rad_0,rho_rad_z,grid_atom)

                DO ir = 2,grid_atom%nr

                   IF (grid_atom%rad(ir) >= atomic_kind%hard_radius) CYCLE

                   DO ia = 1,grid_atom%ng_sphere

                      ! hard part

                      DO idir= 1,3
                         hard_value = 0.0_dp
                         DO iso = 1,harmonics%max_iso_not0
                            hard_value = hard_value + &
                               vh1_rad_h(ir,iso) * harmonics%dslm_dxyz(idir,ia,iso) + &
                               harmonics%slm(ia,iso) * &
                               ( vh1_rad_h(ir - 1,iso) - vh1_rad_h(ir,iso) ) / &
                               ( grid_atom%rad(ir - 1) - grid_atom%rad(ir) ) * &
                               ( harmonics%a(idir,ia) )
                         END DO
                         nablavks_vec_rad_h(idir,1)%r_coef(ir,ia) = hard_value
                      END DO

                      ! soft part

                      DO idir= 1,3
                         soft_value = 0.0_dp
                         DO iso = 1,harmonics%max_iso_not0
                            soft_value = soft_value + &
                               vh1_rad_s(ir,iso) * harmonics%dslm_dxyz(idir,ia,iso) + &
                               harmonics%slm(ia,iso) * &
                               ( vh1_rad_s(ir - 1,iso) - vh1_rad_s(ir,iso) ) / &
                               ( grid_atom%rad(ir - 1) - grid_atom%rad(ir) ) * &
                               ( harmonics%a(idir,ia) )
                         END DO
                         nablavks_vec_rad_s(idir,1)%r_coef(ir,ia) = soft_value
                      END DO

                   END DO ! ia

                END DO ! ir

                DO idir = 1,3
                   nablavks_vec_rad_h(idir,2)%r_coef(:,:) = nablavks_vec_rad_h(idir,1)%r_coef(:,:)
                   nablavks_vec_rad_s(idir,2)%r_coef(:,:) = nablavks_vec_rad_s(idir,1)%r_coef(:,:)
                END DO

             END DO ! iat

             DEALLOCATE(vh1_rad_h,STAT=istat)
             DEALLOCATE(vh1_rad_s,STAT=istat)
             CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

          END IF ! paw_atom

       END DO ! ikind

    END IF ! gapw

    CALL pw_copy(v_hartree_gspace%pw, v_hartree_gtemp%pw, error=error)
    CALL pw_derive(v_hartree_gtemp%pw, (/1,0,0/) , error=error)
    CALL pw_transfer(v_hartree_gtemp%pw, v_hartree_rtemp%pw, error=error)
    CALL pw_copy(v_hartree_rtemp%pw, pwx, error=error)

    CALL pw_copy(v_hartree_gspace%pw, v_hartree_gtemp%pw, error=error)
    CALL pw_derive(v_hartree_gtemp%pw, (/0,1,0/) , error=error)
    CALL pw_transfer(v_hartree_gtemp%pw, v_hartree_rtemp%pw, error=error)
    CALL pw_copy(v_hartree_rtemp%pw, pwy, error=error)

    CALL pw_copy(v_hartree_gspace%pw, v_hartree_gtemp%pw, error=error)
    CALL pw_derive(v_hartree_gtemp%pw, (/0,0,1/) , error=error)
    CALL pw_transfer(v_hartree_gtemp%pw, v_hartree_rtemp%pw, error=error)
    CALL pw_copy(v_hartree_rtemp%pw, pwz, error=error)

    CALL pw_pool_give_back_pw(auxbas_pw_pool,v_hartree_gspace%pw,error=error)
    CALL pw_pool_give_back_pw(auxbas_pw_pool,v_hartree_gtemp%pw,error=error)
    CALL pw_pool_give_back_pw(auxbas_pw_pool,v_hartree_rtemp%pw,error=error)
    CALL pw_pool_give_back_pw(auxbas_pw_pool,rho_tot_gspace%pw,error=error)

    DEALLOCATE(v_hartree_gspace,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "v_hartree_gspace")
    DEALLOCATE(v_hartree_gtemp,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "v_hartree_gtemp")
    DEALLOCATE(v_hartree_rtemp,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "v_hartree_rtemp")
    DEALLOCATE(rho_tot_gspace,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "rho_tot_gspace")

!   -------------------------------------
!   Add Coulomb potential
!   -------------------------------------

    DO ikind = 1,nkind ! loop over atom types

       NULLIFY(atom_list)
       NULLIFY(atomic_kind)
       NULLIFY(grid_atom)
       NULLIFY(harmonics)

       atomic_kind => atomic_kind_set(ikind)

       CALL get_atomic_kind(atomic_kind=atomic_kind,&
                            all_potential=all_potential,&
                            atom_list=atom_list,&
                            grid_atom=grid_atom,&
                            gth_potential=gth_potential,&
                            harmonics=harmonics,&
                            natom=natom,&
                            paw_atom=paw_atom)

       IF (ASSOCIATED(gth_potential)) THEN

          NULLIFY(cexp_ppl)

          CALL get_potential(potential=gth_potential,&
                             zeff=charge,&
                             alpha_ppl=alpha,&
                             nexp_ppl=nexp_ppl,&
                             cexp_ppl=cexp_ppl)

          sqrt_alpha = SQRT(alpha)

          IF (gapw .AND. paw_atom .AND. alpha > gapw_max_alpha) THEN
             make_soft = .TRUE.
          ELSE
             make_soft = .FALSE.
          END IF

          !   -------------------------------------
          !   PW grid part
          !   -------------------------------------

          IF (gth_gspace) THEN

             ALLOCATE(v_coulomb_gspace,STAT=istat)
             IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                              "v_coulomb_gspace",0)
             ALLOCATE(v_coulomb_gtemp,STAT=istat)
             IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                              "v_coulomb_gtemp",0)
             ALLOCATE(v_coulomb_rtemp,STAT=istat)
             IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                              "v_coulomb_rtemp",0)

             CALL pw_pool_create_pw(auxbas_pw_pool,v_coulomb_gspace%pw, &
                                    use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,&
                                    error=error)
             CALL pw_pool_create_pw(auxbas_pw_pool,v_coulomb_gtemp%pw, &
                                    use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,&
                                    error=error)
             CALL pw_pool_create_pw(auxbas_pw_pool,v_coulomb_rtemp%pw,&
                                    use_data=REALDATA3D,in_space=REALSPACE,&
                                    error=error)

             CALL pw_zero(v_coulomb_gspace%pw, error=error)

             DO iat = 1, natom ! natom = # atoms for ikind

                iatom = atom_list(iat)
                ratom = particle_set(iatom)%r

                DO ig = v_coulomb_gspace%pw%pw_grid%first_gne0,v_coulomb_gspace%pw%pw_grid%ngpts_cut_local
                   gtemp = 0.0_dp
                   ! gtemp = - fourpi * charge / cell%deth * &
                   !         EXP ( - v_coulomb_gspace%pw%pw_grid%gsq(ig) / (4.0_dp * alpha) ) &
                   !         / v_coulomb_gspace%pw%pw_grid%gsq(ig)

                   IF (.NOT. make_soft) THEN

                      SELECT CASE (nexp_ppl)
                         CASE(1)
                            gtemp = gtemp + &
                                    (twopi)**(1.5_dp) / ( cell%deth * (2.0_dp * alpha)**(1.5_dp) ) * &
                                    EXP ( - v_coulomb_gspace%pw%pw_grid%gsq(ig) / (4.0_dp * alpha) ) * ( &
                            ! C1
                                    + cexp_ppl(1) &
                                    )
                         CASE(2)
                            gtemp = gtemp + &
                                    (twopi)**(1.5_dp) / ( cell%deth * (2.0_dp * alpha)**(1.5_dp) ) * &
                                    EXP ( - v_coulomb_gspace%pw%pw_grid%gsq(ig) / (4.0_dp * alpha) ) * ( &
                            ! C1
                                    + cexp_ppl(1) &
                            ! C2
                                    + cexp_ppl(2) / (2.0_dp * alpha) * &
                                    ( 3.0_dp - v_coulomb_gspace%pw%pw_grid%gsq(ig)/(2.0_dp * alpha) ) &
                                    )
                         CASE(3)
                            gtemp = gtemp + &
                                    (twopi)**(1.5_dp) / ( cell%deth * (2.0_dp * alpha)**(1.5_dp) ) * &
                                    EXP ( - v_coulomb_gspace%pw%pw_grid%gsq(ig) / (4.0_dp * alpha) ) * ( &
                            ! C1
                                    + cexp_ppl(1) &
                            ! C2
                                    + cexp_ppl(2) / (2.0_dp * alpha) * &
                                    ( 3.0_dp - v_coulomb_gspace%pw%pw_grid%gsq(ig)/(2.0_dp * alpha) ) &
                            ! C3
                                    + cexp_ppl(3) / (2.0_dp * alpha)**2 * &
                                    ( 15.0_dp - 10.0_dp * v_coulomb_gspace%pw%pw_grid%gsq(ig)/(2.0_dp * alpha) &
                                                      + ( v_coulomb_gspace%pw%pw_grid%gsq(ig)/(2.0_dp * alpha) )**2 ) &
                                    )
                         CASE(4)
                            gtemp = gtemp + &
                                    (twopi)**(1.5_dp) / ( cell%deth * (2.0_dp * alpha)**(1.5_dp) ) * &
                                    EXP ( - v_coulomb_gspace%pw%pw_grid%gsq(ig) / (4.0_dp * alpha) ) * ( &
                            ! C1
                                    + cexp_ppl(1) &
                            ! C2
                                    + cexp_ppl(2) / (2.0_dp * alpha) * &
                                    ( 3.0_dp - v_coulomb_gspace%pw%pw_grid%gsq(ig)/(2.0_dp * alpha) ) &
                            ! C3
                                    + cexp_ppl(3) / (2.0_dp * alpha)**2 * &
                                    ( 15.0_dp - 10.0_dp * v_coulomb_gspace%pw%pw_grid%gsq(ig)/(2.0_dp * alpha) &
                                                      + ( v_coulomb_gspace%pw%pw_grid%gsq(ig)/(2.0_dp * alpha) )**2 ) &
                            ! C4
                                    + cexp_ppl(4) / (2.0_dp * alpha)**3 * &
                                    ( 105.0_dp - 105.0_dp * v_coulomb_gspace%pw%pw_grid%gsq(ig)/(2.0_dp * alpha) &
                                                      + 21.0_dp * ( v_coulomb_gspace%pw%pw_grid%gsq(ig)/(2.0_dp * alpha) )**2 &
                                                      - ( v_coulomb_gspace%pw%pw_grid%gsq(ig)/(2.0_dp * alpha) )**3 ) &
                                    )
                      END SELECT

                   END IF

                   arg = DOT_PRODUCT(v_coulomb_gspace%pw%pw_grid%g(:,ig),ratom)

                   gtemp = gtemp * CMPLX(COS(arg),-SIN(arg),KIND=dp)
                   v_coulomb_gspace%pw%cc(ig) = v_coulomb_gspace%pw%cc(ig) + gtemp
                END DO
                IF ( v_coulomb_gspace%pw%pw_grid%have_g0 ) v_coulomb_gspace%pw%cc(1) = 0.0_dp

             END DO

             CALL pw_copy(v_coulomb_gspace%pw, v_coulomb_gtemp%pw,error=error)
             CALL pw_derive(v_coulomb_gtemp%pw, (/1,0,0/),error=error )
             CALL pw_transfer(v_coulomb_gtemp%pw, v_coulomb_rtemp%pw,error=error)
             CALL pw_axpy(v_coulomb_rtemp%pw, pwx,error=error)

             CALL pw_copy(v_coulomb_gspace%pw, v_coulomb_gtemp%pw,error=error)
             CALL pw_derive(v_coulomb_gtemp%pw, (/0,1,0/),error=error )
             CALL pw_transfer(v_coulomb_gtemp%pw, v_coulomb_rtemp%pw,error=error)
             CALL pw_axpy(v_coulomb_rtemp%pw, pwy,error=error)

             CALL pw_copy(v_coulomb_gspace%pw, v_coulomb_gtemp%pw,error=error)
             CALL pw_derive(v_coulomb_gtemp%pw, (/0,0,1/),error=error )
             CALL pw_transfer(v_coulomb_gtemp%pw, v_coulomb_rtemp%pw,error=error)
             CALL pw_axpy(v_coulomb_rtemp%pw, pwz,error=error)

             CALL pw_pool_give_back_pw(auxbas_pw_pool,v_coulomb_gspace%pw,error=error)
             CALL pw_pool_give_back_pw(auxbas_pw_pool,v_coulomb_gtemp%pw,error=error)
             CALL pw_pool_give_back_pw(auxbas_pw_pool,v_coulomb_rtemp%pw,error=error)

             DEALLOCATE(v_coulomb_gspace,STAT=istat)
             IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                              "v_coulomb_gspace")
             DEALLOCATE(v_coulomb_gtemp,STAT=istat)
             IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                              "v_coulomb_gtemp")
             DEALLOCATE(v_coulomb_rtemp,STAT=istat)
             IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                              "v_coulomb_rtemp")
          ELSE

             ! Attic of the atomic parallellisation
             !
             ! bo(2)
             ! bo = get_limit(natom, para_env%num_pe, para_env%mepos)
             ! DO iat =  bo(1),bo(2) ! natom = # atoms for ikind
             ! DO ix = lbound(pwx%cr3d,1), ubound(pwx%cr3d,1)
             ! DO iy = lbound(pwx%cr3d,2), ubound(pwx%cr3d,2)
             ! DO iz = lbound(pwx%cr3d,3), ubound(pwx%cr3d,3)

             bo = pwx%pw_grid%bounds_local

             DO iat =  1, natom ! natom = # atoms for ikind

                iatom = atom_list(iat)
                ratom = particle_set(iatom)%r

                DO ix = bo(1,1),bo(2,1)
                   DO iy = bo(1,2),bo(2,2)
                      DO iz = bo(1,3),bo(2,3)
                          rpoint = (/REAL(ix,dp)*pwx%pw_grid%dr(1),&
                                  REAL(iy,dp)*pwx%pw_grid%dr(2),&
                                  REAL(iz,dp)*pwx%pw_grid%dr(3)/)
                          rpoint = rpoint + roffset
                          rap = rpoint - ratom
                          rap(1)=MODULO(rap(1),cell%hmat(1,1))-cell%hmat(1,1)/2._dp
                          rap(2)=MODULO(rap(2),cell%hmat(2,2))-cell%hmat(2,2)/2._dp
                          rap(3)=MODULO(rap(3),cell%hmat(3,3))-cell%hmat(3,3)/2._dp
                          sqrt_rap = SQRT(DOT_PRODUCT(rap,rap))
                          exp_rap = EXP( - alpha * sqrt_rap**2 )
                          IF (sqrt_rap < 1.e-10_dp ) sqrt_rap = 1.e-10_dp
                          ! d_x

                          pwx%cr3d(ix,iy,iz) = pwx%cr3d(ix,iy,iz) + charge * ( &
                                    - 2.0_dp * sqrt_alpha * EXP( - sqrt_rap**2 * sqrt_alpha**2 ) * rap(1) &
                                    / ( rootpi * sqrt_rap**2 ) &
                                    + erf( sqrt_rap * sqrt_alpha ) * rap(1) &
                                    / sqrt_rap**3 )

                          ! d_y

                          pwy%cr3d(ix,iy,iz) = pwy%cr3d(ix,iy,iz) + charge * ( &
                                    - 2.0_dp * sqrt_alpha * EXP( - sqrt_rap**2 * sqrt_alpha**2 ) * rap(2) &
                                    / ( rootpi * sqrt_rap**2 ) &
                                    + erf( sqrt_rap * sqrt_alpha ) * rap(2) &
                                    / sqrt_rap**3 )

                          ! d_z

                          pwz%cr3d(ix,iy,iz) = pwz%cr3d(ix,iy,iz) + charge * ( &
                                    - 2.0_dp * sqrt_alpha * EXP( - sqrt_rap**2 * sqrt_alpha**2 ) * rap(3) &
                                    / ( rootpi * sqrt_rap**2 ) &
                                    + erf( sqrt_rap * sqrt_alpha ) * rap(3) &
                                    / sqrt_rap**3 )

                          IF (make_soft) CYCLE

                          ! d_x

                          DO iexp = 1,nexp_ppl
                             pwx%cr3d(ix,iy,iz) = pwx%cr3d(ix,iy,iz) + ( &
                                    - 2.0_dp * alpha * rap(1) * exp_rap * &
                                    cexp_ppl(iexp) * ( sqrt_rap**2 )**(iexp - 1) )
                             IF (iexp > 1) THEN
                                pwx%cr3d(ix,iy,iz) = pwx%cr3d(ix,iy,iz) + ( &
                                2.0_dp * exp_rap * cexp_ppl(iexp) * &
                                ( sqrt_rap**2 )**(iexp - 2) * REAL(iexp - 1,dp) * rap(1) )
                             END IF
                          END DO

                          ! d_y

                          DO iexp = 1,nexp_ppl
                             pwy%cr3d(ix,iy,iz) = pwy%cr3d(ix,iy,iz) + ( &
                                    - 2.0_dp * alpha * rap(2) * exp_rap * &
                                    cexp_ppl(iexp) * ( sqrt_rap**2 )**(iexp - 1) )
                             IF (iexp > 1) THEN
                                pwy%cr3d(ix,iy,iz) = pwy%cr3d(ix,iy,iz) + ( &
                                2.0_dp * exp_rap * cexp_ppl(iexp) * &
                                ( sqrt_rap**2 )**(iexp - 2) * REAL(iexp - 1,dp) * rap(2) )
                             END IF
                          END DO

                          ! d_z

                          DO iexp = 1,nexp_ppl
                             pwz%cr3d(ix,iy,iz) = pwz%cr3d(ix,iy,iz) + ( &
                                    - 2.0_dp * alpha * rap(3) * exp_rap * &
                                    cexp_ppl(iexp) * ( sqrt_rap**2 )**(iexp - 1) )
                             IF (iexp > 1) THEN
                                pwz%cr3d(ix,iy,iz) = pwz%cr3d(ix,iy,iz) + ( &
                                2.0_dp * exp_rap * cexp_ppl(iexp) * &
                                ( sqrt_rap**2 )**(iexp - 2) * REAL(iexp - 1,dp) * rap(3) )
                             END IF
                          END DO

                      END DO ! iz
                   END DO ! iy
                END DO ! ix
             END DO ! iat
          END IF ! gth_gspace

          !   -------------------------------------
          !   Atomic grids part
          !   -------------------------------------

          IF (gapw .AND. paw_atom) THEN

             bo_atom = get_limit(natom, para_env%num_pe, para_env%mepos)

             DO iat = bo_atom(1),bo_atom(2) ! natom = # atoms for ikind

                iatom = atom_list(iat)

                nablavks_vec_rad_h => nablavks_atom_set(iatom)%nablavks_vec_rad_h
                nablavks_vec_rad_s => nablavks_atom_set(iatom)%nablavks_vec_rad_s

                DO ir = 1,grid_atom%nr

                   IF (grid_atom%rad(ir) >= atomic_kind%hard_radius) CYCLE

                   exp_rap = EXP( - alpha * grid_atom%rad(ir)**2 )

                   DO ia = 1,grid_atom%ng_sphere

                      DO idir = 1,3
                         hard_value = 0.0_dp
                         hard_value = charge * ( &
                                    - 2.0_dp * sqrt_alpha * EXP( - grid_atom%rad(ir)**2 * sqrt_alpha**2 ) &
                                    * grid_atom%rad(ir)*harmonics%a(idir,ia) &
                                    / ( rootpi * grid_atom%rad(ir)**2 ) &
                                    + erf( grid_atom%rad(ir) * sqrt_alpha ) &
                                    * grid_atom%rad(ir)*harmonics%a(idir,ia) &
                                    / grid_atom%rad(ir)**3 )
                         soft_value = hard_value
                         DO iexp = 1,nexp_ppl
                            hard_value = hard_value + ( &
                                    - 2.0_dp * alpha * grid_atom%rad(ir)*harmonics%a(idir,ia) &
                                    * exp_rap * cexp_ppl(iexp) * ( grid_atom%rad(ir)**2 )**(iexp - 1) )
                            IF (iexp > 1) THEN
                               hard_value = hard_value + ( &
                                    2.0_dp * exp_rap * cexp_ppl(iexp) &
                                    * ( grid_atom%rad(ir)**2 )**(iexp - 2) * REAL(iexp - 1,dp) &
                                    * grid_atom%rad(ir)*harmonics%a(idir,ia) )
                            END IF
                         END DO
                         nablavks_vec_rad_h(idir,1)%r_coef(ir,ia) = &
                            nablavks_vec_rad_h(idir,1)%r_coef(ir,ia) + hard_value
                         IF (make_soft) THEN
                            nablavks_vec_rad_s(idir,1)%r_coef(ir,ia) = &
                               nablavks_vec_rad_s(idir,1)%r_coef(ir,ia) + soft_value
                         ELSE
                            nablavks_vec_rad_s(idir,1)%r_coef(ir,ia) = &
                               nablavks_vec_rad_s(idir,1)%r_coef(ir,ia) + hard_value
                         END IF
                      END DO

                   END DO ! ia
                END DO ! ir

                DO ispin = 2,nspins
                   DO idir = 1,3
                      nablavks_vec_rad_h(idir,ispin)%r_coef(:,:) = nablavks_vec_rad_h(idir,1)%r_coef(:,:)
                      nablavks_vec_rad_s(idir,ispin)%r_coef(:,:) = nablavks_vec_rad_s(idir,1)%r_coef(:,:)
                   END DO
                END DO

             END DO

          END IF

       ELSE IF (ASSOCIATED(all_potential)) THEN

          CALL get_potential(potential=all_potential,&
                             alpha_core_charge=alpha,&
                             zeff=charge)

          sqrt_alpha = SQRT(alpha)

          !   -------------------------------------
          !   Atomic grids part
          !   -------------------------------------

          bo_atom = get_limit(natom, para_env%num_pe, para_env%mepos)

          DO iat = bo_atom(1),bo_atom(2) ! natom = # atoms for ikind

             iatom = atom_list(iat)

             nablavks_vec_rad_h => nablavks_atom_set(iatom)%nablavks_vec_rad_h

             DO ir = 1,grid_atom%nr

                IF (grid_atom%rad(ir) >= atomic_kind%hard_radius) CYCLE

                DO ia = 1,grid_atom%ng_sphere

                   DO idir = 1,3
                      hard_value = 0.0_dp
                      hard_value = charge * ( &
                                 2.0_dp * sqrt_alpha * EXP( - grid_atom%rad(ir)**2 * sqrt_alpha**2 ) &
                                 * grid_atom%rad(ir)*harmonics%a(idir,ia) &
                                 / ( rootpi * grid_atom%rad(ir)**2 ) &
                                 + erfc( grid_atom%rad(ir) * sqrt_alpha ) &
                                 * grid_atom%rad(ir)*harmonics%a(idir,ia) &
                                 / grid_atom%rad(ir)**3 )
                      nablavks_vec_rad_h(idir,1)%r_coef(ir,ia) = &
                         nablavks_vec_rad_h(idir,1)%r_coef(ir,ia) + hard_value
                   END DO

                END DO ! ia
             END DO ! ir

             DO ispin = 2,nspins
                DO idir = 1,3
                   nablavks_vec_rad_h(idir,ispin)%r_coef(:,:) = nablavks_vec_rad_h(idir,1)%r_coef(:,:)
                END DO
             END DO

          END DO

       ELSE
          CYCLE
       END IF

    END DO

    DO idir = 1,3
       CALL pw_copy(nablavks_set(idir,1)%rho%rho_r(1)%pw,nablavks_set(idir,2)%rho%rho_r(1)%pw,&
            error=error)
    END DO

!   -------------------------------------
!   Add V_xc potential
!   -------------------------------------

    ALLOCATE(v_xc_gtemp,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "v_xc_gtemp",0)
    ALLOCATE(v_xc_rtemp,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "v_xc_rtemp",0)

    CALL pw_pool_create_pw(auxbas_pw_pool,v_xc_gtemp%pw, &
                           use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,&
                           error=error)
    CALL pw_pool_create_pw(auxbas_pw_pool,v_xc_rtemp%pw,&
                           use_data=REALDATA3D,in_space=REALSPACE,&
                           error=error)

    xc_section => section_vals_get_subs_vals(input, "PROPERTIES%LINRES%EPR%PRINT%G_TENSOR%XC", error=error)

    IF (gapw_xc) THEN
      CALL qs_vxc_create(qs_env=qs_env, rho_struct=qs_env%rho_xc, xc_section=xc_section,&
           vxc_rho=v_rspace_new, vxc_tau=v_tau_rspace, exc=exc, just_energy=.FALSE., & 
           error=error)
    ELSE
      CALL qs_vxc_create(qs_env=qs_env, rho_struct=qs_env%rho, xc_section=xc_section,&
           vxc_rho=v_rspace_new, vxc_tau=v_tau_rspace, exc=exc, just_energy=.FALSE., & 
           error=error)
    END IF

    IF (ASSOCIATED(v_rspace_new)) THEN

       DO ispin = 1,nspins

          CALL pw_transfer(v_rspace_new(ispin)%pw, v_xc_gtemp%pw, error=error)
          CALL pw_derive(v_xc_gtemp%pw, (/1,0,0/) , error=error)
          CALL pw_transfer(v_xc_gtemp%pw, v_xc_rtemp%pw, error=error)
          CALL pw_axpy(v_xc_rtemp%pw, nablavks_set(1,ispin)%rho%rho_r(1)%pw, error=error)

          CALL pw_transfer(v_rspace_new(ispin)%pw, v_xc_gtemp%pw, error=error)
          CALL pw_derive(v_xc_gtemp%pw, (/0,1,0/) , error=error)
          CALL pw_transfer(v_xc_gtemp%pw, v_xc_rtemp%pw, error=error)
          CALL pw_axpy(v_xc_rtemp%pw, nablavks_set(2,ispin)%rho%rho_r(1)%pw, error=error)

          CALL pw_transfer(v_rspace_new(ispin)%pw, v_xc_gtemp%pw, error=error)
          CALL pw_derive(v_xc_gtemp%pw, (/0,0,1/) , error=error)
          CALL pw_transfer(v_xc_gtemp%pw, v_xc_rtemp%pw, error=error)
          CALL pw_axpy(v_xc_rtemp%pw, nablavks_set(3,ispin)%rho%rho_r(1)%pw, error=error)

          CALL pw_pool_give_back_pw(auxbas_pw_pool,v_rspace_new(ispin)%pw,&
               error=error)

       END DO

       DEALLOCATE(v_rspace_new,stat=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                        "v_rspace_new")
    END IF

    IF (ASSOCIATED(v_tau_rspace)) THEN

       DO ispin = 1,nspins

          CALL pw_transfer(v_tau_rspace(ispin)%pw, v_xc_gtemp%pw, error=error)
          CALL pw_derive(v_xc_gtemp%pw, (/1,0,0/) , error=error)
          CALL pw_transfer(v_xc_gtemp%pw, v_xc_rtemp%pw, error=error)
          CALL pw_axpy(v_xc_rtemp%pw, nablavks_set(1,ispin)%rho%rho_r(1)%pw, error=error)

          CALL pw_transfer(v_tau_rspace(ispin)%pw, v_xc_gtemp%pw, error=error)
          CALL pw_derive(v_xc_gtemp%pw, (/0,1,0/) , error=error)
          CALL pw_transfer(v_xc_gtemp%pw, v_xc_rtemp%pw, error=error)
          CALL pw_axpy(v_xc_rtemp%pw, nablavks_set(2,ispin)%rho%rho_r(1)%pw, error=error)

          CALL pw_transfer(v_tau_rspace(ispin)%pw, v_xc_gtemp%pw, error=error)
          CALL pw_derive(v_xc_gtemp%pw, (/0,0,1/) , error=error)
          CALL pw_transfer(v_xc_gtemp%pw, v_xc_rtemp%pw, error=error)
          CALL pw_axpy(v_xc_rtemp%pw, nablavks_set(3,ispin)%rho%rho_r(1)%pw, error=error)

          CALL pw_pool_give_back_pw(auxbas_pw_pool,v_tau_rspace(ispin)%pw,&
               error=error)

       END DO

       DEALLOCATE(v_tau_rspace,stat=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                        "v_tau_rspace")
    END IF

    CALL pw_pool_give_back_pw(auxbas_pw_pool,v_xc_gtemp%pw,error=error)
    CALL pw_pool_give_back_pw(auxbas_pw_pool,v_xc_rtemp%pw,error=error)

    DEALLOCATE(v_xc_gtemp,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "v_xc_gtemp")
    DEALLOCATE(v_xc_rtemp,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "v_xc_rtemp")

    IF (gapw .OR. gapw_xc) THEN
       CALL calculate_vxc_atom(qs_env=qs_env,energy_only=.FALSE.,&
                               gradient_atom_set=nablavks_atom_set,&
                               error=error)
    END IF

!   -------------------------------------
!   Write Nabla V_KS (local) to cubes
!   -------------------------------------

    IF (BTEST(cp_print_key_should_output(logger%iter_info,lr_section,&
                  "EPR%PRINT%NABLAVKS_CUBES",error=error),cp_p_file)) THEN
       ALLOCATE(wf_r,STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                        "wf_r",0)
       CALL pw_pool_create_pw(auxbas_pw_pool,wf_r%pw,&
               use_data = REALDATA3D,&
               in_space = REALSPACE, error=error)
       DO idir = 1,3
          CALL pw_zero(wf_r%pw, error=error)
          CALL pw_copy(nablavks_set(idir,1)%rho%rho_r(1)%pw,wf_r%pw, error=error) ! RA
          filename="nablavks"
          WRITE(ext,'(a2,I1,a5)')  "_d",idir,".cube"
          unit_nr=cp_print_key_unit_nr(logger,lr_section,"EPR%PRINT%NABLAVKS_CUBES",&
                  extension=TRIM(ext),middle_name=TRIM(filename),&
                  log_filename=.FALSE.,file_position="REWIND",error=error)
          CALL cp_pw_to_cube(wf_r%pw,unit_nr,"NABLA V_KS ",&
                  particles=particles,&
                  stride=section_get_ivals(lr_section,&
                  "EPR%PRINT%NABLAVKS_CUBES%STRIDE",error=error),&
                  error=error)
          CALL cp_print_key_finished_output(unit_nr,logger,lr_section,&
                  "EPR%PRINT%NABLAVKS_CUBES",error=error)
       END DO
       CALL pw_pool_give_back_pw(auxbas_pw_pool,wf_r%pw,error=error)
       DEALLOCATE(wf_r,STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                        "wf_r")
    END IF

    CALL cp_print_key_finished_output(output_unit,logger,lr_section,&
         "PRINT%PROGRAM_RUN_INFO",error=error)

  END SUBROUTINE epr_nablavks

END MODULE qs_linres_epr_nablavks

