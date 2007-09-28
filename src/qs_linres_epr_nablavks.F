!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!


!!****** cp2k/qs_linres_epr_nablavks *
!!
!!   NAME
!!     qs_linres_epr_nablavks
!!
!!   FUNCTION
!!     Calculates Nabla V_KS (local part) on the RS grids
!!       PW grid: PSP G-space (F) and R-space (F, nearest-neighbour only)
!!                AE: Ewald sum (NF)
!!       Atomic grids: (NF)
!!
!!       (N)F: (Not) Functional
!!
!!   NOTE
!!
!!   AUTHOR
!!      RD
!!
!!   MODIFICATION HISTORY
!!      created 06-2007 [RD]
!!
!!   SOURCE
!!****

MODULE qs_linres_epr_nablavks
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE cell_types,                      ONLY: cell_type
  USE coefficient_types,               ONLY: coeff_copy,&
                                             coeff_transform_space,&
                                             coeff_type,&
                                             coeff_zero
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE erf_fn,                          ONLY: erf
  USE external_potential_types,        ONLY: all_potential_type,&
                                             get_potential,&
                                             gth_potential_type
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: fourpi,&
                                             rootpi,&
                                             twopi
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: particle_type
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_poisson_methods,              ONLY: pw_poisson_solve
  USE pw_poisson_types,                ONLY: pw_poisson_type
  USE pw_pool_types,                   ONLY: pw_pool_give_back_coeff,&
                                             pw_pool_init_coeff,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_copy,&
                                             pw_derive,&
                                             pw_sumup,&
                                             pw_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_ks_methods,                   ONLY: calc_rho_tot_gspace,&
                                             qs_vxc_create
  USE qs_linres_types,                 ONLY: epr_env_type,&
                                             get_epr_env
  USE qs_rho_types,                    ONLY: qs_rho_p_type,&
                                             qs_rho_type
  USE realspace_grid_cube,             ONLY: rs_pw_to_cube
  USE termination,                     ONLY: stop_memory
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: epr_nablavks

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_linres_epr_nablavks'

CONTAINS

! *****************************************************************************
!!****** cp2k/qs_linres_epr_nablavks/epr_nablavks
!!
!!   NAME
!!     epr_nablavks
!!
!!   FUNCTION
!!     Evaluates Nabla V_KS on the grids
!!
!!   NOTES
!!
!!   AUTHOR
!!     RD
!!
!!   MODIFICATION HISTORY
!!     06.2006 created [RD]
!!
!****************************************************************************

  SUBROUTINE epr_nablavks(epr_env,qs_env,error)

    TYPE(epr_env_type)                       :: epr_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'epr_nablavks', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=80)                        :: ext, filename
    COMPLEX(KIND=dp)                         :: gtemp
    INTEGER :: iat, iatom, idir, iexp, ig, ikind, ispin, istat, ix, iy, iz, &
      natom, nexp_ppl, nkind, nspins, output_unit, unit_nr
    INTEGER, DIMENSION(2, 3)                 :: bo
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: all_potential_present, gapw, &
                                                gapw_xc, gth_gspace, &
                                                gth_potential_present, ionode
    REAL(KIND=dp)                            :: alpha, arg, charge, ehartree, &
                                                exc, exp_rap, sqrt_alpha, &
                                                sqrt_rap
    REAL(KIND=dp), DIMENSION(3)              :: rap, ratom, roffset, rpoint
    REAL(KIND=dp), DIMENSION(:), POINTER     :: cexp_ppl
    TYPE(all_potential_type), POINTER        :: all_potential
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(coeff_type), DIMENSION(:), POINTER  :: v_rspace_new, v_tau_rspace
    TYPE(coeff_type), POINTER :: rho_tot_gspace, v_coulomb_gspace, &
      v_coulomb_gtemp, v_coulomb_rtemp, v_hartree_gspace, v_hartree_gtemp, &
      v_hartree_rtemp, v_xc_gtemp, v_xc_rtemp, wf_r
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(gth_potential_type), POINTER        :: gth_potential
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_poisson_type), POINTER           :: poisson_env
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(pw_type), POINTER                   :: pwx, pwy, pwz
    TYPE(qs_rho_p_type), DIMENSION(:, :), &
      POINTER                                :: nablavks_set
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(section_vals_type), POINTER         :: epr_section, lr_section

    NULLIFY(auxbas_pw_pool)
    NULLIFY(cell)
    NULLIFY(dft_control)
    NULLIFY(epr_section)
    NULLIFY(logger)
    NULLIFY(lr_section)
    NULLIFY(nablavks_set)
    NULLIFY(para_env)
    NULLIFY(particle_set)
    NULLIFY(particles)
    NULLIFY(pw_env)
    NULLIFY(rho)
    NULLIFY(subsys)
    NULLIFY(v_rspace_new)
    NULLIFY(v_tau_rspace)

    logger => cp_error_get_logger(error)
    lr_section => section_vals_get_subs_vals(qs_env%input,"PROPERTIES%LINRES",error=error)
    ionode = logger%para_env%mepos==logger%para_env%source

    output_unit = cp_print_key_unit_nr(logger,lr_section,"PRINT%PROGRAM_RUN_INFO",&
         extension=".linresLog",error=error)

!   -------------------------------------
!   Set defaults
!   -------------------------------------

    gth_gspace = .FALSE.

!   -------------------------------------
!   Get nablavks arrays
!   -------------------------------------

    CALL get_epr_env(epr_env,nablavks_set=nablavks_set,error=error)

!   -------------------------------------
!   Get grids / atom info
!   -------------------------------------

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    cell=cell,&
                    dft_control=dft_control,&
                    para_env=para_env,&
                    particle_set=particle_set,&
                    pw_env=pw_env,&
                    rho=rho,&
                    subsys=subsys,&
                    error=error)

    CALL pw_env_get(pw_env,auxbas_pw_pool=auxbas_pw_pool,&
                    poisson_env=poisson_env,&
                    error=error)

    CALL cp_subsys_get(subsys,particles=particles,error=error)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             gth_potential_present=gth_potential_present,&
                             all_potential_present=all_potential_present)

    IF (output_unit>0) THEN
        WRITE (UNIT=output_unit,FMT="(/,T2,A,T79,L1)") &
           "epr(debug)|NABLAVKS gth_pot", gth_potential_present
        WRITE (UNIT=output_unit,FMT="(T2,A,T79,L1,/)") &
           "epr(debug)|NABLAVKS all_pot", all_potential_present
    END IF

    gapw    = dft_control%qs_control%gapw
    gapw_xc = dft_control%qs_control%gapw_xc
    nspins  = dft_control%nspins

    NULLIFY(pwx, pwy, pwz)
    pwx => nablavks_set(1,1)%rho%rho_r(1)%pw
    pwy => nablavks_set(2,1)%rho%rho_r(1)%pw
    pwz => nablavks_set(3,1)%rho%rho_r(1)%pw
    roffset = -REAL(MODULO(pwx%pw_grid%npts,2),dp)*pwx%pw_grid%dr/2.0_dp

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

    CALL pw_pool_init_coeff(pool=auxbas_pw_pool,coeff=v_hartree_gspace, &
                            use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,&
                            error=error)
    CALL pw_pool_init_coeff(pool=auxbas_pw_pool,coeff=v_hartree_gtemp, &
                            use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,&
                            error=error)
    CALL pw_pool_init_coeff(pool=auxbas_pw_pool,coeff=v_hartree_rtemp,&
                            use_data=REALDATA3D,in_space=REALSPACE,&
                            error=error)
    CALL pw_pool_init_coeff(pool=auxbas_pw_pool,coeff=rho_tot_gspace,&
                            use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,&
                            error=error)

    CALL calc_rho_tot_gspace(rho_tot_gspace,qs_env,rho,&
                             skip_nuclear_density=.TRUE.,error=error)

    CALL pw_poisson_solve(poisson_env,rho_tot_gspace%pw,ehartree,&
                          v_hartree_gspace%pw,error=error)

    CALL coeff_copy(v_hartree_gspace, v_hartree_gtemp)
    CALL pw_derive(v_hartree_gtemp%pw, (/1,0,0/) )
    CALL coeff_transform_space(v_hartree_gtemp, v_hartree_rtemp)
    CALL pw_copy(v_hartree_rtemp%pw, pwx)

    CALL coeff_copy(v_hartree_gspace, v_hartree_gtemp)
    CALL pw_derive(v_hartree_gtemp%pw, (/0,1,0/) )
    CALL coeff_transform_space(v_hartree_gtemp, v_hartree_rtemp)
    CALL pw_copy(v_hartree_rtemp%pw, pwy)

    CALL coeff_copy(v_hartree_gspace, v_hartree_gtemp)
    CALL pw_derive(v_hartree_gtemp%pw, (/0,0,1/) )
    CALL coeff_transform_space(v_hartree_gtemp, v_hartree_rtemp)
    CALL pw_copy(v_hartree_rtemp%pw, pwz)

    CALL pw_pool_give_back_coeff(auxbas_pw_pool,v_hartree_gspace,error=error)
    CALL pw_pool_give_back_coeff(auxbas_pw_pool,v_hartree_gtemp,error=error)
    CALL pw_pool_give_back_coeff(auxbas_pw_pool,v_hartree_rtemp,error=error)
    CALL pw_pool_give_back_coeff(auxbas_pw_pool,rho_tot_gspace,error=error)

    DEALLOCATE(v_hartree_gspace,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                         "v_hartree_gspace",0)
    DEALLOCATE(v_hartree_gtemp,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                         "v_hartree_gtemp",0)
    DEALLOCATE(v_hartree_rtemp,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                         "v_hartree_rtemp",0)
    DEALLOCATE(rho_tot_gspace,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                         "rho_tot_gspace",0)

!   -------------------------------------
!   Add Coulomb potential
!   -------------------------------------

    nkind = SIZE(atomic_kind_set) ! nkind = number of atom types

    DO ikind = 1,nkind ! loop over atom types
       NULLIFY(atom_list, atomic_kind)
       atomic_kind => atomic_kind_set(ikind)
       IF (gth_potential_present) THEN
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
                               gth_potential=gth_potential,&
                               atom_list=atom_list,natom=natom)
          IF (ASSOCIATED(gth_potential)) THEN
              NULLIFY(cexp_ppl)
              CALL get_potential(potential=gth_potential,&
                                 zeff=charge,&
                                 alpha_ppl=alpha,&
                                 nexp_ppl=nexp_ppl,&
                                 cexp_ppl=cexp_ppl)

              sqrt_alpha = SQRT(alpha)

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
                 
                 CALL pw_pool_init_coeff(pool=auxbas_pw_pool,coeff=v_coulomb_gspace, &
                                         use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,&
                                         error=error)
                 CALL pw_pool_init_coeff(pool=auxbas_pw_pool,coeff=v_coulomb_gtemp, &
                                         use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,&
                                         error=error)
                 CALL pw_pool_init_coeff(pool=auxbas_pw_pool,coeff=v_coulomb_rtemp,&
                                         use_data=REALDATA3D,in_space=REALSPACE,&
                                              error=error)
                 
                 CALL coeff_zero(v_coulomb_gspace)

                 DO iat = 1, natom ! natom = # atoms for ikind

                    iatom = atom_list(iat)
                    ratom = particle_set(iatom)%r

                    DO ig = v_coulomb_gspace%pw%pw_grid%first_gne0,v_coulomb_gspace%pw%pw_grid%ngpts_cut_local
                       gtemp = 0.0_dp
                        gtemp = - fourpi * charge / cell%deth * &
                             EXP ( - v_coulomb_gspace%pw%pw_grid%gsq(ig) / (4.0_dp * alpha) ) &
                                / v_coulomb_gspace%pw%pw_grid%gsq(ig)

                       arg = DOT_PRODUCT(v_coulomb_gspace%pw%pw_grid%g(:,ig),ratom)

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

                       gtemp = gtemp * CMPLX(COS(arg),-SIN(arg),KIND=dp)
                       v_coulomb_gspace%pw%cc(ig) = v_coulomb_gspace%pw%cc(ig) + gtemp 
                    END DO
                    IF ( v_coulomb_gspace%pw%pw_grid%have_g0 ) v_coulomb_gspace%pw%cc(1) = 0.0_dp

                 END DO

                 CALL coeff_copy(v_coulomb_gspace, v_coulomb_gtemp)
                 CALL pw_derive(v_coulomb_gtemp%pw, (/1,0,0/) )
                 CALL coeff_transform_space(v_coulomb_gtemp, v_coulomb_rtemp)
                 CALL pw_sumup(v_coulomb_rtemp%pw, pwx)
             
                 CALL coeff_copy(v_coulomb_gspace, v_coulomb_gtemp)
                 CALL pw_derive(v_coulomb_gtemp%pw, (/0,1,0/) )
                 CALL coeff_transform_space(v_coulomb_gtemp, v_coulomb_rtemp)
                 CALL pw_sumup(v_coulomb_rtemp%pw, pwy)
             
                 CALL coeff_copy(v_coulomb_gspace, v_coulomb_gtemp)
                 CALL pw_derive(v_coulomb_gtemp%pw, (/0,0,1/) )
                 CALL coeff_transform_space(v_coulomb_gtemp, v_coulomb_rtemp)
                 CALL pw_sumup(v_coulomb_rtemp%pw, pwz)
             
                 CALL pw_pool_give_back_coeff(auxbas_pw_pool,v_coulomb_gspace,error=error)
                 CALL pw_pool_give_back_coeff(auxbas_pw_pool,v_coulomb_gtemp,error=error)
                 CALL pw_pool_give_back_coeff(auxbas_pw_pool,v_coulomb_rtemp,error=error)

                 DEALLOCATE(v_coulomb_gspace,STAT=istat)
                 IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                                      "v_coulomb_gspace",0)
                 DEALLOCATE(v_coulomb_gtemp,STAT=istat)
                 IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                                      "v_coulomb_gtemp",0)
                 DEALLOCATE(v_coulomb_rtemp,STAT=istat)
                 IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                                      "v_coulomb_rtemp",0)
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
    
                                    pwy%cr3d(ix,iy,iz) = pwy%cr3d(ix,iy,iz) + charge * ( &
                                            - 2.0_dp * sqrt_alpha * EXP( - sqrt_rap**2 * sqrt_alpha**2 ) * rap(2) &
                                            / ( rootpi * sqrt_rap**2 ) &
                                            + erf( sqrt_rap * sqrt_alpha ) * rap(2) &
                                            / sqrt_rap**3 )
                              
                                 DO iexp = 1,nexp_ppl
                                    pwy%cr3d(ix,iy,iz) = pwy%cr3d(ix,iy,iz) + ( &
                                           - 2.0_dp * alpha * rap(2) * exp_rap * cexp_ppl(iexp) *&
                                           ( sqrt_rap**2 )**(iexp - 1) )
                                    IF (iexp > 1) THEN
                                       pwy%cr3d(ix,iy,iz) = pwy%cr3d(ix,iy,iz) + ( &
                                       2.0_dp * exp_rap * cexp_ppl(iexp) * ( sqrt_rap**2 )**(iexp - 2) *&
                                       REAL(iexp - 1,dp) * rap(2) )
                                    END IF
                                 END DO
   
                              ! d_z
    
                                    pwz%cr3d(ix,iy,iz) = pwz%cr3d(ix,iy,iz) + charge * ( &
                                            - 2.0_dp * sqrt_alpha * EXP( - sqrt_rap**2 * sqrt_alpha**2 ) * rap(3) &
                                            / ( rootpi * sqrt_rap**2 ) &
                                            + erf( sqrt_rap * sqrt_alpha ) * rap(3) &
                                            / sqrt_rap**3 )
   
                                 DO iexp = 1,nexp_ppl
                                    pwz%cr3d(ix,iy,iz) = pwz%cr3d(ix,iy,iz) + ( &
                                           - 2.0_dp * alpha * rap(3) * exp_rap * cexp_ppl(iexp) *&
                                           ( sqrt_rap**2 )**(iexp - 1) )
                                    IF (iexp > 1) THEN
                                       pwz%cr3d(ix,iy,iz) = pwz%cr3d(ix,iy,iz) + ( &
                                       2.0_dp * exp_rap * cexp_ppl(iexp) * ( sqrt_rap**2 )**(iexp - 2) *&
                                       REAL(iexp - 1,dp) * rap(3) )
                                    END IF
                                 END DO
   
                             END DO ! iz
                          END DO ! iy
                       END DO ! ix
   
                 END DO ! iat
              END IF ! gth_gspace
                 
          END IF ! ASSOCIATED(gth_potential)

       ELSE IF (all_potential_present) THEN
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
                               all_potential=all_potential,&
                               atom_list=atom_list,natom=natom)
          IF (ASSOCIATED(all_potential)) THEN
              CALL get_potential(potential=all_potential,&
                                 zeff=charge)
          END IF
       ELSE
          CYCLE
       END IF

    END DO

    DO idir = 1,3
       CALL coeff_copy(nablavks_set(idir,1)%rho%rho_r(1),nablavks_set(idir,2)%rho%rho_r(1))
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

    CALL pw_pool_init_coeff(pool=auxbas_pw_pool,coeff=v_xc_gtemp, &
                            use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,&
                            error=error)
    CALL pw_pool_init_coeff(pool=auxbas_pw_pool,coeff=v_xc_rtemp,&
                            use_data=REALDATA3D,in_space=REALSPACE,&
                            error=error)

    CALL qs_vxc_create(qs_env=qs_env,&
         vxc_rho=v_rspace_new, vxc_tau=v_tau_rspace, exc=exc, &
         just_energy=.FALSE., gapw_xc=gapw_xc, error=error)

    IF (ASSOCIATED(v_rspace_new)) THEN

       DO ispin = 1,nspins

          CALL coeff_transform_space(v_rspace_new(ispin), v_xc_gtemp)
          CALL pw_derive(v_xc_gtemp%pw, (/1,0,0/) )
          CALL coeff_transform_space(v_xc_gtemp, v_xc_rtemp)
          CALL pw_sumup(v_xc_rtemp%pw, nablavks_set(1,ispin)%rho%rho_r(1)%pw)

          CALL coeff_transform_space(v_rspace_new(ispin), v_xc_gtemp)
          CALL pw_derive(v_xc_gtemp%pw, (/0,1,0/) )
          CALL coeff_transform_space(v_xc_gtemp, v_xc_rtemp)
          CALL pw_sumup(v_xc_rtemp%pw, nablavks_set(2,ispin)%rho%rho_r(1)%pw)

          CALL coeff_transform_space(v_rspace_new(ispin), v_xc_gtemp)
          CALL pw_derive(v_xc_gtemp%pw, (/0,0,1/) )
          CALL coeff_transform_space(v_xc_gtemp, v_xc_rtemp)
          CALL pw_sumup(v_xc_rtemp%pw, nablavks_set(3,ispin)%rho%rho_r(1)%pw)

          CALL pw_pool_give_back_coeff(auxbas_pw_pool,v_rspace_new(ispin),&
               error=error)

       END DO

       DEALLOCATE(v_rspace_new,stat=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                            "v_rspace_new",0)
    END IF

    IF (ASSOCIATED(v_tau_rspace)) THEN

       DO ispin = 1,nspins

          CALL coeff_transform_space(v_tau_rspace(ispin), v_xc_gtemp)
          CALL pw_derive(v_xc_gtemp%pw, (/1,0,0/) )
          CALL coeff_transform_space(v_xc_gtemp, v_xc_rtemp)
          CALL pw_sumup(v_xc_rtemp%pw, nablavks_set(1,ispin)%rho%rho_r(1)%pw)

          CALL coeff_transform_space(v_tau_rspace(ispin), v_xc_gtemp)
          CALL pw_derive(v_xc_gtemp%pw, (/0,1,0/) )
          CALL coeff_transform_space(v_xc_gtemp, v_xc_rtemp)
          CALL pw_sumup(v_xc_rtemp%pw, nablavks_set(2,ispin)%rho%rho_r(1)%pw)

          CALL coeff_transform_space(v_tau_rspace(ispin), v_xc_gtemp)
          CALL pw_derive(v_xc_gtemp%pw, (/0,0,1/) )
          CALL coeff_transform_space(v_xc_gtemp, v_xc_rtemp)
          CALL pw_sumup(v_xc_rtemp%pw, nablavks_set(3,ispin)%rho%rho_r(1)%pw)

          CALL pw_pool_give_back_coeff(auxbas_pw_pool,v_tau_rspace(ispin),&
               error=error)

       END DO

       DEALLOCATE(v_tau_rspace,stat=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                            "v_tau_rspace",0)
    END IF

    CALL pw_pool_give_back_coeff(auxbas_pw_pool,v_xc_gtemp,error=error)
    CALL pw_pool_give_back_coeff(auxbas_pw_pool,v_xc_rtemp,error=error)

    DEALLOCATE(v_xc_gtemp,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                         "v_xc_gtemp",0)
    DEALLOCATE(v_xc_rtemp,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                         "v_xc_rtemp",0)

!   -------------------------------------
!   Write Nabla V_KS (local) to cubes
!   -------------------------------------

    epr_section => section_vals_get_subs_vals(qs_env%input,"PROPERTIES%LINRES%EPR",&
                   error=error)
    IF (BTEST(cp_print_key_should_output(logger%iter_info,epr_section,&
                  "PRINT%NABLAVKS_CUBES",error=error),cp_p_file)) THEN
       ALLOCATE(wf_r,STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                            "wf_r",0)
       CALL pw_pool_init_coeff(auxbas_pw_pool,wf_r,&
               use_data = REALDATA3D,&
               in_space = REALSPACE, error=error)
       DO idir = 1,3
          CALL coeff_zero(wf_r)
          CALL coeff_copy(nablavks_set(idir,1)%rho%rho_r(1),wf_r) ! RA
          filename="nablavks"
          WRITE(ext,'(a2,I1,a5)')  "_d",idir,".cube"
          unit_nr=cp_print_key_unit_nr(logger,epr_section,"PRINT%NABLAVKS_CUBES",&
                  extension=TRIM(ext),middle_name=TRIM(filename),&
                  log_filename=.FALSE.,file_position="REWIND",error=error)
          CALL rs_pw_to_cube(wf_r%pw,unit_nr,ionode,"NABLA V_KS ",&
                  particles=particles,&
                  stride=1,&
                  error=error)
          CALL cp_print_key_finished_output(unit_nr,logger,epr_section,&
                  "PRINT%NABLAVKS_CUBES",error=error)
       END DO
       CALL pw_pool_give_back_coeff(auxbas_pw_pool,wf_r,error=error)
       DEALLOCATE(wf_r,STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                            "wf_r",0)
    END IF

    CALL cp_print_key_finished_output(output_unit,logger,lr_section,&
         "PRINT%PROGRAM_RUN_INFO",error=error)

  END SUBROUTINE epr_nablavks

END MODULE  qs_linres_epr_nablavks

