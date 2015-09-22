!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Type definitiona for linear response calculations
!> \author MI
! *****************************************************************************
MODULE qs_linres_types
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cp_array_utils_i,                ONLY: cp_2d_i_p_type
  USE cp_array_utils_r,                ONLY: cp_2d_r_p_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_p_type
  USE cp_fm_types,                     ONLY: cp_fm_p_type
  USE kinds,                           ONLY: dp
  USE qs_grid_atom,                    ONLY: grid_atom_type
  USE qs_harmonics_atom,               ONLY: harmonics_atom_type
  USE qs_kind_types,                   ONLY: get_qs_kind,&
                                             qs_kind_type
  USE qs_loc_types,                    ONLY: qs_loc_env_new_type,&
                                             qs_loc_env_release
  USE qs_rho_atom_types,               ONLY: rho_atom_coeff,&
                                             rho_atom_type
  USE qs_rho_types,                    ONLY: qs_rho_p_type,&
                                             qs_rho_release
  USE realspace_grid_types,            ONLY: realspace_grid_p_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE


  PRIVATE

!****s* qs_linres_types/linres_control_type

! *****************************************************************************
!> \brief General settings for linear response calculations
!> \param property which quantity is to be calculated by LR
!> \param opt_method method to optimize the psi1 by minimization of the second order term of the energy
!> \param preconditioner which kind of preconditioner should be used, if any
!> \param localized_psi 0 : don't use the canonical psi0, but the maximally localized wavefunctions
!> \param do_kernel the kernel is zero if the rho1 is zero as for the magnetic field perturbation
!> \param tolerance convergence criterium for the optimization of the psi1
!> \author MI
! *****************************************************************************
  TYPE linres_control_type
     INTEGER                                   :: ref_count
     INTEGER                                   :: property
     INTEGER                                   :: preconditioner_type
     INTEGER                                   :: restart_every
     REAL(dp)                                  :: energy_gap
     INTEGER                                   :: max_iter
     LOGICAL                                   :: localized_psi0
     LOGICAL                                   :: do_kernel
     LOGICAL                                   :: converged
     LOGICAL                                   :: linres_restart
     LOGICAL                                   :: lr_triplet
     REAL(KIND=dp)                             :: eps
     TYPE(qs_loc_env_new_type), POINTER        :: qs_loc_env
     CHARACTER(LEN=8)                          :: flag
  END TYPE linres_control_type

!****s* qs_linres_types/current_env_type

! *****************************************************************************
!> \param ref_coun t
!> \param full_nmr true if the full correction is calculated
!> \param simplenmr_done , fullnmr_done : flags that indicate what has been
!>                    already calculated: used for restart
!> \param centers_set centers of the maximally localized psi0
!> \param spreads_set spreads of the maximally localized psi0
!> \param p_psi 0      : full matrixes, operator p applied to psi0
!> \param rxp_psi 0    : full matrixes, operator (r-d)xp applied to psi0
!> \param psi 1_p      : response wavefunctions to the perturbation given by
!>                    H1=p (xyz)  applied to psi0
!> \param psi 1_rxp    : response wavefunctions to the perturbation given by
!>                    H1=(r-d_i)xp applied to psi0_i where d_i is the center
!> \param psi 1_D      : response wavefunctions to the perturbation given by
!>                    H1=(d_j-d_i)xp applied to psi0_i where d_i is the center
!>                    and d_j is the center of psi0_j and psi1_D_j is the result
!>                    This operator has to be used in nstate scf calculations,
!>                    one for each psi1_D_j vector
!> \param chemical_shift the tensor for each atom
!> \param chi_tensor the susceptibility tensor
!> \param jrho 1_set   : current density on the global grid, if gapw this is only the soft part
!> \param jrho 1_atom_set : current density on the local atomic grids (only if gapw)
!> \author MI
! *****************************************************************************
  TYPE realspaces_grid_p_type
     TYPE(realspace_grid_p_type), DIMENSION(:), POINTER :: rs
  END TYPE realspaces_grid_p_type

  TYPE current_env_type
     LOGICAL                                     :: full,simple_done(6),simple_converged(6),do_qmmm
     LOGICAL                                     :: use_old_gauge_atom,chi_pbc,do_selected_states,&
                                                    gauge_init, all_pert_op_done 
     LOGICAL, DIMENSION(:,:), POINTER            :: full_done
     INTEGER                                     :: ref_count,nao,nstates(2),gauge,orb_center,nbr_center(2)
     INTEGER, DIMENSION(:    ), POINTER          :: list_cubes,selected_states_on_atom_list
     INTEGER, DIMENSION(:,:,:), POINTER          :: statetrueindex
     CHARACTER(LEN=30)                           :: gauge_name,orb_center_name
     REAL(dp)                                    :: chi_tensor(3,3,2),chi_tensor_loc(3,3,2),gauge_atom_radius
     REAL(dp)                                    :: selected_states_atom_radius
     REAL(dp), DIMENSION(:,:), POINTER           :: basisfun_center
     TYPE(cp_2d_i_p_type),    DIMENSION(:  ), POINTER :: center_list
     TYPE(cp_2d_r_p_type),    DIMENSION(:  ), POINTER :: centers_set
     TYPE(cp_fm_p_type),      DIMENSION(:,:), POINTER :: psi1_p,psi1_rxp,psi1_D,p_psi0,rxp_psi0
     TYPE(jrho_atom_type),    DIMENSION(:  ), POINTER :: jrho1_atom_set
     TYPE(qs_rho_p_type),     DIMENSION(:  ), POINTER :: jrho1_set
     TYPE(realspace_grid_p_type),DIMENSION(:),POINTER :: rs_buf
     TYPE(realspaces_grid_p_type), DIMENSION(:), POINTER :: rs_gauge
     !
     TYPE(cp_fm_p_type)  , DIMENSION(:), POINTER :: psi0_order
  END TYPE current_env_type

!*********************************************************************************************************
! \param type for polarizability calculation using Berry operator

  TYPE polar_env_type
     INTEGER                             :: ref_count
     LOGICAL                             :: do_raman
     REAL(dp), DIMENSION(:,:), POINTER :: polar
     TYPE(cp_fm_p_type), DIMENSION(:,:), POINTER :: psi1_dBerry,dBerry_psi0,mo_derivs
  END TYPE polar_env_type
! ************************************************************************************

  TYPE issc_env_type
     INTEGER                             :: ref_count
     INTEGER                             :: issc_natms
     INTEGER, DIMENSION(:), POINTER      :: issc_on_atom_list
     LOGICAL                             :: interpolate_issc
     LOGICAL                             :: do_fc,do_sd,do_pso,do_dso
     REAL(dp)                            :: issc_gapw_radius,issc_factor,issc_factor_gapw
     REAL(dp), DIMENSION(:,:,:,:,:), POINTER :: issc,issc_loc
     TYPE(cp_fm_p_type),      DIMENSION(:,:), POINTER :: psi1_efg,psi1_pso,efg_psi0,pso_psi0,dso_psi0,psi1_dso!last two not needed
     TYPE(cp_fm_p_type), DIMENSION(:), POINTER :: psi1_fc
     TYPE(cp_fm_p_type), DIMENSION(:), POINTER :: fc_psi0
     TYPE(cp_dbcsr_p_type),DIMENSION(:  ), POINTER :: matrix_efg,matrix_pso,matrix_dso,matrix_fc
  END TYPE issc_env_type

! ************************************************************************************
  TYPE nmr_env_type
     INTEGER                             :: ref_count, n_nics
     INTEGER, DIMENSION(:), POINTER      :: cs_atom_list
     INTEGER, DIMENSION(:), POINTER      :: do_calc_cs_atom
     LOGICAL                             :: do_nics,interpolate_shift
     REAL(dp)                            :: shift_gapw_radius,shift_factor, shift_factor_gapw, chi_factor, &
                                            chi_SI2shiftppm, chi_SI2ppmcgs
     REAL(dp), DIMENSION(:,:  ), POINTER :: r_nics
     REAL(dp), DIMENSION(:,:,:), POINTER :: chemical_shift, chemical_shift_loc, &
                                            chemical_shift_nics_loc, chemical_shift_nics
  END TYPE nmr_env_type

! *****************************************************************************
  TYPE epr_env_type
     INTEGER                                     :: ref_count
     REAL(dp)                                    :: g_free_factor, g_soo_chicorr_factor, g_soo_factor,&
                                                    g_so_factor, g_so_factor_gapw, g_zke_factor, g_zke
     REAL(dp), DIMENSION(:,:), POINTER           :: g_total, g_so, g_soo
     TYPE(qs_rho_p_type),      DIMENSION(:,:), POINTER :: nablavks_set
     TYPE(nablavks_atom_type), DIMENSION(:  ), POINTER :: nablavks_atom_set
     TYPE(qs_rho_p_type),      DIMENSION(:,:), POINTER :: bind_set
     TYPE(rho_atom_coeff),     DIMENSION(:,:), POINTER :: bind_atom_set
     TYPE(rho_atom_type),      DIMENSION(:  ), POINTER :: vks_atom_set
  END TYPE epr_env_type

! *****************************************************************************
  TYPE nablavks_atom_type
     TYPE(rho_atom_coeff), DIMENSION(:,:),&
          POINTER                                :: nablavks_vec_rad_h,&
                                                    nablavks_vec_rad_s
  END TYPE nablavks_atom_type

! *****************************************************************************
  TYPE jrho_atom_p_type
     TYPE(jrho_atom_type), POINTER       :: jrho_atom
  END TYPE jrho_atom_p_type

! *****************************************************************************
  TYPE jrho_atom_type
    TYPE(rho_atom_coeff), DIMENSION(:),&
                               POINTER      :: cjc_h, cjc_s, cjc0_h, cjc0_s
    TYPE(rho_atom_coeff), DIMENSION(:),&
                               POINTER      :: cjc_ii_h, cjc_ii_s
    TYPE(rho_atom_coeff), DIMENSION(:),&
                               POINTER      :: cjc_iii_h, cjc_iii_s
    TYPE(rho_atom_coeff), DIMENSION(:,:),&
                               POINTER      :: jrho_vec_rad_h,&
                                               jrho_vec_rad_s
    TYPE(rho_atom_coeff), DIMENSION(:), &
      POINTER                               :: jrho_h, jrho_s
    TYPE(rho_atom_coeff), DIMENSION(:), &
      POINTER                               :: jrho_a_h, jrho_a_s
    TYPE(rho_atom_coeff), DIMENSION(:), &
      POINTER                               :: jrho_b_h, jrho_b_s
    TYPE(rho_atom_coeff), DIMENSION(:), &
      POINTER                               :: jrho_a_h_ii, jrho_a_s_ii
    TYPE(rho_atom_coeff), DIMENSION(:), &
      POINTER                               :: jrho_b_h_ii, jrho_b_s_ii
    TYPE(rho_atom_coeff), DIMENSION(:), &
      POINTER                               :: jrho_a_h_iii, jrho_a_s_iii
    TYPE(rho_atom_coeff), DIMENSION(:), &
      POINTER                               :: jrho_b_h_iii, jrho_b_s_iii
  END TYPE jrho_atom_type

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_linres_types'

! *** Public data types ***

  PUBLIC :: linres_control_type,&
            nmr_env_type, issc_env_type, jrho_atom_type,&
            epr_env_type,&
            nablavks_atom_type, current_env_type,&
            realspaces_grid_p_type, polar_env_type

! *** Public subroutines ***

  PUBLIC :: allocate_jrho_atom_rad, deallocate_jrho_atom_set, get_nmr_env,&
            get_current_env, allocate_jrho_coeff,&
            init_jrho_atom_set, init_nablavks_atom_set, linres_control_create,&
            linres_control_retain, linres_control_release,&
            set_epr_env, deallocate_nablavks_atom_set,&
            set2zero_jrho_atom_rad, get_epr_env,&
            nmr_env_create, epr_env_create, current_env_create, get_issc_env,&
            set_current_env, issc_env_create, get_polar_env,&
            polar_env_create

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param linres_control ...
! *****************************************************************************
  SUBROUTINE linres_control_create(linres_control)

    TYPE(linres_control_type), POINTER       :: linres_control

    CHARACTER(len=*), PARAMETER :: routineN = 'linres_control_create', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

     failure =.FALSE.

     CPASSERT(.NOT.ASSOCIATED(linres_control))
     ALLOCATE (linres_control)
     linres_control%ref_count=1
     NULLIFY(linres_control%qs_loc_env)
     linres_control%property            = HUGE(0)!is that used?
     linres_control%preconditioner_type = HUGE(0)
     linres_control%restart_every       = HUGE(0)
     linres_control%energy_gap          = HUGE(0.0_dp)
     linres_control%max_iter            = HUGE(0)
     linres_control%localized_psi0      = .FALSE.
     linres_control%converged           = .FALSE.
     linres_control%linres_restart      = .FALSE.
     linres_control%eps                 = HUGE(0.0_dp)
     linres_control%flag                = ""
     linres_control%do_kernel           = .FALSE.
     linres_control%lr_triplet          = .FALSE.

  END SUBROUTINE linres_control_create

! *****************************************************************************
!> \brief ...
!> \param linres_control ...
! *****************************************************************************
  SUBROUTINE linres_control_release(linres_control)

    TYPE(linres_control_type), POINTER       :: linres_control

    CHARACTER(len=*), PARAMETER :: routineN = 'linres_control_release', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

     failure =.FALSE.

     IF (ASSOCIATED(linres_control)) THEN
       CPASSERT(linres_control%ref_count>0)
       linres_control%ref_count=linres_control%ref_count-1
       IF(linres_control%ref_count<1)THEN
         IF(ASSOCIATED(linres_control%qs_loc_env)) THEN
           CALL qs_loc_env_release(linres_control%qs_loc_env)
         END IF
         DEALLOCATE(linres_control)
       END IF
     END IF
     NULLIFY(linres_control)
  END SUBROUTINE linres_control_release

! *****************************************************************************
!> \brief ...
!> \param linres_control ...
! *****************************************************************************
  SUBROUTINE linres_control_retain(linres_control)

    TYPE(linres_control_type), POINTER       :: linres_control

    CHARACTER(len=*), PARAMETER :: routineN = 'linres_control_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

     failure =.FALSE.

     CPASSERT(ASSOCIATED(linres_control))
     CPASSERT(linres_control%ref_count>0)
     linres_control%ref_count=linres_control%ref_count+1

  END SUBROUTINE linres_control_retain

! *****************************************************************************
!> \brief ...
!> \param current_env ...
! *****************************************************************************
  SUBROUTINE current_env_create(current_env)

    TYPE(current_env_type)                   :: current_env

    CHARACTER(len=*), PARAMETER :: routineN = 'current_env_create', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure =.FALSE.

    CPASSERT(current_env%ref_count==0)
    current_env%ref_count = 1
    current_env%nao           = HUGE(1)
    current_env%gauge         = HUGE(1)
    current_env%orb_center    = HUGE(1)
    current_env%nstates(:)    = HUGE(1)
    current_env%nbr_center(:) = HUGE(1)
    current_env%use_old_gauge_atom = .TRUE.
    current_env%chi_pbc = .FALSE.
    current_env%do_selected_states = .FALSE.
    current_env%gauge_init = .FALSE.
    NULLIFY(current_env%full_done)
    NULLIFY(current_env%list_cubes)
    NULLIFY(current_env%statetrueindex)
    NULLIFY(current_env%basisfun_center)
    NULLIFY(current_env%center_list)
    NULLIFY(current_env%centers_set)
    NULLIFY(current_env%psi1_p)
    NULLIFY(current_env%psi1_rxp)
    NULLIFY(current_env%psi1_D)
    NULLIFY(current_env%p_psi0)
    NULLIFY(current_env%rxp_psi0)
    NULLIFY(current_env%jrho1_atom_set)
    NULLIFY(current_env%jrho1_set)
    NULLIFY(current_env%rs_gauge)
    NULLIFY(current_env%rs_buf)
    NULLIFY(current_env%selected_states_on_atom_list)
    NULLIFY(current_env%psi0_order)

  END SUBROUTINE current_env_create
! *****************************************************************************
!> \brief ...
!> \param nmr_env ...
! *****************************************************************************
  SUBROUTINE nmr_env_create(nmr_env)

    TYPE(nmr_env_type)                       :: nmr_env

    CHARACTER(len=*), PARAMETER :: routineN = 'nmr_env_create', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

     failure =.FALSE.

     CPASSERT(nmr_env%ref_count==0)
     nmr_env%ref_count = 1
     NULLIFY(nmr_env%chemical_shift)
     NULLIFY(nmr_env%chemical_shift_loc)
     NULLIFY(nmr_env%chemical_shift_nics_loc)
     NULLIFY(nmr_env%chemical_shift_nics)
     NULLIFY(nmr_env%r_nics)
     NULLIFY(nmr_env%cs_atom_list)
     NULLIFY(nmr_env%do_calc_cs_atom)

  END SUBROUTINE nmr_env_create

! *****************************************************************************
!> \brief ...
!> \param issc_env ...
! *****************************************************************************
  SUBROUTINE issc_env_create(issc_env)

    TYPE(issc_env_type)                      :: issc_env

    CHARACTER(len=*), PARAMETER :: routineN = 'issc_env_create', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

     failure =.FALSE.

     CPASSERT(issc_env%ref_count==0)
     issc_env%ref_count = 1
     NULLIFY(issc_env%issc)
     NULLIFY(issc_env%issc_loc)
     NULLIFY(issc_env%psi1_efg)
     NULLIFY(issc_env%psi1_fc)
     NULLIFY(issc_env%psi1_pso)
     NULLIFY(issc_env%psi1_dso)
     NULLIFY(issc_env%efg_psi0)
     NULLIFY(issc_env%pso_psi0)
     NULLIFY(issc_env%dso_psi0)
     NULLIFY(issc_env%fc_psi0)
     NULLIFY(issc_env%matrix_efg)
     NULLIFY(issc_env%matrix_pso)
     NULLIFY(issc_env%matrix_dso)
     NULLIFY(issc_env%matrix_fc)

   END SUBROUTINE issc_env_create

! *****************************************************************************
!> \brief ...
!> \param epr_env ...
! *****************************************************************************
  SUBROUTINE epr_env_create(epr_env)

    TYPE(epr_env_type)                       :: epr_env

    CHARACTER(len=*), PARAMETER :: routineN = 'epr_env_create', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure =.FALSE.

    CPASSERT(epr_env%ref_count==0)
    epr_env%ref_count = 1
    NULLIFY(epr_env%nablavks_set)
    NULLIFY(epr_env%nablavks_atom_set)
    NULLIFY(epr_env%bind_set)
    NULLIFY(epr_env%bind_atom_set)
    NULLIFY(epr_env%g_total)
    NULLIFY(epr_env%g_so)
    NULLIFY(epr_env%g_soo)
    NULLIFY(epr_env%vks_atom_set)

  END SUBROUTINE epr_env_create

! *****************************************************************************
!> \brief ...
!> \param current_env ...
!> \param simple_done ...
!> \param simple_converged ...
!> \param full_done ...
!> \param ref_count ...
!> \param nao ...
!> \param nstates ...
!> \param gauge ...
!> \param list_cubes ...
!> \param statetrueindex ...
!> \param gauge_name ...
!> \param basisfun_center ...
!> \param nbr_center ...
!> \param center_list ...
!> \param centers_set ...
!> \param psi1_p ...
!> \param psi1_rxp ...
!> \param psi1_D ...
!> \param p_psi0 ...
!> \param rxp_psi0 ...
!> \param jrho1_atom_set ...
!> \param jrho1_set ...
!> \param chi_tensor ...
!> \param chi_tensor_loc ...
!> \param gauge_atom_radius ...
!> \param rs_gauge ...
!> \param use_old_gauge_atom ...
!> \param chi_pbc ...
!> \param psi0_order ...
! *****************************************************************************
  SUBROUTINE get_current_env(current_env,simple_done,simple_converged,full_done,ref_count,nao,&
                             nstates,gauge,list_cubes,statetrueindex,gauge_name,basisfun_center,&
                             nbr_center,center_list,centers_set,psi1_p,psi1_rxp,psi1_D,p_psi0,&
                             rxp_psi0,jrho1_atom_set,jrho1_set,chi_tensor,&
                             chi_tensor_loc,gauge_atom_radius,rs_gauge,use_old_gauge_atom,&
                             chi_pbc,psi0_order)

    TYPE(current_env_type), OPTIONAL         :: current_env
    LOGICAL, OPTIONAL                        :: simple_done(6), &
                                                simple_converged(6)
    LOGICAL, DIMENSION(:, :), OPTIONAL, &
      POINTER                                :: full_done
    INTEGER, OPTIONAL                        :: ref_count, nao, nstates(2), &
                                                gauge
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: list_cubes
    INTEGER, DIMENSION(:, :, :), OPTIONAL, &
      POINTER                                :: statetrueindex
    CHARACTER(LEN=30), OPTIONAL              :: gauge_name
    REAL(dp), DIMENSION(:, :), OPTIONAL, &
      POINTER                                :: basisfun_center
    INTEGER, OPTIONAL                        :: nbr_center(2)
    TYPE(cp_2d_i_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: center_list
    TYPE(cp_2d_r_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: centers_set
    TYPE(cp_fm_p_type), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: psi1_p, psi1_rxp, psi1_D, &
                                                p_psi0, rxp_psi0
    TYPE(jrho_atom_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: jrho1_atom_set
    TYPE(qs_rho_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: jrho1_set
    REAL(dp), INTENT(OUT), OPTIONAL          :: chi_tensor(3,3,2), &
                                                chi_tensor_loc(3,3,2), &
                                                gauge_atom_radius
    TYPE(realspaces_grid_p_type), &
      DIMENSION(:), OPTIONAL, POINTER        :: rs_gauge
    LOGICAL, OPTIONAL                        :: use_old_gauge_atom, chi_pbc
    TYPE(cp_fm_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: psi0_order

    CHARACTER(len=*), PARAMETER :: routineN = 'get_current_env', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

!
!

    failure =.FALSE.

    CPASSERT(current_env%ref_count>0)

    IF(PRESENT(simple_done     )) simple_done(1:6)      = current_env%simple_done(1:6)
    IF(PRESENT(simple_converged)) simple_converged(1:6) = current_env%simple_converged(1:6)
    IF(PRESENT(full_done       )) full_done             => current_env%full_done
    IF(PRESENT(ref_count       )) ref_count             =  current_env%ref_count
    IF(PRESENT(nao             )) nao                   =  current_env%nao
    IF(PRESENT(nstates         )) nstates(1:2)          =  current_env%nstates(1:2)
    IF(PRESENT(gauge           )) gauge                 =  current_env%gauge
    IF(PRESENT(list_cubes      )) list_cubes            => current_env%list_cubes
    IF(PRESENT(statetrueindex  )) statetrueindex        => current_env%statetrueindex
    IF(PRESENT(gauge_name      )) gauge_name            =  current_env%gauge_name
    IF(PRESENT(basisfun_center )) basisfun_center       => current_env%basisfun_center
    IF(PRESENT(nbr_center      )) nbr_center(1:2)       =  current_env%nbr_center(1:2)
    IF(PRESENT(center_list     )) center_list           => current_env%center_list
    IF(PRESENT(centers_set     )) centers_set           => current_env%centers_set
    IF(PRESENT(chi_tensor      )) chi_tensor(:,:,:)     =  current_env%chi_tensor(:,:,:)
    IF(PRESENT(chi_tensor_loc  )) chi_tensor_loc(:,:,:) =  current_env%chi_tensor_loc(:,:,:)
    IF(PRESENT(psi1_p          )) psi1_p                => current_env%psi1_p
    IF(PRESENT(psi1_rxp        )) psi1_rxp              => current_env%psi1_rxp
    IF(PRESENT(psi1_D          )) psi1_D                => current_env%psi1_D
    IF(PRESENT(p_psi0          )) p_psi0                => current_env%p_psi0
    IF(PRESENT(rxp_psi0        )) rxp_psi0              => current_env%rxp_psi0
    IF(PRESENT(jrho1_atom_set  )) jrho1_atom_set        => current_env%jrho1_atom_set
    IF(PRESENT(jrho1_set       )) jrho1_set             => current_env%jrho1_set
    IF(PRESENT(rs_gauge        )) rs_gauge              => current_env%rs_gauge
    IF(PRESENT(psi0_order      )) psi0_order            => current_env%psi0_order
    IF(PRESENT(chi_pbc         )) chi_pbc               =  current_env%chi_pbc
    IF(PRESENT(gauge_atom_radius )) gauge_atom_radius   =  current_env%gauge_atom_radius
    IF(PRESENT(use_old_gauge_atom)) use_old_gauge_atom  =  current_env%use_old_gauge_atom

  END SUBROUTINE get_current_env

! *****************************************************************************
!> \brief ...
!> \param nmr_env ...
!> \param n_nics ...
!> \param cs_atom_list ...
!> \param do_calc_cs_atom ...
!> \param r_nics ...
!> \param chemical_shift ...
!> \param chemical_shift_loc ...
!> \param chemical_shift_nics_loc ...
!> \param chemical_shift_nics ...
!> \param shift_gapw_radius ...
!> \param do_nics ...
!> \param interpolate_shift ...
! *****************************************************************************
  SUBROUTINE get_nmr_env(nmr_env, n_nics, cs_atom_list, do_calc_cs_atom, &
                         r_nics, chemical_shift,chemical_shift_loc, &
                         chemical_shift_nics_loc, chemical_shift_nics, &
                         shift_gapw_radius,do_nics,interpolate_shift)

    TYPE(nmr_env_type)                       :: nmr_env
    INTEGER, INTENT(OUT), OPTIONAL           :: n_nics
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: cs_atom_list, do_calc_cs_atom
    REAL(dp), DIMENSION(:, :), OPTIONAL, &
      POINTER                                :: r_nics
    REAL(dp), DIMENSION(:, :, :), OPTIONAL, &
      POINTER                                :: chemical_shift, &
                                                chemical_shift_loc, &
                                                chemical_shift_nics_loc, &
                                                chemical_shift_nics
    REAL(dp), INTENT(OUT), OPTIONAL          :: shift_gapw_radius
    LOGICAL, INTENT(OUT), OPTIONAL           :: do_nics, interpolate_shift

    CHARACTER(len=*), PARAMETER :: routineN = 'get_nmr_env', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure =.FALSE.

    CPASSERT(nmr_env%ref_count>0)

    IF(PRESENT(n_nics             )) n_nics              =  nmr_env%n_nics
    IF(PRESENT(cs_atom_list       )) cs_atom_list        => nmr_env%cs_atom_list
    IF(PRESENT(do_calc_cs_atom    )) do_calc_cs_atom     => nmr_env%do_calc_cs_atom
    IF(PRESENT(chemical_shift     )) chemical_shift      => nmr_env%chemical_shift
    IF(PRESENT(chemical_shift_loc )) chemical_shift_loc  => nmr_env%chemical_shift_loc
    IF(PRESENT(chemical_shift_nics)) chemical_shift_nics => nmr_env%chemical_shift_nics
    IF(PRESENT(r_nics             )) r_nics              => nmr_env%r_nics
    IF(PRESENT(chemical_shift_nics_loc)) chemical_shift_nics_loc => nmr_env%chemical_shift_nics_loc
    IF(PRESENT(shift_gapw_radius  )) shift_gapw_radius   =  nmr_env%shift_gapw_radius
    IF(PRESENT(do_nics            )) do_nics             =  nmr_env%do_nics
    IF(PRESENT(interpolate_shift  )) interpolate_shift   =  nmr_env%interpolate_shift

  END SUBROUTINE get_nmr_env

! *****************************************************************************
!> \brief ...
!> \param issc_env ...
!> \param issc_on_atom_list ...
!> \param issc_gapw_radius ...
!> \param issc_loc ...
!> \param do_fc ...
!> \param do_sd ...
!> \param do_pso ...
!> \param do_dso ...
!> \param issc ...
!> \param interpolate_issc ...
!> \param psi1_efg ...
!> \param psi1_pso ...
!> \param psi1_dso ...
!> \param psi1_fc ...
!> \param efg_psi0 ...
!> \param pso_psi0 ...
!> \param dso_psi0 ...
!> \param fc_psi0 ...
!> \param matrix_efg ...
!> \param matrix_pso ...
!> \param matrix_dso ...
!> \param matrix_fc ...
! *****************************************************************************
  SUBROUTINE get_issc_env(issc_env,issc_on_atom_list,issc_gapw_radius,issc_loc,&
       do_fc,do_sd,do_pso,do_dso,&
       issc,interpolate_issc,psi1_efg,psi1_pso,psi1_dso,psi1_fc,efg_psi0,pso_psi0,dso_psi0,fc_psi0,&
       matrix_efg,matrix_pso,matrix_dso,matrix_fc)

    TYPE(issc_env_type)                      :: issc_env
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: issc_on_atom_list
    REAL(dp), OPTIONAL                       :: issc_gapw_radius
    REAL(dp), DIMENSION(:, :, :, :, :), &
      OPTIONAL, POINTER                      :: issc_loc
    LOGICAL, OPTIONAL                        :: do_fc, do_sd, do_pso, do_dso
    REAL(dp), DIMENSION(:, :, :, :, :), &
      OPTIONAL, POINTER                      :: issc
    LOGICAL, OPTIONAL                        :: interpolate_issc
    TYPE(cp_fm_p_type), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: psi1_efg, psi1_pso, psi1_dso
    TYPE(cp_fm_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: psi1_fc
    TYPE(cp_fm_p_type), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: efg_psi0, pso_psi0, dso_psi0
    TYPE(cp_fm_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: fc_psi0
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: matrix_efg, matrix_pso, &
                                                matrix_dso, matrix_fc

    CHARACTER(len=*), PARAMETER :: routineN = 'get_issc_env', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure =.FALSE.

    CPASSERT(issc_env%ref_count>0)

    IF(PRESENT(issc_on_atom_list)) issc_on_atom_list => issc_env%issc_on_atom_list
    IF(PRESENT(issc_gapw_radius )) issc_gapw_radius  =  issc_env%issc_gapw_radius
    IF(PRESENT(issc_loc         )) issc_loc          => issc_env%issc_loc
    IF(PRESENT(issc             )) issc              => issc_env%issc
    IF(PRESENT(interpolate_issc )) interpolate_issc  =  issc_env%interpolate_issc
    IF(PRESENT(psi1_efg         )) psi1_efg          => issc_env%psi1_efg
    IF(PRESENT(psi1_pso         )) psi1_pso          => issc_env%psi1_pso
    IF(PRESENT(psi1_dso         )) psi1_dso          => issc_env%psi1_dso
    IF(PRESENT(psi1_fc          )) psi1_fc           => issc_env%psi1_fc
    IF(PRESENT(efg_psi0         )) efg_psi0          => issc_env%efg_psi0
    IF(PRESENT(pso_psi0         )) pso_psi0          => issc_env%pso_psi0
    IF(PRESENT(dso_psi0         )) dso_psi0          => issc_env%dso_psi0
    IF(PRESENT(fc_psi0          )) fc_psi0           => issc_env%fc_psi0
    IF(PRESENT(matrix_efg       )) matrix_efg        => issc_env%matrix_efg
    IF(PRESENT(matrix_pso       )) matrix_pso        => issc_env%matrix_pso
    IF(PRESENT(matrix_fc        )) matrix_fc         => issc_env%matrix_fc
    IF(PRESENT(matrix_dso       )) matrix_dso        => issc_env%matrix_dso
    IF(PRESENT(do_fc            )) do_fc             =  issc_env%do_fc
    IF(PRESENT(do_sd            )) do_sd             =  issc_env%do_sd
    IF(PRESENT(do_pso           )) do_pso            =  issc_env%do_pso
    IF(PRESENT(do_dso           )) do_dso            =  issc_env%do_dso

  END SUBROUTINE get_issc_env

! *****************************************************************************
!> \brief ...
!> \param current_env ...
!> \param jrho1_atom_set ...
!> \param jrho1_set ...
! *****************************************************************************
  SUBROUTINE set_current_env(current_env,jrho1_atom_set,jrho1_set)

    TYPE(current_env_type)                   :: current_env
    TYPE(jrho_atom_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: jrho1_atom_set
    TYPE(qs_rho_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: jrho1_set

    CHARACTER(len=*), PARAMETER :: routineN = 'set_current_env', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: idir
    LOGICAL                                  :: failure

    failure =.FALSE.

    CPASSERT(current_env%ref_count>0)


    IF(PRESENT(jrho1_atom_set)) THEN
       IF(ASSOCIATED(current_env%jrho1_atom_set)) THEN
          CALL deallocate_jrho_atom_set(current_env%jrho1_atom_set)
       ENDIF
       current_env%jrho1_atom_set => jrho1_atom_set
    END IF

    IF(PRESENT(jrho1_set)) THEN
       IF(ASSOCIATED(current_env%jrho1_set)) THEN
          DO idir = 1,3
             CALL qs_rho_release(current_env%jrho1_set(idir)%rho)
          END DO
       END IF
       current_env%jrho1_set => jrho1_set
    END IF

  END SUBROUTINE set_current_env
! *****************************************************************************
!> \brief ...
!> \param nmr_env ...
!> \param shift_factor ...
!> \param chi_factor ...
!> \param chi_SI2shiftppm ...
!> \param chi_SI2ppmcgs ...
! *****************************************************************************
  SUBROUTINE set_nmr_env(nmr_env,shift_factor,chi_factor,chi_SI2shiftppm,chi_SI2ppmcgs)

    TYPE(nmr_env_type)                       :: nmr_env
    REAL(dp), INTENT(IN), OPTIONAL           :: shift_factor, chi_factor, &
                                                chi_SI2shiftppm, chi_SI2ppmcgs

    CHARACTER(len=*), PARAMETER :: routineN = 'set_nmr_env', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure =.FALSE.

    CPASSERT(nmr_env%ref_count>0)

    IF(PRESENT(shift_factor   )) nmr_env%chi_factor      = chi_factor
    IF(PRESENT(shift_factor   )) nmr_env%chi_factor      = chi_factor
    IF(PRESENT(chi_SI2shiftppm)) nmr_env%chi_SI2shiftppm = chi_SI2shiftppm
    IF(PRESENT(chi_SI2ppmcgs  )) nmr_env%chi_SI2ppmcgs   = chi_SI2ppmcgs

  END SUBROUTINE set_nmr_env
! *****************************************************************************
!> \brief ...
!> \param issc_env ...
! *****************************************************************************
  SUBROUTINE set_issc_env(issc_env)

    TYPE(issc_env_type)                      :: issc_env

    CHARACTER(len=*), PARAMETER :: routineN = 'set_issc_env', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure =.FALSE.

    CPASSERT(issc_env%ref_count>0)

  END SUBROUTINE set_issc_env

! *****************************************************************************
!> \brief ...
!> \param epr_env ...
!> \param g_total ...
!> \param g_so ...
!> \param g_soo ...
!> \param nablavks_set ...
!> \param nablavks_atom_set ...
!> \param bind_set ...
!> \param bind_atom_set ...
! *****************************************************************************
  SUBROUTINE get_epr_env(epr_env, g_total, g_so, g_soo, nablavks_set, nablavks_atom_set, &
       bind_set, bind_atom_set)

    TYPE(epr_env_type)                       :: epr_env
    REAL(dp), DIMENSION(:, :), OPTIONAL, &
      POINTER                                :: g_total, g_so, g_soo
    TYPE(qs_rho_p_type), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: nablavks_set
    TYPE(nablavks_atom_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: nablavks_atom_set
    TYPE(qs_rho_p_type), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: bind_set
    TYPE(rho_atom_coeff), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: bind_atom_set

    CHARACTER(len=*), PARAMETER :: routineN = 'get_epr_env', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure =.FALSE.

    CPASSERT(epr_env%ref_count>0)

    IF(PRESENT(g_total)) g_total => epr_env%g_total
    IF(PRESENT(g_so)) g_so => epr_env%g_so
    IF(PRESENT(g_soo)) g_soo => epr_env%g_soo
    IF(PRESENT(nablavks_set)) nablavks_set => epr_env%nablavks_set
    IF(PRESENT(nablavks_atom_set)) nablavks_atom_set => epr_env%nablavks_atom_set
    IF(PRESENT(bind_set)) bind_set => epr_env%bind_set
    IF(PRESENT(bind_atom_set)) bind_atom_set => epr_env%bind_atom_set

  END SUBROUTINE get_epr_env

! *****************************************************************************
!> \brief ...
!> \param epr_env ...
!> \param g_free_factor ...
!> \param g_soo_chicorr_factor ...
!> \param g_soo_factor ...
!> \param g_so_factor ...
!> \param g_so_factor_gapw ...
!> \param g_zke_factor ...
!> \param nablavks_set ...
!> \param nablavks_atom_set ...
! *****************************************************************************
  SUBROUTINE set_epr_env(epr_env,g_free_factor,g_soo_chicorr_factor,&
                         g_soo_factor,g_so_factor,g_so_factor_gapw,&
                         g_zke_factor,nablavks_set,nablavks_atom_set)

    TYPE(epr_env_type)                       :: epr_env
    REAL(dp), INTENT(IN), OPTIONAL :: g_free_factor, g_soo_chicorr_factor, &
      g_soo_factor, g_so_factor, g_so_factor_gapw, g_zke_factor
    TYPE(qs_rho_p_type), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: nablavks_set
    TYPE(nablavks_atom_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: nablavks_atom_set

    CHARACTER(len=*), PARAMETER :: routineN = 'set_epr_env', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: idir, ispin
    LOGICAL                                  :: failure

    failure =.FALSE.

    CPASSERT(epr_env%ref_count>0)

    IF(PRESENT(g_free_factor)) epr_env%g_free_factor=g_free_factor
    IF(PRESENT(g_zke_factor)) epr_env%g_zke_factor=g_zke_factor
    IF(PRESENT(g_so_factor)) epr_env%g_so_factor=g_so_factor
    IF(PRESENT(g_so_factor_gapw)) epr_env%g_so_factor_gapw=g_so_factor_gapw
    IF(PRESENT(g_soo_factor)) epr_env%g_soo_factor=g_soo_factor
    IF(PRESENT(g_soo_chicorr_factor)) epr_env%g_soo_chicorr_factor=g_soo_chicorr_factor

    IF(PRESENT(nablavks_set)) THEN
        IF(ASSOCIATED(epr_env%nablavks_set)) THEN
           DO ispin = 1,2
             DO idir = 1,3
               CALL qs_rho_release(epr_env%nablavks_set(idir,ispin)%rho)
             END DO
           END DO
        END IF
        epr_env%nablavks_set => nablavks_set
    ENDIF

    IF(PRESENT(nablavks_atom_set)) THEN
       IF(ASSOCIATED(epr_env%nablavks_atom_set)) THEN
          CALL deallocate_nablavks_atom_set(epr_env%nablavks_atom_set)
       ENDIF
       epr_env%nablavks_atom_set => nablavks_atom_set
    ENDIF

  END SUBROUTINE set_epr_env

! *****************************************************************************
!> \brief ...
!> \param nablavks_atom_set ...
!> \param natom ...
! *****************************************************************************
  SUBROUTINE allocate_nablavks_atom_set(nablavks_atom_set,natom)

    TYPE(nablavks_atom_type), DIMENSION(:), &
      POINTER                                :: nablavks_atom_set
    INTEGER, INTENT(IN)                      :: natom

    CHARACTER(len=*), PARAMETER :: routineN = 'allocate_nablavks_atom_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iat
    LOGICAL                                  :: failure

    failure = .FALSE.

    ALLOCATE(nablavks_atom_set(natom))

    DO iat = 1,natom
       NULLIFY(nablavks_atom_set(iat)%nablavks_vec_rad_h)
       NULLIFY(nablavks_atom_set(iat)%nablavks_vec_rad_s)
    ENDDO
  END SUBROUTINE allocate_nablavks_atom_set

! *****************************************************************************
!> \brief ...
!> \param nablavks_atom_set ...
! *****************************************************************************
  SUBROUTINE deallocate_nablavks_atom_set(nablavks_atom_set)

    TYPE(nablavks_atom_type), DIMENSION(:), &
      POINTER                                :: nablavks_atom_set

    CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_nablavks_atom_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iat, idir, n, natom
    LOGICAL                                  :: failure

    failure = .FALSE.

    CPASSERT(ASSOCIATED(nablavks_atom_set))
    natom = SIZE(nablavks_atom_set)

    DO iat = 1,natom
       IF(ASSOCIATED(nablavks_atom_set(iat)%nablavks_vec_rad_h)) THEN
          IF(ASSOCIATED(nablavks_atom_set(iat)%nablavks_vec_rad_h(1,1)%r_coef)) THEN
             n = SIZE(nablavks_atom_set(iat)%nablavks_vec_rad_h,2)
             DO i=1,n
                DO idir = 1,3
                   DEALLOCATE(nablavks_atom_set(iat)%nablavks_vec_rad_h(idir,i)%r_coef)
                   DEALLOCATE(nablavks_atom_set(iat)%nablavks_vec_rad_s(idir,i)%r_coef)
                ENDDO
             ENDDO
          ENDIF
          DEALLOCATE(nablavks_atom_set(iat)%nablavks_vec_rad_h)
          DEALLOCATE(nablavks_atom_set(iat)%nablavks_vec_rad_s)
       ENDIF
    ENDDO
    DEALLOCATE(nablavks_atom_set)
  END SUBROUTINE deallocate_nablavks_atom_set

! *****************************************************************************
!> \brief ...
!> \param jrho_atom_set ...
!> \param natom ...
! *****************************************************************************
  SUBROUTINE allocate_jrho_atom_set(jrho_atom_set,natom)

    TYPE(jrho_atom_type), DIMENSION(:), &
      POINTER                                :: jrho_atom_set
    INTEGER, INTENT(IN)                      :: natom

    CHARACTER(len=*), PARAMETER :: routineN = 'allocate_jrho_atom_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iat
    LOGICAL                                  :: failure

    failure = .FALSE.

    ALLOCATE(jrho_atom_set(natom))

    DO iat = 1,natom
       NULLIFY(jrho_atom_set(iat)%cjc0_h)
       NULLIFY(jrho_atom_set(iat)%cjc0_s)
       NULLIFY(jrho_atom_set(iat)%cjc_h)
       NULLIFY(jrho_atom_set(iat)%cjc_s)
       NULLIFY(jrho_atom_set(iat)%cjc_ii_h)
       NULLIFY(jrho_atom_set(iat)%cjc_ii_s)
       NULLIFY(jrho_atom_set(iat)%cjc_iii_h)
       NULLIFY(jrho_atom_set(iat)%cjc_iii_s)
       NULLIFY(jrho_atom_set(iat)%jrho_vec_rad_h)
       NULLIFY(jrho_atom_set(iat)%jrho_vec_rad_s)
       NULLIFY(jrho_atom_set(iat)%jrho_h)
       NULLIFY(jrho_atom_set(iat)%jrho_s)
       NULLIFY(jrho_atom_set(iat)%jrho_a_h)
       NULLIFY(jrho_atom_set(iat)%jrho_a_s)
       NULLIFY(jrho_atom_set(iat)%jrho_b_h)
       NULLIFY(jrho_atom_set(iat)%jrho_b_s)
       NULLIFY(jrho_atom_set(iat)%jrho_a_h_ii)
       NULLIFY(jrho_atom_set(iat)%jrho_a_s_ii)
       NULLIFY(jrho_atom_set(iat)%jrho_b_h_ii)
       NULLIFY(jrho_atom_set(iat)%jrho_b_s_ii)
       NULLIFY(jrho_atom_set(iat)%jrho_a_h_iii)
       NULLIFY(jrho_atom_set(iat)%jrho_a_s_iii)
       NULLIFY(jrho_atom_set(iat)%jrho_b_h_iii)
       NULLIFY(jrho_atom_set(iat)%jrho_b_s_iii)
    ENDDO
  END SUBROUTINE allocate_jrho_atom_set

! *****************************************************************************
!> \brief ...
!> \param jrho_atom_set ...
! *****************************************************************************
  SUBROUTINE deallocate_jrho_atom_set(jrho_atom_set)

    TYPE(jrho_atom_type), DIMENSION(:), &
      POINTER                                :: jrho_atom_set

    CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_jrho_atom_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iat, idir, n, natom

    CPASSERT(ASSOCIATED(jrho_atom_set))
    natom = SIZE(jrho_atom_set)

    DO iat = 1,natom
       IF(ASSOCIATED(jrho_atom_set(iat)%cjc_h)) THEN
          IF(ASSOCIATED(jrho_atom_set(iat)%cjc_h(1)%r_coef)) THEN
             n = SIZE(jrho_atom_set(iat)%cjc_h)
             DO i = 1,n
                !
                ! size = (nsotot,nsotot) replicated
                DEALLOCATE(jrho_atom_set(iat)%cjc0_h(i)%r_coef,&
                           jrho_atom_set(iat)%cjc0_s(i)%r_coef,&
                           jrho_atom_set(iat)%cjc_h(i)%r_coef,&
                           jrho_atom_set(iat)%cjc_s(i)%r_coef,&
                           jrho_atom_set(iat)%cjc_ii_h(i)%r_coef,&
                           jrho_atom_set(iat)%cjc_ii_s(i)%r_coef,&
                           jrho_atom_set(iat)%cjc_iii_h(i)%r_coef,&
                           jrho_atom_set(iat)%cjc_iii_s(i)%r_coef)
             END DO
          END IF
          DEALLOCATE(jrho_atom_set(iat)%cjc0_h,&
                     jrho_atom_set(iat)%cjc0_s,&
                     jrho_atom_set(iat)%cjc_h,&
                     jrho_atom_set(iat)%cjc_s,&
                     jrho_atom_set(iat)%cjc_ii_h,&
                     jrho_atom_set(iat)%cjc_ii_s,&
                     jrho_atom_set(iat)%cjc_iii_h,&
                     jrho_atom_set(iat)%cjc_iii_s)
       END IF

       IF(ASSOCIATED(jrho_atom_set(iat)%jrho_a_h)) THEN
          IF(ASSOCIATED(jrho_atom_set(iat)%jrho_a_h(1)%r_coef)) THEN
             n = SIZE(jrho_atom_set(iat)%jrho_a_h)
             DO i = 1,n
                !
                ! size = (nr,max_iso_not0) distributed
                DEALLOCATE(jrho_atom_set(iat)%jrho_h(i)%r_coef,&
                           jrho_atom_set(iat)%jrho_s(i)%r_coef,&
                           jrho_atom_set(iat)%jrho_a_h(i)%r_coef,&
                           jrho_atom_set(iat)%jrho_a_s(i)%r_coef,&
                           jrho_atom_set(iat)%jrho_b_h(i)%r_coef,&
                           jrho_atom_set(iat)%jrho_b_s(i)%r_coef,&
                           jrho_atom_set(iat)%jrho_a_h_ii(i)%r_coef,&
                           jrho_atom_set(iat)%jrho_a_s_ii(i)%r_coef,&
                           jrho_atom_set(iat)%jrho_b_h_ii(i)%r_coef,&
                           jrho_atom_set(iat)%jrho_b_s_ii(i)%r_coef,&
                           jrho_atom_set(iat)%jrho_a_h_iii(i)%r_coef,&
                           jrho_atom_set(iat)%jrho_a_s_iii(i)%r_coef,&
                           jrho_atom_set(iat)%jrho_b_h_iii(i)%r_coef,&
                           jrho_atom_set(iat)%jrho_b_s_iii(i)%r_coef)
             END DO
          END IF
          DEALLOCATE(jrho_atom_set(iat)%jrho_h,&
                     jrho_atom_set(iat)%jrho_s,&
                     jrho_atom_set(iat)%jrho_a_h,&
                     jrho_atom_set(iat)%jrho_a_s,&
                     jrho_atom_set(iat)%jrho_b_h,&
                     jrho_atom_set(iat)%jrho_b_s,&
                     jrho_atom_set(iat)%jrho_a_h_ii,&
                     jrho_atom_set(iat)%jrho_a_s_ii,&
                     jrho_atom_set(iat)%jrho_b_h_ii,&
                     jrho_atom_set(iat)%jrho_b_s_ii,&
                     jrho_atom_set(iat)%jrho_a_h_iii,&
                     jrho_atom_set(iat)%jrho_a_s_iii,&
                     jrho_atom_set(iat)%jrho_b_h_iii,&
                     jrho_atom_set(iat)%jrho_b_s_iii)
       END IF

       IF(ASSOCIATED(jrho_atom_set(iat)%jrho_vec_rad_h)) THEN
          IF(ASSOCIATED(jrho_atom_set(iat)%jrho_vec_rad_h(1,1)%r_coef)) THEN
             n = SIZE(jrho_atom_set(iat)%jrho_vec_rad_h,2)
             DO i = 1,n
                DO idir = 1,3
                   !
                   ! size =(nr,na) distributed
                   DEALLOCATE(jrho_atom_set(iat)%jrho_vec_rad_h(idir,i)%r_coef,&
                              jrho_atom_set(iat)%jrho_vec_rad_s(idir,i)%r_coef)
                END DO
             END DO
          ENDIF
          DEALLOCATE(jrho_atom_set(iat)%jrho_vec_rad_h,&
                     jrho_atom_set(iat)%jrho_vec_rad_s)
       END IF
    END DO
    DEALLOCATE(jrho_atom_set)

  END SUBROUTINE deallocate_jrho_atom_set

! *****************************************************************************
!> \brief ...
!> \param jrho1_atom ...
!> \param ispin ...
!> \param nr ...
!> \param na ...
!> \param max_iso_not0 ...
! *****************************************************************************
  SUBROUTINE  allocate_jrho_atom_rad(jrho1_atom,ispin,nr,na,max_iso_not0)

    TYPE(jrho_atom_type), POINTER            :: jrho1_atom
    INTEGER, INTENT(IN)                      :: ispin, nr, na, max_iso_not0

    CHARACTER(len=*), PARAMETER :: routineN = 'allocate_jrho_atom_rad', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, idir
    LOGICAL                                  :: failure

    CALL timeset(routineN,handle)

    failure = .FALSE.

    CPASSERT(ASSOCIATED(jrho1_atom))


    DO idir = 1,3
       ALLOCATE(jrho1_atom%jrho_vec_rad_h(idir,ispin)%r_coef(nr,na),&
                jrho1_atom%jrho_vec_rad_s(idir,ispin)%r_coef(nr,na))
       jrho1_atom%jrho_vec_rad_h(idir,ispin)%r_coef = 0.0_dp
       jrho1_atom%jrho_vec_rad_s(idir,ispin)%r_coef = 0.0_dp
    ENDDO

    ALLOCATE(jrho1_atom%jrho_h(ispin)%r_coef(nr,max_iso_not0),&
             jrho1_atom%jrho_s(ispin)%r_coef(nr,max_iso_not0),&
             jrho1_atom%jrho_a_h(ispin)%r_coef(nr,max_iso_not0),&
             jrho1_atom%jrho_a_s(ispin)%r_coef(nr,max_iso_not0),&
             jrho1_atom%jrho_b_h(ispin)%r_coef(nr,max_iso_not0),&
             jrho1_atom%jrho_b_s(ispin)%r_coef(nr,max_iso_not0),&
             jrho1_atom%jrho_a_h_ii(ispin)%r_coef(nr,max_iso_not0),&
             jrho1_atom%jrho_a_s_ii(ispin)%r_coef(nr,max_iso_not0),&
             jrho1_atom%jrho_b_h_ii(ispin)%r_coef(nr,max_iso_not0),&
             jrho1_atom%jrho_b_s_ii(ispin)%r_coef(nr,max_iso_not0),&
             jrho1_atom%jrho_a_h_iii(ispin)%r_coef(nr,max_iso_not0),&
             jrho1_atom%jrho_a_s_iii(ispin)%r_coef(nr,max_iso_not0),&
             jrho1_atom%jrho_b_h_iii(ispin)%r_coef(nr,max_iso_not0),&
             jrho1_atom%jrho_b_s_iii(ispin)%r_coef(nr,max_iso_not0))
    !
    jrho1_atom%jrho_h(ispin)%r_coef = 0.0_dp
    jrho1_atom%jrho_s(ispin)%r_coef = 0.0_dp
    jrho1_atom%jrho_a_h(ispin)%r_coef = 0.0_dp
    jrho1_atom%jrho_a_s(ispin)%r_coef = 0.0_dp
    jrho1_atom%jrho_b_h(ispin)%r_coef = 0.0_dp
    jrho1_atom%jrho_b_s(ispin)%r_coef = 0.0_dp
    jrho1_atom%jrho_a_h_ii(ispin)%r_coef = 0.0_dp
    jrho1_atom%jrho_a_s_ii(ispin)%r_coef = 0.0_dp
    jrho1_atom%jrho_b_h_ii(ispin)%r_coef = 0.0_dp
    jrho1_atom%jrho_b_s_ii(ispin)%r_coef = 0.0_dp
    jrho1_atom%jrho_a_h_iii(ispin)%r_coef = 0.0_dp
    jrho1_atom%jrho_a_s_iii(ispin)%r_coef = 0.0_dp
    jrho1_atom%jrho_b_h_iii(ispin)%r_coef = 0.0_dp
    jrho1_atom%jrho_b_s_iii(ispin)%r_coef = 0.0_dp
    CALL timestop(handle)

  END SUBROUTINE allocate_jrho_atom_rad

! *****************************************************************************
!> \brief ...
!> \param jrho1_atom ...
!> \param ispin ...
! *****************************************************************************
  SUBROUTINE set2zero_jrho_atom_rad(jrho1_atom,ispin)
    !
    TYPE(jrho_atom_type), POINTER            :: jrho1_atom
    INTEGER, INTENT(IN)                      :: ispin

    CHARACTER(len=*), PARAMETER :: routineN = 'set2zero_jrho_atom_rad', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    !
    CPASSERT(ASSOCIATED(jrho1_atom))
    !
    jrho1_atom%jrho_h(ispin)%r_coef = 0.0_dp
    jrho1_atom%jrho_s(ispin)%r_coef = 0.0_dp
    !
    jrho1_atom%jrho_a_h(ispin)%r_coef = 0.0_dp
    jrho1_atom%jrho_a_s(ispin)%r_coef = 0.0_dp
    jrho1_atom%jrho_b_h(ispin)%r_coef = 0.0_dp
    jrho1_atom%jrho_b_s(ispin)%r_coef = 0.0_dp
    !
    jrho1_atom%jrho_a_h_ii(ispin)%r_coef = 0.0_dp
    jrho1_atom%jrho_a_s_ii(ispin)%r_coef = 0.0_dp
    jrho1_atom%jrho_b_h_ii(ispin)%r_coef = 0.0_dp
    jrho1_atom%jrho_b_s_ii(ispin)%r_coef = 0.0_dp
    !
    jrho1_atom%jrho_a_h_iii(ispin)%r_coef = 0.0_dp
    jrho1_atom%jrho_a_s_iii(ispin)%r_coef = 0.0_dp
    jrho1_atom%jrho_b_h_iii(ispin)%r_coef = 0.0_dp
    jrho1_atom%jrho_b_s_iii(ispin)%r_coef = 0.0_dp
    !
  END SUBROUTINE set2zero_jrho_atom_rad

! *****************************************************************************

! *****************************************************************************
!> \brief ...
!> \param jrho1_atom_set ...
!> \param iatom ...
!> \param nsotot ...
! *****************************************************************************
  SUBROUTINE allocate_jrho_coeff(jrho1_atom_set,iatom,nsotot)

    TYPE(jrho_atom_type), DIMENSION(:), &
      POINTER                                :: jrho1_atom_set
    INTEGER, INTENT(IN)                      :: iatom, nsotot

    CHARACTER(len=*), PARAMETER :: routineN = 'allocate_jrho_coeff', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i
    LOGICAL                                  :: failure

    CALL timeset(routineN,handle)
    failure = .FALSE.
    CPASSERT(ASSOCIATED(jrho1_atom_set))
    DO i = 1,SIZE(jrho1_atom_set(iatom)%cjc0_h,1)
       ALLOCATE(jrho1_atom_set(iatom)%cjc0_h(i)%r_coef(nsotot,nsotot),&
                jrho1_atom_set(iatom)%cjc0_s(i)%r_coef(nsotot,nsotot),&
                jrho1_atom_set(iatom)%cjc_h(i)%r_coef(nsotot,nsotot),&
                jrho1_atom_set(iatom)%cjc_s(i)%r_coef(nsotot,nsotot),&
                jrho1_atom_set(iatom)%cjc_ii_h(i)%r_coef(nsotot,nsotot),&
                jrho1_atom_set(iatom)%cjc_ii_s(i)%r_coef(nsotot,nsotot),&
                jrho1_atom_set(iatom)%cjc_iii_h(i)%r_coef(nsotot,nsotot),&
                jrho1_atom_set(iatom)%cjc_iii_s(i)%r_coef(nsotot,nsotot))
       jrho1_atom_set(iatom)%cjc0_h(i)%r_coef = 0.0_dp
       jrho1_atom_set(iatom)%cjc0_s(i)%r_coef = 0.0_dp
       jrho1_atom_set(iatom)%cjc_h(i)%r_coef = 0.0_dp
       jrho1_atom_set(iatom)%cjc_s(i)%r_coef = 0.0_dp
       jrho1_atom_set(iatom)%cjc_ii_h(i)%r_coef = 0.0_dp
       jrho1_atom_set(iatom)%cjc_ii_s(i)%r_coef = 0.0_dp
       jrho1_atom_set(iatom)%cjc_iii_h(i)%r_coef = 0.0_dp
       jrho1_atom_set(iatom)%cjc_iii_s(i)%r_coef = 0.0_dp
    ENDDO
    CALL timestop(handle)
  END SUBROUTINE allocate_jrho_coeff

! *****************************************************************************

! *****************************************************************************
!> \brief ...
!> \param jrho1_atom_set ...
!> \param iatom ...
! *****************************************************************************
  SUBROUTINE deallocate_jrho_coeff(jrho1_atom_set,iatom)

    TYPE(jrho_atom_type), DIMENSION(:), &
      POINTER                                :: jrho1_atom_set
    INTEGER, INTENT(IN)                      :: iatom

    CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_jrho_coeff', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i
    LOGICAL                                  :: failure

    CALL timeset(routineN,handle)
    failure = .FALSE.
    CPASSERT(ASSOCIATED(jrho1_atom_set))
    DO i = 1,SIZE(jrho1_atom_set(iatom)%cjc0_h,1)
       DEALLOCATE(jrho1_atom_set(iatom)%cjc0_h(i)%r_coef,&
                  jrho1_atom_set(iatom)%cjc0_s(i)%r_coef,&
                  jrho1_atom_set(iatom)%cjc_h(i)%r_coef,&
                  jrho1_atom_set(iatom)%cjc_s(i)%r_coef,&
                  jrho1_atom_set(iatom)%cjc_ii_h(i)%r_coef,&
                  jrho1_atom_set(iatom)%cjc_ii_s(i)%r_coef,&
                  jrho1_atom_set(iatom)%cjc_iii_h(i)%r_coef,&
                  jrho1_atom_set(iatom)%cjc_iii_s(i)%r_coef)
    ENDDO
    CALL timestop(handle)
  END SUBROUTINE deallocate_jrho_coeff

! *****************************************************************************

! *****************************************************************************
!> \brief ...
!> \param jrho1_atom_set ...
!> \param iatom ...
!> \param cjc_h ...
!> \param cjc_s ...
!> \param cjc_ii_h ...
!> \param cjc_ii_s ...
!> \param cjc_iii_h ...
!> \param cjc_iii_s ...
!> \param jrho_vec_rad_h ...
!> \param jrho_vec_rad_s ...
! *****************************************************************************
  SUBROUTINE get_jrho_atom(jrho1_atom_set,iatom,cjc_h,cjc_s,cjc_ii_h,cjc_ii_s,&
             cjc_iii_h,cjc_iii_s,jrho_vec_rad_h,jrho_vec_rad_s)

    TYPE(jrho_atom_type), DIMENSION(:), &
      POINTER                                :: jrho1_atom_set
    INTEGER, INTENT(IN)                      :: iatom
    TYPE(rho_atom_coeff), DIMENSION(:), &
      OPTIONAL, POINTER                      :: cjc_h, cjc_s, cjc_ii_h, &
                                                cjc_ii_s, cjc_iii_h, cjc_iii_s
    TYPE(rho_atom_coeff), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: jrho_vec_rad_h, jrho_vec_rad_s

    CHARACTER(len=*), PARAMETER :: routineN = 'get_jrho_atom', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.

    CPASSERT(ASSOCIATED(jrho1_atom_set))

    IF(PRESENT(cjc_h         )) cjc_h          => jrho1_atom_set(iatom)%cjc_h
    IF(PRESENT(cjc_s         )) cjc_s          => jrho1_atom_set(iatom)%cjc_s
    IF(PRESENT(cjc_ii_h      )) cjc_ii_h       => jrho1_atom_set(iatom)%cjc_ii_h
    IF(PRESENT(cjc_ii_s      )) cjc_ii_s       => jrho1_atom_set(iatom)%cjc_ii_s
    IF(PRESENT(cjc_iii_h     )) cjc_iii_h      => jrho1_atom_set(iatom)%cjc_iii_h
    IF(PRESENT(cjc_iii_s     )) cjc_iii_s      => jrho1_atom_set(iatom)%cjc_iii_s
    IF(PRESENT(jrho_vec_rad_h)) jrho_vec_rad_h => jrho1_atom_set(iatom)%jrho_vec_rad_h
    IF(PRESENT(jrho_vec_rad_s)) jrho_vec_rad_s => jrho1_atom_set(iatom)%jrho_vec_rad_s

  END SUBROUTINE get_jrho_atom

! *****************************************************************************
!> \brief ...
!> \param jrho1_atom_set ...
!> \param atomic_kind_set ...
!> \param nspins ...
! *****************************************************************************
  SUBROUTINE init_jrho_atom_set(jrho1_atom_set,atomic_kind_set,nspins)
    TYPE(jrho_atom_type), DIMENSION(:), &
      POINTER                                :: jrho1_atom_set
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    INTEGER, INTENT(IN)                      :: nspins

    CHARACTER(len=*), PARAMETER :: routineN = 'init_jrho_atom_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iat, iatom, idir, &
                                                ikind, ispin, nat, natom, &
                                                nkind
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: failure

    CALL timeset(routineN,handle)

    failure = .FALSE.

    CPASSERT(ASSOCIATED(atomic_kind_set))

    IF(ASSOCIATED(jrho1_atom_set)) THEN
       CALL deallocate_jrho_atom_set(jrho1_atom_set)
    END IF

    CALL get_atomic_kind_set(atomic_kind_set, natom = natom)

    CALL allocate_jrho_atom_set(jrho1_atom_set,natom)

    nkind = SIZE(atomic_kind_set)

    DO ikind = 1,nkind

       CALL get_atomic_kind(atomic_kind_set(ikind), atom_list=atom_list, natom=nat)

       DO iat = 1,nat
          iatom = atom_list(iat)

          !*** allocate the radial density for each LM,for each atom ***
          ALLOCATE(jrho1_atom_set(iatom)%jrho_vec_rad_h(3,nspins),&
                   jrho1_atom_set(iatom)%jrho_vec_rad_s(3,nspins),&
                   jrho1_atom_set(iatom)%jrho_h(nspins),&
                   jrho1_atom_set(iatom)%jrho_s(nspins),&
                   jrho1_atom_set(iatom)%jrho_a_h(nspins),&
                   jrho1_atom_set(iatom)%jrho_a_s(nspins),&
                   jrho1_atom_set(iatom)%jrho_b_h(nspins),&
                   jrho1_atom_set(iatom)%jrho_b_s(nspins),&
                   jrho1_atom_set(iatom)%jrho_a_h_ii(nspins),&
                   jrho1_atom_set(iatom)%jrho_a_s_ii(nspins),&
                   jrho1_atom_set(iatom)%jrho_b_s_ii(nspins),&
                   jrho1_atom_set(iatom)%jrho_b_h_ii(nspins),&
                   jrho1_atom_set(iatom)%jrho_a_h_iii(nspins),&
                   jrho1_atom_set(iatom)%jrho_a_s_iii(nspins),&
                   jrho1_atom_set(iatom)%jrho_b_s_iii(nspins),&
                   jrho1_atom_set(iatom)%jrho_b_h_iii(nspins),&
                   jrho1_atom_set(iatom)%cjc0_h(nspins),&
                   jrho1_atom_set(iatom)%cjc0_s(nspins),&
                   jrho1_atom_set(iatom)%cjc_h(nspins),&
                   jrho1_atom_set(iatom)%cjc_s(nspins),&
                   jrho1_atom_set(iatom)%cjc_ii_h(nspins),&
                   jrho1_atom_set(iatom)%cjc_ii_s(nspins),&
                   jrho1_atom_set(iatom)%cjc_iii_h(nspins),&
                   jrho1_atom_set(iatom)%cjc_iii_s(nspins))

          DO ispin = 1,nspins
             DO idir = 1,3
                NULLIFY(jrho1_atom_set(iatom)%jrho_vec_rad_h(idir,ispin)%r_coef)
                NULLIFY(jrho1_atom_set(iatom)%jrho_vec_rad_s(idir,ispin)%r_coef)
             END DO
             NULLIFY(jrho1_atom_set(iatom)%jrho_h(ispin)%r_coef)
             NULLIFY(jrho1_atom_set(iatom)%jrho_s(ispin)%r_coef)
             NULLIFY(jrho1_atom_set(iatom)%jrho_a_h(ispin)%r_coef)
             NULLIFY(jrho1_atom_set(iatom)%jrho_a_s(ispin)%r_coef)
             NULLIFY(jrho1_atom_set(iatom)%jrho_b_h(ispin)%r_coef)
             NULLIFY(jrho1_atom_set(iatom)%jrho_b_s(ispin)%r_coef)
             NULLIFY(jrho1_atom_set(iatom)%jrho_a_h_ii(ispin)%r_coef)
             NULLIFY(jrho1_atom_set(iatom)%jrho_a_s_ii(ispin)%r_coef)
             NULLIFY(jrho1_atom_set(iatom)%jrho_b_h_ii(ispin)%r_coef)
             NULLIFY(jrho1_atom_set(iatom)%jrho_b_s_ii(ispin)%r_coef)
             NULLIFY(jrho1_atom_set(iatom)%jrho_a_h_iii(ispin)%r_coef)
             NULLIFY(jrho1_atom_set(iatom)%jrho_a_s_iii(ispin)%r_coef)
             NULLIFY(jrho1_atom_set(iatom)%jrho_b_h_iii(ispin)%r_coef)
             NULLIFY(jrho1_atom_set(iatom)%jrho_b_s_iii(ispin)%r_coef)
             NULLIFY(jrho1_atom_set(iatom)%cjc0_h(ispin)%r_coef)
             NULLIFY(jrho1_atom_set(iatom)%cjc0_s(ispin)%r_coef)
             NULLIFY(jrho1_atom_set(iatom)%cjc_h(ispin)%r_coef)
             NULLIFY(jrho1_atom_set(iatom)%cjc_s(ispin)%r_coef)
             NULLIFY(jrho1_atom_set(iatom)%cjc_ii_h(ispin)%r_coef)
             NULLIFY(jrho1_atom_set(iatom)%cjc_ii_s(ispin)%r_coef)
             NULLIFY(jrho1_atom_set(iatom)%cjc_iii_h(ispin)%r_coef)
             NULLIFY(jrho1_atom_set(iatom)%cjc_iii_s(ispin)%r_coef)

          ENDDO ! ispin

       END DO  ! iat

    END DO  ! ikind

    CALL timestop(handle)

  END SUBROUTINE init_jrho_atom_set

! *****************************************************************************
!> \brief ...
!> \param nablavks_atom_set ...
!> \param atomic_kind_set ...
!> \param qs_kind_set ...
!> \param nspins ...
! *****************************************************************************
  SUBROUTINE init_nablavks_atom_set(nablavks_atom_set,atomic_kind_set,qs_kind_set,nspins)

    TYPE(nablavks_atom_type), DIMENSION(:), &
      POINTER                                :: nablavks_atom_set
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    INTEGER, INTENT(IN)                      :: nspins

    CHARACTER(len=*), PARAMETER :: routineN = 'init_nablavks_atom_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iat, iatom, idir, &
                                                ikind, ispin, max_iso_not0, &
                                                maxso, na, nat, natom, nkind, &
                                                nr, nset, nsotot
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: failure
    TYPE(grid_atom_type), POINTER            :: grid_atom
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(harmonics_atom_type), POINTER       :: harmonics

    CALL timeset(routineN,handle)

    failure = .FALSE.

    CPASSERT(ASSOCIATED(qs_kind_set))

    IF(ASSOCIATED(nablavks_atom_set)) THEN
      CALL deallocate_nablavks_atom_set(nablavks_atom_set)
    END IF

    CALL get_atomic_kind_set(atomic_kind_set, natom=natom)

    CALL allocate_nablavks_atom_set(nablavks_atom_set,natom)

    nkind = SIZE(atomic_kind_set)

    DO ikind = 1,nkind
       CALL get_atomic_kind(atomic_kind_set(ikind), atom_list=atom_list,natom=nat)
       CALL get_qs_kind(qs_kind_set(ikind),&
                        basis_set=orb_basis_set, &
                        harmonics=harmonics,&
                        grid_atom=grid_atom)

       na = grid_atom%ng_sphere
       nr = grid_atom%nr

       CALL get_gto_basis_set(gto_basis_set=orb_basis_set, &
                              maxso=maxso, nset=nset)
       nsotot = maxso * nset
       max_iso_not0 = harmonics%max_iso_not0
       DO iat = 1,nat
          iatom = atom_list(iat)
          !*** allocate the radial density for each LM,for each atom ***

            ALLOCATE (nablavks_atom_set(iatom)%nablavks_vec_rad_h(3,nspins))
            ALLOCATE (nablavks_atom_set(iatom)%nablavks_vec_rad_s(3,nspins))
            DO ispin = 1,nspins
              DO idir = 1,3
                NULLIFY(nablavks_atom_set(iatom)%nablavks_vec_rad_h(idir,ispin)%r_coef)
                NULLIFY(nablavks_atom_set(iatom)%nablavks_vec_rad_s(idir,ispin)%r_coef)
                ALLOCATE(nablavks_atom_set(iatom)%nablavks_vec_rad_h(idir,ispin)%r_coef(nr,na))
                ALLOCATE(nablavks_atom_set(iatom)%nablavks_vec_rad_s(idir,ispin)%r_coef(nr,na))
              END DO
            END DO  ! ispin
       END DO  ! iat

    END DO  ! ikind

    CALL timestop(handle)

  END SUBROUTINE init_nablavks_atom_set
! *****************************************************************************
!> \brief ...
!> \param polar_env ...
! *****************************************************************************
  SUBROUTINE polar_env_create(polar_env)

    TYPE(polar_env_type)                     :: polar_env

    CHARACTER(len=*), PARAMETER :: routineN = 'polar_env_create', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

     failure =.FALSE.

     CPASSERT(polar_env%ref_count==0)
     polar_env%ref_count = 1
!        polar_env%do_raman=.FALSE.
     NULLIFY(polar_env%polar)
     NULLIFY(polar_env%psi1_dBerry)
     NULLIFY(polar_env%dBerry_psi0)
     NULLIFY(polar_env%mo_derivs)

   END SUBROUTINE polar_env_create

! *****************************************************************************
!> \brief ...
!> \param polar_env ...
!> \param do_raman ...
!> \param dBerry_psi0 ...
!> \param polar ...
!> \param psi1_dBerry ...
!> \param mo_derivs ...
! *****************************************************************************
  SUBROUTINE get_polar_env(polar_env,do_raman,dBerry_psi0, polar, &
       psi1_dBerry,mo_derivs)

    TYPE(polar_env_type)                     :: polar_env
    LOGICAL, OPTIONAL                        :: do_raman
    TYPE(cp_fm_p_type), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: dBerry_psi0
    REAL(dp), DIMENSION(:, :), OPTIONAL, &
      POINTER                                :: polar
    TYPE(cp_fm_p_type), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: psi1_dBerry, mo_derivs

    CHARACTER(len=*), PARAMETER :: routineN = 'get_polar_env', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure =.FALSE.

    CPASSERT(polar_env%ref_count>0)

    IF(PRESENT(polar))             polar => polar_env%polar
    IF(PRESENT(psi1_dBerry )) psi1_dBerry => polar_env%psi1_dBerry
    IF(PRESENT(dBerry_psi0)) dBerry_psi0          => polar_env%dBerry_psi0
    IF(PRESENT(mo_derivs            )) mo_derivs             => polar_env%mo_derivs
    IF(PRESENT(do_raman           )) do_raman            =  polar_env%do_raman
  END SUBROUTINE get_polar_env

! *****************************************************************************

! *****************************************************************************
!> \brief ...
!> \param polar_env ...
! *****************************************************************************
  SUBROUTINE set_polar_env(polar_env)

    TYPE(polar_env_type)                     :: polar_env

    CHARACTER(len=*), PARAMETER :: routineN = 'set_polar_env', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure =.FALSE.

    CPASSERT(polar_env%ref_count>0)

  END SUBROUTINE set_polar_env
! *****************************************************************************



END MODULE qs_linres_types

