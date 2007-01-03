!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!****** cp2k/qs_linres_types [1.0] *
!!
!!   NAME
!!     qs_linres_types
!!
!!   FUNCTION
!!     Type definitiona for linear response calculations
!!
!!   AUTHOR
!!     MI 
!!
!!   MODIFICATION HISTORY
!!     
!!
!!   SOURCE
!******************************************************************************

MODULE qs_linres_types

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cp_array_r_utils,                ONLY: cp_2d_r_p_type
  USE cp_fm_pool_types,                ONLY: cp_fm_pool_p_type
  USE cp_fm_types,                     ONLY: cp_fm_p_type
  USE kinds,                           ONLY: dp
  USE qs_harmonics_atom,               ONLY: harmonics_atom_type
  USE qs_loc_control,                  ONLY: localized_wfn_control_create,&
                                             localized_wfn_control_release,&
                                             localized_wfn_control_type
  USE qs_rho_atom_types,               ONLY: rho_atom_coeff
  USE qs_rho_types,                    ONLY: qs_rho_p_type,&
                                             qs_rho_release
  USE sparse_matrix_types,             ONLY: real_matrix_p_type,&
                                             real_matrix_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

!****s* qs_linres_types/linres_control_type
!!
!! NAME
!!   linres_control_type
!!
!! FUNCTION
!!   General settings for linear response calculations
!!
!! ATTRIBUTES
!!   - property : which quantity is to be calculated by LR
!!   - opt_method : method to optimize the psi1 by minimization of the second order term of the energy
!!   - preconditioner : which kind of preconditioner should be used, if any
!!   - localized_psi0 : don't use the canonical psi0, but the maximally localized wavefunctions
!!   - do_kernel  : the kernel is zero if the rho1 is zero as for the magnetic field perturbation 
!!   - tolerance : convergence criterium for the optimization of the psi1
!!
!! NOTES
!!   
!! AUTHOR
!!   MI 
!!
!!***
  TYPE linres_control_type
     INTEGER                                   :: ref_count
     INTEGER                                   :: property
     INTEGER                                   :: preconditioner_type
     REAL(dp)                                  :: energy_gap
     INTEGER                                   :: opt_method
     INTEGER                                   :: ls_method
     REAL(dp)                                  :: ds_min
     INTEGER                                   :: max_scf
     LOGICAL                                   :: localized_psi0
     LOGICAL                                   :: do_kernel
     LOGICAL                                   :: converged
     LOGICAL                                   :: check_ortho
     REAL(KIND=dp)                             :: eps_scf
     TYPE(localized_wfn_control_type),POINTER  :: localized_wfn_control
     CHARACTER(LEN=8)                          :: flag
  END TYPE linres_control_type

!****s* qs_linres_types/nmr_env_type
!!
!! NAME
!!   nmr_env_type
!!
!! FUNCTION
!!
!! ATTRIBUTES
!!   - ref_count
!!   - full_nmr    : true if the full correction is calculated
!!   - simplenmr_done, fullnmr_done : flags that indicate what has been
!!                   already calculated: used for restart
!!   - centers_set : centers of the maximally localized psi0
!!   - spreads_set : spreads of the maximally localized psi0 
!!   - op_p_ao     : sparse matrixes for the px, py and pz operator 
!!                   in the contracted basis set representation 
!!   - op_rmd_ao   : sparse matrixes used to construct the (r-d)xp operator
!!                   in the contracted basis set representation
!!   - p_psi0      : full matrixes, operator p applied to psi0
!!   - rxp_psi0    : full matrixes, operator (r-d)xp applied to psi0
!!   - psi1_p      : response wavefunctions to the perturbation given by 
!!                   H1=p (xyz)  applied to psi0
!!   - psi1_rxp    : response wavefunctions to the perturbation given by 
!!                   H1=(r-d_i)xp applied to psi0_i where d_i is the center
!!   - psi1_D      : response wavefunctions to the perturbation given by 
!!                   H1=(d_j-d_i)xp applied to psi0_i where d_i is the center
!!                   and d_j is the center of psi0_j and psi1_D_j is the result
!!                   This operator has to be used in nstate scf calculations,
!!                   one for each psi1_D_j vector
!!   - chemical_shift: the tensor for each atom
!!   - chi_tensor  : the susceptibility tensor
!!   - jp1_ao :      current density matrices. 
!!   - jp2_ao :      density matrices used to compute the r-dependent terms
!!   - jrho1_set   : current density on the global grid, if gapw this is only the soft part
!!   - jrho1_atom_set : current density on the local atomic grids (only if gapw)
!!
!! NOTES
!!   
!! AUTHOR
!!   MI 
!!***

  TYPE nmr_env_type
     INTEGER                                     :: ref_count, n_nics, nao, nstates(2)
     INTEGER, DIMENSION(:), POINTER              :: list_cubes
     INTEGER, DIMENSION(:), POINTER              :: cs_atom_list
     INTEGER, DIMENSION(:), POINTER              :: do_calc_cs_atom
     INTEGER, DIMENSION(:,:), POINTER            :: statetrueindex
     LOGICAL                                     :: do_nics,&
                                                    full_nmr,simplenmr_done(6),&
                                                    simplenmr_converged(6),&
                                                    store_current,restart_nmr
     LOGICAL, DIMENSION(:,:), POINTER            :: fullnmr_done
     REAL(dp)                                    :: shift_factor, chi_factor
     REAL(dp)                                    :: chi_SI2shiftppm, chi_SI2ppmcgs
     REAL(dp)                                    :: chi_tensor(3,3)
     REAL(dp)                                    :: chi_tensor_loc(3,3)
     REAL(dp), DIMENSION(:,:), POINTER           :: basisfun_center
     REAL(dp), DIMENSION(:,:), POINTER           :: r_nics
     REAL(dp), DIMENSION(:,:,:), POINTER         :: chemical_shift
     REAL(dp), DIMENSION(:,:,:), POINTER         :: chemical_shift_loc
     REAL(dp), DIMENSION(:,:,:), POINTER         :: chemical_shift_loc_nics
     REAL(dp), DIMENSION(:,:,:), POINTER         :: chemical_shift_nics
     TYPE(cp_2d_r_p_type), DIMENSION(:),POINTER  :: centers_set
     TYPE(cp_fm_pool_p_type), DIMENSION(:), &
       POINTER                                   :: ao_mo_fm_pools
     TYPE(cp_fm_p_type), DIMENSION(:,:), POINTER :: psi1_p, psi1_rxp, psi1_D
     TYPE(cp_fm_p_type), DIMENSION(:,:), POINTER :: p_psi0, rxp_psi0
     TYPE(jrho_atom_type),DIMENSION(:),&
       POINTER                                   :: jrho1_atom_set
     TYPE(qs_rho_p_type),DIMENSION(:,:), POINTER :: jrho1_set
     TYPE(real_matrix_p_type), DIMENSION(:), &
          POINTER                                :: op_p_ao, op_rmd_ao
     TYPE(real_matrix_p_type), DIMENSION(:), &
          POINTER                                :: jp2_ao
     TYPE(real_matrix_type), POINTER             :: jp1_ao
  END TYPE nmr_env_type

  TYPE jrho_atom_p_type
     TYPE(jrho_atom_type), POINTER       :: jrho_atom
  END TYPE jrho_atom_p_type 

  TYPE jrho_atom_type
    TYPE(rho_atom_coeff), DIMENSION(:),&
                               POINTER      :: cjc_h, cjc_s
    TYPE(rho_atom_coeff), DIMENSION(:),&
                               POINTER      :: cjc_ii_h, cjc_ii_s
    TYPE(rho_atom_coeff), DIMENSION(:),&
                               POINTER      :: cjc_iii_h, cjc_iii_s
    TYPE(rho_atom_coeff), DIMENSION(:),&
                               POINTER      :: jrho_rad_h,&
                                               jrho_rad_s
    TYPE(rho_atom_coeff), DIMENSION(:,:,:),&
                               POINTER      :: jrho_vec_rad_h,&
                                               jrho_vec_rad_s
    TYPE(rho_atom_coeff), DIMENSION(:), &
      POINTER                               :: jrho_a_h, jrho_a_s
    TYPE(rho_atom_coeff), DIMENSION(:), &
      POINTER                               :: jrho_b_h, jrho_b_s
    TYPE(rho_atom_coeff), DIMENSION(:), &
      POINTER                               :: jrho_a_h_ii, jrho_a_s_ii
    TYPE(rho_atom_coeff), DIMENSION(:), &
      POINTER                               :: jrho_b_h_ii, jrho_b_s_ii
    TYPE(rho_atom_coeff), DIMENSION(:), &
      POINTER                               :: jrho_c_h_ii, jrho_c_s_ii
    TYPE(rho_atom_coeff), DIMENSION(:), &
      POINTER                               :: jrho_a_h_iii, jrho_a_s_iii
    TYPE(rho_atom_coeff), DIMENSION(:), &
      POINTER                               :: jrho_b_h_iii, jrho_b_s_iii
    TYPE(rho_atom_coeff), DIMENSION(:), &
      POINTER                               :: jrho_c_h_iii, jrho_c_s_iii
  END TYPE jrho_atom_type



  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_linres_types'

! *** Public data types ***

  PUBLIC ::  allocate_jrho_atom_rad, deallocate_jrho_atom_set, get_nmr_env, &
             init_jrho_atom_set, linres_control_type, &
             nmr_env_type, jrho_atom_type, set_nmr_env, set2zero_jrho_atom_rad


! *** Public subroutines ***

  PUBLIC :: linres_control_create, linres_control_retain, linres_control_release,&
            nmr_env_create

!!***
! *****************************************************************************

CONTAINS  

  SUBROUTINE linres_control_create(linres_control,error)
 
    TYPE(linres_control_type), POINTER       :: linres_control
    TYPE(cp_error_type), INTENT(inout)  &
                                             :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'linres_control_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat
    LOGICAL                                  :: failure

     failure =.FALSE.

     CPPrecondition(.NOT.ASSOCIATED(linres_control),cp_failure_level,routineP,error,failure)
     IF (.NOT. failure) THEN
       ALLOCATE (linres_control,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       linres_control%ref_count=1
       NULLIFY(linres_control%localized_wfn_control)
       CALL localized_wfn_control_create(linres_control%localized_wfn_control, error=error)

     END IF

  END SUBROUTINE linres_control_create

! *****************************************************************************
  SUBROUTINE linres_control_release(linres_control,error)
 
    TYPE(linres_control_type), POINTER       :: linres_control
    TYPE(cp_error_type), INTENT(inout)  &
                                             :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'linres_control_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat
    LOGICAL                                  :: failure

     failure =.FALSE.
     
     IF (ASSOCIATED(linres_control)) THEN
       CPPostcondition(linres_control%ref_count>0,cp_failure_level,routineP,error,failure)
       linres_control%ref_count=linres_control%ref_count-1
       IF(linres_control%ref_count<1)THEN
         CALL localized_wfn_control_release(linres_control%localized_wfn_control, error=error)
         DEALLOCATE(linres_control,STAT=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       END IF
     END IF
     NULLIFY(linres_control)
  END SUBROUTINE linres_control_release

! *****************************************************************************
  SUBROUTINE linres_control_retain(linres_control,error)
 
    TYPE(linres_control_type), POINTER       :: linres_control
    TYPE(cp_error_type), INTENT(inout)  &
                                             :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'linres_control_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

     failure =.FALSE.

     CPPrecondition(ASSOCIATED(linres_control),cp_failure_level,routineP,error,failure)
     IF (.NOT. failure) THEN
       CPPostcondition(linres_control%ref_count>0,cp_failure_level,routineP,error,failure)
       linres_control%ref_count=linres_control%ref_count+1
     END IF
     
  END SUBROUTINE linres_control_retain

! *****************************************************************************

  SUBROUTINE nmr_env_create(nmr_env,error)

    TYPE(nmr_env_type)                       :: nmr_env
    TYPE(cp_error_type), INTENT(inout)  &
                                             :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'nmr_env_create', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

     failure =.FALSE.

     CPPrecondition(nmr_env%ref_count==0, cp_failure_level,routineP,error,failure)
     IF(.NOT. failure) THEN
       nmr_env%ref_count = 1
       NULLIFY(nmr_env%centers_set)
       NULLIFY(nmr_env%ao_mo_fm_pools)
       NULLIFY(nmr_env%op_p_ao)
       NULLIFY(nmr_env%op_rmd_ao)
       NULLIFY(nmr_env%psi1_p)
       NULLIFY(nmr_env%psi1_rxp)
       NULLIFY(nmr_env%psi1_D)
       NULLIFY(nmr_env%jp1_ao)
       NULLIFY(nmr_env%jp2_ao)
       NULLIFY(nmr_env%p_psi0)
       NULLIFY(nmr_env%rxp_psi0)
       NULLIFY(nmr_env%jrho1_set)
       NULLIFY(nmr_env%jrho1_atom_set)
       NULLIFY(nmr_env%fullnmr_done)
       NULLIFY(nmr_env%basisfun_center)
       NULLIFY(nmr_env%chemical_shift)
       NULLIFY(nmr_env%chemical_shift_loc)
       NULLIFY(nmr_env%chemical_shift_loc_nics)
       NULLIFY(nmr_env%chemical_shift_nics)
       NULLIFY(nmr_env%r_nics)
       NULLIFY(nmr_env%statetrueindex)
       NULLIFY(nmr_env%list_cubes)
       NULLIFY(nmr_env%cs_atom_list)
       NULLIFY(nmr_env%do_calc_cs_atom)
     END IF 


  END SUBROUTINE nmr_env_create
       

! *****************************************************************************

  SUBROUTINE get_nmr_env(nmr_env, n_nics, nao, nstates,&
             chi_tensor,list_cubes, cs_atom_list, do_calc_cs_atom, statetrueindex, &
             chi_tensor_loc, r_nics, chemical_shift,&
             chemical_shift_loc, chemical_shift_loc_nics, chemical_shift_nics,&
             basisfun_center, centers_set,&
             psi1_p, psi1_rxp, psi1_D, p_psi0, rxp_psi0,&
             jrho1_atom_set, jrho1_set, op_p_ao, op_rmd_ao, &
             jp1_ao, jp2_ao, error)

    TYPE(nmr_env_type)                       :: nmr_env
    INTEGER, INTENT(OUT), OPTIONAL           :: n_nics, nao, nstates(2)
    REAL(dp), INTENT(out), OPTIONAL          :: chi_tensor(3,3)
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: list_cubes, cs_atom_list, &
                                                do_calc_cs_atom
    INTEGER, DIMENSION(:, :), OPTIONAL, &
      POINTER                                :: statetrueindex
    REAL(dp), INTENT(out), OPTIONAL          :: chi_tensor_loc(3,3)
    REAL(dp), DIMENSION(:, :), OPTIONAL, &
      POINTER                                :: r_nics
    REAL(dp), DIMENSION(:, :, :), OPTIONAL, &
      POINTER                                :: chemical_shift, &
                                                chemical_shift_loc, &
                                                chemical_shift_loc_nics, &
                                                chemical_shift_nics
    REAL(dp), DIMENSION(:, :), OPTIONAL, &
      POINTER                                :: basisfun_center
    TYPE(cp_2d_r_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: centers_set
    TYPE(cp_fm_p_type), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: psi1_p, psi1_rxp, psi1_D, &
                                                p_psi0, rxp_psi0
    TYPE(jrho_atom_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: jrho1_atom_set
    TYPE(qs_rho_p_type), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: jrho1_set
    TYPE(real_matrix_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: op_p_ao, op_rmd_ao
    TYPE(real_matrix_type), OPTIONAL, &
      POINTER                                :: jp1_ao
    TYPE(real_matrix_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: jp2_ao
    TYPE(cp_error_type), INTENT(inout)  &
                                             :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_nmr_env', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

     failure =.FALSE.

     CPPrecondition(nmr_env%ref_count>0, cp_failure_level,routineP,error,failure)

     IF(PRESENT(n_nics)) n_nics=nmr_env%n_nics
     IF(PRESENT(nao)) nao=nmr_env%nao
     IF(PRESENT(nstates)) nstates(1:2)=nmr_env%nstates(1:2)
     IF(PRESENT(list_cubes)) list_cubes => nmr_env%list_cubes
     IF(PRESENT(cs_atom_list)) cs_atom_list => nmr_env%cs_atom_list
     IF(PRESENT(do_calc_cs_atom)) do_calc_cs_atom => nmr_env%do_calc_cs_atom
     IF(PRESENT(statetrueindex)) statetrueindex => nmr_env%statetrueindex
     IF(PRESENT(centers_set)) centers_set => nmr_env%centers_set
     IF(PRESENT(chi_tensor)) chi_tensor = nmr_env%chi_tensor
     IF(PRESENT(chi_tensor_loc)) chi_tensor_loc = nmr_env%chi_tensor_loc
     IF(PRESENT(chemical_shift)) chemical_shift => nmr_env%chemical_shift
     IF(PRESENT(r_nics)) r_nics => nmr_env%r_nics
     IF(PRESENT(basisfun_center)) basisfun_center => nmr_env%basisfun_center
     IF(PRESENT(centers_set)) centers_set => nmr_env%centers_set
     IF(PRESENT(psi1_p)) psi1_p => nmr_env%psi1_p
     IF(PRESENT(psi1_rxp)) psi1_rxp => nmr_env%psi1_rxp
     IF(PRESENT(psi1_D)) psi1_D => nmr_env%psi1_D
     IF(PRESENT(p_psi0)) p_psi0 => nmr_env%p_psi0
     IF(PRESENT(rxp_psi0)) rxp_psi0 => nmr_env%rxp_psi0
     IF(PRESENT(jrho1_atom_set)) jrho1_atom_set => nmr_env%jrho1_atom_set
     IF(PRESENT(jrho1_set)) jrho1_set => nmr_env%jrho1_set
     IF(PRESENT(op_rmd_ao)) op_rmd_ao => nmr_env%op_rmd_ao
     IF(PRESENT(op_p_ao)) op_p_ao => nmr_env%op_p_ao
     IF(PRESENT(jp2_ao)) jp2_ao => nmr_env%jp2_ao
     IF(PRESENT(jp1_ao)) jp1_ao => nmr_env%jp1_ao
     
  END SUBROUTINE get_nmr_env

  SUBROUTINE set_nmr_env(nmr_env,shift_factor,chi_factor,chi_SI2shiftppm,chi_SI2ppmcgs,&
             chi_tensor,chi_tensor_loc,jrho1_atom_set,jrho1_set,error)

    TYPE(nmr_env_type)                       :: nmr_env
    REAL(dp), INTENT(IN), OPTIONAL :: shift_factor, chi_factor, &
      chi_SI2shiftppm, chi_SI2ppmcgs, chi_tensor(3,3), chi_tensor_loc(3,3)
    TYPE(jrho_atom_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: jrho1_atom_set
    TYPE(qs_rho_p_type), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: jrho1_set
    TYPE(cp_error_type), INTENT(inout)  &
                                             :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'set_nmr_env', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iB, idir
    LOGICAL                                  :: failure

     failure =.FALSE.

     CPPrecondition(nmr_env%ref_count>0, cp_failure_level,routineP,error,failure)

     IF(PRESENT(shift_factor)) nmr_env%chi_factor=chi_factor
     IF(PRESENT(shift_factor)) nmr_env%chi_factor=chi_factor
     IF(PRESENT(chi_SI2shiftppm)) nmr_env%chi_SI2shiftppm=chi_SI2shiftppm
     IF(PRESENT(chi_SI2ppmcgs))   nmr_env%chi_SI2ppmcgs=chi_SI2ppmcgs
     IF(PRESENT(chi_tensor))   THEN
       DO iB=1,3
       DO idir=1,3
         nmr_env%chi_tensor(idir,iB) = chi_tensor(idir,iB)
       END DO 
       END DO 
     END IF
     IF(PRESENT(chi_tensor_loc))   THEN
       DO iB=1,3
       DO idir=1,3
         nmr_env%chi_tensor_loc(idir,iB) = chi_tensor_loc(idir,iB)
       END DO 
       END DO 
     END IF

     IF(PRESENT(jrho1_atom_set)) THEN
       IF(ASSOCIATED(nmr_env%jrho1_atom_set)) THEN
         CALL deallocate_jrho_atom_set(nmr_env%jrho1_atom_set,error=error)
       END IF
       nmr_env%jrho1_atom_set => jrho1_atom_set
     END IF

     IF(PRESENT(jrho1_set)) THEN
       IF(ASSOCIATED(nmr_env%jrho1_set)) THEN
          DO iB = 1,3
            DO idir = 1,3
              CALL qs_rho_release(nmr_env%jrho1_set(idir,iB)%rho,error=error)
            END DO
          END DO
        END IF
       nmr_env%jrho1_set => jrho1_set
     END IF

  END SUBROUTINE set_nmr_env

! *****************************************************************************
  SUBROUTINE allocate_jrho_atom_set(jrho_atom_set,natom,error)

    TYPE(jrho_atom_type), DIMENSION(:), &
      POINTER                                :: jrho_atom_set
    INTEGER, INTENT(IN)                      :: natom
    TYPE(cp_error_type), INTENT(inout)  &
                                             :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_jrho_atom_set'

    INTEGER                                  :: iat, istat
    LOGICAL                                  :: failure

     failure = .FALSE.

     ALLOCATE(jrho_atom_set(natom), STAT=istat)
     CPPrecondition(istat==0, cp_failure_level,routineN,error,failure)
     
     DO iat = 1,natom

       NULLIFY(jrho_atom_set(iat)%cjc_h)
       NULLIFY(jrho_atom_set(iat)%cjc_s)
       NULLIFY(jrho_atom_set(iat)%cjc_ii_h)
       NULLIFY(jrho_atom_set(iat)%cjc_ii_s)
       NULLIFY(jrho_atom_set(iat)%cjc_iii_h)
       NULLIFY(jrho_atom_set(iat)%cjc_iii_s)
       NULLIFY(jrho_atom_set(iat)%jrho_rad_h)
       NULLIFY(jrho_atom_set(iat)%jrho_rad_s)
       NULLIFY(jrho_atom_set(iat)%jrho_vec_rad_h)
       NULLIFY(jrho_atom_set(iat)%jrho_vec_rad_s)
       NULLIFY(jrho_atom_set(iat)%jrho_a_h)
       NULLIFY(jrho_atom_set(iat)%jrho_a_s)
       NULLIFY(jrho_atom_set(iat)%jrho_b_h)
       NULLIFY(jrho_atom_set(iat)%jrho_b_s)
       NULLIFY(jrho_atom_set(iat)%jrho_a_h_ii)
       NULLIFY(jrho_atom_set(iat)%jrho_a_s_ii)
       NULLIFY(jrho_atom_set(iat)%jrho_b_h_ii)
       NULLIFY(jrho_atom_set(iat)%jrho_b_s_ii)
       NULLIFY(jrho_atom_set(iat)%jrho_c_h_ii)
       NULLIFY(jrho_atom_set(iat)%jrho_c_s_ii)
       NULLIFY(jrho_atom_set(iat)%jrho_a_h_iii)
       NULLIFY(jrho_atom_set(iat)%jrho_a_s_iii)
       NULLIFY(jrho_atom_set(iat)%jrho_b_h_iii)
       NULLIFY(jrho_atom_set(iat)%jrho_b_s_iii)
       NULLIFY(jrho_atom_set(iat)%jrho_c_h_iii)
       NULLIFY(jrho_atom_set(iat)%jrho_c_s_iii)

     END DO

  END SUBROUTINE allocate_jrho_atom_set

! *****************************************************************************
  SUBROUTINE deallocate_jrho_atom_set(jrho_atom_set,error)

    TYPE(jrho_atom_type), DIMENSION(:), &
      POINTER                                :: jrho_atom_set
    TYPE(cp_error_type), INTENT(inout)  &
                                             :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'deallocate_jrho_atom_set'

    INTEGER                                  :: i, iat, iB, idir, istat, n, &
                                                natom
    LOGICAL                                  :: failure

     CPPrecondition(ASSOCIATED(jrho_atom_set),cp_failure_level,routineN,error,failure)
     natom = SIZE(jrho_atom_set)

     DO iat = 1,natom

       IF(ASSOCIATED(jrho_atom_set(iat)%cjc_h)) THEN
         IF(ASSOCIATED(jrho_atom_set(iat)%cjc_h(1)%r_coef)) THEN
           n = SIZE(jrho_atom_set(iat)%cjc_h)
           DO i = 1,n
             DEALLOCATE(jrho_atom_set(iat)%cjc_h(i)%r_coef,STAT=istat)
             CPPrecondition(istat==0, cp_failure_level,routineN,error,failure)
             DEALLOCATE(jrho_atom_set(iat)%cjc_s(i)%r_coef,STAT=istat)
             CPPrecondition(istat==0, cp_failure_level,routineN,error,failure)
             DEALLOCATE(jrho_atom_set(iat)%cjc_ii_h(i)%r_coef,STAT=istat)
             CPPrecondition(istat==0, cp_failure_level,routineN,error,failure)
             DEALLOCATE(jrho_atom_set(iat)%cjc_ii_s(i)%r_coef,STAT=istat)
             CPPrecondition(istat==0, cp_failure_level,routineN,error,failure)
             DEALLOCATE(jrho_atom_set(iat)%cjc_iii_h(i)%r_coef,STAT=istat)
             CPPrecondition(istat==0, cp_failure_level,routineN,error,failure)
             DEALLOCATE(jrho_atom_set(iat)%cjc_iii_s(i)%r_coef,STAT=istat)
             CPPrecondition(istat==0, cp_failure_level,routineN,error,failure)
           END DO
         END IF
         DEALLOCATE(jrho_atom_set(iat)%cjc_h,STAT=istat)
       CPPrecondition(istat==0, cp_failure_level,routineN,error,failure)
         DEALLOCATE(jrho_atom_set(iat)%cjc_s,STAT=istat)
       CPPrecondition(istat==0, cp_failure_level,routineN,error,failure)
         DEALLOCATE(jrho_atom_set(iat)%cjc_ii_h,STAT=istat)
       CPPrecondition(istat==0, cp_failure_level,routineN,error,failure)
         DEALLOCATE(jrho_atom_set(iat)%cjc_ii_s,STAT=istat)
       CPPrecondition(istat==0, cp_failure_level,routineN,error,failure)
         DEALLOCATE(jrho_atom_set(iat)%cjc_iii_h,STAT=istat)
       CPPrecondition(istat==0, cp_failure_level,routineN,error,failure)
         DEALLOCATE(jrho_atom_set(iat)%cjc_iii_s,STAT=istat)
       CPPrecondition(istat==0, cp_failure_level,routineN,error,failure)
       END IF

       IF(ASSOCIATED(jrho_atom_set(iat)%jrho_a_h)) THEN
         IF(ASSOCIATED(jrho_atom_set(iat)%jrho_a_h(1)%r_coef)) THEN
           n = SIZE(jrho_atom_set(iat)%jrho_a_h)
           DO i = 1,n
              DEALLOCATE(jrho_atom_set(iat)%jrho_a_h(i)%r_coef)
              DEALLOCATE(jrho_atom_set(iat)%jrho_a_s(i)%r_coef)
              DEALLOCATE(jrho_atom_set(iat)%jrho_b_h(i)%r_coef)
              DEALLOCATE(jrho_atom_set(iat)%jrho_b_s(i)%r_coef)

              DEALLOCATE(jrho_atom_set(iat)%jrho_a_h_ii(i)%r_coef)
              DEALLOCATE(jrho_atom_set(iat)%jrho_a_s_ii(i)%r_coef)
              DEALLOCATE(jrho_atom_set(iat)%jrho_b_h_ii(i)%r_coef)
              DEALLOCATE(jrho_atom_set(iat)%jrho_b_s_ii(i)%r_coef)
              DEALLOCATE(jrho_atom_set(iat)%jrho_c_h_ii(i)%r_coef)
              DEALLOCATE(jrho_atom_set(iat)%jrho_c_s_ii(i)%r_coef)

              DEALLOCATE(jrho_atom_set(iat)%jrho_a_h_iii(i)%r_coef)
              DEALLOCATE(jrho_atom_set(iat)%jrho_a_s_iii(i)%r_coef)
              DEALLOCATE(jrho_atom_set(iat)%jrho_b_h_iii(i)%r_coef)
              DEALLOCATE(jrho_atom_set(iat)%jrho_b_s_iii(i)%r_coef)
              DEALLOCATE(jrho_atom_set(iat)%jrho_c_h_iii(i)%r_coef)
              DEALLOCATE(jrho_atom_set(iat)%jrho_c_s_iii(i)%r_coef)
           END DO
         END IF
         DEALLOCATE(jrho_atom_set(iat)%jrho_a_h)
         DEALLOCATE(jrho_atom_set(iat)%jrho_a_s)
         DEALLOCATE(jrho_atom_set(iat)%jrho_b_h)
         DEALLOCATE(jrho_atom_set(iat)%jrho_b_s)
         DEALLOCATE(jrho_atom_set(iat)%jrho_a_h_ii)
         DEALLOCATE(jrho_atom_set(iat)%jrho_a_s_ii)
         DEALLOCATE(jrho_atom_set(iat)%jrho_b_h_ii)
         DEALLOCATE(jrho_atom_set(iat)%jrho_b_s_ii)
         DEALLOCATE(jrho_atom_set(iat)%jrho_c_h_ii)
         DEALLOCATE(jrho_atom_set(iat)%jrho_c_s_ii)
         DEALLOCATE(jrho_atom_set(iat)%jrho_a_h_iii)
         DEALLOCATE(jrho_atom_set(iat)%jrho_a_s_iii)
         DEALLOCATE(jrho_atom_set(iat)%jrho_b_h_iii)
         DEALLOCATE(jrho_atom_set(iat)%jrho_b_s_iii)
         DEALLOCATE(jrho_atom_set(iat)%jrho_c_h_iii)
         DEALLOCATE(jrho_atom_set(iat)%jrho_c_s_iii)
       END IF


       IF(ASSOCIATED(jrho_atom_set(iat)%jrho_vec_rad_h)) THEN
         NULLIFY(jrho_atom_set(iat)%jrho_rad_h,jrho_atom_set(iat)%jrho_rad_s)
         n = SIZE(jrho_atom_set(iat)%jrho_vec_rad_h,3)
         DO i = 1,n
           DO iB=1,3
             DO idir = 1,3
               DEALLOCATE(jrho_atom_set(iat)%jrho_vec_rad_h(idir,iB,i)%r_coef,STAT=istat)
               CPPrecondition(istat==0, cp_failure_level,routineN,error,failure)
               DEALLOCATE(jrho_atom_set(iat)%jrho_vec_rad_s(idir,iB,i)%r_coef,STAT=istat)
               CPPrecondition(istat==0, cp_failure_level,routineN,error,failure)
             END DO
           END DO
         END DO
         DEALLOCATE(jrho_atom_set(iat)%jrho_vec_rad_h,STAT=istat)
         CPPrecondition(istat==0, cp_failure_level,routineN,error,failure)
         DEALLOCATE(jrho_atom_set(iat)%jrho_vec_rad_s,STAT=istat)
         CPPrecondition(istat==0, cp_failure_level,routineN,error,failure)
       ELSE
         n = SIZE(jrho_atom_set(iat)%jrho_rad_h,1)
         DO i = 1,n
            DEALLOCATE(jrho_atom_set(iat)%jrho_rad_h(i)%r_coef,STAT=istat)
            CPPrecondition(istat==0, cp_failure_level,routineN,error,failure)
            DEALLOCATE(jrho_atom_set(iat)%jrho_rad_s(i)%r_coef,STAT=istat)
            CPPrecondition(istat==0, cp_failure_level,routineN,error,failure)
         END DO
         DEALLOCATE(jrho_atom_set(iat)%jrho_rad_h,STAT=istat)
         CPPrecondition(istat==0, cp_failure_level,routineN,error,failure)
         DEALLOCATE(jrho_atom_set(iat)%jrho_rad_s,STAT=istat)
         CPPrecondition(istat==0, cp_failure_level,routineN,error,failure)
       END IF
     END DO
     DEALLOCATE(jrho_atom_set, STAT=istat)
     CPPrecondition(istat==0, cp_failure_level,routineN,error,failure)

  END SUBROUTINE deallocate_jrho_atom_set

  SUBROUTINE  allocate_jrho_atom_rad(jrho1_atom,ispin,nr,max_iso_not0,error)
  
    TYPE(jrho_atom_type), POINTER            :: jrho1_atom
    INTEGER, INTENT(IN)                      :: ispin, nr, max_iso_not0
    TYPE(cp_error_type), INTENT(inout)  &
                                             :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'allocate_jrho_atom_rad', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, istat
    LOGICAL                                  :: failure

    CALL timeset("calculate_jrho_atom_rad","I"," ",handle)
    CPPrecondition(ASSOCIATED(jrho1_atom),cp_failure_level,routineP,error,failure)
    IF(.NOT.failure) THEN
         ALLOCATE (jrho1_atom%jrho_a_h(ispin)%r_coef(nr,max_iso_not0),STAT=istat)
         CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
         jrho1_atom%jrho_a_h(ispin)%r_coef = 0.0_dp
         ALLOCATE (jrho1_atom%jrho_a_s(ispin)%r_coef(nr,max_iso_not0),STAT=istat)
         CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
         jrho1_atom%jrho_a_s(ispin)%r_coef = 0.0_dp

         ALLOCATE (jrho1_atom%jrho_b_h(ispin)%r_coef(nr,max_iso_not0),STAT=istat)
         CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
         jrho1_atom%jrho_b_h(ispin)%r_coef = 0.0_dp
         ALLOCATE (jrho1_atom%jrho_b_s(ispin)%r_coef(nr,max_iso_not0),STAT=istat)
         CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
         jrho1_atom%jrho_b_s(ispin)%r_coef = 0.0_dp

         ALLOCATE (jrho1_atom%jrho_a_h_ii(ispin)%r_coef(nr,max_iso_not0),STAT=istat)
         CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
         jrho1_atom%jrho_a_h_ii(ispin)%r_coef = 0.0_dp
         ALLOCATE (jrho1_atom%jrho_a_s_ii(ispin)%r_coef(nr,max_iso_not0),STAT=istat)
         CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
         jrho1_atom%jrho_a_s_ii(ispin)%r_coef = 0.0_dp

         ALLOCATE (jrho1_atom%jrho_b_h_ii(ispin)%r_coef(nr,max_iso_not0),STAT=istat)
         CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
         jrho1_atom%jrho_b_h_ii(ispin)%r_coef = 0.0_dp
         ALLOCATE (jrho1_atom%jrho_b_s_ii(ispin)%r_coef(nr,max_iso_not0),STAT=istat)
         CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
         jrho1_atom%jrho_b_s_ii(ispin)%r_coef = 0.0_dp

         ALLOCATE (jrho1_atom%jrho_c_h_ii(ispin)%r_coef(nr,max_iso_not0),STAT=istat)
         CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
         jrho1_atom%jrho_c_h_ii(ispin)%r_coef = 0.0_dp
         ALLOCATE (jrho1_atom%jrho_c_s_ii(ispin)%r_coef(nr,max_iso_not0),STAT=istat)
         CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
         jrho1_atom%jrho_c_s_ii(ispin)%r_coef = 0.0_dp


         ALLOCATE (jrho1_atom%jrho_a_h_iii(ispin)%r_coef(nr,max_iso_not0),STAT=istat)
         CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
         jrho1_atom%jrho_a_h_iii(ispin)%r_coef = 0.0_dp
         ALLOCATE (jrho1_atom%jrho_a_s_iii(ispin)%r_coef(nr,max_iso_not0),STAT=istat)
         CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
         jrho1_atom%jrho_a_s_iii(ispin)%r_coef = 0.0_dp

         ALLOCATE (jrho1_atom%jrho_b_h_iii(ispin)%r_coef(nr,max_iso_not0),STAT=istat)
         CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
         jrho1_atom%jrho_b_h_iii(ispin)%r_coef = 0.0_dp
         ALLOCATE (jrho1_atom%jrho_b_s_iii(ispin)%r_coef(nr,max_iso_not0),STAT=istat)
         CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
         jrho1_atom%jrho_b_s_iii(ispin)%r_coef = 0.0_dp

         ALLOCATE (jrho1_atom%jrho_c_h_iii(ispin)%r_coef(nr,max_iso_not0),STAT=istat)
         CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
         jrho1_atom%jrho_c_h_iii(ispin)%r_coef = 0.0_dp
         ALLOCATE (jrho1_atom%jrho_c_s_iii(ispin)%r_coef(nr,max_iso_not0),STAT=istat)
         CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
         jrho1_atom%jrho_c_s_iii(ispin)%r_coef = 0.0_dp

    END IF
    CALL timestop(0.0_dp,handle)

  END SUBROUTINE allocate_jrho_atom_rad

! *****************************************************************************

  SUBROUTINE set2zero_jrho_atom_rad(jrho1_atom,ispin,error)
  
    TYPE(jrho_atom_type), POINTER            :: jrho1_atom
    INTEGER, INTENT(IN)                      :: ispin
    TYPE(cp_error_type), INTENT(inout)  &
                                             :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'set2zero_jrho_atom_rad', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    CPPrecondition(ASSOCIATED(jrho1_atom),cp_failure_level,routineP,error,failure)
    IF(.NOT.failure) THEN

      jrho1_atom%jrho_a_h(ispin)%r_coef = 0.0_dp
      jrho1_atom%jrho_a_s(ispin)%r_coef = 0.0_dp
      jrho1_atom%jrho_b_h(ispin)%r_coef = 0.0_dp
      jrho1_atom%jrho_b_s(ispin)%r_coef = 0.0_dp

      jrho1_atom%jrho_a_h_ii(ispin)%r_coef = 0.0_dp
      jrho1_atom%jrho_a_s_ii(ispin)%r_coef = 0.0_dp
      jrho1_atom%jrho_b_h_ii(ispin)%r_coef = 0.0_dp
      jrho1_atom%jrho_b_s_ii(ispin)%r_coef = 0.0_dp
      jrho1_atom%jrho_c_h_ii(ispin)%r_coef = 0.0_dp
      jrho1_atom%jrho_c_s_ii(ispin)%r_coef = 0.0_dp

      jrho1_atom%jrho_a_h_iii(ispin)%r_coef = 0.0_dp
      jrho1_atom%jrho_a_s_iii(ispin)%r_coef = 0.0_dp
      jrho1_atom%jrho_b_h_iii(ispin)%r_coef = 0.0_dp
      jrho1_atom%jrho_b_s_iii(ispin)%r_coef = 0.0_dp
      jrho1_atom%jrho_c_h_iii(ispin)%r_coef = 0.0_dp
      jrho1_atom%jrho_c_s_iii(ispin)%r_coef = 0.0_dp

    END IF

  END SUBROUTINE set2zero_jrho_atom_rad

! *****************************************************************************


  SUBROUTINE get_jrho_atom(jrho1_atom_set,iatom,cjc_h,cjc_s,cjc_ii_h,cjc_ii_s,&
             cjc_iii_h,cjc_iii_s,jrho_rad_h,jrho_rad_s,jrho_vec_rad_h,jrho_vec_rad_s,error)

    TYPE(jrho_atom_type), DIMENSION(:), &
      POINTER                                :: jrho1_atom_set
    INTEGER, INTENT(IN)                      :: iatom
    TYPE(rho_atom_coeff), DIMENSION(:), &
      OPTIONAL, POINTER                      :: cjc_h, cjc_s, cjc_ii_h, &
                                                cjc_ii_s, cjc_iii_h, &
                                                cjc_iii_s, jrho_rad_h, &
                                                jrho_rad_s
    TYPE(rho_atom_coeff), &
      DIMENSION(:, :, :), OPTIONAL, POINTER  :: jrho_vec_rad_h, jrho_vec_rad_s
    TYPE(cp_error_type), INTENT(inout)  &
                                             :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_jrho_atom', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

     CPPrecondition(ASSOCIATED(jrho1_atom_set), cp_failure_level,routineN,error,failure)

     IF(PRESENT(cjc_h)) cjc_h => jrho1_atom_set(iatom)%cjc_h
     IF(PRESENT(cjc_s)) cjc_s => jrho1_atom_set(iatom)%cjc_s
     IF(PRESENT(cjc_ii_h)) cjc_ii_h => jrho1_atom_set(iatom)%cjc_ii_h
     IF(PRESENT(cjc_ii_s)) cjc_ii_s => jrho1_atom_set(iatom)%cjc_ii_s
     IF(PRESENT(cjc_iii_h)) cjc_iii_h => jrho1_atom_set(iatom)%cjc_iii_h
     IF(PRESENT(cjc_iii_s)) cjc_iii_s => jrho1_atom_set(iatom)%cjc_iii_s
     IF(PRESENT(jrho_rad_h)) jrho_rad_h => jrho1_atom_set(iatom)%jrho_rad_h
     IF(PRESENT(jrho_rad_s)) jrho_rad_s => jrho1_atom_set(iatom)%jrho_rad_s
     IF(PRESENT(jrho_vec_rad_h)) jrho_vec_rad_h => jrho1_atom_set(iatom)%jrho_vec_rad_h
     IF(PRESENT(jrho_vec_rad_s)) jrho_vec_rad_s => jrho1_atom_set(iatom)%jrho_vec_rad_s

  END SUBROUTINE get_jrho_atom

! *****************************************************************************

  SUBROUTINE init_jrho_atom_set(jrho1_atom_set,atomic_kind_set,nspins,store_current,error)

    TYPE(jrho_atom_type), DIMENSION(:), &
      POINTER                                :: jrho1_atom_set
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    INTEGER, INTENT(IN)                      :: nspins
    LOGICAL, INTENT(IN)                      :: store_current
    TYPE(cp_error_type), INTENT(inout)  &
                                             :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'init_jrho_atom_set'

    INTEGER :: handle, iat, iatom, iB, idir, ikind, ispin, istat, &
      max_iso_not0, maxso, na, nat, natom, nkind, nr, nset, nsotot
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: failure, paw_atom
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(harmonics_atom_type), POINTER       :: harmonics

    CALL timeset("init_jrho_atom_set","I"," ",handle)

    failure = .FALSE.

    CPPrecondition(ASSOCIATED(atomic_kind_set),cp_failure_level,routineN,error,failure)

    IF(ASSOCIATED(jrho1_atom_set)) THEN
      CALL deallocate_jrho_atom_set(jrho1_atom_set,error=error)
    END IF

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             natom = natom)

    CALL allocate_jrho_atom_set(jrho1_atom_set,natom,error=error)

    nkind = SIZE(atomic_kind_set)

    DO ikind = 1,nkind

       atomic_kind => atomic_kind_set(ikind)
       CALL get_atomic_kind(atomic_kind=atomic_kind,&
                            orb_basis_set=orb_basis_set, &
                            atom_list=atom_list,natom=nat, &
                            paw_atom=paw_atom,&
                            harmonics=harmonics,&
                            ngrid_rad=nr,&
                            ngrid_ang=na)

       CALL get_gto_basis_set(gto_basis_set=orb_basis_set, &
                              maxso=maxso, nset=nset)
       nsotot = maxso * nset 
       max_iso_not0 = harmonics%max_iso_not0
       DO iat = 1,nat
          iatom = atom_list(iat) 
          !*** allocate the radial density for each LM,for each atom ***

          IF(store_current) THEN
            ALLOCATE (jrho1_atom_set(iatom)%jrho_vec_rad_h(3,3,nspins),STAT=istat)
            CPPrecondition(istat==0,cp_failure_level,routineN,error,failure)
            ALLOCATE (jrho1_atom_set(iatom)%jrho_vec_rad_s(3,3,nspins),STAT=istat)
            CPPrecondition(istat==0,cp_failure_level,routineN,error,failure)
            DO ispin = 1,nspins
              DO iB = 1,3
                DO idir = 1,3
                  NULLIFY(jrho1_atom_set(iatom)%jrho_vec_rad_h(idir,iB,ispin)%r_coef)
                  NULLIFY(jrho1_atom_set(iatom)%jrho_vec_rad_s(idir,iB,ispin)%r_coef)
                  ALLOCATE(jrho1_atom_set(iatom)%jrho_vec_rad_h(idir,iB,ispin)%r_coef(nr,na),STAT=istat)
                  CPPrecondition(istat==0,cp_failure_level,routineN,error,failure)
                  ALLOCATE(jrho1_atom_set(iatom)%jrho_vec_rad_s(idir,iB,ispin)%r_coef(nr,na),STAT=istat)
                  CPPrecondition(istat==0,cp_failure_level,routineN,error,failure)
                END DO
              END DO
            END DO  ! ispin
          ELSE
            ALLOCATE (jrho1_atom_set(iatom)%jrho_rad_h(nspins),STAT=istat)
            CPPrecondition(istat==0,cp_failure_level,routineN,error,failure)
            ALLOCATE (jrho1_atom_set(iatom)%jrho_rad_s(nspins),STAT=istat)
            CPPrecondition(istat==0,cp_failure_level,routineN,error,failure)
 
            DO ispin = 1,nspins
               NULLIFY(jrho1_atom_set(iatom)%jrho_rad_h(ispin)%r_coef)
               NULLIFY(jrho1_atom_set(iatom)%jrho_rad_s(ispin)%r_coef)

               ALLOCATE(jrho1_atom_set(iatom)%jrho_rad_h(ispin)%r_coef(nr,na),STAT=istat)
              CPPrecondition(istat==0,cp_failure_level,routineN,error,failure)
               jrho1_atom_set(iatom)%jrho_rad_h(ispin)%r_coef = 0.0_dp
               ALLOCATE(jrho1_atom_set(iatom)%jrho_rad_s(ispin)%r_coef(nr,na),STAT=istat)
              CPPrecondition(istat==0,cp_failure_level,routineN,error,failure)
               jrho1_atom_set(iatom)%jrho_rad_s(ispin)%r_coef = 0.0_dp
            END DO  ! ispin 
          END IF

          ALLOCATE (jrho1_atom_set(iatom)%jrho_a_h(nspins),STAT=istat)
         CPPrecondition(istat==0,cp_failure_level,routineN,error,failure)
          ALLOCATE (jrho1_atom_set(iatom)%jrho_a_s(nspins),STAT=istat)
         CPPrecondition(istat==0,cp_failure_level,routineN,error,failure)
          ALLOCATE (jrho1_atom_set(iatom)%jrho_b_h(nspins),STAT=istat)
         CPPrecondition(istat==0,cp_failure_level,routineN,error,failure)
          ALLOCATE (jrho1_atom_set(iatom)%jrho_b_s(nspins),STAT=istat)
         CPPrecondition(istat==0,cp_failure_level,routineN,error,failure)

          ALLOCATE (jrho1_atom_set(iatom)%jrho_a_h_ii(nspins),STAT=istat)
         CPPrecondition(istat==0,cp_failure_level,routineN,error,failure)
          ALLOCATE (jrho1_atom_set(iatom)%jrho_a_s_ii(nspins),STAT=istat)
         CPPrecondition(istat==0,cp_failure_level,routineN,error,failure)
          ALLOCATE (jrho1_atom_set(iatom)%jrho_b_s_ii(nspins),STAT=istat)
         CPPrecondition(istat==0,cp_failure_level,routineN,error,failure)
          ALLOCATE (jrho1_atom_set(iatom)%jrho_b_h_ii(nspins),STAT=istat)
         CPPrecondition(istat==0,cp_failure_level,routineN,error,failure)
          ALLOCATE (jrho1_atom_set(iatom)%jrho_c_s_ii(nspins),STAT=istat)
         CPPrecondition(istat==0,cp_failure_level,routineN,error,failure)
          ALLOCATE (jrho1_atom_set(iatom)%jrho_c_s_ii(nspins),STAT=istat)
         CPPrecondition(istat==0,cp_failure_level,routineN,error,failure)

          ALLOCATE (jrho1_atom_set(iatom)%jrho_a_h_iii(nspins),STAT=istat)
         CPPrecondition(istat==0,cp_failure_level,routineN,error,failure)
          ALLOCATE (jrho1_atom_set(iatom)%jrho_a_s_iii(nspins),STAT=istat)
         CPPrecondition(istat==0,cp_failure_level,routineN,error,failure)
          ALLOCATE (jrho1_atom_set(iatom)%jrho_b_s_iii(nspins),STAT=istat)
         CPPrecondition(istat==0,cp_failure_level,routineN,error,failure)
          ALLOCATE (jrho1_atom_set(iatom)%jrho_b_h_iii(nspins),STAT=istat)
         CPPrecondition(istat==0,cp_failure_level,routineN,error,failure)
          ALLOCATE (jrho1_atom_set(iatom)%jrho_c_s_iii(nspins),STAT=istat)
         CPPrecondition(istat==0,cp_failure_level,routineN,error,failure)
          ALLOCATE (jrho1_atom_set(iatom)%jrho_c_s_iii(nspins),STAT=istat)
         CPPrecondition(istat==0,cp_failure_level,routineN,error,failure)

          ALLOCATE (jrho1_atom_set(iatom)%cjc_h(nspins), STAT=istat) 
         CPPrecondition(istat==0,cp_failure_level,routineN,error,failure)
          ALLOCATE (jrho1_atom_set(iatom)%cjc_s(nspins), STAT=istat)
         CPPrecondition(istat==0,cp_failure_level,routineN,error,failure)
          ALLOCATE (jrho1_atom_set(iatom)%cjc_ii_h(nspins), STAT=istat)
         CPPrecondition(istat==0,cp_failure_level,routineN,error,failure)
          ALLOCATE (jrho1_atom_set(iatom)%cjc_ii_s(nspins), STAT=istat)
         CPPrecondition(istat==0,cp_failure_level,routineN,error,failure)
          ALLOCATE (jrho1_atom_set(iatom)%cjc_iii_h(nspins), STAT=istat)
         CPPrecondition(istat==0,cp_failure_level,routineN,error,failure)
          ALLOCATE (jrho1_atom_set(iatom)%cjc_iii_s(nspins), STAT=istat)
         CPPrecondition(istat==0,cp_failure_level,routineN,error,failure)

          IF (paw_atom) THEN 

            DO ispin = 1,nspins

              NULLIFY(jrho1_atom_set(iatom)%jrho_a_h(ispin)%r_coef)
              NULLIFY(jrho1_atom_set(iatom)%jrho_a_s(ispin)%r_coef)
              NULLIFY(jrho1_atom_set(iatom)%jrho_b_h(ispin)%r_coef)
              NULLIFY(jrho1_atom_set(iatom)%jrho_b_s(ispin)%r_coef)

              NULLIFY(jrho1_atom_set(iatom)%jrho_a_h_ii(ispin)%r_coef)
              NULLIFY(jrho1_atom_set(iatom)%jrho_a_s_ii(ispin)%r_coef)
              NULLIFY(jrho1_atom_set(iatom)%jrho_b_h_ii(ispin)%r_coef)
              NULLIFY(jrho1_atom_set(iatom)%jrho_b_s_ii(ispin)%r_coef)
              NULLIFY(jrho1_atom_set(iatom)%jrho_c_h_ii(ispin)%r_coef)
              NULLIFY(jrho1_atom_set(iatom)%jrho_c_s_ii(ispin)%r_coef)

              NULLIFY(jrho1_atom_set(iatom)%jrho_a_h_iii(ispin)%r_coef)
              NULLIFY(jrho1_atom_set(iatom)%jrho_a_s_iii(ispin)%r_coef)
              NULLIFY(jrho1_atom_set(iatom)%jrho_b_h_iii(ispin)%r_coef)
              NULLIFY(jrho1_atom_set(iatom)%jrho_b_s_iii(ispin)%r_coef)
              NULLIFY(jrho1_atom_set(iatom)%jrho_c_h_iii(ispin)%r_coef)
              NULLIFY(jrho1_atom_set(iatom)%jrho_c_s_iii(ispin)%r_coef)

              NULLIFY(jrho1_atom_set(iatom)%cjc_h(ispin)%r_coef)
              ALLOCATE(jrho1_atom_set(iatom)%cjc_h(ispin)%r_coef(1:nsotot,1:nsotot),STAT=istat)
              CPPrecondition(istat==0,cp_failure_level,routineN,error,failure) 
              NULLIFY(jrho1_atom_set(iatom)%cjc_s(ispin)%r_coef)
              ALLOCATE(jrho1_atom_set(iatom)%cjc_s(ispin)%r_coef(1:nsotot,1:nsotot),STAT=istat)
              CPPrecondition(istat==0,cp_failure_level,routineN,error,failure) 
              NULLIFY(jrho1_atom_set(iatom)%cjc_ii_h(ispin)%r_coef)
              ALLOCATE(jrho1_atom_set(iatom)%cjc_ii_h(ispin)%r_coef(1:nsotot,1:nsotot),STAT=istat)
              CPPrecondition(istat==0,cp_failure_level,routineN,error,failure) 
              NULLIFY(jrho1_atom_set(iatom)%cjc_ii_s(ispin)%r_coef)
              ALLOCATE(jrho1_atom_set(iatom)%cjc_ii_s(ispin)%r_coef(1:nsotot,1:nsotot),STAT=istat)
              CPPrecondition(istat==0,cp_failure_level,routineN,error,failure) 
              NULLIFY(jrho1_atom_set(iatom)%cjc_iii_h(ispin)%r_coef)
              ALLOCATE(jrho1_atom_set(iatom)%cjc_iii_h(ispin)%r_coef(1:nsotot,1:nsotot),STAT=istat)
              CPPrecondition(istat==0,cp_failure_level,routineN,error,failure) 
              NULLIFY(jrho1_atom_set(iatom)%cjc_iii_s(ispin)%r_coef)
              ALLOCATE(jrho1_atom_set(iatom)%cjc_iii_s(ispin)%r_coef(1:nsotot,1:nsotot),STAT=istat)
              CPPrecondition(istat==0,cp_failure_level,routineN,error,failure) 


            END DO  ! ispin
          ELSE
             DO ispin = 1,nspins
                NULLIFY(jrho1_atom_set(iatom)%cjc_h(ispin)%r_coef)
                NULLIFY(jrho1_atom_set(iatom)%cjc_s(ispin)%r_coef)
                NULLIFY(jrho1_atom_set(iatom)%cjc_ii_h(ispin)%r_coef)
                NULLIFY(jrho1_atom_set(iatom)%cjc_ii_s(ispin)%r_coef)
                NULLIFY(jrho1_atom_set(iatom)%cjc_iii_h(ispin)%r_coef)
                NULLIFY(jrho1_atom_set(iatom)%cjc_iii_s(ispin)%r_coef)
             END DO  ! ispin

          END IF  ! paw_atom

       END DO  ! iat

    END DO  ! ikind

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE init_jrho_atom_set

END MODULE qs_linres_types

