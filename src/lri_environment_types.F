!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief contains the types and subroutines for dealing with the lri_env
!>        lri : local resolution of the identity
!> \author JGH August 2012
! *****************************************************************************
MODULE lri_environment_types
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             set_atomic_kind
  USE basis_set_types,                 ONLY: allocate_gto_basis_set,&
                                             deallocate_gto_basis_set,&
                                             get_gto_basis_set,&
                                             gto_basis_set_p_type,&
                                             gto_basis_set_type,&
                                             init_orb_basis_set,&
                                             write_gto_basis_set
  USE cp_control_types,                ONLY: lripaw_control_type
  USE cp_dbcsr_operations,             ONLY: cp_dbcsr_deallocate_matrix_set
  USE cp_dbcsr_types,                  ONLY: cp_dbcsr_p_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: dp
  USE lebedev,                         ONLY: deallocate_lebedev_grids,&
                                             get_number_of_lebedev_grid,&
                                             init_lebedev_grids,&
                                             lebedev_grid
  USE mathconstants,                   ONLY: dfac,&
                                             pi,&
                                             rootpi,&
                                             twopi
  USE mathlib,                         ONLY: invert_matrix
  USE orbital_pointers,                ONLY: indco,&
                                             indso,&
                                             init_orbital_pointers,&
                                             nco,&
                                             ncoset,&
                                             nso,&
                                             nsoset
  USE orbital_symbols,                 ONLY: cgf_symbol,&
                                             sgf_symbol
  USE qs_grid_atom,                    ONLY: allocate_grid_atom,&
                                             create_grid_atom,&
                                             deallocate_grid_atom,&
                                             grid_atom_type
  USE qs_interactions,                 ONLY: init_interaction_radii_orb_basis
  USE qs_neighbor_list_types,          ONLY: deallocate_neighbor_list_set,&
                                             neighbor_list_set_p_type
  USE qs_util,                         ONLY: exp_radius
  USE soft_basis_set,                  ONLY: create_soft_basis
  USE spherical_harmonics,             ONLY: drvy_lm,&
                                             y_lm
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

! *****************************************************************************
  TYPE lri_rhoab_type
     ! number of basis functions (a)
     INTEGER                                                 :: nba
     ! number of basis functions (b)
     INTEGER                                                 :: nbb
     ! number of basis functions (ai)
     INTEGER                                                 :: nfa
     ! number of basis functions (bi)
     INTEGER                                                 :: nfb
     ! expansion coeffs for RI density
     REAL(KIND=dp), DIMENSION(:), POINTER                    :: avec
     ! projection coeffs for RI density: SUM_ab (ab,i)*Pab
     REAL(KIND=dp), DIMENSION(:), POINTER                    :: tvec
     ! Lagrange parameter
     REAL(KIND=dp)                                           :: lambda
     ! Charge of pair density
     REAL(KIND=dp)                                           :: charge
  END TYPE lri_rhoab_type

  TYPE lri_int_type
     ! number of basis functions (a)
     INTEGER                                                 :: nba
     ! number of basis functions (b)
     INTEGER                                                 :: nbb
     ! number of basis functions (ai)
     INTEGER                                                 :: nfa
     ! number of basis functions (bi)
     INTEGER                                                 :: nfb
     ! integrals (a,b,ai)
     REAL(KIND=dp), DIMENSION(:,:,:), POINTER                :: abaint
     ! integrals (a,b,bi)
     REAL(KIND=dp), DIMENSION(:,:,:), POINTER                :: abbint
     ! integrals (a,b)
     REAL(KIND=dp), DIMENSION(:,:), POINTER                  :: soo
     ! integrals (ai,bi)
     REAL(KIND=dp), DIMENSION(:,:), POINTER                  :: sab
     ! integrals -SUM_i (a,b,ai) sn(i)/nsn
     REAL(KIND=dp), DIMENSION(:,:), POINTER                  :: qab
     ! inverse of integrals (ai,bi)
     REAL(KIND=dp), DIMENSION(:,:), POINTER                  :: sinv
     ! sinv * (ai)
     REAL(KIND=dp), DIMENSION(:), POINTER                    :: sn
     ! (ai) * sinv * (ai)
     REAL(KIND=dp)                                           :: nsn
  END TYPE lri_int_type

  TYPE lri_paw_type
     INTEGER                                                 :: nb
     REAL(KIND=dp), DIMENSION(:), POINTER                    :: fhh, fhs
     REAL(KIND=dp), DIMENSION(:), POINTER                    :: fxh, fxs
     INTEGER                                                 :: na
     REAL(KIND=dp), DIMENSION(:), POINTER                    :: flhh, flhs
     REAL(KIND=dp), DIMENSION(:), POINTER                    :: flxh, flxs
  END TYPE lri_paw_type

  TYPE lri_int_node
     INTEGER                                                 :: nnode
     TYPE(lri_int_type), DIMENSION(:), POINTER               :: lri_int
     TYPE(lri_rhoab_type), DIMENSION(:), POINTER             :: lri_rhoab
     TYPE(lri_paw_type), POINTER                             :: lri_paw
     REAL(KIND=dp), DIMENSION(:), POINTER                    :: acoef
     REAL(KIND=dp), DIMENSION(:), POINTER                    :: fcoef
  END TYPE lri_int_node

  TYPE lri_int_atom
     INTEGER                                                 :: natom
     TYPE(lri_int_node), DIMENSION(:), POINTER               :: lri_node
  END TYPE lri_int_atom

  TYPE lri_int_list
     INTEGER                                                 :: nkind
     TYPE(lri_int_atom), DIMENSION(:), POINTER               :: lri_atom
  END TYPE lri_int_list

  TYPE lri_int_p_list
     TYPE(lri_int_list),POINTER                              :: int_list
  END TYPE lri_int_p_list
! *****************************************************************************
  TYPE lri_bas_norm
     REAL(KIND=dp), DIMENSION(:), POINTER                    :: norm
  END TYPE lri_bas_norm
  TYPE lri_bas_overlap
     REAL(KIND=dp), DIMENSION(:,:), POINTER                  :: ovlp
  END TYPE lri_bas_overlap
  TYPE harmonics_info_type
     INTEGER                                                 :: lmax
     INTEGER                                                 :: maxs
     INTEGER                                                 :: ngrid
     REAL(dp), DIMENSION(:,:),    POINTER                    :: coord
     REAL(dp), DIMENSION(:,:),    POINTER                    :: slm
     REAL(dp), DIMENSION(:,:,:),  POINTER                    :: dslm
  END TYPE harmonics_info_type
  !
  TYPE lri_atom_info
     ! atomic grids for PAW
     TYPE(grid_atom_type), POINTER                           :: grid_atom
     ! spherical harmoincs 
     TYPE(harmonics_info_type), POINTER                      :: harmonics
  END TYPE lri_atom_info
! *****************************************************************************
  TYPE lri_environment_type
     INTEGER :: id_nr, ref_count, in_use
     TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER       :: orb_basis
     TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER       :: ri_basis
     TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER       :: ri_basis_soft
     TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER       :: paw_basis
     TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER       :: proj_basis
     ! info on atoms
     TYPE(lri_atom_info), DIMENSION(:), POINTER              :: atom
     ! orb_basis neighborlist, used for RI projection integrals
     TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER   :: soo_list
     ! neighborlist for <ri_basis_soft|proj_basis> integrals
     TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER   :: sip_list
     ! Overlap matrix ri basis/paw projectors
     TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER            :: rip_matrix
     ! local RI projection integrals
     TYPE(lri_int_list), POINTER                             :: lri_ints
     ! norm of ri basis
     TYPE(lri_bas_norm), DIMENSION(:), POINTER               :: bas_norm
     ! self overlap of ri basis
     TYPE(lri_bas_overlap), DIMENSION(:), POINTER            :: bas_ovlp
  END TYPE lri_environment_type

  TYPE lri_environment_p_type
     TYPE(lri_environment_type),POINTER :: lri_env
  END TYPE lri_environment_p_type
! *****************************************************************************
  TYPE lri_density_type
     INTEGER :: id_nr, ref_count, in_use
     INTEGER                                                 :: nspin
     ! pair density expansion (nspin)
     TYPE(lri_int_p_list), DIMENSION(:),POINTER              :: lri_rhos
     ! coefficients of RI expansion and gradients (nspin)
     TYPE(lri_int_p_list), DIMENSION(:),POINTER              :: lri_coef
     ! coefficients of PAW expansion and gradients (nspin)
     TYPE(lri_int_p_list), DIMENSION(:),POINTER              :: lri_paw
  END TYPE lri_density_type

  TYPE lri_density_p_type
     TYPE(lri_density_type),POINTER :: lri_density
  END TYPE lri_density_p_type
! *****************************************************************************

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'lri_environment_types'
  INTEGER, SAVE, PRIVATE :: last_lri_env_id=0
  INTEGER, SAVE, PRIVATE :: last_lri_density_id=0

  PUBLIC :: lri_environment_type, lri_environment_p_type, &
            lri_env_create, lri_env_release, &
            lri_int_type, lri_int_node, lri_int_atom, lri_int_list, &
            deallocate_lri_ints, &
            lri_atom_init, lri_basis_init, write_lripaw_basis_sets
  PUBLIC :: lri_density_type, lri_density_p_type, &
            lri_density_create, lri_density_release, lri_rhoab_type

! *****************************************************************************

CONTAINS

! *****************************************************************************
!> \brief creates and initializes an lri_env
!> \param lri_env the lri_environment you want to create
!>
!> \author JGH
! *****************************************************************************
  SUBROUTINE lri_env_create ( lri_env, error )

    TYPE(lri_environment_type), POINTER      :: lri_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'lri_env_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat

    ALLOCATE ( lri_env, stat=istat )
    CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)

    last_lri_env_id=last_lri_env_id+1
    lri_env%id_nr=last_lri_env_id
    lri_env%ref_count=1
    lri_env%in_use=0
  
    NULLIFY(lri_env%orb_basis)
    NULLIFY(lri_env%ri_basis)
    NULLIFY(lri_env%ri_basis_soft)
    NULLIFY(lri_env%paw_basis)
    NULLIFY(lri_env%proj_basis)

    NULLIFY(lri_env%atom)

    NULLIFY(lri_env%soo_list)
    NULLIFY(lri_env%sip_list)

    NULLIFY(lri_env%rip_matrix)
    NULLIFY(lri_env%lri_ints)

    NULLIFY(lri_env%bas_norm)
    NULLIFY(lri_env%bas_ovlp)

  END SUBROUTINE lri_env_create

! *****************************************************************************
!> \brief releases the given lri_env
!> \param lri_env the lri environment to release
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!>
!> \author JGH
! *****************************************************************************
  SUBROUTINE lri_env_release(lri_env, error)
    TYPE(lri_environment_type), POINTER      :: lri_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'lri_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ikind, nkind, stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(lri_env)) THEN
       lri_env%ref_count=0

       ! deallocate basis sets
       IF(ASSOCIATED(lri_env%orb_basis)) THEN
          nkind = SIZE(lri_env%orb_basis)
          DO ikind=1,nkind
             CALL deallocate_gto_basis_set(lri_env%orb_basis(ikind)%gto_basis_set,error)
          END DO
          DEALLOCATE(lri_env%orb_basis,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
       IF(ASSOCIATED(lri_env%ri_basis)) THEN
          nkind = SIZE(lri_env%ri_basis)
          DO ikind=1,nkind
             CALL deallocate_gto_basis_set(lri_env%ri_basis(ikind)%gto_basis_set,error)
          END DO
          DEALLOCATE(lri_env%ri_basis,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
       IF(ASSOCIATED(lri_env%ri_basis_soft)) THEN
          nkind = SIZE(lri_env%ri_basis_soft)
          DO ikind=1,nkind
             CALL deallocate_gto_basis_set(lri_env%ri_basis_soft(ikind)%gto_basis_set,error)
          END DO
          DEALLOCATE(lri_env%ri_basis_soft,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
       IF(ASSOCIATED(lri_env%paw_basis)) THEN
          nkind = SIZE(lri_env%paw_basis)
          DO ikind=1,nkind
             CALL deallocate_gto_basis_set(lri_env%paw_basis(ikind)%gto_basis_set,error)
          END DO
          DEALLOCATE(lri_env%paw_basis,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
       IF(ASSOCIATED(lri_env%proj_basis)) THEN
          nkind = SIZE(lri_env%proj_basis)
          DO ikind=1,nkind
             CALL deallocate_gto_basis_set(lri_env%proj_basis(ikind)%gto_basis_set,error)
          END DO
          DEALLOCATE(lri_env%proj_basis,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF

       IF(ASSOCIATED(lri_env%atom)) THEN
          nkind = SIZE(lri_env%atom)
          DO ikind=1,nkind
             CALL deallocate_grid_atom(lri_env%atom(ikind)%grid_atom,error)
             CALL deallocate_harmonics(lri_env%atom(ikind)%harmonics,error)
          END DO
          DEALLOCATE(lri_env%atom,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
  
       IF (ASSOCIATED(lri_env%soo_list)) THEN
          DO i=1,SIZE(lri_env%soo_list)
             CALL deallocate_neighbor_list_set(lri_env%soo_list(i)%neighbor_list_set)
          END DO
          DEALLOCATE(lri_env%soo_list,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
       IF (ASSOCIATED(lri_env%sip_list)) THEN
          DO i=1,SIZE(lri_env%sip_list)
             CALL deallocate_neighbor_list_set(lri_env%sip_list(i)%neighbor_list_set)
          END DO
          DEALLOCATE(lri_env%sip_list,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
       IF (ASSOCIATED(lri_env%rip_matrix)) THEN
          CALL cp_dbcsr_deallocate_matrix_set(lri_env%rip_matrix,error)
       END IF
       IF (ASSOCIATED(lri_env%lri_ints)) THEN
          CALL deallocate_lri_ints(lri_env%lri_ints,error)
       END IF
       IF (ASSOCIATED(lri_env%bas_norm)) THEN
          DO i=1,SIZE(lri_env%bas_norm)
             DEALLOCATE(lri_env%bas_norm(i)%norm,stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END DO
          DEALLOCATE(lri_env%bas_norm,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
       IF (ASSOCIATED(lri_env%bas_ovlp)) THEN
          DO i=1,SIZE(lri_env%bas_ovlp)
             DEALLOCATE(lri_env%bas_ovlp(i)%ovlp,stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END DO
          DEALLOCATE(lri_env%bas_ovlp,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF

       DEALLOCATE(lri_env,stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    END IF
    NULLIFY(lri_env)

  END SUBROUTINE lri_env_release

! *****************************************************************************
!> \brief creates and initializes an lri_env
!> \param lri_env the lri_environment you want to create
!>
!> \author JGH
! *****************************************************************************
  SUBROUTINE lri_density_create ( lri_density, error )

    TYPE(lri_density_type), POINTER          :: lri_density
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'lri_density_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat

    ALLOCATE ( lri_density, stat=istat )
    CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)

    last_lri_density_id=last_lri_density_id+1
    lri_density%id_nr=last_lri_density_id
    lri_density%ref_count=1
    lri_density%in_use=0

    lri_density%nspin=0
  
    NULLIFY(lri_density%lri_rhos)
    NULLIFY(lri_density%lri_coef)
    NULLIFY(lri_density%lri_paw)

  END SUBROUTINE lri_density_create

! *****************************************************************************
!> \brief releases the given lri_env
!> \param lri_env the lri environment to release
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!>
!> \author JGH
! *****************************************************************************
  SUBROUTINE lri_density_release(lri_density, error)
    TYPE(lri_density_type), POINTER          :: lri_density
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'lri_density_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(lri_density)) THEN
       lri_density%ref_count=0

       IF (ASSOCIATED(lri_density%lri_rhos)) THEN
          DO i=1,SIZE(lri_density%lri_rhos)
             CALL deallocate_lri_ints(lri_density%lri_rhos(i)%int_list,error)
          END DO
          DEALLOCATE(lri_density%lri_rhos,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF

       IF (ASSOCIATED(lri_density%lri_coef)) THEN
          DO i=1,SIZE(lri_density%lri_coef)
             CALL deallocate_lri_ints(lri_density%lri_coef(i)%int_list,error)
          END DO
          DEALLOCATE(lri_density%lri_coef,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF

       IF (ASSOCIATED(lri_density%lri_paw)) THEN
          DO i=1,SIZE(lri_density%lri_paw)
             CALL deallocate_lri_ints(lri_density%lri_paw(i)%int_list,error)
          END DO
          DEALLOCATE(lri_density%lri_paw,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF

       DEALLOCATE(lri_density,stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    END IF
    NULLIFY(lri_density)

  END SUBROUTINE lri_density_release

! *****************************************************************************

  SUBROUTINE lri_basis_init(lri_env,lri_control,atomic_kind_set,error)
    TYPE(lri_environment_type), POINTER      :: lri_env
    TYPE(lripaw_control_type), POINTER       :: lri_control
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'lri_basis_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikind, nkind, stat
    LOGICAL                                  :: failure, paw_atom
    REAL(KIND=dp)                            :: eps_fit, eps_pgf, rc
    REAL(KIND=dp), DIMENSION(:), POINTER     :: ri_norm
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: ri_ovlp
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(gto_basis_set_type), POINTER        :: orb_basis, paw_basis, &
                                                proj_basis, ri_basis, ri_soft

    failure=.FALSE.

    IF (ASSOCIATED(lri_env)) THEN
       eps_fit = lri_control%eps_fit
       IF (ASSOCIATED(lri_env%orb_basis)) THEN
          CPPreconditionNoFail(ASSOCIATED(lri_env%ri_basis),cp_warning_level,routineP,error)
          nkind = SIZE(lri_env%orb_basis)
          ALLOCATE (lri_env%ri_basis_soft(nkind),STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE (lri_env%paw_basis(nkind),STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE (lri_env%proj_basis(nkind),STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE (lri_env%bas_norm(nkind),STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE (lri_env%bas_ovlp(nkind),STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
          DO ikind = 1, nkind
             atomic_kind => atomic_kind_set(ikind)
             NULLIFY(orb_basis,ri_basis,ri_soft,paw_basis,proj_basis)
             NULLIFY(lri_env%ri_basis_soft(ikind)%gto_basis_set)
             NULLIFY(lri_env%paw_basis(ikind)%gto_basis_set)
             NULLIFY(lri_env%proj_basis(ikind)%gto_basis_set)
             orb_basis => lri_env%orb_basis(ikind)%gto_basis_set
             IF (ASSOCIATED(orb_basis)) THEN
                ri_basis => lri_env%ri_basis(ikind)%gto_basis_set
                CPPreconditionNoFail(ASSOCIATED(ri_basis),cp_warning_level,routineP,error)
                CALL get_atomic_kind(atomic_kind=atomic_kind,hard_radius=rc)
                ! soft lri basis
                CALL allocate_gto_basis_set(ri_soft,error)
                CALL create_soft_basis(ri_basis,ri_soft,eps_fit,rc,paw_atom,&
                                       .TRUE.,.FALSE.,error)
                eps_pgf = eps_fit/100._dp
                CALL init_interaction_radii_orb_basis(ri_soft,eps_pgf,error=error)
                CPPreconditionNoFail(paw_atom,cp_warning_level,routineP,error)
                lri_env%ri_basis_soft(ikind)%gto_basis_set => ri_soft
                CALL set_atomic_kind(atomic_kind=atomic_kind,paw_atom=paw_atom)
                ! paw and corresponding proj basis
                CALL allocate_gto_basis_set(paw_basis,error)
                CALL allocate_gto_basis_set(proj_basis,error)
                CALL build_projector_basis(proj_basis,paw_basis,rc,eps_fit,error)
                lri_env%paw_basis(ikind)%gto_basis_set => paw_basis
                lri_env%proj_basis(ikind)%gto_basis_set => proj_basis
                ! norm and self-overlap of RI basis
                NULLIFY(ri_norm)
                CALL basis_norm(ri_basis,ri_norm,error)
                lri_env%bas_norm(ikind)%norm => ri_norm
                NULLIFY(ri_ovlp)
                CALL basis_ovlp(ri_basis,ri_ovlp,error)
                lri_env%bas_ovlp(ikind)%ovlp => ri_ovlp
             END IF
          END DO
       END IF
    END IF

  END SUBROUTINE lri_basis_init

! *****************************************************************************

  SUBROUTINE basis_norm(basis,norm,error)

    TYPE(gto_basis_set_type), POINTER        :: basis
    REAL(dp), DIMENSION(:), POINTER          :: norm
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'basis_norm', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ipgf, iset, isgf, ishell, l, &
                                                nbas, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: aa, cc, pp

    nbas = basis%nsgf
    ALLOCATE (norm(nbas),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    norm = 0._dp

    DO iset=1,basis%nset
       DO ishell=1,basis%nshell(iset)
          l = basis%l(ishell,iset)
          IF (l/=0) CYCLE
          DO isgf=basis%first_sgf(ishell,iset),basis%last_sgf(ishell,iset)
             DO ipgf=1,basis%npgf(iset)
                cc = basis%gcc(ipgf,ishell,iset)
                aa = basis%zet(ipgf,iset)
                pp = (0.5_dp*pi)**0.25_dp / aa**0.75_dp
                norm(isgf) = norm(isgf) + cc*pp
             END DO
          END DO
       END DO
    END DO

  END SUBROUTINE basis_norm

! *****************************************************************************

  SUBROUTINE basis_ovlp(basis,ovlp,error)

    TYPE(gto_basis_set_type), POINTER        :: basis
    REAL(dp), DIMENSION(:, :), POINTER       :: ovlp
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'basis_ovlp', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ipgf, iset, isgf, ishell, &
                                                jpgf, jset, jsgf, jshell, l, &
                                                li, lj, nbas, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: aai, aaj, cci, ccj, el, il, &
                                                oo, ppi, ppj, ppl

    nbas = basis%nsgf
    ALLOCATE (ovlp(nbas,nbas),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ovlp = 0._dp

    DO iset=1,basis%nset
       DO ishell=1,basis%nshell(iset)
          li = basis%l(ishell,iset)
          DO jset=1,basis%nset
             DO jshell=1,basis%nshell(jset)
                lj = basis%l(jshell,jset)
                IF (li == lj) THEN
                   l = li
                   ppl = 2._dp**(l+2)/SQRT(dfac(2*l+1))/twopi**0.25_dp
                   el = (2._dp*l+3._dp)/4._dp
                   il = dfac(2*l+1)/2._dp**(l+1)*SQRT(pi)
                   DO isgf=basis%first_sgf(ishell,iset),basis%last_sgf(ishell,iset)
                      DO jsgf=basis%first_sgf(jshell,jset),basis%last_sgf(jshell,jset)
                         DO ipgf=1,basis%npgf(iset)
                            cci = basis%gcc(ipgf,ishell,iset)
                            aai = basis%zet(ipgf,iset)
                            ppi = cci*ppl*aai**el
                            DO jpgf=1,basis%npgf(jset)
                               ccj = basis%gcc(jpgf,jshell,jset)
                               aaj = basis%zet(jpgf,jset)
                               ppj = ccj*ppl*aaj**el
                               oo = 0.5_dp*(aai+aaj)**(2._dp*el) 
                               ovlp(isgf,jsgf) = ovlp(isgf,jsgf) + ppi*ppj*il*oo
                            END DO
                         END DO
                      END DO
                   END DO
                END IF
             END DO
          END DO
       END DO
    END DO

  END SUBROUTINE basis_ovlp
! *****************************************************************************
!> \brief   initialize the projector-type set data set.
!> \version 1.0
! *****************************************************************************
  SUBROUTINE build_projector_basis(prj_basis,paw_basis,rc,eps_fit,error)

    TYPE(gto_basis_set_type), POINTER        :: prj_basis, paw_basis
    REAL(dp), INTENT(IN)                     :: rc, eps_fit
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'build_projector_basis', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ico, ipgf, is, iset, ishell, &
                                                jpgf, l, lshell, m, maxco, &
                                                maxl, maxpgf, maxshell, na, &
                                                ncgf, np, nset, nsgf, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: npgf, nshell
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: ll, nn
    LOGICAL                                  :: failure
    REAL(dp)                                 :: eps_pgf, myerr, prefac, &
                                                radius, za, zetmin, zp
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: sinv, smat, tmat

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(prj_basis),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(paw_basis),cp_failure_level,routineP,error,failure)

    radius = rc

    paw_basis%name = "Expansion_basis"
    prj_basis%name = "Projector_basis"
    maxl = 4 ! fix for soft expansion

    ! each l quantum number is a set
    nset = maxl + 1
    ALLOCATE (paw_basis%lmax(nset),paw_basis%lmin(nset),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (prj_basis%lmax(nset),prj_basis%lmin(nset),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    DO iset=1,nset
       paw_basis%lmax(iset) = iset - 1
       paw_basis%lmin(iset) = iset - 1
       prj_basis%lmax(iset) = iset - 1
       prj_basis%lmin(iset) = iset - 1
    END DO

    ALLOCATE (npgf(nset),nshell(nset),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    DO l=0,maxl
       npgf(l+1) = MAX(6 - l/2,1)
    END DO
    nshell = npgf
    maxshell = MAXVAL(nshell)
    maxpgf = MAXVAL(npgf)

    ALLOCATE (nn(maxshell,nset),ll(maxshell,nset),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    nn = 0
    ll = 0
    DO iset=1,nset
       DO is=1,nshell(iset)
         ll(is,iset) = paw_basis%lmax(iset)
         nn(is,iset) = paw_basis%lmax(iset) + is
       END DO
    END DO

    ALLOCATE (paw_basis%npgf(nset),paw_basis%nshell(nset),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (paw_basis%l(maxshell,nset),paw_basis%n(maxshell,nset),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    paw_basis%nset = nset
    paw_basis%npgf = npgf
    paw_basis%nshell = nshell
    paw_basis%l = ll
    paw_basis%n = nn
    paw_basis%norm_type = 1
    paw_basis%kind_radius = 0._dp
    ALLOCATE (paw_basis%set_radius(nset),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    paw_basis%set_radius = 0._dp
    ALLOCATE (paw_basis%pgf_radius(maxpgf,nset),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    paw_basis%pgf_radius = 0._dp

    ALLOCATE (prj_basis%npgf(nset),prj_basis%nshell(nset),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (prj_basis%l(maxshell,nset),prj_basis%n(maxshell,nset),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    prj_basis%nset = nset
    prj_basis%npgf = npgf
    prj_basis%nshell = nshell
    prj_basis%l = ll
    prj_basis%n = nn
    prj_basis%norm_type = 1
    prj_basis%kind_radius = 0._dp
    ALLOCATE (prj_basis%set_radius(nset),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    prj_basis%set_radius = 0._dp
    ALLOCATE (prj_basis%pgf_radius(maxpgf,nset),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    prj_basis%pgf_radius = 0._dp

    ALLOCATE (paw_basis%ncgf_set(nset),paw_basis%nsgf_set(nset),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (paw_basis%first_cgf(maxshell,nset),paw_basis%last_cgf(maxshell,nset),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (paw_basis%first_sgf(maxshell,nset),paw_basis%last_sgf(maxshell,nset),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ncgf = 0
    nsgf = 0
    maxco = 0
    DO iset=1,nset
       paw_basis%ncgf_set(iset) = 0
       paw_basis%nsgf_set(iset) = 0
       DO ishell=1,nshell(iset)
          lshell = paw_basis%l(ishell,iset)
          paw_basis%first_cgf(ishell,iset) = ncgf + 1
          ncgf = ncgf + nco(lshell)
          paw_basis%last_cgf(ishell,iset) = ncgf
          paw_basis%ncgf_set(iset) = paw_basis%ncgf_set(iset) + nco(lshell)
          paw_basis%first_sgf(ishell,iset) = nsgf + 1
          nsgf = nsgf + nso(lshell)
          paw_basis%last_sgf(ishell,iset) = nsgf
          paw_basis%nsgf_set(iset) = paw_basis%nsgf_set(iset) + nso(lshell)
       END DO
       maxco = MAX(maxco,npgf(iset)*ncoset(paw_basis%lmax(iset)))
    END DO
    paw_basis%ncgf = ncgf
    paw_basis%nsgf = nsgf
    ALLOCATE (paw_basis%norm_cgf(ncgf),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    paw_basis%norm_cgf = 0._dp
    ALLOCATE (paw_basis%zet(maxpgf,nset),paw_basis%cphi(maxco,ncgf),&
              paw_basis%sphi(maxco,nsgf),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (paw_basis%gcc(maxpgf,maxshell,nset),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE (paw_basis%cgf_symbol(ncgf),paw_basis%sgf_symbol(nsgf),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (paw_basis%lx(ncgf),paw_basis%ly(ncgf),paw_basis%lz(ncgf),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (paw_basis%m(nsgf),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ncgf = 0
    nsgf = 0
    DO iset=1,nset
       DO ishell=1,nshell(iset)
          lshell = paw_basis%l(ishell,iset)
          DO ico=ncoset(lshell-1)+1,ncoset(lshell)
             ncgf = ncgf + 1
             paw_basis%lx(ncgf) = indco(1,ico)
             paw_basis%ly(ncgf) = indco(2,ico)
             paw_basis%lz(ncgf) = indco(3,ico)
             paw_basis%cgf_symbol(ncgf) =&
                  cgf_symbol(nn(ishell,iset),(/paw_basis%lx(ncgf),&
                  paw_basis%ly(ncgf),paw_basis%lz(ncgf)/))
          END DO
          DO m=-lshell,lshell
             nsgf = nsgf + 1
             paw_basis%m(nsgf) = m
             paw_basis%sgf_symbol(nsgf) = sgf_symbol(nn(ishell,iset),lshell,m)
          END DO
       END DO
    END DO

    ALLOCATE (prj_basis%ncgf_set(nset),prj_basis%nsgf_set(nset),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (prj_basis%first_cgf(maxshell,nset),prj_basis%last_cgf(maxshell,nset),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (prj_basis%first_sgf(maxshell,nset),prj_basis%last_sgf(maxshell,nset),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ncgf = 0
    nsgf = 0
    maxco = 0
    DO iset=1,nset
       prj_basis%ncgf_set(iset) = 0
       prj_basis%nsgf_set(iset) = 0
       DO ishell=1,nshell(iset)
          lshell = prj_basis%l(ishell,iset)
          prj_basis%first_cgf(ishell,iset) = ncgf + 1
          ncgf = ncgf + nco(lshell)
          prj_basis%last_cgf(ishell,iset) = ncgf
          prj_basis%ncgf_set(iset) = prj_basis%ncgf_set(iset) + nco(lshell)
          prj_basis%first_sgf(ishell,iset) = nsgf + 1
          nsgf = nsgf + nso(lshell)
          prj_basis%last_sgf(ishell,iset) = nsgf
          prj_basis%nsgf_set(iset) = prj_basis%nsgf_set(iset) + nso(lshell)
       END DO
       maxco = MAX(maxco,npgf(iset)*ncoset(prj_basis%lmax(iset)))
    END DO
    prj_basis%ncgf = ncgf
    prj_basis%nsgf = nsgf
    ALLOCATE (prj_basis%norm_cgf(ncgf),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    prj_basis%norm_cgf = 0._dp
    ALLOCATE (prj_basis%zet(maxpgf,nset),prj_basis%cphi(maxco,ncgf),&
              prj_basis%sphi(maxco,nsgf),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (prj_basis%gcc(maxpgf,maxshell,nset),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE (prj_basis%cgf_symbol(ncgf),prj_basis%sgf_symbol(nsgf),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (prj_basis%lx(ncgf),prj_basis%ly(ncgf),prj_basis%lz(ncgf),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (prj_basis%m(nsgf),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ncgf = 0
    nsgf = 0
    DO iset=1,nset
       DO ishell=1,nshell(iset)
          lshell = prj_basis%l(ishell,iset)
          DO ico=ncoset(lshell-1)+1,ncoset(lshell)
             ncgf = ncgf + 1
             prj_basis%lx(ncgf) = indco(1,ico)
             prj_basis%ly(ncgf) = indco(2,ico)
             prj_basis%lz(ncgf) = indco(3,ico)
             prj_basis%cgf_symbol(ncgf) =&
                  cgf_symbol(nn(ishell,iset),(/prj_basis%lx(ncgf),&
                  prj_basis%ly(ncgf),prj_basis%lz(ncgf)/))
          END DO
          DO m=-lshell,lshell
             nsgf = nsgf + 1
             prj_basis%m(nsgf) = m
             prj_basis%sgf_symbol(nsgf) = sgf_symbol(nn(ishell,iset),lshell,m)
          END DO
       END DO
    END DO

    ! set up exponents of atom basis
    paw_basis%gcc = 0._dp
    ! unique set for all soft basis sets
    DO iset=1,paw_basis%nset
       DO ipgf=1,paw_basis%npgf(iset)
          l = iset+ipgf-2
          paw_basis%zet(ipgf,iset) = 0.1_dp*2.0_dp**l
          paw_basis%gcc(ipgf,ipgf,iset) = 1._dp
       END DO
    END DO
    CALL init_orb_basis_set(paw_basis,error)

    ! set up exponents of projector basis
    prj_basis%gcc = 0._dp
    DO iset=1,prj_basis%nset
       zetmin = minexp(prj_basis%lmax(iset),radius,eps_fit)
       DO ipgf=1,prj_basis%npgf(iset)
          prj_basis%zet(ipgf,iset) = zetmin*2.0_dp**(ipgf-1)
       END DO
    END DO

    ! calculate overlap matrix <basis|projector>
    ALLOCATE (smat(maxpgf,maxpgf),sinv(maxpgf,maxpgf),tmat(maxpgf,maxpgf),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    DO iset=1,paw_basis%nset
       na = paw_basis%npgf(iset)
       np = prj_basis%npgf(iset)
       CPPrecondition(na==np,cp_failure_level,routineP,error,failure)
       l = paw_basis%lmax(iset)
       prefac = 0.5_dp**(l + 2)*rootpi*dfac(2*l + 1)
       smat = 0._dp
       DO ipgf=1,na
          za = paw_basis%zet(ipgf,iset)
          DO jpgf=1,np
             zp = prj_basis%zet(jpgf,iset)
             smat(ipgf,jpgf) = prefac/(za + zp)**(REAL(l,dp)+1.5_dp)
          END DO
       END DO
       ! contraction of atom basis (most likely not the case)
       sinv = 0._dp
       sinv(1:na,1:np) = MATMUL(TRANSPOSE(paw_basis%gcc(1:na,1:na,iset)),smat(1:na,1:np))
       ! invert matrix
       CALL invert_matrix(sinv(1:na,1:np),tmat(1:na,1:np),myerr,"T",error=error)
       ! projector contraction
       prj_basis%gcc(1:np,1:np,iset) = tmat(1:na,1:np)
    END DO
    DEALLOCATE (smat,sinv,tmat,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    ! setup transformation matrices for projector basis
    CALL init_orb_basis_set(prj_basis,error)

    ! setup interaction radii
    eps_pgf = eps_fit/100._dp
    CALL init_interaction_radii_orb_basis(paw_basis,eps_pgf,error=error)
    CALL init_interaction_radii_orb_basis(prj_basis,eps_pgf,error=error)

    DEALLOCATE (npgf,nshell,nn,ll,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE  build_projector_basis

  FUNCTION minexp(lmax,radius,eps_fit) RESULT(zet)
    INTEGER                                  :: lmax
    REAL(dp)                                 :: radius, eps_fit, zet

    REAL(dp)                                 :: rad

    zet = 1._dp
    DO
       rad = exp_radius(lmax,zet,eps_fit,1.0_dp)
       IF(rad < radius) EXIT
       zet = zet + 0.5_dp
    END DO

  END FUNCTION minexp

! *****************************************************************************

  SUBROUTINE deallocate_lri_ints (lri_ints,error)

    TYPE(lri_int_list), POINTER              :: lri_ints
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_lri_ints', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iatom, ijkind, inode, istat, &
                                                natom, nkind, nnode
    LOGICAL                                  :: failure

    CPPrecondition(ASSOCIATED(lri_ints),cp_failure_level,routineP,error,failure)
    nkind = lri_ints%nkind

    IF(nkind > 0) THEN
       DO ijkind = 1,SIZE(lri_ints%lri_atom)
          natom = lri_ints%lri_atom(ijkind)%natom
          IF(natom > 0) THEN
             DO iatom=1,natom
                nnode = lri_ints%lri_atom(ijkind)%lri_node(iatom)%nnode
                IF(nnode > 0) THEN
                   IF(ASSOCIATED(lri_ints%lri_atom(ijkind)%lri_node(iatom)%lri_int)) THEN
                      DO inode = 1,nnode
                         DEALLOCATE (lri_ints%lri_atom(ijkind)%lri_node(iatom)%lri_int(inode)%abaint,&
                                     lri_ints%lri_atom(ijkind)%lri_node(iatom)%lri_int(inode)%abbint,&
                                     lri_ints%lri_atom(ijkind)%lri_node(iatom)%lri_int(inode)%soo,&
                                     lri_ints%lri_atom(ijkind)%lri_node(iatom)%lri_int(inode)%sab,&
                                     lri_ints%lri_atom(ijkind)%lri_node(iatom)%lri_int(inode)%qab,&
                                     lri_ints%lri_atom(ijkind)%lri_node(iatom)%lri_int(inode)%sinv,&
                                     lri_ints%lri_atom(ijkind)%lri_node(iatom)%lri_int(inode)%sn,&
                                     STAT=istat)
                         CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
                      END DO
                      DEALLOCATE (lri_ints%lri_atom(ijkind)%lri_node(iatom)%lri_int,STAT=istat)
                      CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
                   END IF
                   IF(ASSOCIATED(lri_ints%lri_atom(ijkind)%lri_node(iatom)%lri_rhoab)) THEN
                      DO inode = 1,nnode
                         DEALLOCATE (lri_ints%lri_atom(ijkind)%lri_node(iatom)%lri_rhoab(inode)%avec,&
                                     lri_ints%lri_atom(ijkind)%lri_node(iatom)%lri_rhoab(inode)%tvec,&
                                     STAT=istat)
                         CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
                      END DO
                      DEALLOCATE (lri_ints%lri_atom(ijkind)%lri_node(iatom)%lri_rhoab,STAT=istat)
                      CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
                   END IF
                   IF(ASSOCIATED(lri_ints%lri_atom(ijkind)%lri_node(iatom)%acoef)) THEN
                      DEALLOCATE (lri_ints%lri_atom(ijkind)%lri_node(iatom)%acoef,STAT=istat)
                      CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
                   END IF
                   IF(ASSOCIATED(lri_ints%lri_atom(ijkind)%lri_node(iatom)%fcoef)) THEN
                      DEALLOCATE (lri_ints%lri_atom(ijkind)%lri_node(iatom)%fcoef,STAT=istat)
                      CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
                   END IF
                   IF(ASSOCIATED(lri_ints%lri_atom(ijkind)%lri_node(iatom)%lri_paw)) THEN
                      IF(ASSOCIATED(lri_ints%lri_atom(ijkind)%lri_node(iatom)%lri_paw%fhh)) THEN
                        DEALLOCATE (lri_ints%lri_atom(ijkind)%lri_node(iatom)%lri_paw%fhh,STAT=istat)
                        CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
                      END IF
                      IF(ASSOCIATED(lri_ints%lri_atom(ijkind)%lri_node(iatom)%lri_paw%fhs)) THEN
                        DEALLOCATE (lri_ints%lri_atom(ijkind)%lri_node(iatom)%lri_paw%fhs,STAT=istat)
                        CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
                      END IF
                      IF(ASSOCIATED(lri_ints%lri_atom(ijkind)%lri_node(iatom)%lri_paw%fxh)) THEN
                        DEALLOCATE (lri_ints%lri_atom(ijkind)%lri_node(iatom)%lri_paw%fxh,STAT=istat)
                        CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
                      END IF
                      IF(ASSOCIATED(lri_ints%lri_atom(ijkind)%lri_node(iatom)%lri_paw%fxs)) THEN
                        DEALLOCATE (lri_ints%lri_atom(ijkind)%lri_node(iatom)%lri_paw%fxs,STAT=istat)
                        CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
                      END IF
                      IF(ASSOCIATED(lri_ints%lri_atom(ijkind)%lri_node(iatom)%lri_paw%flhh)) THEN
                        DEALLOCATE (lri_ints%lri_atom(ijkind)%lri_node(iatom)%lri_paw%flhh,STAT=istat)
                        CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
                      END IF
                      IF(ASSOCIATED(lri_ints%lri_atom(ijkind)%lri_node(iatom)%lri_paw%flhs)) THEN
                        DEALLOCATE (lri_ints%lri_atom(ijkind)%lri_node(iatom)%lri_paw%flhs,STAT=istat)
                        CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
                      END IF
                      IF(ASSOCIATED(lri_ints%lri_atom(ijkind)%lri_node(iatom)%lri_paw%flxh)) THEN
                        DEALLOCATE (lri_ints%lri_atom(ijkind)%lri_node(iatom)%lri_paw%flxh,STAT=istat)
                        CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
                      END IF
                      IF(ASSOCIATED(lri_ints%lri_atom(ijkind)%lri_node(iatom)%lri_paw%flxs)) THEN
                        DEALLOCATE (lri_ints%lri_atom(ijkind)%lri_node(iatom)%lri_paw%flxs,STAT=istat)
                        CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
                      END IF
                      DEALLOCATE (lri_ints%lri_atom(ijkind)%lri_node(iatom)%lri_paw,STAT=istat)
                      CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
                   END IF
                END IF
             END DO
             DEALLOCATE (lri_ints%lri_atom(ijkind)%lri_node,STAT=istat)
             CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
           END IF
       END DO
       DEALLOCATE (lri_ints%lri_atom,STAT=istat)
       CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    END IF
    DEALLOCATE (lri_ints,STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE deallocate_lri_ints

! *****************************************************************************
  SUBROUTINE lri_atom_init (lri_env,lri_control,atomic_kind_set,error)

    TYPE(lri_environment_type), POINTER      :: lri_env
    TYPE(lripaw_control_type), POINTER       :: lri_control
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'lri_atom_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikind, istat, la, ll, llmax, &
                                                lpaw, lri, na, nkind, nr
    LOGICAL                                  :: failure
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(grid_atom_type), POINTER            :: grid_atom
    TYPE(gto_basis_set_type), POINTER        :: paw_basis, ri_basis
    TYPE(harmonics_info_type), POINTER       :: harmonics

    failure = .FALSE.

    nkind = SIZE(lri_env%atom)

    CALL init_lebedev_grids()

    DO ikind = 1, nkind
       atomic_kind => atomic_kind_set(ikind)
       CALL get_atomic_kind(atomic_kind=atomic_kind,ngrid_rad=nr,ngrid_ang=na)
       ri_basis => lri_env%ri_basis(ikind)%gto_basis_set
       paw_basis => lri_env%paw_basis(ikind)%gto_basis_set
       CALL get_gto_basis_set(gto_basis_set=ri_basis,maxl=lri)
       CALL get_gto_basis_set(gto_basis_set=paw_basis,maxl=lpaw)
       llmax = MAX(lri,lpaw)
       ll = get_number_of_lebedev_grid(n=na)
       na = lebedev_grid(ll)%n
       la = lebedev_grid(ll)%l

       IF (llmax > la) THEN
          WRITE(6,'(/,72("*"))')
          WRITE(6,'(T2,A,T66,I4)') &
             "WARNING: the lebedev grid is built for angular momentum l up to ",la,&
             "         the max l of spherical harmonics is larger, l_max = ", llmax,&
             "         good integration is guaranteed only for l <= ",la
          WRITE(6,'(72("*"),/)')
       END IF

       CALL init_orbital_pointers(la)

       ! initialize radial and Lebedev grids
       NULLIFY(grid_atom)
       CALL allocate_grid_atom(grid_atom,error)
       grid_atom%ng_sphere = na
       grid_atom%nr = nr
       CALL create_grid_atom(grid_atom,nr,na,llmax,ll,lri_control%quadrature)
       lri_env%atom(ikind)%grid_atom => grid_atom

       ! initialize spherical harmonics
       ALLOCATE (harmonics,STAT=istat)
       CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
       CALL create_harmonics(harmonics,ll,error)
       lri_env%atom(ikind)%harmonics => harmonics
    END DO

    CALL deallocate_lebedev_grids()

  END SUBROUTINE lri_atom_init

! *****************************************************************************

  SUBROUTINE create_harmonics(harmonics,lleb,error)

    TYPE(harmonics_info_type), POINTER       :: harmonics
    INTEGER, INTENT(IN)                      :: lleb
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_harmonics', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: is, istat, l, la, m, ms, na
    LOGICAL                                  :: failure

    failure = .FALSE.

    na = lebedev_grid(lleb)%n
    la = lebedev_grid(lleb)%l
    ms = nsoset(la)

    harmonics%lmax   = la
    harmonics%maxs   = ms
    harmonics%ngrid  = na

    ALLOCATE (harmonics%coord(3,na),&
              harmonics%slm(na,ms),&
              harmonics%dslm(3,na,ms),STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)

    harmonics%coord(1:3,1:na) = lebedev_grid(lleb)%r(1:3,1:na)

    DO is=1,ms
       l = indso(1,is)
       m = indso(2,is)
       CALL y_lm(harmonics%coord,harmonics%slm(1:na,is),l,m)
       CALL drvy_lm(harmonics%coord,harmonics%dslm(1:3,1:na,is),l,m)
    END DO

  END SUBROUTINE create_harmonics

  SUBROUTINE deallocate_harmonics(harmonics,error)

    TYPE(harmonics_info_type), POINTER       :: harmonics
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_harmonics', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat
    LOGICAL                                  :: failure

    failure = .FALSE.

    IF(ASSOCIATED(harmonics%coord)) THEN
       DEALLOCATE (harmonics%coord,STAT=istat)
       CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    END IF
    IF(ASSOCIATED(harmonics%slm)) THEN
       DEALLOCATE (harmonics%slm,STAT=istat)
       CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    END IF
    IF(ASSOCIATED(harmonics%dslm)) THEN
       DEALLOCATE (harmonics%dslm,STAT=istat)
       CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    END IF

    DEALLOCATE (harmonics,STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE deallocate_harmonics

! *****************************************************************************
!> \brief Write all the GTO basis sets in use for the LRIPAW method to output
! *****************************************************************************
  SUBROUTINE write_lripaw_basis_sets(lri_env,atomic_kind_set,subsys_section,error)

    TYPE(lri_environment_type), POINTER      :: lri_env
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'write_lripaw_basis_sets', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikind, nkind, output_unit
    LOGICAL                                  :: failure
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_logger_type), POINTER            :: logger

    failure = .FALSE.
    NULLIFY (logger)
    logger => cp_error_get_logger(error)
    output_unit = cp_print_key_unit_nr(logger,subsys_section,&
         "PRINT%KINDS/BASIS_SET",&
         extension=".Log",error=error)
    IF (output_unit > 0) THEN
       IF (ASSOCIATED(lri_env)) THEN
          WRITE (UNIT=output_unit,FMT="(/,/,T2,A)")&
               "LRIPAW BASIS SET INFORMATION"
          nkind = SIZE(atomic_kind_set)
          DO ikind=1,nkind
             atomic_kind => atomic_kind_set(ikind)
             WRITE (UNIT=output_unit,FMT="(/,T2,I2,A)")&
                  atomic_kind%kind_number,". Atomic kind: "//TRIM(atomic_kind%name)
             CALL write_gto_basis_set(lri_env%ri_basis(ikind)%gto_basis_set,&
                  output_unit,"LRI Basis Set",error)
             CALL write_gto_basis_set(lri_env%ri_basis_soft(ikind)%gto_basis_set,&
                  output_unit,"LRI Soft Basis Set",error)
             CALL write_gto_basis_set(lri_env%paw_basis(ikind)%gto_basis_set,&
                  output_unit,"PAW Expansion Basis Set",error)
             CALL write_gto_basis_set(lri_env%proj_basis(ikind)%gto_basis_set,&
                  output_unit,"PAW Projector Basis Set",error)
          END DO
       ELSE
          CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
       END IF
    END IF

    CALL cp_print_key_finished_output(output_unit,logger,subsys_section,&
         "PRINT%KINDS/BASIS_SET",error=error)

  END SUBROUTINE write_lripaw_basis_sets
! *****************************************************************************

END MODULE lri_environment_types

