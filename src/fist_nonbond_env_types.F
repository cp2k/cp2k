!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      none
!> \author HAF
! *****************************************************************************
MODULE fist_nonbond_env_types
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE cell_types,                      ONLY: cell_release,&
                                             cell_type
  USE f77_blas
  USE fist_neighbor_list_types,        ONLY: fist_neighbor_deallocate,&
                                             fist_neighbor_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE pair_potential_types,            ONLY: pair_potential_pp_release,&
                                             pair_potential_pp_type,&
                                             siepmann_type,&
                                             tersoff_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'fist_nonbond_env_types'
  INTEGER, PRIVATE, SAVE :: last_fist_nonbond_env_id_nr=0
  PUBLIC :: fist_nonbond_env_type, fist_nonbond_env_set, &
            fist_nonbond_env_get, fist_nonbond_env_create, &
            fist_nonbond_env_retain, fist_nonbond_env_release, pos_type, eam_type, &
            quip_data_type

! *****************************************************************************
  TYPE pos_type
    REAL (KIND=dp) :: r(3)
  END TYPE

  TYPE eam_type
    REAL (KIND=dp) :: f_embed, rho
  END TYPE

  TYPE quip_data_type
    INTEGER, POINTER :: use_indices(:)
    REAL (KIND=dp), POINTER :: force(:,:)
    REAL (KIND=dp) :: virial(3,3)
  END TYPE

! *****************************************************************************
  TYPE fist_nonbond_env_type
    INTEGER :: ref_count, id_nr, natom_types
    INTEGER :: counter, last_update, num_update
    LOGICAL :: do_nonbonded, shift_cutoff
    CHARACTER ( len = default_string_length )   :: unit_type
    REAL (KIND=dp)                              :: lup, aup, ei_scale14, &
                                                   vdw_scale14
    REAL (KIND=dp)                              :: long_range_correction
    REAL (KIND=dp), DIMENSION ( :, : ), POINTER :: rlist_cut, rlist_lowsq, ij_kind_full_fac
    REAL (KIND=dp), DIMENSION ( : ), POINTER    :: charges
    TYPE(fist_neighbor_type), POINTER           :: nonbonded
    TYPE(pair_potential_pp_type), POINTER       :: potparm14, potparm
    TYPE(cell_type), POINTER                    :: cell_last_update
    TYPE(pos_type), DIMENSION(:), POINTER       :: r_last_update, r_last_update_pbc
    TYPE(pos_type), DIMENSION(:), POINTER       :: rshell_last_update_pbc
    TYPE(pos_type), DIMENSION(:), POINTER       :: rcore_last_update_pbc
    TYPE(eam_type), DIMENSION(:), POINTER       :: eam_data
    TYPE(quip_data_type), POINTER               :: quip_data
  END TYPE fist_nonbond_env_type

CONTAINS

! *****************************************************************************
!> \brief sets a fist_nonbond_env
!> \param fist_nonbond_env the object to create
!> \param potparm14 ...
!> \param potparm ...
!> \param nonbonded ...
!> \param rlist_cut ...
!> \param rlist_lowsq ...
!> \param aup ...
!> \param lup ...
!> \param ei_scale14 ...
!> \param vdw_scale14 ...
!> \param shift_cutoff ...
!> \param r_last_update ...
!> \param r_last_update_pbc ...
!> \param rshell_last_update_pbc ...
!> \param rcore_last_update_pbc ...
!> \param cell_last_update ...
!> \param num_update ...
!> \param last_update ...
!> \param counter ...
!> \param natom_types ...
!> \param long_range_correction ...
!> \param ij_kind_full_fac ...
!> \param eam_data ...
!> \param quip_data ...
!> \param charges ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      12.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE fist_nonbond_env_get(fist_nonbond_env, potparm14, potparm, &
       nonbonded, rlist_cut, rlist_lowsq, aup, lup, ei_scale14, vdw_scale14, &
       shift_cutoff, r_last_update, r_last_update_pbc, rshell_last_update_pbc, &
       rcore_last_update_pbc, cell_last_update, num_update, last_update, &
       counter, natom_types, long_range_correction, ij_kind_full_fac, eam_data, &
       quip_data, charges, error)

    TYPE(fist_nonbond_env_type), POINTER     :: fist_nonbond_env
    TYPE(pair_potential_pp_type), OPTIONAL, &
      POINTER                                :: potparm14, potparm
    TYPE(fist_neighbor_type), OPTIONAL, &
      POINTER                                :: nonbonded
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: rlist_cut, rlist_lowsq
    REAL(KIND=dp), OPTIONAL                  :: aup, lup, ei_scale14, &
                                                vdw_scale14
    LOGICAL, INTENT(OUT), OPTIONAL           :: shift_cutoff
    TYPE(pos_type), DIMENSION(:), OPTIONAL, &
      POINTER                                :: r_last_update, &
                                                r_last_update_pbc, &
                                                rshell_last_update_pbc, &
                                                rcore_last_update_pbc
    TYPE(cell_type), OPTIONAL, POINTER       :: cell_last_update
    INTEGER, OPTIONAL                        :: num_update, last_update, &
                                                counter, natom_types
    REAL(KIND=dp), OPTIONAL                  :: long_range_correction
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: ij_kind_full_fac
    TYPE(eam_type), DIMENSION(:), OPTIONAL, &
      POINTER                                :: eam_data
    TYPE(quip_data_type), OPTIONAL, POINTER  :: quip_data
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: charges
    TYPE(cp_error_type), INTENT(inout)       :: error

    IF ( PRESENT ( charges ) ) charges => fist_nonbond_env%charges
    IF ( PRESENT ( potparm14 ) ) potparm14 => fist_nonbond_env%potparm14
    IF ( PRESENT ( eam_data ) ) eam_data => fist_nonbond_env%eam_data
    IF ( PRESENT ( quip_data ) ) quip_data => fist_nonbond_env%quip_data
    IF ( PRESENT ( potparm ) ) potparm => fist_nonbond_env%potparm
    IF ( PRESENT ( rlist_cut ) ) rlist_cut => fist_nonbond_env%rlist_cut
    IF ( PRESENT ( rlist_lowsq ) ) rlist_lowsq => fist_nonbond_env%rlist_lowsq
    IF ( PRESENT ( ij_kind_full_fac ) ) ij_kind_full_fac => fist_nonbond_env%ij_kind_full_fac
    IF ( PRESENT ( nonbonded ) ) nonbonded => fist_nonbond_env%nonbonded
    IF ( PRESENT ( r_last_update ) ) &
         r_last_update => fist_nonbond_env%r_last_update
    IF ( PRESENT ( r_last_update_pbc ) ) &
         r_last_update_pbc => fist_nonbond_env%r_last_update_pbc
    IF ( PRESENT ( rshell_last_update_pbc ) ) &
         rshell_last_update_pbc => fist_nonbond_env%rshell_last_update_pbc
    IF ( PRESENT ( rcore_last_update_pbc ) ) &
         rcore_last_update_pbc => fist_nonbond_env%rcore_last_update_pbc
    IF ( PRESENT ( cell_last_update ) ) &
         cell_last_update => fist_nonbond_env%cell_last_update
    IF ( PRESENT ( lup ) ) lup = fist_nonbond_env%lup
    IF ( PRESENT ( aup ) ) aup = fist_nonbond_env%aup
    IF ( PRESENT ( ei_scale14 ) ) ei_scale14 = fist_nonbond_env%ei_scale14
    IF ( PRESENT ( vdw_scale14 ) ) vdw_scale14 = fist_nonbond_env%vdw_scale14
    IF ( PRESENT ( shift_cutoff ) ) &
       shift_cutoff = fist_nonbond_env%shift_cutoff
    IF ( PRESENT ( natom_types ) ) natom_types = fist_nonbond_env%natom_types
    IF ( PRESENT ( counter ) ) counter = fist_nonbond_env%counter
    IF ( PRESENT ( last_update ) ) last_update = fist_nonbond_env%last_update
    IF ( PRESENT ( num_update ) ) num_update = fist_nonbond_env%num_update
    IF ( PRESENT ( long_range_correction ) ) &
         long_range_correction = fist_nonbond_env%long_range_correction
  END SUBROUTINE fist_nonbond_env_get

! *****************************************************************************
!> \brief sets a fist_nonbond_env
!> \param fist_nonbond_env the object to create
!> \param potparm14 ...
!> \param potparm ...
!> \param rlist_cut ...
!> \param rlist_lowsq ...
!> \param nonbonded ...
!> \param aup ...
!> \param lup ...
!> \param ei_scale14 ...
!> \param vdw_scale14 ...
!> \param shift_cutoff ...
!> \param r_last_update ...
!> \param r_last_update_pbc ...
!> \param rshell_last_update_pbc ...
!> \param rcore_last_update_pbc ...
!> \param cell_last_update ...
!> \param num_update ...
!> \param last_update ...
!> \param counter ...
!> \param natom_types ...
!> \param long_range_correction ...
!> \param eam_data ...
!> \param quip_data ...
!> \param charges ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      12.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE fist_nonbond_env_set(fist_nonbond_env, potparm14, potparm, &
       rlist_cut, rlist_lowsq, nonbonded, aup, lup, ei_scale14, vdw_scale14, &
       shift_cutoff, r_last_update, r_last_update_pbc, rshell_last_update_pbc, &
       rcore_last_update_pbc, cell_last_update, num_update, last_update, &
       counter, natom_types, long_range_correction, eam_data, quip_data, charges, error)

    TYPE(fist_nonbond_env_type), POINTER     :: fist_nonbond_env
    TYPE(pair_potential_pp_type), OPTIONAL, &
      POINTER                                :: potparm14, potparm
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: rlist_cut, rlist_lowsq
    TYPE(fist_neighbor_type), OPTIONAL, &
      POINTER                                :: nonbonded
    REAL(KIND=dp), OPTIONAL                  :: aup, lup, ei_scale14, &
                                                vdw_scale14
    LOGICAL, INTENT(IN), OPTIONAL            :: shift_cutoff
    TYPE(pos_type), DIMENSION(:), OPTIONAL, &
      POINTER                                :: r_last_update, &
                                                r_last_update_pbc, &
                                                rshell_last_update_pbc, &
                                                rcore_last_update_pbc
    TYPE(cell_type), OPTIONAL, POINTER       :: cell_last_update
    INTEGER, OPTIONAL                        :: num_update, last_update, &
                                                counter, natom_types
    REAL(KIND=dp), OPTIONAL                  :: long_range_correction
    TYPE(eam_type), DIMENSION(:), OPTIONAL, &
      POINTER                                :: eam_data
    TYPE(quip_data_type), OPTIONAL, POINTER  :: quip_data
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: charges
    TYPE(cp_error_type), INTENT(inout)       :: error

    IF ( PRESENT ( potparm14 ) ) fist_nonbond_env%potparm14 => potparm14
    IF ( PRESENT ( eam_data ) ) fist_nonbond_env%eam_data => eam_data
    IF ( PRESENT ( quip_data ) ) fist_nonbond_env%quip_data => quip_data
    IF ( PRESENT ( potparm ) ) fist_nonbond_env%potparm => potparm
    IF ( PRESENT ( rlist_cut ) ) fist_nonbond_env%rlist_cut => rlist_cut
    IF ( PRESENT ( charges ) ) fist_nonbond_env%charges => charges
    IF ( PRESENT ( rlist_lowsq ) ) fist_nonbond_env%rlist_lowsq => rlist_lowsq
    IF ( PRESENT ( nonbonded ) ) fist_nonbond_env%nonbonded => nonbonded
    IF ( PRESENT ( r_last_update ) ) &
         fist_nonbond_env%r_last_update => r_last_update
    IF ( PRESENT ( r_last_update_pbc ) ) &
         fist_nonbond_env%r_last_update_pbc => r_last_update_pbc
    IF ( PRESENT ( rshell_last_update_pbc ) ) &
         fist_nonbond_env%rshell_last_update_pbc => rshell_last_update_pbc
    IF ( PRESENT ( rcore_last_update_pbc ) ) &
         fist_nonbond_env%rcore_last_update_pbc => rcore_last_update_pbc
    IF ( PRESENT ( cell_last_update ) ) &
         fist_nonbond_env%cell_last_update => cell_last_update
    IF ( PRESENT ( lup ) ) fist_nonbond_env%lup = lup
    IF ( PRESENT ( aup ) ) fist_nonbond_env%aup = aup
    IF ( PRESENT ( ei_scale14 ) ) fist_nonbond_env%ei_scale14 = ei_scale14
    IF ( PRESENT ( vdw_scale14 ) ) fist_nonbond_env%vdw_scale14 = vdw_scale14
    IF ( PRESENT ( shift_cutoff ) ) &
       fist_nonbond_env%shift_cutoff = shift_cutoff
    IF ( PRESENT ( natom_types ) ) fist_nonbond_env%natom_types = natom_types
    IF ( PRESENT ( counter ) ) fist_nonbond_env%counter = counter
    IF ( PRESENT ( last_update ) ) fist_nonbond_env%last_update = last_update
    IF ( PRESENT ( num_update ) ) fist_nonbond_env%num_update = num_update
    IF ( PRESENT ( long_range_correction ) ) &
         fist_nonbond_env%long_range_correction = long_range_correction
  END SUBROUTINE fist_nonbond_env_set

! *****************************************************************************
!> \brief allocates and intitializes a fist_nonbond_env
!> \param fist_nonbond_env the object to create
!> \param atomic_kind_set ...
!> \param potparm14 ...
!> \param potparm ...
!> \param do_nonbonded ...
!> \param verlet_skin ...
!> \param ewald_rcut ...
!> \param ei_scale14 ...
!> \param vdw_scale14 ...
!> \param shift_cutoff ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      12.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE fist_nonbond_env_create(fist_nonbond_env, atomic_kind_set, &
       potparm14, potparm, do_nonbonded, verlet_skin, ewald_rcut, ei_scale14, &
       vdw_scale14, shift_cutoff, error)
    TYPE(fist_nonbond_env_type), POINTER     :: fist_nonbond_env
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(pair_potential_pp_type), OPTIONAL, &
      POINTER                                :: potparm14, potparm
    LOGICAL, INTENT(IN)                      :: do_nonbonded
    REAL(KIND=dp), INTENT(IN)                :: verlet_skin, ewald_rcut, &
                                                ei_scale14, vdw_scale14
    LOGICAL, INTENT(IN)                      :: shift_cutoff
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fist_nonbond_env_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.
    ALLOCATE(fist_nonbond_env, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    NULLIFY(fist_nonbond_env%potparm14)
    NULLIFY(fist_nonbond_env%potparm)
    NULLIFY(fist_nonbond_env%rlist_cut)
    NULLIFY(fist_nonbond_env%rlist_lowsq)
    NULLIFY(fist_nonbond_env%ij_kind_full_fac)
    NULLIFY(fist_nonbond_env%nonbonded)
    NULLIFY(fist_nonbond_env%cell_last_update)
    NULLIFY(fist_nonbond_env%r_last_update)
    NULLIFY(fist_nonbond_env%r_last_update_pbc)
    NULLIFY(fist_nonbond_env%rshell_last_update_pbc)
    NULLIFY(fist_nonbond_env%rcore_last_update_pbc)
    NULLIFY(fist_nonbond_env%eam_data)
    NULLIFY(fist_nonbond_env%quip_data)
    NULLIFY(fist_nonbond_env%charges)
    CALL init_fist_nonbond_env(fist_nonbond_env, atomic_kind_set, potparm14, &
         potparm, do_nonbonded, verlet_skin, ewald_rcut, ei_scale14, vdw_scale14, &
         shift_cutoff, error)
  END SUBROUTINE fist_nonbond_env_create

! *****************************************************************************
!> \brief Purpose: Initialise the FIST nonbond environment.
!> \param fist_nonbond_env the object to create
!> \param atomic_kind_set ...
!> \param potparm14 ...
!> \param potparm ...
!> \param do_nonbonded ...
!> \param verlet_skin ...
!> \param ewald_rcut ...
!> \param ei_scale14 ...
!> \param vdw_scale14 ...
!> \param shift_cutoff ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE init_fist_nonbond_env(fist_nonbond_env, atomic_kind_set, &
       potparm14, potparm, do_nonbonded, verlet_skin, ewald_rcut, ei_scale14, &
       vdw_scale14, shift_cutoff, error)

    TYPE(fist_nonbond_env_type), POINTER     :: fist_nonbond_env
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(pair_potential_pp_type), OPTIONAL, &
      POINTER                                :: potparm14, potparm
    LOGICAL, INTENT(IN)                      :: do_nonbonded
    REAL(KIND=dp), INTENT(IN)                :: verlet_skin, ewald_rcut, &
                                                ei_scale14, vdw_scale14
    LOGICAL, INTENT(IN)                      :: shift_cutoff
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'init_fist_nonbond_env', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: idim, jdim, natom_types, stat
    LOGICAL                                  :: check, failure, use_potparm, &
                                                use_potparm14
    REAL(KIND=dp)                            :: fac, rcut, rlow

    use_potparm14 = PRESENT(potparm14)
    IF (use_potparm14) use_potparm14 = use_potparm14.OR.ASSOCIATED(potparm14)
    use_potparm = PRESENT(potparm)
    IF (use_potparm) use_potparm = use_potparm.OR.ASSOCIATED(potparm)
    NULLIFY (fist_nonbond_env%nonbonded)
    NULLIFY (fist_nonbond_env%r_last_update)
    NULLIFY (fist_nonbond_env%r_last_update_pbc)
    NULLIFY (fist_nonbond_env%rshell_last_update_pbc)
    NULLIFY (fist_nonbond_env%rcore_last_update_pbc)
    NULLIFY (fist_nonbond_env%cell_last_update)
    NULLIFY (fist_nonbond_env%rlist_cut)
    NULLIFY (fist_nonbond_env%rlist_lowsq)
    NULLIFY (fist_nonbond_env%ij_kind_full_fac)
    fist_nonbond_env%unit_type = "ANGSTROM"
    fist_nonbond_env%do_nonbonded = do_nonbonded
    fist_nonbond_env%lup = 0
    fist_nonbond_env%aup = 0
    fist_nonbond_env%ei_scale14 = ei_scale14
    fist_nonbond_env%vdw_scale14 = vdw_scale14
    fist_nonbond_env%shift_cutoff = shift_cutoff
    fist_nonbond_env%counter = 0
    fist_nonbond_env%last_update = 0
    fist_nonbond_env%num_update = 0
    fist_nonbond_env%long_range_correction = 0
    fist_nonbond_env%ref_count=1
    fist_nonbond_env%id_nr=last_fist_nonbond_env_id_nr
    last_fist_nonbond_env_id_nr=last_fist_nonbond_env_id_nr+1
    IF (do_nonbonded) THEN
       natom_types = 1
       ! Determine size of kind arrays
       natom_types = SIZE(atomic_kind_set)
       IF (use_potparm14) THEN
          check = (SIZE(potparm14%pot,1)==natom_types)
          CPPostcondition(check,cp_failure_level,routineP,error,failure)
       END IF
       IF (use_potparm) THEN
          check = (SIZE(potparm%pot,1)==natom_types)
          CPPostcondition(check,cp_failure_level,routineP,error,failure)
       END IF
       ALLOCATE ( fist_nonbond_env%rlist_cut(natom_types, natom_types), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE ( fist_nonbond_env%rlist_lowsq(natom_types, natom_types), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE ( fist_nonbond_env%ij_kind_full_fac(natom_types, natom_types), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       fist_nonbond_env%ij_kind_full_fac = 1.0_dp
       DO idim = 1, natom_types
          DO jdim = idim, natom_types
             IF ((use_potparm) .OR. (use_potparm14)) THEN
                IF (use_potparm) THEN
                    rcut = SQRT(potparm%pot(idim,jdim)%pot%rcutsq)
                    fac = potparm%pot(idim,jdim)%pot%spl_f%rscale(1)
                    rlow = fac / (potparm%pot(idim,jdim)%pot%pair_spline_data(1)%spline_data%xn)
                ELSE
                    rcut = SQRT(potparm14%pot(idim,jdim)%pot%rcutsq)
                    fac = potparm14%pot(idim,jdim)%pot%spl_f%rscale(1)
                    rlow = fac / (potparm14%pot(idim,jdim)%pot%pair_spline_data(1)%spline_data%xn)
                END IF
                ! Warning: rlist_rcut should only be used by the neighbor list
                ! algorithm. It is not the cutoff for the evaluation of the
                ! interactions because rlist_rcut includes the Verlet skin.
                rcut = MAX(rcut, ewald_rcut) + verlet_skin
                fist_nonbond_env%rlist_cut(idim,jdim) = rcut
                fist_nonbond_env%rlist_cut(jdim,idim) = rcut
                rlow  = rlow*(1.06_dp)**2 ! 1.06_dp in order to have 1/2 Emax_spline
                fist_nonbond_env%rlist_lowsq(idim,jdim) = rlow
                fist_nonbond_env%rlist_lowsq(jdim,idim) = rlow
                ! In case of manybody potential the neighbor list will be full.
                ! This means that for each atom pair (a,b) of the current types,
                ! atom a is in the neighbor list of b and b is in the neighbor
                ! list of a. ij_kind_full_fac is used to correct for the double
                ! counting in the conventional pair potentials cause by this
                ! situation.
                IF (ANY(potparm%pot(idim,jdim)%pot%type==tersoff_type)) THEN
                   ! TODO: what if 14 is not of tersoff type while the normal
                   ! nonbond is? (or the reverse). We'd better impose
                   ! consistency.
                   fist_nonbond_env%ij_kind_full_fac(idim,jdim) = 0.5_dp
                   fist_nonbond_env%ij_kind_full_fac(idim,jdim) = 0.5_dp
                END IF
                IF (ANY(potparm%pot(idim,jdim)%pot%type==siepmann_type)) THEN
                   ! TODO:see tersoff_type
                   fist_nonbond_env%ij_kind_full_fac(idim,jdim) = 0.5_dp
                   fist_nonbond_env%ij_kind_full_fac(idim,jdim) = 0.5_dp
                END IF
             ELSE
                ! In case we don't use potparm for initialization let's account
                ! only for the real-space part of the Ewald sum.
                fist_nonbond_env%rlist_cut(idim,jdim) = ewald_rcut
                fist_nonbond_env%rlist_cut(jdim,idim) = ewald_rcut
                fist_nonbond_env%rlist_lowsq(idim,jdim) = 0.0_dp
                fist_nonbond_env%rlist_lowsq(jdim,idim) = 0.0_dp
             END IF
          END DO
       END DO
       IF (use_potparm14) fist_nonbond_env%potparm14 => potparm14
       IF (use_potparm) fist_nonbond_env%potparm => potparm
       fist_nonbond_env%natom_types = natom_types
    ELSE
       NULLIFY(fist_nonbond_env%potparm)
       NULLIFY(fist_nonbond_env%potparm14)
    END IF
  END SUBROUTINE init_fist_nonbond_env

! *****************************************************************************
!> \brief retains the given fist_nonbond_env (see doc/ReferenceCounting.html)
!> \param fist_nonbond_env the object to retain
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      12.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE fist_nonbond_env_retain(fist_nonbond_env,error)
    TYPE(fist_nonbond_env_type), POINTER     :: fist_nonbond_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fist_nonbond_env_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(fist_nonbond_env),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CPPrecondition(fist_nonbond_env%ref_count>0,cp_failure_level,routineP,error,failure)
       fist_nonbond_env%ref_count=fist_nonbond_env%ref_count+1
    END IF
  END SUBROUTINE fist_nonbond_env_retain

! *****************************************************************************
!> \brief releases the given fist_nonbond_env (see doc/ReferenceCounting.html)
!> \param fist_nonbond_env the object to release
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      12.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE fist_nonbond_env_release(fist_nonbond_env,error)
    TYPE(fist_nonbond_env_type), POINTER     :: fist_nonbond_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fist_nonbond_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.
    IF (ASSOCIATED(fist_nonbond_env)) THEN
       CPPrecondition(fist_nonbond_env%ref_count>0,cp_failure_level,routineP,error,failure)
       fist_nonbond_env%ref_count=fist_nonbond_env%ref_count-1
       IF (fist_nonbond_env%ref_count<1) THEN
          IF (ASSOCIATED(fist_nonbond_env%nonbonded)) THEN
             CALL fist_neighbor_deallocate (fist_nonbond_env%nonbonded, error )
          END IF
          ! Release potparm
          CALL pair_potential_pp_release ( fist_nonbond_env%potparm, error )
          ! Release potparm14
          CALL pair_potential_pp_release ( fist_nonbond_env%potparm14, error )
          IF (ASSOCIATED(fist_nonbond_env%r_last_update)) THEN
             DEALLOCATE (fist_nonbond_env%r_last_update,STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ENDIF
          IF (ASSOCIATED(fist_nonbond_env%r_last_update_pbc)) THEN
             DEALLOCATE (fist_nonbond_env%r_last_update_pbc,STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ENDIF
          IF (ASSOCIATED(fist_nonbond_env%charges)) THEN
             DEALLOCATE (fist_nonbond_env%charges,STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ENDIF
          IF (ASSOCIATED(fist_nonbond_env%eam_data)) THEN
             DEALLOCATE (fist_nonbond_env%eam_data,STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ENDIF
          IF (ASSOCIATED(fist_nonbond_env%quip_data)) THEN
             IF (ASSOCIATED(fist_nonbond_env%quip_data%force)) THEN
                DEALLOCATE (fist_nonbond_env%quip_data%force,STAT=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             ENDIF
             IF (ASSOCIATED(fist_nonbond_env%quip_data%use_indices)) THEN
                DEALLOCATE (fist_nonbond_env%quip_data%use_indices,STAT=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             ENDIF
             DEALLOCATE (fist_nonbond_env%quip_data,STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ENDIF
          IF (ASSOCIATED(fist_nonbond_env%rshell_last_update_pbc)) THEN
             DEALLOCATE (fist_nonbond_env%rshell_last_update_pbc,STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ENDIF
          IF (ASSOCIATED(fist_nonbond_env%rcore_last_update_pbc)) THEN
             DEALLOCATE (fist_nonbond_env%rcore_last_update_pbc,STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ENDIF
          IF (ASSOCIATED(fist_nonbond_env%cell_last_update)) THEN
             CALL cell_release(fist_nonbond_env%cell_last_update,error)
          ENDIF
          IF (ASSOCIATED(fist_nonbond_env%ij_kind_full_fac)) THEN
             DEALLOCATE (fist_nonbond_env%ij_kind_full_fac,STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ENDIF
          IF (ASSOCIATED(fist_nonbond_env%rlist_cut)) THEN
             DEALLOCATE (fist_nonbond_env%rlist_cut,STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ENDIF
          IF (ASSOCIATED(fist_nonbond_env%rlist_lowsq)) THEN
             DEALLOCATE (fist_nonbond_env%rlist_lowsq,STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ENDIF
          DEALLOCATE(fist_nonbond_env, STAT=stat)
          CPPrecondition(stat==0,cp_warning_level,routineP,error,failure)
       END IF
    END IF
  END SUBROUTINE fist_nonbond_env_release

END MODULE fist_nonbond_env_types
