!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2010  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      September 2005 - Introduced the Born-Mayer-Huggins-Fumi-Tosi  Potential (BMHTF)
!>      2006 - Major rewriting of the routines.. Linear scaling setup of splines
!>      2007 - Teodoro Laino - University of Zurich - Multiple potential 
!>             Major rewriting nr.2
!> \author CJM
! *****************************************************************************
MODULE pair_potential

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE ewald_environment_types,         ONLY: ewald_env_get,&
                                             ewald_environment_type
  USE f77_blas
  USE fparser,                         ONLY: finalizef,&
                                             initf,&
                                             parsef
  USE kinds,                           ONLY: default_path_length,&
                                             dp
  USE memory_utilities,                ONLY: reallocate
  USE pair_potential_types,            ONLY: &
       b4_type, bm_type, compare_pot, ea_type, ft_type, gp_type, gw_type, &
       ip_type, list_pot, lj_charmm_type, lj_type, multi_type, nn_type, &
       pair_potential_pp_type, pair_potential_single_type, &
       potential_single_allocation, tersoff_type, wl_type
  USE pair_potential_util,             ONLY: ener_coul,&
                                             ener_pot,&
                                             zbl_matching_polinomial
  USE shell_potential_types,           ONLY: get_shell,&
                                             shell_kind_type
  USE splines_methods,                 ONLY: init_spline,&
                                             init_splinexy,&
                                             potential_s
  USE splines_types,                   ONLY: spline_data_p_type,&
                                             spline_data_type,&
                                             spline_env_create,&
                                             spline_environment_type,&
                                             spline_factor_create,&
                                             spline_factor_release,&
                                             spline_factor_type
  USE string_table,                    ONLY: str2id
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: sort
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pair_potential'
  REAL(KIND=dp), PARAMETER, PRIVATE    :: MIN_HICUT_VALUE=1.0E-15_dp,&
                                          DEFAULT_HICUT_VALUE=1.0E3_dp
  INTEGER, PARAMETER, PRIVATE          :: MAX_POINTS=2000000

  PUBLIC :: spline_nonbond_control,&
            get_nonbond_storage

CONTAINS

! *****************************************************************************
!> \brief Initialize genpot
!> \par History
!>      Teo 2007.06 - Zurich University
! *****************************************************************************
  SUBROUTINE init_genpot(potparm, ntype)
    TYPE(pair_potential_pp_type), POINTER    :: potparm
    INTEGER, INTENT(IN)                      :: ntype

    INTEGER                                  :: i, j, k, ngp
    TYPE(pair_potential_single_type), &
      POINTER                                :: pot

    NULLIFY(pot)
    ngp    = 0
    ! Prescreen for general potential type
    DO i = 1, ntype ! i:  first  atom type
       DO j = 1, i  ! j:  second atom type
          pot => potparm%pot(i,j)%pot
          ngp = ngp + COUNT(pot%type==gp_type)
       END DO
    END DO
    CALL initf(ngp)
    ngp = 0
    DO i = 1, ntype ! i:  first  atom type
       DO j = 1, i  ! j:  second atom type
          pot => potparm%pot(i,j)%pot
          DO k = 1, SIZE(pot%type)
             IF (pot%type(k)==gp_type) THEN 
                ngp = ngp + 1
                pot%set(k)%gp%myid = ngp
                CALL parsef(ngp,TRIM(pot%set(k)%gp%potential),pot%set(k)%gp%parameters)
             END IF
          END DO
       END DO
    END DO
  END SUBROUTINE init_genpot

! *****************************************************************************
!> \brief creates the splines for the potentials
!> \par History
!>      Teo 2006.05 : Improved speed and accuracy. Linear scaling of the setup
! *****************************************************************************
  SUBROUTINE spline_nonbond_control ( spline_env, potparm, atomic_kind_set, &
       ewald_env, ei_scale14, vdw_scale14, eps_spline, max_energy, rlow_nb,&
       emax_spline, iw, iw2, iw3, do_14, do_zbl, scale_cutoff, spline_c_coul, &
       spline_cs_coul, spline_sc_coul, error)

    TYPE(spline_environment_type), POINTER   :: spline_env
    TYPE(pair_potential_pp_type), POINTER    :: potparm
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(ewald_environment_type), OPTIONAL, &
      POINTER                                :: ewald_env
    REAL(KIND=dp), OPTIONAL                  :: ei_scale14, vdw_scale14
    REAL(KIND=dp), INTENT(IN)                :: eps_spline, max_energy, &
                                                rlow_nb, emax_spline
    INTEGER, INTENT(IN)                      :: iw, iw2, iw3
    LOGICAL, INTENT(IN), OPTIONAL            :: do_14, do_zbl
    LOGICAL, INTENT(IN)                      :: scale_cutoff
    TYPE(spline_environment_type), &
      OPTIONAL, POINTER                      :: spline_c_coul, &
                                                spline_cs_coul, spline_sc_coul
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'spline_nonbond_control', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ewald_type, handle, i, j, k, &
                                                n, ncount, npoints, ntype
    LOGICAL :: do_i_shell_model, do_j_shell_model, do_multipoles, &
      do_shell_model, failure, found_locut, my_do_14, my_do_zbl
    REAL(KIND=dp) :: alpha, c_cutoff_e, cs_cutoff_e, cutoff_e, cutoff_n, &
      e_cutoff_coul, ec_cc, ec_cs, ec_sc, energy_cutoff, fac, hicut, hicut0, &
      locut, my_ei_scale14, my_vdw_scale14, sc_cutoff_e
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(pair_potential_single_type), &
      POINTER                                :: pot

    failure= .FALSE.
    n      = 0
    ncount = 0
    my_ei_scale14 = 1.0_dp ; IF (PRESENT(ei_scale14))   my_ei_scale14 = ei_scale14
    my_vdw_scale14= 1.0_dp ; IF (PRESENT(vdw_scale14)) my_vdw_scale14 = vdw_scale14
    my_do_14= .FALSE.      ; IF (PRESENT(do_14))             my_do_14 = do_14
    my_do_zbl= .FALSE.     ; IF (PRESENT(do_zbl))           my_do_zbl = do_zbl
    ntype = SIZE ( atomic_kind_set )
    CALL timeset(routineN,handle)
    IF (iw3>0) THEN
       WRITE(iw3,'(T2,A,I0,A,I0,A)') 'SPLINE_INFO| Generating: ',&
            (ntype*(ntype+1))/2,' splines due to ',ntype,' different atomic kinds'
    ENDIF
    CALL init_genpot(potparm, ntype)
    ! Real computation of splines
    CALL ewald_env_get(ewald_env,ewald_type=ewald_type,alpha=alpha,&
                       do_multipoles=do_multipoles,error=error)
    IF (scale_cutoff)  fac = my_ei_scale14
    IF (do_multipoles) fac = 0.0_dp
    DO i = 1, ntype
       ! First Atom
       atomic_kind=> atomic_kind_set(i)
       CALL get_atom_info(atomic_kind, my_do_14, do_shell_model=do_i_shell_model, error=error)
       DO j = 1, i
          ! Second Atom
          pot => potparm%pot(i,j)%pot
          IF (iw3>0 .AND. iw<=0) THEN
             IF (MOD(i*(i-1)/2+j,MAX(1,(ntype*(ntype+1))/(2*10)))==0) WRITE(iw3,'(T2,A3,I0)') '...',i*(i-1)/2+j
          ENDIF
          atomic_kind=> atomic_kind_set(j)
          CALL get_atom_info(atomic_kind, my_do_14, do_shell_model=do_j_shell_model, error=error)
          do_shell_model   = do_i_shell_model .OR. do_j_shell_model

          ! Setup of Exclusion Types
          pot%only_mb = .TRUE.
          pot%only_pp = .TRUE.
          DO k = 1, SIZE(pot%type)
             SELECT CASE(pot%type(k))
             CASE (lj_type,lj_charmm_type,wl_type,gw_type,ft_type,ip_type,b4_type,bm_type,gp_type,ea_type)
                pot%only_mb = .FALSE.
             CASE (tersoff_type)
                pot%only_pp = .FALSE.
             CASE (nn_type)
                ! Do nothing..
             CASE DEFAULT
                ! Never reach this point
                CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
             END SELECT
             ! Special case for EAM
             SELECT CASE(pot%type(k))
             CASE (ea_type)
                pot%only_pp = .FALSE.
             END SELECT
          END DO

          ! Starting SetUp of splines
          IF (.NOT.pot%undef) CYCLE
          ncount = ncount + 1
          n = spline_env % spltab(i,j)
          locut = rlow_nb
          hicut0 = SQRT (pot%rcutsq)
          IF (ABS(hicut0) <= MIN_HICUT_VALUE) hicut0 = DEFAULT_HICUT_VALUE
          hicut = hicut0 / SQRT (pot%spl_f%rcutsq_f)

          energy_cutoff = pot%spl_f%cutoff_n
          e_cutoff_coul = pot%spl_f%cutoff_e
          IF (do_shell_model) THEN
             ec_cc = pot%spl_c%cutoff_e
             ec_sc = pot%spl_sc%cutoff_e
             ec_cs = pot%spl_cs%cutoff_e
          END IF

          ! Find the real locut according emax_spline
          CALL get_spline_cutoff(hicut, locut, found_locut, do_shell_model, pot, my_vdw_scale14,&
          my_ei_scale14, my_do_14, ewald_type, do_multipoles, my_do_zbl, alpha, e_cutoff_coul, &
               ec_cc, ec_sc, ec_cs, energy_cutoff, emax_spline)
          locut=MAX(locut*SQRT(pot%spl_f%rcutsq_f),rlow_nb)

          ! Real Generation of the Spline
          CALL generate_spline_low(spline_env%spl_pp(n)%spl_p, npoints, locut, hicut, my_vdw_scale14,&
               my_ei_scale14, my_do_14, ewald_type, do_multipoles, alpha, pot%e_fac, scale_cutoff,&
               eps_spline, iw, iw2, i, j, n, ncount, max_energy, e_cutoff_coul, pot, energy_cutoff,&
               label='', found_locut=found_locut, do_nonbonded=.TRUE.,hicut0=hicut0, &
               do_zbl=my_do_zbl, error=error)

          ! This is Specific for Shell-Model (Only Coulomb terms in the splines)
          IF(do_shell_model) THEN
             CPPostcondition(PRESENT(spline_c_coul),cp_failure_level,routineP,error,failure)
             CPPostcondition(PRESENT(spline_cs_coul),cp_failure_level,routineP,error,failure)
             CPPostcondition(PRESENT(spline_sc_coul),cp_failure_level,routineP,error,failure) 

             ! Real Generation of the Spline Core-Core
             CALL generate_spline_low(spline_c_coul%spl_pp(n)%spl_p, npoints, locut, hicut,&
                  0.0_dp, 1.0_dp, .FALSE., ewald_type, do_multipoles, alpha, pot%e_fcc, &
                  .FALSE., eps_spline, iw, iw2, i, j, n, ncount, max_energy, ec_cc,&
                  pot, label='Core-Core Coulomb term in Shell-Model',&
                  found_locut=found_locut,do_nonbonded=.FALSE.,hicut0=hicut0,do_zbl=my_do_zbl,&
                  error=error)

             ! Both are shell
             IF (do_i_shell_model.AND.do_j_shell_model) THEN
                ! Real Generation of the Spline Core-Shell
                CALL generate_spline_low(spline_cs_coul%spl_pp(n)%spl_p, npoints, locut, hicut,&
                     0.0_dp, 1.0_dp, .FALSE., ewald_type, do_multipoles, alpha, pot%e_fcs, &
                     .FALSE., eps_spline, iw, iw2, i, j, n, ncount, max_energy, ec_cs,&
                     pot,label='Core-Shell Coulomb term in Shell-Model',&
                     found_locut=found_locut,do_nonbonded=.FALSE.,hicut0=hicut0,do_zbl=my_do_zbl,&
                     error=error)

                ! Real Generation of the Spline Shell-Core
                CALL generate_spline_low(spline_sc_coul%spl_pp(n)%spl_p, npoints, locut, hicut,&
                     0.0_dp, 1.0_dp, .FALSE., ewald_type, do_multipoles, alpha, pot%e_fsc, &
                     .FALSE., eps_spline, iw, iw2, i, j, n, ncount, max_energy, ec_sc,&
                     pot,label='Shell-Core Coulomb term in Shell-Model',&
                     found_locut=found_locut,do_nonbonded=.FALSE.,hicut0=hicut0,do_zbl=my_do_zbl,&
                     error=error)
             END IF
          END IF
          pot%undef = .FALSE.
          ! Unique Spline working only for a pure LJ potential..
          IF (SIZE(pot%type)==1) THEN
             IF (ANY(potential_single_allocation==pot%type(1))) THEN
                ! Restoring the proper values for the generating spline pot
                IF ((pot%type(1)==lj_type).OR.(pot%type(1)==lj_charmm_type)) THEN
                   pot%set(1)%lj%sigma6  = pot%set(1)%lj%sigma6*pot%spl_f%rscale(2)**3
                   pot%set(1)%lj%sigma12 = pot%set(1)%lj%sigma6**2
                   pot%set(1)%lj%epsilon = pot%set(1)%lj%epsilon*pot%spl_f%fscale(2)
                END IF
             END IF
          END IF
          ! Correct Cutoff...
          IF (scale_cutoff) THEN
             cutoff_n = ener_pot(pot, hicut0, 0.0_dp, my_vdw_scale14)
             cutoff_e = ener_coul(pot, 0.0_dp, fac*pot%e_fac*pot%spl_f%fscale(1), hicut0, ewald_type, alpha, &
                        my_do_14, zbl = my_do_zbl)
             pot%spl_f%cutoff_n = pot%spl_f%cutoff_n*pot%spl_f%fscale(2) - cutoff_n
             pot%spl_f%cutoff_e = pot%spl_f%cutoff_e*pot%spl_f%fscale(1) - cutoff_e
             pot%spl_f%cutoff   = pot%spl_f%cutoff_n + pot%spl_f%cutoff_e
             IF(do_shell_model) THEN
                c_cutoff_e  = ener_coul(pot,0.0_dp, fac*pot%e_fcc*pot%spl_c%fscale(1), hicut0, ewald_type, alpha, &
                    my_do_14, zbl=my_do_zbl )
                sc_cutoff_e = ener_coul(pot,0.0_dp, fac*pot%e_fsc*pot%spl_sc%fscale(1),hicut0, ewald_type, alpha, &
                    my_do_14, zbl=my_do_zbl )
                cs_cutoff_e = ener_coul(pot,0.0_dp, fac*pot%e_fcs*pot%spl_cs%fscale(1),hicut0, ewald_type, alpha, &
                    my_do_14, zbl=my_do_zbl )
                pot%spl_c%cutoff_n  = 0.0_dp
                pot%spl_sc%cutoff_n = 0.0_dp
                pot%spl_cs%cutoff_n = 0.0_dp
                pot%spl_c%cutoff_e  = pot%spl_c%cutoff_e  * pot%spl_c%fscale(1)  - c_cutoff_e
                pot%spl_sc%cutoff_e = pot%spl_sc%cutoff_e * pot%spl_sc%fscale(1) - sc_cutoff_e
                pot%spl_cs%cutoff_e = pot%spl_cs%cutoff_e * pot%spl_cs%fscale(1) - cs_cutoff_e
                pot%spl_c%cutoff    = pot%spl_c%cutoff_n  + pot%spl_c%cutoff_e
                pot%spl_cs%cutoff   = pot%spl_cs%cutoff_n + pot%spl_cs%cutoff_e
                pot%spl_sc%cutoff   = pot%spl_sc%cutoff_n + pot%spl_sc%cutoff_e
             END IF
          END IF
       END DO
    END DO
    CALL finalizef()
    IF (iw3>0) THEN
       WRITE(iw3,'()')
    ENDIF
    IF (iw>0) THEN
       WRITE(iw,'("SPLINE_INFO|",5X,A,I7)')"Number of pair potential Splines allocated:",&
            MAXVAL(spline_env % spltab)
    END IF
    IF (iw3>0) THEN
       WRITE(iw3,'(T2,A,I7)') 'SPLINE_INFO| Number of unique splines computed:',&
            MAXVAL(spline_env % spltab)
       WRITE(iw3,'(T2,A,I0)') 'SPLINE_INFO| Done'
    ENDIF
    CALL timestop(handle)
  END SUBROUTINE spline_nonbond_control

! *****************************************************************************
!> \brief Gather atom info for the evaluation of the splines
!> \par History
!>      Splitting in order to make some season cleaning..
!> \author Teodoro Laino [tlaino] 2007.06
! *****************************************************************************
  SUBROUTINE get_atom_info(atomic_kind, my_do_14, qi, qci, qsi, z, rc, do_shell_model, error)
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    LOGICAL, INTENT(IN)                      :: my_do_14
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: qi, qci, qsi, z, rc
    LOGICAL, INTENT(OUT), OPTIONAL           :: do_shell_model
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_atom_info', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: number
    LOGICAL                                  :: is_i_shell, my_do_shell_model
    REAL(KIND=dp)                            :: qlci, qli, qlsi, rcov
    TYPE(shell_kind_type), POINTER           :: shell_i

    NULLIFY(shell_i)
    CALL get_atomic_kind(atomic_kind=atomic_kind,qeff=qli,shell_active=is_i_shell,shell=shell_i,z=number, rcov=rcov)
    my_do_shell_model = is_i_shell .AND.(.NOT. my_do_14)
    IF(my_do_shell_model) THEN
       CALL get_shell(shell=shell_i, charge_core=qlci, charge_shell=qlsi, error=error) 
    ELSE
       qlsi = qli
       qlci = qli
    END IF
    IF (PRESENT(qi))  qi  = qli
    IF (PRESENT(qci)) qci = qlci
    IF (PRESENT(qsi)) qsi = qlsi
    IF (PRESENT(z)) z = REAL(number)
    IF (PRESENT(rc)) rc = rcov
    IF (PRESENT(do_shell_model)) do_shell_model = my_do_shell_model
  END SUBROUTINE get_atom_info

! *****************************************************************************
!> \brief Finds the cutoff for the generation of the spline
!>      In a two pass approach, first with low resolution, refine in a second iteration
!> \par History
!>      Splitting in order to make some season cleaning..
!> \author Teodoro Laino [tlaino] 2007.06
! *****************************************************************************
  SUBROUTINE  get_spline_cutoff(hicut, locut, found_locut, do_shell_model, pot,&
  my_vdw_scale14, my_ei_scale14, my_do_14, ewald_type, do_multipoles, my_do_zbl, alpha,&
       e_cutoff_coul, ec_cc, ec_sc, ec_cs, energy_cutoff, emax_spline)

    REAL(KIND=dp), INTENT(IN)                :: hicut
    REAL(KIND=dp), INTENT(INOUT)             :: locut
    LOGICAL, INTENT(OUT)                     :: found_locut
    LOGICAL, INTENT(IN)                      :: do_shell_model
    TYPE(pair_potential_single_type), &
      OPTIONAL, POINTER                      :: pot
    REAL(KIND=dp), INTENT(IN)                :: my_vdw_scale14, my_ei_scale14
    LOGICAL, INTENT(IN)                      :: my_do_14
    INTEGER, INTENT(IN)                      :: ewald_type
    LOGICAL, INTENT(IN)                      :: do_multipoles, my_do_zbl
    REAL(KIND=dp), INTENT(IN)                :: alpha, e_cutoff_coul, ec_cc, &
                                                ec_sc, ec_cs, energy_cutoff, &
                                                emax_spline

    CHARACTER(len=*), PARAMETER :: routineN = 'get_spline_cutoff', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ilevel, jx
    REAL(KIND=dp)                            :: dx2, e, f, fac, locut_found, x

    dx2 = (hicut-locut)
    x   = hicut
    locut_found = locut
    found_locut = .FALSE.
    fac = 1.0_dp
    IF (do_multipoles) fac = 0.0_dp
    IF (do_shell_model) THEN
       DO ilevel=1,2
          dx2 = dx2 / 100.0_dp
          DO jx = 1, 100
             e = ener_pot( pot, x , energy_cutoff, my_vdw_scale14 )
             f = fac * pot%e_fac
             e = e + ener_coul ( pot, e_cutoff_coul, f, x, ewald_type, alpha, .FALSE., zbl=my_do_zbl )
             f = fac * pot%e_fsc
             e = e + ener_coul ( pot,        ec_sc, f, x, ewald_type, alpha, .FALSE., zbl=my_do_zbl )
             f = fac * pot%e_fcs
             e = e + ener_coul ( pot,        ec_cs, f, x, ewald_type, alpha, .FALSE., zbl=my_do_zbl )
             f = fac * pot%e_fcc
             e = e + ener_coul ( pot,        ec_cc, f, x, ewald_type, alpha, .FALSE., zbl=my_do_zbl )  
             IF (ABS(e) > emax_spline) THEN
                locut_found = x
                found_locut = .TRUE.
                EXIT
             END IF
             x = x - dx2
          END DO
          x= x + dx2
       ENDDO
    ELSE
       DO ilevel=1,2
          dx2 = dx2 / 100.0_dp
          DO jx = 1, 100
             e = ener_pot(pot, x , energy_cutoff, my_vdw_scale14 )
             f =  fac * my_ei_scale14 * pot%e_fac
             e = e + ener_coul (pot,e_cutoff_coul, f, x, ewald_type, alpha, my_do_14, zbl=my_do_zbl )
             IF (ABS(e) > emax_spline) THEN
                locut_found = x
                found_locut = .TRUE.
                EXIT
             END IF
             x = x - dx2
          END DO
          x= x + dx2
       ENDDO
    END IF
    locut = locut_found

  END SUBROUTINE get_spline_cutoff

! *****************************************************************************
!> \brief Real Generation of spline.. 
!> \par History
!>      Splitting in order to make some season cleaning..
!> \author Teodoro Laino [tlaino] 2007.06
! *****************************************************************************
  SUBROUTINE generate_spline_low(spl_p, npoints, locut, hicut, my_vdw_scale14,&
       my_ei_scale14, my_do_14, ewald_type, do_multipoles, alpha, e_fac,&
       scale_cutoff, eps_spline, iw, iw2, i, j, n, ncount, max_energy, e_cutoff_coul,&
       pot, energy_cutoff, label, found_locut, do_nonbonded, hicut0, do_zbl, error)

    TYPE(spline_data_p_type), DIMENSION(:), &
      POINTER                                :: spl_p
    INTEGER, INTENT(INOUT)                   :: npoints
    REAL(KIND=dp), INTENT(IN)                :: locut, hicut, my_vdw_scale14, &
                                                my_ei_scale14
    LOGICAL, INTENT(IN)                      :: my_do_14
    INTEGER, INTENT(IN)                      :: ewald_type
    LOGICAL, INTENT(IN)                      :: do_multipoles
    REAL(KIND=dp), INTENT(IN)                :: alpha, e_fac
    LOGICAL, INTENT(IN)                      :: scale_cutoff
    REAL(KIND=dp), INTENT(IN)                :: eps_spline
    INTEGER, INTENT(IN)                      :: iw, iw2, i, j, n, ncount
    REAL(KIND=dp), INTENT(IN)                :: max_energy, e_cutoff_coul
    TYPE(pair_potential_single_type), &
      POINTER                                :: pot
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: energy_cutoff
    CHARACTER(LEN=*), INTENT(IN)             :: label
    LOGICAL, INTENT(IN)                      :: found_locut, do_nonbonded
    REAL(KIND=dp), INTENT(IN)                :: hicut0
    LOGICAL, INTENT(IN)                      :: do_zbl
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'generate_spline_low', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: message
    INTEGER                                  :: jx, mfac
    LOGICAL                                  :: check, failure
    REAL(KIND=dp)                            :: diffmax, dx2, e, f, myfac1, &
                                                myfac2, x, x2, xdum, xdum1, &
                                                xsav
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(spline_data_type), POINTER          :: spline_data
    TYPE(spline_factor_type), POINTER        :: spl_f

    failure = .FALSE.
    NULLIFY(logger, spl_f)
    logger => cp_error_get_logger(error)
    check = PRESENT(energy_cutoff)
    CALL spline_factor_create(spl_f, error)
    IF (do_nonbonded) THEN
       CPPostcondition(check,cp_failure_level,routineP,error,failure)
    END IF
    mfac    =  5
    npoints = 20
    spline_data => spl_p(1)%spline_data
    IF (.NOT. found_locut) npoints = 2
    f = my_ei_scale14 * e_fac
    IF (do_multipoles) f = 0.0_dp
    DO WHILE (.TRUE.)
       CALL init_splinexy ( spline_data, npoints+1)
       dx2 = (1.0_dp/locut**2-1.0_dp/hicut**2)/REAL(npoints,KIND=dp)
       x2  = 1.0_dp/hicut**2
       spline_data % x1 = x2
       DO jx = 1, npoints + 1
          ! jx: loop over 1/distance**2
          x = SQRT ( 1.0_dp/x2 )
          e = ener_coul(pot,e_cutoff_coul, f, x, ewald_type, alpha, my_do_14, do_zbl)
          IF (do_nonbonded) THEN
             e =  e + ener_pot(pot, x , energy_cutoff, my_vdw_scale14 )
          END IF
          spline_data % y(jx) = e
          x2 = x2 + dx2
       END DO
       IF (my_do_14) THEN
          myfac1 = my_ei_scale14 * e_fac
          myfac2 = 0.0_dp
          IF (scale_cutoff) THEN
             myfac2 = my_ei_scale14 * e_fac / hicut0
          END IF
          CALL init_spline(spline_data,dx=dx2,&
                           ei_scale14=myfac1,&
                           ei_cutoff14=myfac2,&
                           dual=.FALSE.)
       ELSE
          CALL init_spline(spline_data,dx=dx2,dual=.FALSE.)
       END IF
       ! This is the check for required accuracy on spline setup
       dx2 = (hicut - locut)/REAL(mfac*npoints+1,KIND=dp)
       x2  = locut + dx2
       diffmax = 0.0_dp
       xsav = hicut
       DO jx = 1, mfac*npoints
          x = x2
          e = ener_coul(pot,e_cutoff_coul, f, x, ewald_type, alpha, my_do_14, do_zbl )
          IF (do_nonbonded) THEN
             e =  e + ener_pot(pot, x , energy_cutoff, my_vdw_scale14 )
          END IF
          IF (ABS(e)<max_energy) THEN
             xdum1   = ABS(e-potential_s(spl_p,x*x,xdum,spl_f,logger))
             diffmax = MAX(diffmax,xdum1)
             xsav = MIN(x,xsav)
          END IF
          x2 = x2 + dx2
          IF (x2>hicut) EXIT
       END DO
       IF (npoints > MAX_POINTS) THEN
          WRITE(message,'(A,I8,A,G12.6,A)') "SPLINE_INFO| Number of points: ",npoints,&
               " obtained accuracy ",diffmax,". MM SPLINE: no convergence on required"//&
               " accuracy (adjust EPS_SPLINE and rerun)"
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed, routineP,&
               TRIM(message)//" "//&
CPSourceFileRef,&
               only_ionode=.TRUE.)
       END IF
       IF (diffmax > eps_spline ) THEN
          npoints   = CEILING(1.2_dp*REAL(npoints,KIND=dp))
       ELSE
          EXIT
       END IF
    END DO
    ! Print Spline info
    IF (iw>0) THEN
       WRITE(iw,*)
       WRITE(iw,'("SPLINE_INFO|",5X,A,I7,A)')"Spline number:",i*(i-1)/2+j,TRIM(label)
       WRITE(iw,'("SPLINE_INFO|",5X,2(A,I7,2X),A,2I7)')"Spline number",ncount,&
            "Unique number",n," Kinds involved:",i,j
       WRITE(iw,'("SPLINE_INFO|",5X,A,I7)')"Number of points used in the Splines:",npoints
       WRITE(iw,'("SPLINE_INFO|",5X,A,F18.12)')"Accuracy requested:",eps_spline
       WRITE(iw,'("SPLINE_INFO|",5X,A,F18.12)')"Accuracy achieved: ",diffmax
       WRITE(iw,'("SPLINE_INFO|",5X,A,2F12.6)')"RMIN and RMAX  used for spline (bohr):               ",&
            locut,hicut
       WRITE(iw,'("SPLINE_INFO|",5X,A,2F12.6)')"RMIN and RMAX used to achieve spline accuracy (bohr):",&
            xsav,hicut
       dx2 = (hicut - locut)/REAL(npoints+1,KIND=dp)
       x  = locut + dx2
       WRITE(iw,'("SPLINE_INFO|",5X,A,F15.9)') "Spline value at RMIN (Hartree):",&
            potential_s(spl_p,x*x,xdum,spl_f,logger)
       WRITE(iw,'("SPLINE_INFO|",5X,A,F15.9)') "Spline value at RMAX (Hartree):",&
            potential_s(spl_p,hicut*hicut,xdum,spl_f,logger)
       WRITE(iw,'("SPLINE_INFO|",5X,A,F15.9)') "Electrostatic cutoff:          ",e_cutoff_coul
       IF (do_nonbonded) THEN
          WRITE(iw,'("SPLINE_INFO|",5X,A,F15.9)') "Non-bonded cutoff:             ",energy_cutoff
       END IF
    END IF
    ! Print spline data on file if requested
    IF (iw2>0) THEN
       dx2 = (hicut - locut)/REAL(5000+1,KIND=dp)
       x2  = locut + dx2
       WRITE(1000+n,*)"KINDS: ",i,j,"  SPLINE NUMBER:",ncount
       DO jx = 1, 5000
          x = x2
          e =  ener_coul(pot,e_cutoff_coul, f, x, ewald_type, alpha, my_do_14, do_zbl )
          IF (do_nonbonded) THEN
             e =  e + ener_pot(pot, x , energy_cutoff, my_vdw_scale14 )
          END IF
          WRITE(1000+n,*) x, e, potential_s(spl_p,x*x,xdum,spl_f,logger)
          x2 = x2 + dx2
       END DO
       WRITE(1000+n,*)
    END IF
    ! In case of potential allowing a cumulative allocation...
    IF (ANY(potential_single_allocation==pot%type(1))) THEN
       ! Setup spline for coulomb potential
       spline_data => spl_p(1)%spline_data
       CALL init_splinexy (spline_data, npoints+1)
       dx2 = (1.0_dp/locut**2-1.0_dp/hicut**2)/REAL(npoints,KIND=dp)
       x2  = 1.0_dp/hicut**2
       spline_data % x1 = x2
       DO jx = 1, npoints + 1
          ! jx: loop over 1/distance**2
          x = SQRT ( 1.0_dp/x2 )
          e = ener_coul(pot,e_cutoff_coul, f, x, ewald_type, alpha, my_do_14, do_zbl )
          spline_data % y(jx) = e
          x2 = x2 + dx2
       END DO
       IF (my_do_14) THEN
          myfac1 = my_ei_scale14 * e_fac
          myfac2 = 0.0_dp
          IF (scale_cutoff) THEN
             myfac2 = my_ei_scale14 * e_fac / hicut0
          END IF
          CALL init_spline(spline_data,dx=dx2,&
                           ei_scale14=myfac1,&
                           ei_cutoff14=myfac2,&
                           dual=.TRUE.)
       ELSE
          CALL init_spline(spline_data,dx=dx2,dual=.TRUE.)
       END IF
       ! Setup spline for the non-bonded potential
       spline_data => spl_p(2)%spline_data
       IF (do_nonbonded) THEN
          CALL init_splinexy (spline_data, npoints+1)
          dx2 = (1.0_dp/locut**2-1.0_dp/hicut**2)/REAL(npoints,KIND=dp)
          x2  = 1.0_dp/hicut**2
          spline_data % x1 = x2
          DO jx = 1, npoints + 1
             ! jx: loop over 1/distance**2
             x = SQRT ( 1.0_dp/x2 )
             e = ener_pot(pot, x , energy_cutoff, my_vdw_scale14 )
             spline_data % y(jx) = e
             x2 = x2 + dx2
          END DO
          IF (my_do_14) THEN
             myfac1 = my_ei_scale14 * e_fac
             myfac2 = 0.0_dp
             IF (scale_cutoff) THEN
                myfac2 = my_ei_scale14 * e_fac / hicut0
             END IF
             CALL init_spline(spline_data,dx=dx2,&
                              ei_scale14=myfac1,&
                              ei_cutoff14=myfac2,&
                              dual=.TRUE.)
          ELSE
             CALL init_spline(spline_data,dx=dx2,dual=.TRUE.)
          END IF
       END IF
    END IF
    CALL spline_factor_release(spl_f, error)
  END SUBROUTINE generate_spline_low

! *****************************************************************************
!> \brief Prescreening of the effective bonds evaluations. linear scaling algorithm
!> \author Teodoro Laino [tlaino] 2006.05
! *****************************************************************************
  SUBROUTINE get_nonbond_storage ( spline_env, potparm, atomic_kind_set, &
       spline_c_coul_env, spline_cs_coul_env, spline_sc_coul_env,&
       ei_scale14, vdw_scale14, do_14, do_zbl, scale_cutoff, ewald_env, error)

    TYPE(spline_environment_type), POINTER   :: spline_env
    TYPE(pair_potential_pp_type), POINTER    :: potparm
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(spline_environment_type), &
      OPTIONAL, POINTER                      :: spline_c_coul_env, &
                                                spline_cs_coul_env, &
                                                spline_sc_coul_env
    REAL(KIND=dp), OPTIONAL                  :: ei_scale14, vdw_scale14
    LOGICAL, INTENT(IN), OPTIONAL            :: do_14, do_zbl
    LOGICAL, INTENT(IN)                      :: scale_cutoff
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_nonbond_storage', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, i, idim, iend, istart, j, k, locij, n, ndim, nk, &
      ntype, nunique, nvar, pot_target, stat, tmpij(2), tmpij0(2)
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: Iwork1, Iwork2, my_index
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: tmp_index
    INTEGER, DIMENSION(:), POINTER           :: array_double_allocation
    LOGICAL :: at_least_one, check, do_i_shell_model, do_j_shell_model, &
      failure, my_do_14, my_do_zbl, shell_present
    REAL(KIND=dp)                            :: my_ei_scale14, &
                                                my_vdw_scale14, qci, qcj, qi, &
                                                qj, qsi, qsj
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: Cwork, Rwork, wtmp
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: pot_par
    TYPE(atomic_kind_type), POINTER          :: atomic_kind

    failure     = .FALSE.
    CALL timeset(routineN,handle)
    ntype = SIZE ( atomic_kind_set )
    my_ei_scale14 = 1.0_dp ; IF (PRESENT(ei_scale14))  my_ei_scale14  = ei_scale14
    my_vdw_scale14= 1.0_dp ; IF (PRESENT(vdw_scale14)) my_vdw_scale14 = vdw_scale14
    my_do_14= .FALSE.      ; IF (PRESENT(do_14))             my_do_14 = do_14
    my_do_zbl= .FALSE.     ; IF (PRESENT(do_zbl))           my_do_zbl = do_zbl
    DO i = 1, ntype
       DO j = 1, i
          potparm%pot(i,j)%pot%undef = .FALSE.
       END DO
    END DO
    ALLOCATE(tmp_index(ntype,ntype),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(array_double_allocation(0),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    !
    nunique = 0
    tmp_index = HUGE(0)
    DO pot_target = MINVAL(list_pot), MAXVAL(list_pot)
       ndim  = 0
       DO i = 1, ntype
          DO j = 1, i
             IF (SIZE(potparm%pot(i,j)%pot%type)/=1) CYCLE
             IF (potparm%pot(i,j)%pot%type(1)==pot_target) THEN
                tmp_index(i,j) = 1
                tmp_index(j,i) = 1
                ndim = ndim + 1
             END IF
          END DO
       END DO
       IF (ndim==0) CYCLE ! No potential of this kind found
       nvar = 0
       IF (my_ei_scale14 /=0.0_dp) nvar = nvar + 6
       SELECT CASE(pot_target)
       CASE (lj_type,lj_charmm_type)
          IF (my_vdw_scale14/=0.0_dp) nvar = 3 + nvar
       CASE (wl_type)
          nvar = 3 + nvar
       CASE (gw_type)
          nvar = 5 + nvar
       CASE (ea_type)
          nvar = 4 + nvar
       CASE (ft_type)
          nvar = 4 + nvar
       CASE (ip_type)
          nvar = 3 + nvar
       CASE (b4_type)
          nvar = 6 + nvar
       CASE (bm_type)
          nvar = 9 + nvar
       CASE (gp_type)
          nvar = 2 + nvar
       CASE (tersoff_type)
          nvar = 13 + nvar
       CASE (nn_type)
          nvar = nvar
       CASE DEFAULT
          CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
       END SELECT
       ! Setup a table of the indexes..
       ALLOCATE(my_index(ndim),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       n  = 0
       nk = 0
       DO i = 1, ntype
          DO j = 1, i
             n = n + 1
             IF (SIZE(potparm%pot(i,j)%pot%type)/=1) CYCLE
             IF (potparm%pot(i,j)%pot%type(1)==pot_target) THEN
                nk = nk + 1
                my_index(nk) = n
             END IF
          END DO
       END DO
       IF (nvar /=0) THEN
          ALLOCATE(pot_par(ndim,nvar),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          n  = 0
          nk = 0
          DO i = 1, ntype
             atomic_kind=> atomic_kind_set(i)
             CALL get_atom_info(atomic_kind, my_do_14, qi, qci, qsi, do_shell_model=do_i_shell_model, error=error)
             DO j = 1, i
                atomic_kind=> atomic_kind_set(j)
                CALL get_atom_info(atomic_kind, my_do_14, qj, qcj, qsj, do_shell_model=do_j_shell_model, error=error)
                n = n + 1
                IF (SIZE(potparm%pot(i,j)%pot%type)/=1) CYCLE
                IF (potparm%pot(i,j)%pot%type(1)==pot_target) THEN
                   nk = nk + 1
                   my_index(nk) = n
                   SELECT CASE(pot_target)
                   CASE (lj_type,lj_charmm_type)
                      pot_par(nk,1)  = potparm%pot(i,j)%pot%set(1)%lj%epsilon
                      pot_par(nk,2)  = potparm%pot(i,j)%pot%set(1)%lj%sigma6
                      pot_par(nk,3)  = potparm%pot(i,j)%pot%set(1)%lj%sigma12
                   CASE(gp_type)
                      pot_par(nk,1)  = str2id(potparm%pot(i,j)%pot%set(1)%gp%potential)
                      pot_par(nk,2)  = str2id(potparm%pot(i,j)%pot%set(1)%gp%variables)
                   CASE (wl_type)
                      pot_par(nk,1)  = potparm%pot(i,j)%pot%set(1)%willis%a
                      pot_par(nk,2)  = potparm%pot(i,j)%pot%set(1)%willis%b
                      pot_par(nk,3)  = potparm%pot(i,j)%pot%set(1)%willis%c
                   CASE (gw_type)
                      pot_par(nk,1)  = potparm%pot(i,j)%pot%set(1)%goodwin%vr0
                      pot_par(nk,2)  = potparm%pot(i,j)%pot%set(1)%goodwin%m
                      pot_par(nk,3)  = potparm%pot(i,j)%pot%set(1)%goodwin%mc
                      pot_par(nk,4)  = potparm%pot(i,j)%pot%set(1)%goodwin%d
                      pot_par(nk,5)  = potparm%pot(i,j)%pot%set(1)%goodwin%dc
                   CASE (ea_type)
                      pot_par(nk,1)  = potparm%pot(i,j)%pot%set(1)%eam%drar
                      pot_par(nk,2)  = potparm%pot(i,j)%pot%set(1)%eam%drhoar
                      pot_par(nk,3)  = potparm%pot(i,j)%pot%set(1)%eam%acutal
                      pot_par(nk,4)  = potparm%pot(i,j)%pot%set(1)%eam%npoints
                   CASE (ft_type)
                      pot_par(nk,1)  = potparm%pot(i,j)%pot%set(1)%ft%A
                      pot_par(nk,2)  = potparm%pot(i,j)%pot%set(1)%ft%B
                      pot_par(nk,3)  = potparm%pot(i,j)%pot%set(1)%ft%C
                      pot_par(nk,4)  = potparm%pot(i,j)%pot%set(1)%ft%D
                   CASE (ip_type)
                      pot_par(nk,1)  = potparm%pot(i,j)%pot%set(1)%ipbv%rcore
                      pot_par(nk,2)  = potparm%pot(i,j)%pot%set(1)%ipbv%m
                      pot_par(nk,3)  = potparm%pot(i,j)%pot%set(1)%ipbv%b
                   CASE (b4_type)
                      pot_par(nk,1)  = potparm%pot(i,j)%pot%set(1)%buck4r%a
                      pot_par(nk,2)  = potparm%pot(i,j)%pot%set(1)%buck4r%b
                      pot_par(nk,3)  = potparm%pot(i,j)%pot%set(1)%buck4r%c
                      pot_par(nk,4)  = potparm%pot(i,j)%pot%set(1)%buck4r%r1
                      pot_par(nk,5)  = potparm%pot(i,j)%pot%set(1)%buck4r%r2
                      pot_par(nk,6)  = potparm%pot(i,j)%pot%set(1)%buck4r%r3
                   CASE (bm_type)
                      pot_par(nk,1)  = potparm%pot(i,j)%pot%set(1)%buckmo%f0
                      pot_par(nk,2)  = potparm%pot(i,j)%pot%set(1)%buckmo%a1
                      pot_par(nk,3)  = potparm%pot(i,j)%pot%set(1)%buckmo%a2
                      pot_par(nk,4)  = potparm%pot(i,j)%pot%set(1)%buckmo%b1
                      pot_par(nk,5)  = potparm%pot(i,j)%pot%set(1)%buckmo%b2
                      pot_par(nk,6)  = potparm%pot(i,j)%pot%set(1)%buckmo%c
                      pot_par(nk,7)  = potparm%pot(i,j)%pot%set(1)%buckmo%d
                      pot_par(nk,8)  = potparm%pot(i,j)%pot%set(1)%buckmo%r0
                      pot_par(nk,9)  = potparm%pot(i,j)%pot%set(1)%buckmo%beta
                   CASE (tersoff_type)
                      pot_par(nk,1)  = potparm%pot(i,j)%pot%set(1)%tersoff%A
                      pot_par(nk,2)  = potparm%pot(i,j)%pot%set(1)%tersoff%B
                      pot_par(nk,3)  = potparm%pot(i,j)%pot%set(1)%tersoff%lambda1
                      pot_par(nk,4)  = potparm%pot(i,j)%pot%set(1)%tersoff%lambda2
                      pot_par(nk,5)  = potparm%pot(i,j)%pot%set(1)%tersoff%alpha
                      pot_par(nk,6)  = potparm%pot(i,j)%pot%set(1)%tersoff%beta
                      pot_par(nk,7)  = potparm%pot(i,j)%pot%set(1)%tersoff%n
                      pot_par(nk,8)  = potparm%pot(i,j)%pot%set(1)%tersoff%c
                      pot_par(nk,9)  = potparm%pot(i,j)%pot%set(1)%tersoff%d
                      pot_par(nk,10) = potparm%pot(i,j)%pot%set(1)%tersoff%h
                      pot_par(nk,11) = potparm%pot(i,j)%pot%set(1)%tersoff%lambda3
                      pot_par(nk,12) = potparm%pot(i,j)%pot%set(1)%tersoff%bigR
                      pot_par(nk,13) = potparm%pot(i,j)%pot%set(1)%tersoff%bigD
                   CASE (nn_type)
                      ! no checks
                   CASE DEFAULT
                      CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
                   END SELECT
                   IF (my_ei_scale14 /=0.0_dp) THEN
                      pot_par(nk,nvar-5)=qci
                      pot_par(nk,nvar-4)=qcj
                      pot_par(nk,nvar-3)=qsi
                      pot_par(nk,nvar-2)=qsj
                      pot_par(nk,nvar-1)=qi
                      pot_par(nk,nvar  )=qj
                   END IF
                   IF (ANY(potential_single_allocation==pot_target))  pot_par(nk,:) = REAL(pot_target, KIND=dp)
                END IF
             END DO
          END DO
          ! Main Sorting Loop
          ALLOCATE(Rwork(ndim),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(Iwork1(ndim),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(Iwork2(ndim),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(wtmp(nvar),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          CALL sort(pot_par(:,1),ndim,Iwork1)
          ! Sort all the other components of the potential
          DO k = 2, nvar
             Rwork = pot_par(:,k)
             DO i = 1, ndim
                pot_par(i,k)=Rwork(Iwork1(i))
             END DO
          END DO
          Iwork2  = my_index
          DO i = 1, ndim
             my_index(i) = Iwork2(Iwork1(i))
          END DO
          ! Iterative sorting
          DO k = 2, nvar
             wtmp(1:k-1) = pot_par(1,1:k-1)
             istart   = 1
             at_least_one = .FALSE.
             DO j = 1, ndim
                Rwork(j)  = pot_par(j,k)
                IF (ALL(pot_par(j,1:k-1)==wtmp(1:k-1))) CYCLE
                iend   = j - 1
                wtmp(1:k-1) = pot_par(j,1:k-1)
                ! If the ordered array has no two same consecutive elements
                ! does not make any sense to proceed ordering the others
                ! related parameters..
                idim = iend-istart+1
                CALL sort(Rwork(istart:iend),idim,Iwork1(istart:iend))
                Iwork1(istart:iend)=Iwork1(istart:iend)-1+istart
                IF (idim/=1) at_least_one = .TRUE.
                istart = j
             END DO
             iend = ndim
             idim = iend-istart+1
             CALL sort(Rwork(istart:iend),idim,Iwork1(istart:iend))
             Iwork1(istart:iend)=Iwork1(istart:iend)-1+istart
             IF (idim/=1) at_least_one = .TRUE.
             pot_par(:,k) = Rwork
             IF (.NOT.at_least_one) EXIT
             ! Sort other components
             DO j = k+1, nvar
                Rwork = pot_par(:,j)
                DO i = 1, ndim
                   pot_par(i,j)=Rwork(Iwork1(i))
                END DO
             END DO
             Iwork2  = my_index
             DO i = 1, ndim
                my_index(i) = Iwork2(Iwork1(i))
             END DO
          END DO
          DEALLOCATE(wtmp,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(Iwork1,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(Iwork2,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(Rwork,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          !
          ! Let's determine the number of unique potentials and tag them
          !
          ALLOCATE(Cwork(nvar),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          Cwork = pot_par(1,:)
          locij = my_index(1)
          CALL get_indexes(locij,ntype,tmpij0)
          istart = 1
          DO j = 1, ndim
             ! Special cases for EAM and IPBV
             locij = my_index(j)
             CALL get_indexes(locij,ntype,tmpij)
             SELECT CASE(pot_target)
             CASE(ea_type,ip_type)
                ! check the array components
                CALL compare_pot(potparm%pot(tmpij(1) ,tmpij(2) )%pot,&
                                 potparm%pot(tmpij0(1),tmpij0(2))%pot,&
                                 check, error)
             CASE(gp_type)
                check = .TRUE.
                IF  (ASSOCIATED(potparm%pot(tmpij(1) ,tmpij(2) )%pot%set(1)%gp%parameters).AND.&
                     ASSOCIATED(potparm%pot(tmpij0(1),tmpij0(2))%pot%set(1)%gp%parameters)) THEN
                   IF  (SIZE(potparm%pot(tmpij(1) ,tmpij(2) )%pot%set(1)%gp%parameters)==&
                        SIZE(potparm%pot(tmpij0(1),tmpij0(2))%pot%set(1)%gp%parameters)) THEN
                      IF  (ANY(potparm%pot(tmpij(1) ,tmpij(2) )%pot%set(1)%gp%parameters /= &
                               potparm%pot(tmpij0(1),tmpij0(2))%pot%set(1)%gp%parameters)) check = .FALSE.
                   END IF
                END IF
                IF  (ASSOCIATED(potparm%pot(tmpij(1) ,tmpij(2) )%pot%set(1)%gp%values).AND.&
                     ASSOCIATED(potparm%pot(tmpij0(1),tmpij0(2))%pot%set(1)%gp%values)) THEN
                   IF  (SIZE(potparm%pot(tmpij(1) ,tmpij(2) )%pot%set(1)%gp%values)==&
                        SIZE(potparm%pot(tmpij0(1),tmpij0(2))%pot%set(1)%gp%values)) THEN
                      IF (ANY(potparm%pot(tmpij(1) ,tmpij(2) )%pot%set(1)%gp%values /= &
                              potparm%pot(tmpij0(1),tmpij0(2))%pot%set(1)%gp%values))     check = .FALSE.
                   END IF
                END IF
             CASE default
                check = .TRUE.
             END SELECT
             IF (ALL(Cwork==pot_par(j,:)).AND.check) CYCLE
             Cwork = pot_par(j,:)
             nunique = nunique + 1
             iend    = j - 1
             CALL set_potparm_index(potparm, my_index(istart:iend), pot_target, ntype, tmpij,&
                  atomic_kind_set, my_do_14, scale_cutoff, ewald_env, my_ei_scale14,&
                  my_vdw_scale14, my_do_zbl, error)
             ! Single allocation potential
             IF (ANY(potential_single_allocation==potparm%pot(tmpij(2),tmpij(1))%pot%type(1))) THEN
                CALL reallocate(array_double_allocation, 1, SIZE(array_double_allocation)+1)
                array_double_allocation(SIZE(array_double_allocation)) = nunique
             END IF
             !
             DO i = istart, iend
                locij = my_index(i)
                CALL get_indexes(locij,ntype,tmpij)
                tmp_index(tmpij(1),tmpij(2))=nunique
                tmp_index(tmpij(2),tmpij(1))=nunique
             END DO
             istart = j
             locij = my_index(j)
             CALL get_indexes(locij,ntype,tmpij0)
          END DO
          nunique = nunique + 1
          iend    = ndim
          CALL set_potparm_index(potparm, my_index(istart:iend), pot_target, ntype, tmpij,&
               atomic_kind_set, my_do_14, scale_cutoff, ewald_env, my_ei_scale14,&
               my_vdw_scale14, my_do_zbl, error)
          ! Single allocation potential
          IF (ANY(potential_single_allocation==potparm%pot(tmpij(2),tmpij(1))%pot%type(1))) THEN
             CALL reallocate(array_double_allocation, 1, SIZE(array_double_allocation)+1)
             array_double_allocation(SIZE(array_double_allocation)) = nunique
          END IF
          DO i = istart, iend
             locij = my_index(i)
             CALL get_indexes(locij,ntype,tmpij)
             tmp_index(tmpij(1),tmpij(2))=nunique
             tmp_index(tmpij(2),tmpij(1))=nunique
          END DO
          DEALLOCATE(Cwork,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(pot_par,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ELSE
          nunique = nunique + 1
          CALL set_potparm_index(potparm, my_index, pot_target, ntype, tmpij,&
               atomic_kind_set, my_do_14, scale_cutoff, ewald_env, my_ei_scale14,&
               my_vdw_scale14, my_do_zbl, error)
       END IF
       DEALLOCATE(my_index,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END DO
    ! Multiple defined potential
    n     = 0
    DO i = 1, ntype
       DO j = 1, i
          n  = n + 1
          IF (SIZE(potparm%pot(i,j)%pot%type)==1) CYCLE
          nunique = nunique + 1
          tmp_index(i,j) = nunique
          tmp_index(j,i) = nunique
          !
          CALL set_potparm_index(potparm, (/n/), multi_type, ntype, tmpij,&
               atomic_kind_set, my_do_14, scale_cutoff, ewald_env, my_ei_scale14,&
               my_vdw_scale14, my_do_zbl, error)         
       END DO
    END DO
    ! Concluding the postprocess..
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, shell_present=shell_present)
    CALL spline_env_create ( spline_env, ntype, nunique, array_double_allocation, error=error )
    spline_env%spltab=tmp_index
    IF(shell_present) THEN
      IF(PRESENT(spline_c_coul_env)) THEN
        CALL spline_env_create ( spline_c_coul_env, ntype, nunique, array_double_allocation, error=error )
        spline_c_coul_env%spltab=tmp_index
      END IF
      IF(PRESENT(spline_cs_coul_env)) THEN
        CALL spline_env_create ( spline_cs_coul_env, ntype, nunique, array_double_allocation, error=error )
        spline_cs_coul_env%spltab=tmp_index
      END IF
      IF(PRESENT(spline_sc_coul_env)) THEN
        CALL spline_env_create ( spline_sc_coul_env, ntype, nunique, array_double_allocation, error=error )
        spline_sc_coul_env%spltab=tmp_index
      END IF
    END IF
    DEALLOCATE(tmp_index,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(array_double_allocation,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL timestop(handle)
  END SUBROUTINE get_nonbond_storage

! *****************************************************************************
!> \brief Trivial for non LJ potential.. gives back in the case of LJ
!>      the potparm with the smallest sigma..
!> \author Teodoro Laino [tlaino] 2007.06
! *****************************************************************************
  SUBROUTINE set_potparm_index(potparm, my_index, pot_target, ntype, tmpij_out,&
       atomic_kind_set, my_do_14, scale_cutoff, ewald_env, my_ei_scale14, my_vdw_scale14,&
       my_do_zbl, error)
    TYPE(pair_potential_pp_type), POINTER    :: potparm
    INTEGER, INTENT(IN)                      :: my_index(:), pot_target, ntype
    INTEGER, INTENT(OUT)                     :: tmpij_out(2)
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    LOGICAL, INTENT(IN)                      :: my_do_14, scale_cutoff
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    REAL(KIND=dp), INTENT(IN)                :: my_ei_scale14, my_vdw_scale14
    LOGICAL, INTENT(IN)                      :: my_do_zbl
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'set_potparm_index', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ewald_type, i, min_val, &
                                                nvalues, stat, tmpij(2), value
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: wrk
    LOGICAL                                  :: check, do_i_shell_model, &
                                                do_j_shell_model, &
                                                do_multipoles, &
                                                do_shell_model, failure
    REAL(KIND=dp) :: alpha, fac, hicut0, l_elec, l_elec_c, l_elec_cs, &
      l_elec_sc, l_epsilon, l_sigma6, m_elec, m_elec_c, m_elec_cs, m_elec_sc, &
      m_epsilon, m_sigma6, min_sigma6, qci, qcj, qi, qj, qsi, qsj, rci, rcj, &
      zi, zj
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: sigma6
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(pair_potential_single_type), &
      POINTER                                :: pot, pot_ref

    NULLIFY(pot, pot_ref)
    failure = .FALSE.
    do_shell_model = .FALSE.
    nvalues = SIZE(my_index)
    IF ((pot_target==lj_type).OR.(pot_target==lj_charmm_type)) THEN
       ALLOCATE(sigma6(nvalues),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(wrk(nvalues),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       min_sigma6= HUGE(0.0_dp)
       m_epsilon =-HUGE(0.0_dp)
       m_elec    = HUGE(0.0_dp)
       m_elec_c  = HUGE(0.0_dp)
       m_elec_cs = HUGE(0.0_dp)
       m_elec_sc = HUGE(0.0_dp)
       DO i = 1, nvalues
          value = my_index(i)
          CALL get_indexes(value,ntype,tmpij)
          pot => potparm%pot(tmpij(1),tmpij(2))%pot
          ! Preliminary check.. 
          check = SIZE(pot%type)==1
          CPPostcondition(check,cp_failure_level,routineP,error,failure)

          atomic_kind=> atomic_kind_set(tmpij(1))
          CALL get_atom_info(atomic_kind, my_do_14, qi, qci, qsi, do_shell_model=do_i_shell_model, error=error)
          atomic_kind=> atomic_kind_set(tmpij(2))
          CALL get_atom_info(atomic_kind, my_do_14, qj, qcj, qsj, do_shell_model=do_j_shell_model, error=error)
          do_shell_model   = do_shell_model .OR. do_i_shell_model .OR. do_j_shell_model
          sigma6(i) = pot%set(1)%lj%sigma6
          l_epsilon = pot%set(1)%lj%epsilon
          l_elec    = ABS(qsi * qsj)
          l_elec_c  = ABS(qci * qcj)
          l_elec_cs = ABS(qci * qsj)
          l_elec_sc = ABS(qsi * qcj)
          IF (sigma6(i) /= 0.0_dp) min_sigma6 =  MIN(min_sigma6, sigma6(i))
          IF (sigma6(i) == 0.0_dp) sigma6(i)  = -HUGE(0.0_dp)
          IF (l_epsilon /= 0.0_dp) m_epsilon  =  MAX(m_epsilon, l_epsilon)
          IF (l_elec    /= 0.0_dp) m_elec     =  MIN(m_elec   , l_elec   )
          IF (l_elec_c  /= 0.0_dp) m_elec_c   =  MIN(m_elec_c , l_elec_c )
          IF (l_elec_cs /= 0.0_dp) m_elec_cs  =  MIN(m_elec_cs, l_elec_cs)
          IF (l_elec_sc /= 0.0_dp) m_elec_sc  =  MIN(m_elec_sc, l_elec_sc)
       END DO
       CALL sort(sigma6,nvalues,wrk)
       min_val  = my_index(wrk(nvalues))
       m_sigma6 = sigma6(nvalues)
       ! In case there are only zeros.. let's consider them properly..
       IF (m_sigma6   == -HUGE(0.0_dp)) m_sigma6   = 1.0_dp
       IF (m_epsilon  == -HUGE(0.0_dp)) m_epsilon  = 0.0_dp
       IF (min_sigma6 ==  HUGE(0.0_dp)) min_sigma6 = 0.0_dp
       IF (m_elec     ==  HUGE(0.0_dp)) m_elec     = 0.0_dp
       IF (m_elec_c   ==  HUGE(0.0_dp)) m_elec_c   = 0.0_dp
       IF (m_elec_cs  ==  HUGE(0.0_dp)) m_elec_cs  = 0.0_dp
       IF (m_elec_sc  ==  HUGE(0.0_dp)) m_elec_sc  = 0.0_dp
       DEALLOCATE(sigma6,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(wrk,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ELSE
       DO i = 1, nvalues
          value = my_index(i)
          CALL get_indexes(value,ntype,tmpij)
          atomic_kind=> atomic_kind_set(tmpij(1))
          CALL get_atom_info(atomic_kind, my_do_14, do_shell_model=do_i_shell_model, error=error)
          atomic_kind=> atomic_kind_set(tmpij(2))
          CALL get_atom_info(atomic_kind, my_do_14, do_shell_model=do_j_shell_model, error=error)
          do_shell_model   = do_shell_model .OR. do_i_shell_model .OR. do_j_shell_model
       END DO
       min_val = MINVAL(my_index(:))
    END IF
    CALL get_indexes(min_val,ntype,tmpij)
    tmpij_out = tmpij
    pot => potparm%pot(tmpij(1),tmpij(2))%pot
    pot%undef = .TRUE.
    IF (scale_cutoff) THEN 
       CALL ewald_env_get(ewald_env,ewald_type=ewald_type,alpha=alpha,&
                          do_multipoles=do_multipoles,error=error)
       fac = my_ei_scale14
       IF (do_multipoles) fac = 0.0_dp
       hicut0 = SQRT(pot%rcutsq)
       IF (ABS(hicut0) <= MIN_HICUT_VALUE) hicut0 = DEFAULT_HICUT_VALUE
    END IF
    IF (pot_target==gp_type) CALL init_genpot(potparm, ntype)
    DO i = 1, nvalues
       value = my_index(i)
       CALL get_indexes(value,ntype,tmpij)
       pot => potparm%pot(tmpij(1),tmpij(2))%pot
       CALL spline_factor_create(pot%spl_f, error)
       IF (do_shell_model) THEN
          CALL spline_factor_create(pot%spl_c,  error)
          CALL spline_factor_create(pot%spl_sc, error)
          CALL spline_factor_create(pot%spl_cs, error)
       END IF
       atomic_kind=> atomic_kind_set(tmpij(1))
       CALL get_atom_info(atomic_kind, my_do_14, qi, qci, qsi, z=zi, rc=rci, do_shell_model=do_i_shell_model, error=error)
       atomic_kind=> atomic_kind_set(tmpij(2))
       CALL get_atom_info(atomic_kind, my_do_14, qj, qcj, qsj, z=zj, rc=rcj, do_shell_model=do_j_shell_model, error=error)
       l_elec    = qsi * qsj
       l_elec_c  = qci * qcj
       l_elec_cs = qci * qsj
       l_elec_sc = qsi * qcj

       IF (ANY(potential_single_allocation==pot_target)) THEN
          check = SIZE(pot%type)==1
          CPPostcondition(check,cp_failure_level,routineP,error,failure)
          ! Undef potential.. this will be used to compute the splines..
          ! Coulomb Potential
          pot%e_fac = -m_elec
          pot%z1 = zi
          pot%z2 = zj

          pot%spl_f%rscale(1) = 1.0_dp 
          pot%spl_f%fscale(1) = 0.0_dp
          IF (l_elec /=0.0_dp) THEN
             pot%spl_f%fscale(1) = -l_elec/m_elec
          END IF
          ! Non bonded potential
          IF ((pot_target==lj_type).OR.(pot_target==lj_charmm_type)) THEN
             l_sigma6  = pot%set(1)%lj%sigma6
             l_epsilon = pot%set(1)%lj%epsilon
             ! Undef potential.. this will be used to compute the splines..
             IF (pot%undef) THEN
                pot%set(1)%lj%sigma6  = m_sigma6
                pot%set(1)%lj%sigma12 = m_sigma6**2
                pot%set(1)%lj%epsilon = m_epsilon
             END IF
             pot%spl_f%rscale(2) = 1.0_dp
             pot%spl_f%fscale(2) = 0.0_dp
             IF (l_sigma6*l_epsilon/=0.0_dp) THEN
                pot%spl_f%rcutsq_f  = (min_sigma6/m_sigma6)**(1.0_dp/3.0_dp)
                pot%spl_f%rscale(2) = (l_sigma6/m_sigma6)**(1.0_dp/3.0_dp)
                pot%spl_f%fscale(2) =  l_epsilon/m_epsilon
             END IF
          END IF
          IF (do_shell_model) THEN
             ! Undef potential.. this will be used to compute the splines..
             pot%e_fcc = -m_elec_c
             pot%e_fcs = -m_elec_cs
             pot%e_fsc = -m_elec_sc

             ! Core-Core
             pot%spl_c%rscale(1)  = 1.0_dp
             pot%spl_c%fscale(1)  = 0.0_dp
             IF (l_elec_c /=0.0_dp) THEN
                pot%spl_c%fscale(1)  = -l_elec_c/m_elec_c
             END IF
             pot%spl_c%rscale(2)  = 0.0_dp
             pot%spl_c%fscale(2)  = 0.0_dp
             ! Shell-Core
             pot%spl_sc%rscale(1) = 1.0_dp
             pot%spl_sc%fscale(1) = 0.0_dp
             IF (l_elec_sc /=0.0_dp) THEN
                pot%spl_sc%fscale(1) = -l_elec_sc/m_elec_sc
             END IF
             pot%spl_sc%rscale(2) = 0.0_dp
             pot%spl_sc%fscale(2) = 0.0_dp
             ! Core-Shell
             pot%spl_cs%rscale(1) = 1.0_dp
             pot%spl_cs%fscale(1) = 0.0_dp
             IF (l_elec_cs /=0.0_dp) THEN
                pot%spl_cs%fscale(1) = -l_elec_cs/m_elec_cs
             END IF
             pot%spl_cs%rscale(2) = 0.0_dp
             pot%spl_cs%fscale(2) = 0.0_dp
          END IF
          IF(my_do_zbl) THEN
             CALL zbl_matching_polinomial(pot, rci, rcj)
          END IF
       ELSE
          pot%spl_f%rcutsq_f  = 1.0_dp
          pot%spl_f%rscale    = 1.0_dp
          pot%spl_f%fscale    = 1.0_dp
          pot%e_fac = l_elec
          pot%e_fcc = l_elec_c
          pot%e_fcs = l_elec_cs
          pot%e_fsc = l_elec_sc
          pot%z1 = zi
          pot%z2 = zj
          IF (do_shell_model) THEN
             pot%spl_c%rscale  = 1.0_dp
             pot%spl_c%fscale  = 1.0_dp
             
             pot%spl_sc%rscale = 1.0_dp
             pot%spl_sc%fscale = 1.0_dp
             
             pot%spl_cs%rscale = 1.0_dp
             pot%spl_cs%fscale = 1.0_dp
          END IF
          IF (my_do_zbl) THEN
             CALL zbl_matching_polinomial(pot, rci, rcj)
          END IF
       END IF
       ! Derivative factors
       pot%spl_f%dscale  = pot%spl_f%fscale  / pot%spl_f%rscale
       IF (do_shell_model) THEN
          pot%spl_c%dscale  = pot%spl_c%fscale  / pot%spl_c%rscale
          pot%spl_sc%dscale = pot%spl_sc%fscale / pot%spl_sc%rscale
          pot%spl_cs%dscale = pot%spl_cs%fscale / pot%spl_cs%rscale
       END IF
       ! Cutoff for the potentials on splines
       IF (scale_cutoff) THEN
          ! Cutoff NonBonded
          pot%spl_f%cutoff_n = ener_pot(pot, hicut0, 0.0_dp, my_vdw_scale14)
          ! Cutoff Coulomb
          pot%spl_f%cutoff_e = ener_coul(pot,0.0_dp, fac*l_elec, hicut0, ewald_type, alpha, my_do_14, zbl=my_do_zbl )
          IF(do_shell_model) THEN
             pot%spl_c%cutoff_e  = ener_coul(pot,0.0_dp, fac*l_elec_c,  hicut0, ewald_type, alpha, my_do_14,&
                zbl=my_do_zbl )
             pot%spl_sc%cutoff_e = ener_coul(pot,0.0_dp, fac*l_elec_sc, hicut0, ewald_type, alpha, my_do_14,&
                zbl=my_do_zbl )
             pot%spl_cs%cutoff_e = ener_coul(pot,0.0_dp, fac*l_elec_cs, hicut0, ewald_type, alpha, my_do_14,&
                zbl=my_do_zbl )
          END IF
       END IF
    END DO
    ! Handle the cutoff
    IF (scale_cutoff) THEN
       pot_ref => potparm%pot(tmpij_out(1),tmpij_out(2))%pot
       DO i = 1, nvalues
          value = my_index(i)
          CALL get_indexes(value,ntype,tmpij)
          pot => potparm%pot(tmpij(1),tmpij(2))%pot
          IF (value==min_val) CYCLE
          ! Cutoff NonBonded
          pot%spl_f%cutoff_n = pot_ref%spl_f%cutoff_n*pot%spl_f%fscale(2) - pot%spl_f%cutoff_n
          ! Cutoff Coulomb
          pot%spl_f%cutoff_e = pot_ref%spl_f%cutoff_e*pot%spl_f%fscale(1) - pot%spl_f%cutoff_e 
          pot%spl_f%cutoff   = pot%spl_f%cutoff_n + pot%spl_f%cutoff_e
          IF(do_shell_model) THEN
             pot%spl_c%cutoff_n  = 0.0_dp
             pot%spl_sc%cutoff_n = 0.0_dp
             pot%spl_cs%cutoff_n = 0.0_dp
             pot%spl_c%cutoff_e  = pot_ref%spl_c%cutoff_e  * pot%spl_c%fscale(1)  - pot%spl_c%cutoff_e 
             pot%spl_sc%cutoff_e = pot_ref%spl_sc%cutoff_e * pot%spl_sc%fscale(1) - pot%spl_sc%cutoff_e
             pot%spl_cs%cutoff_e = pot_ref%spl_cs%cutoff_e * pot%spl_cs%fscale(1) - pot%spl_cs%cutoff_e
             pot%spl_c%cutoff    = pot%spl_c%cutoff_n  + pot%spl_c%cutoff_e
             pot%spl_cs%cutoff   = pot%spl_cs%cutoff_n + pot%spl_cs%cutoff_e
             pot%spl_sc%cutoff   = pot%spl_sc%cutoff_n + pot%spl_sc%cutoff_e
          END IF
       END DO
    END IF
    IF (pot_target==gp_type) CALL finalizef()
  END SUBROUTINE set_potparm_index

! *****************************************************************************
!> \brief Gives back the indices of the matrix w.r.t. the collective array index
!> \author Teodoro Laino [tlaino] 2006.05
! *****************************************************************************
  SUBROUTINE get_indexes(Inind, ndim, ij)
    INTEGER, INTENT(IN)                      :: Inind, ndim
    INTEGER, DIMENSION(2), INTENT(OUT)       :: ij

    INTEGER                                  :: i, tmp

    tmp = 0
    ij  = HUGE(0)
    DO i = 1, ndim
       tmp = tmp + i
       IF (tmp >= Inind) THEN
          ij(1) = i
          ij(2) = Inind - tmp + i
          EXIT
       END IF
    END DO
  END SUBROUTINE get_indexes

END MODULE pair_potential

