!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
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
                                             get_atomic_kind
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE fparser,                         ONLY: finalizef,&
                                             initf,&
                                             parsef
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE pair_potential_types,            ONLY: &
       b4_type, bm_type, compare_pot, ea_type, ft_type, ftd_type, gp_type, &
       gw_type, ip_type, list_pot, lj_charmm_type, lj_type, multi_type, &
       nn_type, pair_potential_pp_type, pair_potential_single_type, &
       potential_single_allocation, quip_type, siepmann_type, tersoff_type, &
       wl_type
  USE pair_potential_util,             ONLY: ener_pot,&
                                             ener_zbl,&
                                             zbl_matching_polinomial
  USE physcon,                         ONLY: bohr,&
                                             evolt,&
                                             kjmol
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
#include "./common/cp_common_uses.f90"

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
!> \param potparm ...
!> \param ntype ...
!> \par History
!>      Teo 2007.06 - Zurich University
! *****************************************************************************
  SUBROUTINE init_genpot(potparm, ntype)
    TYPE(pair_potential_pp_type), POINTER    :: potparm
    INTEGER, INTENT(IN)                      :: ntype

    CHARACTER(len=*), PARAMETER :: routineN = 'init_genpot', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, j, k, ngp
    TYPE(pair_potential_single_type), &
      POINTER                                :: pot

    CALL timeset(routineN,handle)

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
    CALL timestop(handle)

  END SUBROUTINE init_genpot

! *****************************************************************************
!> \brief creates the splines for the potentials
!> \param spline_env ...
!> \param potparm ...
!> \param atomic_kind_set ...
!> \param eps_spline ...
!> \param max_energy ...
!> \param rlow_nb ...
!> \param emax_spline ...
!> \param npoints ...
!> \param iw ...
!> \param iw2 ...
!> \param iw3 ...
!> \param do_zbl ...
!> \param shift_cutoff ...
!> \param nonbonded_type ...
!> \par History
!>      Teo 2006.05 : Improved speed and accuracy. Linear scaling of the setup
! *****************************************************************************
  SUBROUTINE spline_nonbond_control(spline_env, potparm, atomic_kind_set, &
       eps_spline, max_energy, rlow_nb, emax_spline, npoints, iw, iw2, iw3, do_zbl, &
       shift_cutoff, nonbonded_type)

    TYPE(spline_environment_type), POINTER   :: spline_env
    TYPE(pair_potential_pp_type), POINTER    :: potparm
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    REAL(KIND=dp), INTENT(IN)                :: eps_spline, max_energy, &
                                                rlow_nb, emax_spline
    INTEGER, INTENT(IN)                      :: npoints, iw, iw2, iw3
    LOGICAL, INTENT(IN)                      :: do_zbl, shift_cutoff
    CHARACTER(LEN=*), INTENT(IN)             :: nonbonded_type

    CHARACTER(len=*), PARAMETER :: routineN = 'spline_nonbond_control', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ip, j, k, n, &
                                                ncount, npoints_spline, ntype
    LOGICAL                                  :: failure, found_locut
    REAL(KIND=dp)                            :: energy_cutoff, hicut, hicut0, &
                                                locut
    TYPE(pair_potential_single_type), &
      POINTER                                :: pot

    failure= .FALSE.
    n      = 0
    ncount = 0

    ntype = SIZE ( atomic_kind_set )
    CALL timeset(routineN,handle)
    IF (iw3 > 0) THEN
       WRITE (iw3,"(/,T2,A,I0,A,I0,A)")&
        "SPLINE_INFO| Generating ",(ntype*(ntype+1))/2," splines for "//&
        TRIM(ADJUSTL(nonbonded_type))//" interactions "
       WRITE (iw3,"(T2,A,I0,A)")&
        "             Due to ",ntype," different atomic kinds"
    END IF
    CALL init_genpot(potparm, ntype)
    ! Real computation of splines
    ip = 0
    DO i = 1, ntype
       DO j = 1, i
          pot => potparm%pot(i,j)%pot
          IF (iw3>0 .AND. iw<=0) THEN
             IF (MOD(i*(i-1)/2+j,MAX(1,(ntype*(ntype+1))/(2*10)))==0) THEN
                 WRITE(UNIT=iw3,ADVANCE="NO",FMT='(2X,A3,I0)') '...',i*(i-1)/2+j
                 ip = ip+1
                 IF (ip >= 11) THEN
                   WRITE(iw3,*)
                   ip = 0
                 END IF
             END IF
          ENDIF
          ! Setup of Exclusion Types
          pot%no_pp = .TRUE.
          pot%no_mb = .TRUE.
          DO k = 1, SIZE(pot%type)
             SELECT CASE(pot%type(k))
             CASE (lj_type,lj_charmm_type,wl_type,gw_type,ft_type, ftd_type, ip_type,&
                   b4_type,bm_type,gp_type,ea_type,quip_type)
                pot%no_pp = .FALSE.
             CASE (tersoff_type)
                pot%no_mb = .FALSE.
             CASE (siepmann_type)
                pot%no_mb = .FALSE.
             CASE (nn_type)
                ! Do nothing..
             CASE DEFAULT
                ! Never reach this point
                CPASSERT(.FALSE.)
             END SELECT
             ! Special case for EAM
             SELECT CASE(pot%type(k))
             CASE (ea_type,quip_type)
                pot%no_mb = .FALSE.
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

          energy_cutoff = pot%spl_f%cutoff

          ! Find the real locut according emax_spline
          CALL get_spline_cutoff(hicut, locut, found_locut, pot, do_zbl, &
               energy_cutoff, emax_spline)
          locut=MAX(locut*SQRT(pot%spl_f%rcutsq_f),rlow_nb)

          ! Real Generation of the Spline
          npoints_spline=npoints
          CALL generate_spline_low(spline_env%spl_pp(n)%spl_p, npoints_spline, locut, &
               hicut, eps_spline, iw, iw2, i, j, n, ncount, max_energy, pot, &
               energy_cutoff, found_locut, do_zbl, atomic_kind_set, &
               nonbonded_type)

          pot%undef = .FALSE.
          ! Unique Spline working only for a pure LJ potential..
          IF (SIZE(pot%type)==1) THEN
             IF (ANY(potential_single_allocation==pot%type(1))) THEN
                ! Restoring the proper values for the generating spline pot
                IF ((pot%type(1)==lj_type).OR.(pot%type(1)==lj_charmm_type)) THEN
                   pot%set(1)%lj%sigma6  = pot%set(1)%lj%sigma6*pot%spl_f%rscale(1)**3
                   pot%set(1)%lj%sigma12 = pot%set(1)%lj%sigma6**2
                   pot%set(1)%lj%epsilon = pot%set(1)%lj%epsilon*pot%spl_f%fscale(1)
                END IF
             END IF
          END IF
          ! Correct Cutoff...
          IF (shift_cutoff) THEN
             pot%spl_f%cutoff = pot%spl_f%cutoff*pot%spl_f%fscale(1) - &
                ener_pot(pot, hicut0, 0.0_dp)
          END IF
       END DO
    END DO
    CALL finalizef()

    IF (iw > 0) THEN
       WRITE(iw,'(/,T2,A,I0)')&
        "SPLINE_INFO| Number of pair potential splines allocated:   ",MAXVAL(spline_env%spltab)
    END IF
    IF (iw3 > 0) THEN
       WRITE (iw3,'(/,T2,A,I0,/,T2,A)')&
        "SPLINE_INFO| Number of unique splines computed:            ",MAXVAL(spline_env%spltab),&
        "SPLINE_INFO| Done"
    END IF

    CALL timestop(handle)

  END SUBROUTINE spline_nonbond_control

! *****************************************************************************
!> \brief Finds the cutoff for the generation of the spline
!>      In a two pass approach, first with low resolution, refine in a second iteration
!> \param hicut ...
!> \param locut ...
!> \param found_locut ...
!> \param pot ...
!> \param do_zbl ...
!> \param energy_cutoff ...
!> \param emax_spline ...
!> \par History
!>      Splitting in order to make some season cleaning..
!> \author Teodoro Laino [tlaino] 2007.06
! *****************************************************************************
  SUBROUTINE get_spline_cutoff(hicut, locut, found_locut, pot, do_zbl, &
       energy_cutoff, emax_spline)

    REAL(KIND=dp), INTENT(IN)                :: hicut
    REAL(KIND=dp), INTENT(INOUT)             :: locut
    LOGICAL, INTENT(OUT)                     :: found_locut
    TYPE(pair_potential_single_type), &
      OPTIONAL, POINTER                      :: pot
    LOGICAL, INTENT(IN)                      :: do_zbl
    REAL(KIND=dp), INTENT(IN)                :: energy_cutoff, emax_spline

    CHARACTER(len=*), PARAMETER :: routineN = 'get_spline_cutoff', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ilevel, jx
    REAL(KIND=dp)                            :: dx2, e, locut_found, x

    dx2 = (hicut-locut)
    x = hicut
    locut_found = locut
    found_locut = .FALSE.
    DO ilevel=1,2
       dx2 = dx2 / 100.0_dp
       DO jx = 1, 100
          e = ener_pot(pot, x, energy_cutoff)
          IF (do_zbl) THEN
             e = e + ener_zbl(pot, x)
          END IF
          IF (ABS(e) > emax_spline) THEN
             locut_found = x
             found_locut = .TRUE.
             EXIT
          END IF
          x = x - dx2
       END DO
       x= x + dx2
    ENDDO
    locut = locut_found

  END SUBROUTINE get_spline_cutoff

! *****************************************************************************
!> \brief Real Generation of spline..
!> \param spl_p ...
!> \param npoints ...
!> \param locut ...
!> \param hicut ...
!> \param eps_spline ...
!> \param iw ...
!> \param iw2 ...
!> \param i ...
!> \param j ...
!> \param n ...
!> \param ncount ...
!> \param max_energy ...
!> \param pot ...
!> \param energy_cutoff ...
!> \param found_locut ...
!> \param do_zbl ...
!> \param atomic_kind_set ...
!> \param nonbonded_type ...
!> \par History
!>      Splitting in order to make some season cleaning..
!> \author Teodoro Laino [tlaino] 2007.06
! *****************************************************************************
  SUBROUTINE generate_spline_low(spl_p, npoints, locut, hicut,  eps_spline, &
       iw, iw2, i, j, n, ncount, max_energy, pot, energy_cutoff, &
       found_locut, do_zbl, atomic_kind_set, nonbonded_type)

    TYPE(spline_data_p_type), DIMENSION(:), &
      POINTER                                :: spl_p
    INTEGER, INTENT(INOUT)                   :: npoints
    REAL(KIND=dp), INTENT(IN)                :: locut, hicut, eps_spline
    INTEGER, INTENT(IN)                      :: iw, iw2, i, j, n, ncount
    REAL(KIND=dp), INTENT(IN)                :: max_energy
    TYPE(pair_potential_single_type), &
      POINTER                                :: pot
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: energy_cutoff
    LOGICAL, INTENT(IN)                      :: found_locut, do_zbl
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    CHARACTER(LEN=*), INTENT(IN)             :: nonbonded_type

    CHARACTER(len=*), PARAMETER :: routineN = 'generate_spline_low', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=2*default_string_length)   :: message, tmp
    CHARACTER(LEN=default_path_length)       :: file_name
    INTEGER                                  :: ix, jx, mfac, nppa, nx, &
                                                unit_number
    LOGICAL                                  :: fixed_spline_points
    REAL(KIND=dp)                            :: df, dg, dh, diffmax, dx, dx2, &
                                                e, e_spline, f, g, h, r, &
                                                rcut, x, x2, xdum, xdum1, xsav
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(spline_data_type), POINTER          :: spline_data
    TYPE(spline_factor_type), POINTER        :: spl_f

    NULLIFY (logger,spl_f)
    logger => cp_get_default_logger()

    CALL spline_factor_create(spl_f)
    mfac    =  5
    IF (npoints>0) THEN
       fixed_spline_points=.TRUE.
    ELSE
       fixed_spline_points=.FALSE.
       npoints = 20
       IF (.NOT. found_locut) npoints = 2
    ENDIF
    spline_data => spl_p(1)%spline_data
    DO WHILE (.TRUE.)
       CALL init_splinexy ( spline_data, npoints+1)
       dx2 = (1.0_dp/locut**2-1.0_dp/hicut**2)/REAL(npoints,KIND=dp)
       x2  = 1.0_dp/hicut**2
       spline_data % x1 = x2
       DO jx = 1, npoints + 1
          ! jx: loop over 1/distance**2
          x = SQRT ( 1.0_dp/x2 )
          e = ener_pot(pot, x, energy_cutoff)
          IF (do_zbl) THEN
            e = e + ener_zbl(pot, x)
          END IF
          spline_data % y(jx) = e
          x2 = x2 + dx2
       END DO
       CALL init_spline(spline_data,dx=dx2)
       ! This is the check for required accuracy on spline setup
       dx2 = (hicut - locut)/REAL(mfac*npoints+1,KIND=dp)
       x2  = locut + dx2
       diffmax = -1.0_dp
       xsav = hicut
       ! if a fixed number of points is requested, no check on its error
       IF (fixed_spline_points) EXIT
       DO jx = 1, mfac*npoints
          x = x2
          e = ener_pot(pot, x, energy_cutoff)
          IF (do_zbl) THEN
            e = e + ener_zbl(pot, x)
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
       ! accuracy is poor or we have found no points below max_energy, refine mesh
       IF (diffmax > eps_spline .OR. diffmax < 0.0_dp ) THEN
          npoints   = CEILING(1.2_dp*REAL(npoints,KIND=dp))
       ELSE
          EXIT
       END IF
    END DO
    ! Print spline info to STDOUT if requested
    IF (iw > 0) THEN
       WRITE (UNIT=iw,&
              FMT="(/,A,I0,/,A,I0,/,A,I0,1X,I0,/,A,/,A,I0,2(/,A,ES13.6),2(/,A,2ES13.6))")&
        " SPLINE_INFO| Spline number:                                ",ncount,&
        " SPLINE_INFO| Unique spline number:                         ",n,&
        " SPLINE_INFO| Atomic kind numbers:                          ",i,j,&
        " SPLINE_INFO| Atomic kind names:                            "//TRIM(ADJUSTL(atomic_kind_set(i)%name))//" "//&
                                                                        TRIM(ADJUSTL(atomic_kind_set(j)%name)),&
        " SPLINE_INFO| Number of spline points:                      ",npoints,&
        " SPLINE_INFO| Requested accuracy [Hartree]:                ",eps_spline,&
        " SPLINE_INFO| Achieved accuracy [Hartree]:                 ",diffmax,&
        " SPLINE_INFO| Spline range [bohr]:                         ",locut,hicut,&
        " SPLINE_INFO| Spline range used to achieve accuracy [bohr]:",xsav,hicut
       dx2 = (hicut - locut)/REAL(npoints+1,KIND=dp)
       x  = locut + dx2
       WRITE (UNIT=iw,FMT='(A,ES17.9)')&
        " SPLINE_INFO| Spline value at RMIN [Hartree]:             ",potential_s(spl_p,x*x,xdum,spl_f,logger),&
        " SPLINE_INFO| Spline value at RMAX [Hartree]:             ",potential_s(spl_p,hicut*hicut,xdum,spl_f,logger),&
        " SPLINE_INFO| Non-bonded energy cutoff [Hartree]:         ",energy_cutoff
    END IF
    ! Print spline data on file if requested
    IF (iw2 > 0) THEN
       ! Set increment to 200 points per Angstrom
       nppa = 200
       dx = bohr/REAL(nppa,KIND=dp)
       nx = NINT(hicut/dx)
       file_name = ""
       tmp = ADJUSTL(cp_to_string(n))
       WRITE (UNIT=file_name,FMT="(A,I0,A)")&
        TRIM(ADJUSTL(nonbonded_type))//"_SPLINE_"//TRIM(tmp)//"_"//&
        TRIM(ADJUSTL(atomic_kind_set(i)%name))//"_"//&
        TRIM(ADJUSTL(atomic_kind_set(j)%name))
       CALL open_file(file_name=file_name,&
                      file_status="UNKNOWN",&
                      file_form="FORMATTED",&
                      file_action="WRITE",&
                      unit_number=unit_number)
       WRITE (UNIT=unit_number,&
              FMT="(2(A,I0,/),A,I0,1X,I0,/,A,/,A,I0,2(/,A,ES13.6),2(/,A,2ES13.6),/,A,ES13.6,/,A,I0,A,/,A)")&
        "# Spline number:                                    ",ncount,&
        "# Unique spline number:                             ",n,&
        "# Atomic kind numbers:                              ",i,j,&
        "# Atomic kind names:                                "//TRIM(ADJUSTL(atomic_kind_set(i)%name))//" "//&
                                                                TRIM(ADJUSTL(atomic_kind_set(j)%name)),&
        "# Number of spline points:                          ",npoints,&
        "# Requested accuracy [eV]:                         ",eps_spline*evolt,&
        "# Achieved accuracy [eV]:                          ",diffmax*evolt,&
        "# Spline range [Angstrom]:                         ",locut/bohr,hicut/bohr,&
        "# Spline range used to achieve accuracy [Angstrom]:",xsav/bohr,hicut/bohr,&
        "# Non-bonded energy cutoff [eV]:                   ",energy_cutoff*evolt,&
        "# Test spline using ",nppa," points per Angstrom:",&
        "#     Abscissa [Angstrom]              Energy [eV]      Splined energy [eV] Derivative [eV/Angstrom]"//&
        "      |Energy error| [eV]"
       x = 0.0_dp
       DO jx=0,nx
          IF (x > hicut) x = hicut
          IF (x > locut) THEN
             e = ener_pot(pot,x,energy_cutoff)
             IF (do_zbl) e = e + ener_zbl(pot,x)
             e_spline = potential_s(spl_p,x*x,xdum,spl_f,logger)
             WRITE (UNIT=unit_number,FMT="(5ES25.12)")&
              x/bohr,e*evolt,e_spline*evolt,-bohr*x*xdum*evolt,ABS((e - e_spline)*evolt)
          END IF
          x = x + dx
       END DO
       CALL close_file(unit_number=unit_number)
       !MK Write table.xvf for GROMACS 4.5.5
       WRITE (UNIT=file_name,FMT="(A,I0,A)")&
        "table_"//&
        TRIM(ADJUSTL(atomic_kind_set(i)%name))//"_"//&
        TRIM(ADJUSTL(atomic_kind_set(j)%name))//".xvg"
       CALL open_file(file_name=file_name,&
                      file_status="UNKNOWN",&
                      file_form="FORMATTED",&
                      file_action="WRITE",&
                      unit_number=unit_number)
       ! Recommended increment for dp is 0.0005 nm = 0.005 Angstrom
       ! which are 200 points/Angstrom
       rcut = 0.1_dp*hicut/bohr 
       x = 0.0_dp
       DO jx=0,nx
          IF (x > hicut) x = hicut
          r = 0.1_dp*x/bohr ! Convert bohr to nm
          IF (x <= locut) THEN
             WRITE (UNIT=unit_number,FMT="(7ES25.12)")&
              r,(0.0_dp,ix=1,6)
          ELSE
             e_spline = potential_s(spl_p,x*x,xdum,spl_f,logger)
             f  =  1.0_dp/r
             df = -1.0_dp/r**2
             g  = -1.0_dp/r**6 + 1.0_dp/rcut**6
             dg =  6.0_dp/r**7
             h  =  e_spline*kjmol
             dh = -10.0_dp*bohr*x*xdum*kjmol
             WRITE (UNIT=unit_number,FMT="(7ES25.12)")&
              r,f,-df,& ! r, f(r), -f'(r) => probably not used
              g,-dg,&   !    g(r), -g'(r) => not used, if C = 0
              h,-dh     !    h(r), -h'(r) => used, if A = 1
          END IF
          x = x + dx
       END DO
       CALL close_file(unit_number=unit_number)
    END IF

    CALL spline_factor_release(spl_f)

  END SUBROUTINE generate_spline_low

! *****************************************************************************
!> \brief Prescreening of the effective bonds evaluations. linear scaling algorithm
!> \param spline_env ...
!> \param potparm ...
!> \param atomic_kind_set ...
!> \param do_zbl ...
!> \param shift_cutoff ...
!> \author Teodoro Laino [tlaino] 2006.05
! *****************************************************************************
  SUBROUTINE get_nonbond_storage(spline_env, potparm, atomic_kind_set, do_zbl, &
       shift_cutoff)

    TYPE(spline_environment_type), POINTER   :: spline_env
    TYPE(pair_potential_pp_type), POINTER    :: potparm
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    LOGICAL, INTENT(IN)                      :: do_zbl, shift_cutoff

    CHARACTER(len=*), PARAMETER :: routineN = 'get_nonbond_storage', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, i, idim, iend, istart, j, k, locij, n, ndim, nk, &
      ntype, nunique, nvar, pot_target, tmpij(2), tmpij0(2)
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: Iwork1, Iwork2, my_index
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: tmp_index
    LOGICAL                                  :: at_least_one, check, failure
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: Cwork, Rwork, wtmp
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: pot_par

    failure     = .FALSE.
    CALL timeset(routineN,handle)

    ntype = SIZE(atomic_kind_set)
    DO i = 1, ntype
       DO j = 1, i
          potparm%pot(i,j)%pot%undef = .FALSE.
       END DO
    END DO
    ALLOCATE(tmp_index(ntype,ntype))
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
       SELECT CASE(pot_target)
       CASE (lj_type,lj_charmm_type)
          nvar = 3 + nvar
       CASE (wl_type)
          nvar = 3 + nvar
       CASE (gw_type)
          nvar = 5 + nvar
       CASE (ea_type)
          nvar = 4 + nvar
       CASE (quip_type)
          nvar = 1 + nvar
       CASE (ft_type)
          nvar = 4 + nvar
       CASE (ftd_type)
          nvar = 5 + nvar
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
       CASE (siepmann_type)
          nvar = 5 + nvar
       CASE (nn_type)
          nvar = nvar
       CASE DEFAULT
          CPASSERT(.FALSE.)
       END SELECT
       ! Setup a table of the indexes..
       ALLOCATE(my_index(ndim))
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
          ALLOCATE(pot_par(ndim,nvar))
          n  = 0
          nk = 0
          DO i = 1, ntype
             DO j = 1, i
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
                   CASE (quip_type)
                      pot_par(nk,1)  = str2id(&
                                        TRIM(potparm%pot(i,j)%pot%set(1)%quip%quip_file_name)// &
                                        TRIM(potparm%pot(i,j)%pot%set(1)%quip%init_args)//&
                                        TRIM(potparm%pot(i,j)%pot%set(1)%quip%calc_args))
                   CASE (ft_type)
                      pot_par(nk,1)  = potparm%pot(i,j)%pot%set(1)%ft%A
                      pot_par(nk,2)  = potparm%pot(i,j)%pot%set(1)%ft%B
                      pot_par(nk,3)  = potparm%pot(i,j)%pot%set(1)%ft%C
                      pot_par(nk,4)  = potparm%pot(i,j)%pot%set(1)%ft%D
                   CASE (ftd_type)
                      pot_par(nk,1)  = potparm%pot(i,j)%pot%set(1)%ftd%A
                      pot_par(nk,2)  = potparm%pot(i,j)%pot%set(1)%ftd%B
                      pot_par(nk,3)  = potparm%pot(i,j)%pot%set(1)%ftd%C
                      pot_par(nk,4)  = potparm%pot(i,j)%pot%set(1)%ftd%D
                      pot_par(nk,5)  = potparm%pot(i,j)%pot%set(1)%ftd%BD
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
                   CASE (siepmann_type)
                      pot_par(nk,1)  = potparm%pot(i,j)%pot%set(1)%siepmann%B
                      pot_par(nk,2)  = potparm%pot(i,j)%pot%set(1)%siepmann%D
                      pot_par(nk,3)  = potparm%pot(i,j)%pot%set(1)%siepmann%E
                      pot_par(nk,4)  = potparm%pot(i,j)%pot%set(1)%siepmann%F
                      pot_par(nk,5)  = potparm%pot(i,j)%pot%set(1)%siepmann%beta
                   CASE (nn_type)
                      ! no checks
                   CASE DEFAULT
                      CPASSERT(.FALSE.)
                   END SELECT
                   IF (ANY(potential_single_allocation==pot_target)) THEN
                      pot_par(nk,:) = REAL(pot_target, KIND=dp)
                   END IF
                END IF
             END DO
          END DO
          ! Main Sorting Loop
          ALLOCATE(Rwork(ndim))
          ALLOCATE(Iwork1(ndim))
          ALLOCATE(Iwork2(ndim))
          ALLOCATE(wtmp(nvar))
          CALL sort(pot_par(:,1),ndim,Iwork1)
          ! Sort all the other components of the potential
          DO k = 2, nvar
             Rwork(:) = pot_par(:,k)
             DO i = 1, ndim
                pot_par(i,k)=Rwork(Iwork1(i))
             END DO
          END DO
          Iwork2(:) = my_index
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
                Rwork(:) = pot_par(:,j)
                DO i = 1, ndim
                   pot_par(i,j)=Rwork(Iwork1(i))
                END DO
             END DO
             Iwork2(:) = my_index
             DO i = 1, ndim
                my_index(i) = Iwork2(Iwork1(i))
             END DO
          END DO
          DEALLOCATE(wtmp)
          DEALLOCATE(Iwork1)
          DEALLOCATE(Iwork2)
          DEALLOCATE(Rwork)
          !
          ! Let's determine the number of unique potentials and tag them
          !
          ALLOCATE(Cwork(nvar))
          Cwork(:) = pot_par(1,:)
          locij = my_index(1)
          CALL get_indexes(locij,ntype,tmpij0)
          istart = 1
          DO j = 1, ndim
             ! Special cases for EAM and IPBV
             locij = my_index(j)
             CALL get_indexes(locij,ntype,tmpij)
             SELECT CASE(pot_target)
             !NB should do something about QUIP here?
             CASE(ea_type,ip_type)
                ! check the array components
                CALL compare_pot(potparm%pot(tmpij(1) ,tmpij(2) )%pot,&
                                 potparm%pot(tmpij0(1),tmpij0(2))%pot,&
                                 check)
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
             Cwork(:) = pot_par(j,:)
             nunique = nunique + 1
             iend    = j - 1
             CALL set_potparm_index(potparm, my_index(istart:iend), pot_target, &
                  ntype, tmpij, atomic_kind_set, shift_cutoff, do_zbl)
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
          CALL set_potparm_index(potparm, my_index(istart:iend), pot_target, &
               ntype, tmpij, atomic_kind_set, shift_cutoff, do_zbl)
          DO i = istart, iend
             locij = my_index(i)
             CALL get_indexes(locij,ntype,tmpij)
             tmp_index(tmpij(1),tmpij(2))=nunique
             tmp_index(tmpij(2),tmpij(1))=nunique
          END DO
          DEALLOCATE(Cwork)
          DEALLOCATE(pot_par)
       ELSE
          nunique = nunique + 1
          CALL set_potparm_index(potparm, my_index, pot_target, ntype, tmpij,&
               atomic_kind_set, shift_cutoff, do_zbl)
       END IF
       DEALLOCATE(my_index)
    END DO
    ! Multiple defined potential
    n = 0
    DO i = 1, ntype
       DO j = 1, i
          n  = n + 1
          IF (SIZE(potparm%pot(i,j)%pot%type)==1) CYCLE
          nunique = nunique + 1
          tmp_index(i,j) = nunique
          tmp_index(j,i) = nunique
          !
          CALL set_potparm_index(potparm, (/n/), multi_type, ntype, tmpij, &
               atomic_kind_set, shift_cutoff, do_zbl)
       END DO
    END DO
    ! Concluding the postprocess..
    CALL spline_env_create(spline_env, ntype, nunique)
    spline_env%spltab=tmp_index
    DEALLOCATE(tmp_index)
    CALL timestop(handle)
  END SUBROUTINE get_nonbond_storage

! *****************************************************************************
!> \brief Trivial for non LJ potential.. gives back in the case of LJ
!>      the potparm with the smallest sigma..
!> \param potparm ...
!> \param my_index ...
!> \param pot_target ...
!> \param ntype ...
!> \param tmpij_out ...
!> \param atomic_kind_set ...
!> \param shift_cutoff ...
!> \param do_zbl ...
!> \author Teodoro Laino [tlaino] 2007.06
! *****************************************************************************
  SUBROUTINE set_potparm_index(potparm, my_index, pot_target, ntype, tmpij_out,&
       atomic_kind_set, shift_cutoff, do_zbl)

    TYPE(pair_potential_pp_type), POINTER    :: potparm
    INTEGER, INTENT(IN)                      :: my_index(:), pot_target, ntype
    INTEGER, INTENT(OUT)                     :: tmpij_out(2)
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    LOGICAL, INTENT(IN)                      :: shift_cutoff, do_zbl

    CHARACTER(len=*), PARAMETER :: routineN = 'set_potparm_index', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, min_val, nvalues, &
                                                tmpij(2), value, zi, zj
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: wrk
    LOGICAL                                  :: check, failure
    REAL(KIND=dp)                            :: hicut0, l_epsilon, l_sigma6, &
                                                m_epsilon, m_sigma6, &
                                                min_sigma6, rcovi, rcovj
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: sigma6
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(pair_potential_single_type), &
      POINTER                                :: pot, pot_ref

    CALL timeset(routineN,handle)

    NULLIFY(pot, pot_ref)
    failure = .FALSE.
    nvalues = SIZE(my_index)
    IF ((pot_target==lj_type).OR.(pot_target==lj_charmm_type)) THEN
       ALLOCATE(sigma6(nvalues))
       ALLOCATE(wrk(nvalues))
       min_sigma6= HUGE(0.0_dp)
       m_epsilon =-HUGE(0.0_dp)
       DO i = 1, nvalues
          value = my_index(i)
          CALL get_indexes(value,ntype,tmpij)
          pot => potparm%pot(tmpij(1),tmpij(2))%pot
          ! Preliminary check..
          check = SIZE(pot%type)==1
          CPASSERT(check)

          sigma6(i) = pot%set(1)%lj%sigma6
          l_epsilon = pot%set(1)%lj%epsilon
          IF (sigma6(i) /= 0.0_dp) min_sigma6 =  MIN(min_sigma6, sigma6(i))
          IF (sigma6(i) == 0.0_dp) sigma6(i)  = -HUGE(0.0_dp)
          IF (l_epsilon /= 0.0_dp) m_epsilon  =  MAX(m_epsilon, l_epsilon)
       END DO
       CALL sort(sigma6,nvalues,wrk)
       min_val  = my_index(wrk(nvalues))
       m_sigma6 = sigma6(nvalues)
       ! In case there are only zeros.. let's consider them properly..
       IF (m_sigma6   == -HUGE(0.0_dp)) m_sigma6   = 1.0_dp
       IF (m_epsilon  == -HUGE(0.0_dp)) m_epsilon  = 0.0_dp
       IF (min_sigma6 ==  HUGE(0.0_dp)) min_sigma6 = 0.0_dp
       DEALLOCATE(sigma6)
       DEALLOCATE(wrk)
    ELSE
       min_val = MINVAL(my_index(:))
    END IF
    CALL get_indexes(min_val,ntype,tmpij)
    tmpij_out = tmpij
    pot => potparm%pot(tmpij(1),tmpij(2))%pot
    pot%undef = .TRUE.
    IF (shift_cutoff) THEN
       hicut0 = SQRT(pot%rcutsq)
       IF (ABS(hicut0) <= MIN_HICUT_VALUE) hicut0 = DEFAULT_HICUT_VALUE
    END IF
    CALL init_genpot(potparm, ntype)

    DO i = 1, nvalues
          value = my_index(i)
          CALL get_indexes(value,ntype,tmpij)
          pot => potparm%pot(tmpij(1),tmpij(2))%pot
          CALL spline_factor_create(pot%spl_f)
          pot%spl_f%rcutsq_f  = 1.0_dp
          pot%spl_f%rscale    = 1.0_dp
          pot%spl_f%fscale    = 1.0_dp
    ENDDO

    IF(ANY(potential_single_allocation==pot_target)) THEN
       DO i = 1, nvalues
          value = my_index(i)
          CALL get_indexes(value,ntype,tmpij)
          pot => potparm%pot(tmpij(1),tmpij(2))%pot

          check = SIZE(pot%type)==1
          CPASSERT(check)
          ! Undef potential.. this will be used to compute the splines..
          IF ((pot_target==lj_type).OR.(pot_target==lj_charmm_type)) THEN
             l_sigma6 = pot%set(1)%lj%sigma6
             l_epsilon = pot%set(1)%lj%epsilon
             ! Undef potential.. this will be used to compute the splines..
             IF (pot%undef) THEN
                pot%set(1)%lj%sigma6  = m_sigma6
                pot%set(1)%lj%sigma12 = m_sigma6**2
                pot%set(1)%lj%epsilon = m_epsilon
             END IF
             pot%spl_f%rscale(1) = 1.0_dp
             pot%spl_f%fscale(1) = 0.0_dp
             IF (l_sigma6*l_epsilon/=0.0_dp) THEN
                pot%spl_f%rcutsq_f  = (min_sigma6/m_sigma6)**(1.0_dp/3.0_dp)
                pot%spl_f%rscale(1) = (l_sigma6/m_sigma6)**(1.0_dp/3.0_dp)
                pot%spl_f%fscale(1) = l_epsilon/m_epsilon
             END IF
          END IF
       ENDDO
    END IF

    DO i = 1, nvalues
       value = my_index(i)
       CALL get_indexes(value,ntype,tmpij)
       pot => potparm%pot(tmpij(1),tmpij(2))%pot

       IF (do_zbl) THEN
          atomic_kind => atomic_kind_set(tmpij(1))
          CALL get_atomic_kind(atomic_kind, rcov=rcovi, z=zi)
          atomic_kind => atomic_kind_set(tmpij(2))
          CALL get_atomic_kind(atomic_kind, rcov=rcovj, z=zj)
          CALL zbl_matching_polinomial(pot, rcovi, rcovj, REAL(zi,KIND=dp), &
               REAL(zj,KIND=dp))
       END IF
       ! Derivative factors
       pot%spl_f%dscale  = pot%spl_f%fscale  / pot%spl_f%rscale
       ! Cutoff for the potentials on splines
       IF (shift_cutoff) THEN
          ! Cutoff NonBonded
          pot%spl_f%cutoff = ener_pot(pot, hicut0, 0.0_dp)
       END IF
    END DO

    ! Handle the cutoff
    IF (shift_cutoff) THEN
       pot_ref => potparm%pot(tmpij_out(1),tmpij_out(2))%pot
       DO i = 1, nvalues
          value = my_index(i)
          CALL get_indexes(value,ntype,tmpij)
          pot => potparm%pot(tmpij(1),tmpij(2))%pot
          IF (value==min_val) CYCLE
          ! Cutoff NonBonded
          pot%spl_f%cutoff = pot_ref%spl_f%cutoff*pot%spl_f%fscale(1) - pot%spl_f%cutoff
       END DO
    END IF
    CALL finalizef()

    CALL timestop(handle)

  END SUBROUTINE set_potparm_index

! *****************************************************************************
!> \brief Gives back the indices of the matrix w.r.t. the collective array index
!> \param Inind ...
!> \param ndim ...
!> \param ij ...
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

