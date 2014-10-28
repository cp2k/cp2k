!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
MODULE qmmm_per_elpot

! *****************************************************************************
!> \brief Setting up the potential for QM/MM periodic boundary conditions calculations
!> \par History
!>      07.2005 created [tlaino]
!> \author Teodoro Laino
  USE ao_util,                         ONLY: exp_radius
  USE cell_types,                      ONLY: cell_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE ewald_spline_util,               ONLY: Setup_Ewald_Spline
  USE f77_blas
  USE input_constants,                 ONLY: do_qmmm_coulomb,&
                                             do_qmmm_gauss,&
                                             do_qmmm_pcharge,&
                                             do_qmmm_swave
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi
  USE qmmm_gaussian_types,             ONLY: qmmm_gaussian_p_type,&
                                             qmmm_gaussian_type
  USE qmmm_types,                      ONLY: qmmm_per_pot_p_type,&
                                             qmmm_per_pot_type,&
                                             qmmm_pot_p_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qmmm_per_elpot'
  PUBLIC :: qmmm_per_potential_init

CONTAINS

! *****************************************************************************
!> \brief Initialize the QMMM potential stored on vector,
!>      according the qmmm_coupl_type
!> \param qmmm_coupl_type ...
!> \param per_potentials ...
!> \param potentials ...
!> \param pgfs ...
!> \param qm_cell_small ...
!> \param mm_cell ...
!> \param para_env ...
!> \param compatibility ...
!> \param qmmm_periodic ...
!> \param print_section ...
!> \param eps_mm_rspace ...
!> \param maxchrg ...
!> \param ncp ...
!> \param ncpl ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      09.2004 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE qmmm_per_potential_init(qmmm_coupl_type, per_potentials, potentials,&
    pgfs, qm_cell_small, mm_cell, para_env, compatibility, qmmm_periodic, print_section,&
    eps_mm_rspace, maxchrg, ncp, ncpl, error)
    INTEGER, INTENT(IN)                      :: qmmm_coupl_type
    TYPE(qmmm_per_pot_p_type), &
      DIMENSION(:), POINTER                  :: per_potentials
    TYPE(qmmm_pot_p_type), DIMENSION(:), &
      POINTER                                :: potentials
    TYPE(qmmm_gaussian_p_type), &
      DIMENSION(:), POINTER                  :: pgfs
    TYPE(cell_type), POINTER                 :: qm_cell_small, mm_cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    LOGICAL, INTENT(IN)                      :: compatibility
    TYPE(section_vals_type), POINTER         :: qmmm_periodic, print_section
    REAL(KIND=dp), INTENT(IN)                :: eps_mm_rspace, maxchrg
    INTEGER, INTENT(IN)                      :: ncp(3), ncpl(3)
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_per_potential_init', &
      routineP = moduleN//':'//routineN

    INTEGER :: I, idim, ig, ig_start, iw, ix, iy, iz, K, Kmax(3), &
      n_rep_real(3), n_rep_real_val, ncoarsel, ncoarset, Ndim, output_unit, &
      stat
    INTEGER, DIMENSION(:), POINTER           :: mm_atom_index
    LOGICAL                                  :: failure
    REAL(KIND=dp) :: Ak, alpha, box(3), Fac(3), fs, g, g2, Gk, Gmax, &
      mymaxradius, npl, npt, Prefactor, rc, rc2, Rmax, tmp, vec(3), vol
    REAL(KIND=dp), DIMENSION(:), POINTER     :: gx, gy, gz, Lg
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(qmmm_gaussian_type), POINTER        :: pgf

    NULLIFY(Lg, gx, gy, gz)
    ncoarset = PRODUCT(ncp)
    ncoarsel = PRODUCT(ncpl)
    logger => cp_error_get_logger(error)
    output_unit= cp_logger_get_default_io_unit(logger)
    Rmax = SQRT(mm_cell%hmat(1,1)**2+&
                mm_cell%hmat(2,2)**2+&
                mm_cell%hmat(3,3)**2 )
    CALL section_vals_val_get(qmmm_periodic,"GMAX",r_val=Gmax,error=error)
    CALL section_vals_val_get(qmmm_periodic,"REPLICA",i_val=n_rep_real_val,error=error)
    fac = 2.0e0_dp*Pi/(/mm_cell%hmat(1,1),mm_cell%hmat(2,2),mm_cell%hmat(3,3)/)
    Kmax    = CEILING(Gmax/Fac)
    Vol     = mm_cell%hmat(1,1)*&
              mm_cell%hmat(2,2)*&
              mm_cell%hmat(3,3)
    Ndim    = (Kmax(1)+1)*(2*Kmax(2)+1)*(2*Kmax(3)+1)
    ig_start = 1
    n_rep_real = n_rep_real_val
    IF (compatibility.AND.(qmmm_coupl_type==do_qmmm_gauss)) ig_start = 2

    CPPrecondition(.NOT.ASSOCIATED(per_potentials),cp_failure_level,routineP,error,failure)
    ALLOCATE(per_potentials(SIZE(pgfs)),stat=stat)
    CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
    CPPrecondition(SIZE(pgfs)==SIZE(potentials),cp_failure_level,routineP,error,failure)
    Potential_Type: DO K = 1, SIZE(pgfs)

       rc = pgfs(K)%pgf%Elp_Radius
       ALLOCATE(per_potentials(K)%Pot, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       SELECT CASE(qmmm_coupl_type)
       CASE(do_qmmm_coulomb,do_qmmm_pcharge)
          ! Not yet implemented for this case
          CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
       CASE(do_qmmm_gauss,do_qmmm_swave)
          ALLOCATE(Lg(Ndim), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(gx(Ndim), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(gy(Ndim), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(gz(Ndim), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END SELECT

       LG = 0.0_dp
       gx = 0.0_dp
       gy = 0.0_dp
       gz = 0.0_dp

       SELECT CASE(qmmm_coupl_type)
       CASE(do_qmmm_coulomb,do_qmmm_pcharge)
          ! Not yet implemented for this case
          CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
       CASE(do_qmmm_gauss,do_qmmm_swave)
          pgf => pgfs(K)%pgf
          idim = 0
          DO ix = 0, kmax(1)
             DO iy = -kmax(2), kmax(2)
                DO iz = -kmax(3), kmax(3)
                   idim = idim + 1
                   IF (ix == 0.AND.iy == 0.AND.iz == 0) THEN
                      DO Ig = ig_start, pgf%number_of_gaussians
                         Gk = pgf%Gk(Ig)
                         Ak = pgf%Ak(Ig)*Pi**(3.0_dp/2.0_dp)*Gk**3.0_dp
                         LG (idim) = LG (idim) - Ak
                      END DO
                   ELSE
                      fs = 2.0_dp; IF (ix==0) fs = 1.0_dp
                      vec =  fac*(/REAL(ix,KIND=dp),REAL(iy,KIND=dp),REAL(iz,KIND=dp)/)
                      g2  =  DOT_PRODUCT(vec,vec)
                      rc2 =  rc*rc
                      g   =  SQRT(g2)
                      IF     (qmmm_coupl_type==do_qmmm_gauss) THEN
                         LG (idim) = 4.0_dp*Pi/g2 * EXP(-(g2*rc2)/4.0_dp)
                      ELSEIF (qmmm_coupl_type==do_qmmm_swave) THEN
                         tmp = 4.0_dp/rc2
                         LG (idim) = 4.0_dp*Pi*tmp**2 / (g2*(g2+tmp)**2)
                      END IF
                      DO Ig = ig_start, pgf%number_of_gaussians
                         Gk = pgf%Gk(Ig)
                         Ak = pgf%Ak(Ig)*Pi**(3.0_dp/2.0_dp)*Gk**3.0_dp
                         LG (idim) = LG (idim) - Ak * EXP(-(g*Gk)**2.0_dp/4.0_dp)
                      END DO
                   ENDIF
                   LG (idim) = fs * LG (idim) * 1.0_dp/Vol
                   gx (idim) = fac(1)*REAL(ix,KIND=dp)
                   gy (idim) = fac(2)*REAL(iy,KIND=dp)
                   gz (idim) = fac(3)*REAL(iz,KIND=dp)
                END DO
             END DO
          END DO

          IF (ALL(n_rep_real == -1)) THEN
             mymaxradius = 0.0_dp
             DO I = 1, pgf%number_of_gaussians
                IF ( pgf%Gk(I) /= 0.0_dp) THEN
                   alpha = 1.0_dp / pgf%Gk(I)
                   alpha = alpha*alpha
                   Prefactor = pgf%Ak(I)*maxchrg
                   mymaxradius=MAX(mymaxradius,exp_radius( 0, alpha, eps_mm_rspace, Prefactor))
                END IF
             END DO
             box (1) = (qm_cell_small%hmat(1,1) - mm_cell%hmat(1,1))/2.0_dp
             box (2) = (qm_cell_small%hmat(2,2) - mm_cell%hmat(2,2))/2.0_dp
             box (3) = (qm_cell_small%hmat(3,3) - mm_cell%hmat(3,3))/2.0_dp
             IF (ANY(box > 0.0_dp)) THEN
                CPPostcondition(.FALSE.,cp_fatal_level,routineP,error,failure)
             END IF
             n_rep_real(1) = CEILING(( box(1) + mymaxradius ) / mm_cell%hmat(1,1))
             n_rep_real(2) = CEILING(( box(2) + mymaxradius ) / mm_cell%hmat(2,2))
             n_rep_real(3) = CEILING(( box(3) + mymaxradius ) / mm_cell%hmat(3,3))
          END IF

       CASE DEFAULT
          DEALLOCATE(per_potentials(K)%Pot)
          CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
          NULLIFY(per_potentials(K)%Pot)
          IF (output_unit>0) WRITE(output_unit,'(A)')" QMMM Periodic Potential - not Initialized!"
          CYCLE Potential_Type
       END SELECT

       NULLIFY(mm_atom_index)
       ALLOCATE(mm_atom_index(SIZE(potentials(K)%pot%mm_atom_index)),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,Failure)
       mm_atom_index = potentials(K)%pot%mm_atom_index

       NULLIFY(per_potentials(K)%Pot%LG, per_potentials(K)%Pot%mm_atom_index,&
               per_potentials(K)%Pot%gx, per_potentials(K)%Pot%gy, per_potentials(K)%Pot%gz)
       CALL qmmm_per_pot_type_create(per_potentials(K)%Pot, LG=LG, gx=gx, gy=gy, gz=gz,&
                                     Gmax=Gmax, Kmax=Kmax, n_rep_real=n_rep_real,&
                                     Fac=Fac, mm_atom_index=mm_atom_index, &
                                     mm_cell=mm_cell, para_env=para_env,&
                                     qmmm_per_section=qmmm_periodic, print_section=print_section,&
                                     error=error)

       iw=cp_print_key_unit_nr(logger,print_section,"PERIODIC_INFO",&
                                         extension=".log",error=error)
       IF ( iw > 0) THEN
             npt = REAL(ncoarset,KIND=dp)*REAL(ndim,KIND=dp)*REAL(SIZE(mm_atom_index),KIND=dp)
             npl = REAL(ncoarsel,KIND=dp)*REAL(ndim,KIND=dp)*REAL(SIZE(mm_atom_index),KIND=dp)
             WRITE (UNIT=iw,FMT="(/,T2,A)") REPEAT("-",79)
             WRITE (UNIT=iw,FMT="(T2,A,T20,A,T80,A)")"-","QMMM PERIODIC BOUNDARY CONDITION INFO","-"
             WRITE (UNIT=iw,FMT="(T2,A)") REPEAT("-",79)
             WRITE (UNIT=iw,FMT="(T2,A,T10,A,F15.6,T50,A,3I5,T80,A)")"-","RADIUS  =",rc,"REPLICA =",n_rep_real,"-"
             WRITE (UNIT=iw,FMT="(T2,A,T10,A,F15.6,T50,A,I15,T80,A)")"-","MINGVAL =",MINVAL(ABS(Lg)),&
                  "GPOINTS =",ndim,"-"
             WRITE (UNIT=iw,FMT="(T2,A,T10,A,3I5,T50,A,3I5,T80,A)")  "-","NCOARSL =",ncpl,&
                  "NCOARST =",ncp,"-"
             WRITE (UNIT=iw,FMT="(T2,A,T10,A,F15.0,T50,A,F15.0,T80,A)")  "-","NFLOP-L ~",npl,&
                  "NFLOP-T ~",npt,"-"
             WRITE (UNIT=iw,FMT="(T2,A)") REPEAT("-",79)
       END IF
       CALL cp_print_key_finished_output(iw,logger,print_section,&
               "PERIODIC_INFO", error=error)

    END DO Potential_Type

  END SUBROUTINE qmmm_per_potential_init

! *****************************************************************************
!> \brief Creates the qmmm_pot_type structure
!> \param Pot ...
!> \param LG ...
!> \param gx ...
!> \param gy ...
!> \param gz ...
!> \param GMax ...
!> \param Kmax ...
!> \param n_rep_real ...
!> \param Fac ...
!> \param mm_atom_index ...
!> \param mm_cell ...
!> \param para_env ...
!> \param qmmm_per_section ...
!> \param print_section ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      09.2004 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE qmmm_per_pot_type_create(Pot, LG, gx, gy, gz, GMax, Kmax, n_rep_real,&
       Fac, mm_atom_index, mm_cell, para_env, qmmm_per_section, print_section,&
       error)
    TYPE(qmmm_per_pot_type), POINTER         :: Pot
    REAL(KIND=dp), DIMENSION(:), POINTER     :: LG, gx, gy, gz
    REAL(KIND=dp), INTENT(IN)                :: Gmax
    INTEGER, INTENT(IN)                      :: Kmax(3), n_rep_real(3)
    REAL(KIND=dp), INTENT(IN)                :: Fac(3)
    INTEGER, DIMENSION(:), POINTER           :: mm_atom_index
    TYPE(cell_type), POINTER                 :: mm_cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: qmmm_per_section, &
                                                print_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_per_pot_type_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: npts(3)
    INTEGER, DIMENSION(:), POINTER           :: ngrids
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: hmat(3,3)
    TYPE(section_vals_type), POINTER         :: grid_print_section

    failure = .FALSE.
    IF (.NOT.failure) THEN
       Pot%LG            => LG
       Pot%gx            => gx
       Pot%gy            => gy
       Pot%gz            => gz
       Pot%mm_atom_index => mm_atom_index
       Pot%Gmax          =  Gmax
       Pot%Kmax          =  Kmax
       Pot%n_rep_real    =  n_rep_real
       Pot%Fac           =  Fac
       !
       ! Setting Up Fit Procedure
       !
       NULLIFY(Pot%pw_grid)
       NULLIFY(Pot%pw_pool)
       NULLIFY(Pot%TabLR, ngrids)
       CALL section_vals_val_get(qmmm_per_section,"ngrids", i_vals=ngrids, error=error)
       npts =  ngrids
       hmat =  mm_cell%hmat

       grid_print_section=> section_vals_get_subs_vals(print_section,"GRID_INFORMATION",error=error)
       CALL Setup_Ewald_Spline(pw_grid=Pot%pw_grid, pw_pool=Pot%pw_pool, coeff=Pot%TabLR,&
            LG=LG, gx=gx, gy=gy, gz=gz, hmat=hmat, npts=npts, param_section=qmmm_per_section,&
            tag="qmmm",para_env=para_env, print_section=grid_print_section,error=error)

    END IF

  END SUBROUTINE qmmm_per_pot_type_create

END MODULE qmmm_per_elpot

