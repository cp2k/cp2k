!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
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
  USE ewald_environment_types,         ONLY: ewald_env_create,&
                                             ewald_env_get,&
                                             ewald_env_set,&
                                             ewald_environment_type,&
                                             read_ewald_section
  USE ewald_pw_types,                  ONLY: ewald_pw_create,&
                                             ewald_pw_type
  USE ewald_spline_util,               ONLY: Setup_Ewald_Spline
  USE input_constants,                 ONLY: do_qmmm_coulomb,&
                                             do_qmmm_gauss,&
                                             do_qmmm_pcharge,&
                                             do_qmmm_swave
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi
  USE pw_poisson_types,                ONLY: do_ewald_ewald,&
                                             do_ewald_none,&
                                             do_ewald_pme,&
                                             do_ewald_spme
  USE qmmm_gaussian_types,             ONLY: qmmm_gaussian_p_type,&
                                             qmmm_gaussian_type
  USE qmmm_types_low,                  ONLY: qmmm_per_pot_p_type,&
                                             qmmm_per_pot_type,&
                                             qmmm_pot_p_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qmmm_per_elpot'
  PUBLIC :: qmmm_per_potential_init, qmmm_ewald_potential_init

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
!> \par History
!>      09.2004 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE qmmm_per_potential_init(qmmm_coupl_type, per_potentials, potentials,&
    pgfs, qm_cell_small, mm_cell, para_env, compatibility, qmmm_periodic, print_section,&
    eps_mm_rspace, maxchrg, ncp, ncpl)
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

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_per_potential_init', &
      routineP = moduleN//':'//routineN

    INTEGER :: I, idim, ig, ig_start, iw, ix, iy, iz, K, Kmax(3), &
      n_rep_real(3), n_rep_real_val, ncoarsel, ncoarset, Ndim, output_unit, &
      stat
    INTEGER, DIMENSION(:), POINTER           :: mm_atom_index
    REAL(KIND=dp) :: Ak, alpha, box(3), Fac(3), fs, g, g2, Gk, Gmax, &
      mymaxradius, npl, npt, Prefactor, rc, rc2, Rmax, tmp, vec(3), vol
    REAL(KIND=dp), DIMENSION(:), POINTER     :: gx, gy, gz, Lg
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(qmmm_gaussian_type), POINTER        :: pgf

    NULLIFY(Lg, gx, gy, gz)
    ncoarset = PRODUCT(ncp)
    ncoarsel = PRODUCT(ncpl)
    logger => cp_get_default_logger()
    output_unit= cp_logger_get_default_io_unit(logger)
    Rmax = SQRT(mm_cell%hmat(1,1)**2+&
                mm_cell%hmat(2,2)**2+&
                mm_cell%hmat(3,3)**2 )
    CALL section_vals_val_get(qmmm_periodic,"GMAX",r_val=Gmax)
    CALL section_vals_val_get(qmmm_periodic,"REPLICA",i_val=n_rep_real_val)
    fac = 2.0e0_dp*Pi/(/mm_cell%hmat(1,1),mm_cell%hmat(2,2),mm_cell%hmat(3,3)/)
    Kmax    = CEILING(Gmax/Fac)
    Vol     = mm_cell%hmat(1,1)*&
              mm_cell%hmat(2,2)*&
              mm_cell%hmat(3,3)
    Ndim    = (Kmax(1)+1)*(2*Kmax(2)+1)*(2*Kmax(3)+1)
    ig_start = 1
    n_rep_real = n_rep_real_val
    IF (compatibility.AND.(qmmm_coupl_type==do_qmmm_gauss)) ig_start = 2

    CPASSERT(.NOT.ASSOCIATED(per_potentials))
    ALLOCATE(per_potentials(SIZE(pgfs)))
    CPASSERT(SIZE(pgfs)==SIZE(potentials))
    Potential_Type: DO K = 1, SIZE(pgfs)

       rc = pgfs(K)%pgf%Elp_Radius
       ALLOCATE(per_potentials(K)%Pot)
       SELECT CASE(qmmm_coupl_type)
       CASE(do_qmmm_coulomb,do_qmmm_pcharge)
          ! Not yet implemented for this case
          CPASSERT(.FALSE.)
       CASE(do_qmmm_gauss,do_qmmm_swave)
          ALLOCATE(Lg(Ndim))
          ALLOCATE(gx(Ndim))
          ALLOCATE(gy(Ndim))
          ALLOCATE(gz(Ndim))
       END SELECT

       LG = 0.0_dp
       gx = 0.0_dp
       gy = 0.0_dp
       gz = 0.0_dp

       SELECT CASE(qmmm_coupl_type)
       CASE(do_qmmm_coulomb,do_qmmm_pcharge)
          ! Not yet implemented for this case
          CPASSERT(.FALSE.)
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
                CPASSERT(.FALSE.)
             END IF
             n_rep_real(1) = CEILING(( box(1) + mymaxradius ) / mm_cell%hmat(1,1))
             n_rep_real(2) = CEILING(( box(2) + mymaxradius ) / mm_cell%hmat(2,2))
             n_rep_real(3) = CEILING(( box(3) + mymaxradius ) / mm_cell%hmat(3,3))
          END IF

       CASE DEFAULT
          DEALLOCATE(per_potentials(K)%Pot)
          CPASSERT(stat==0)
          NULLIFY(per_potentials(K)%Pot)
          IF (output_unit>0) WRITE(output_unit,'(A)')" QMMM Periodic Potential - not Initialized!"
          CYCLE Potential_Type
       END SELECT

       NULLIFY(mm_atom_index)
       ALLOCATE(mm_atom_index(SIZE(potentials(K)%pot%mm_atom_index)))
       mm_atom_index = potentials(K)%pot%mm_atom_index

       NULLIFY(per_potentials(K)%Pot%LG, per_potentials(K)%Pot%mm_atom_index,&
               per_potentials(K)%Pot%gx, per_potentials(K)%Pot%gy, per_potentials(K)%Pot%gz)
       CALL qmmm_per_pot_type_create(per_potentials(K)%Pot, LG=LG, gx=gx, gy=gy, gz=gz,&
                                     Gmax=Gmax, Kmax=Kmax, n_rep_real=n_rep_real,&
                                     Fac=Fac, mm_atom_index=mm_atom_index, &
                                     mm_cell=mm_cell, para_env=para_env,&
                                     qmmm_per_section=qmmm_periodic, print_section=print_section)

       iw=cp_print_key_unit_nr(logger,print_section,"PERIODIC_INFO",&
                                         extension=".log")
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
               "PERIODIC_INFO")

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
!> \par History
!>      09.2004 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE qmmm_per_pot_type_create(Pot, LG, gx, gy, gz, GMax, Kmax, n_rep_real,&
       Fac, mm_atom_index, mm_cell, para_env, qmmm_per_section, print_section)
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

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_per_pot_type_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: npts(3)
    INTEGER, DIMENSION(:), POINTER           :: ngrids
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: hmat(3,3)
    TYPE(section_vals_type), POINTER         :: grid_print_section

    failure = .FALSE.
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
    CALL section_vals_val_get(qmmm_per_section,"ngrids", i_vals=ngrids)
    npts =  ngrids
    hmat =  mm_cell%hmat

    grid_print_section=> section_vals_get_subs_vals(print_section,"GRID_INFORMATION")
    CALL Setup_Ewald_Spline(pw_grid=Pot%pw_grid, pw_pool=Pot%pw_pool, coeff=Pot%TabLR,&
         LG=LG, gx=gx, gy=gy, gz=gz, hmat=hmat, npts=npts, param_section=qmmm_per_section,&
         tag="qmmm",para_env=para_env, print_section=grid_print_section)


  END SUBROUTINE qmmm_per_pot_type_create

! *****************************************************************************
!> \brief Initialize the QMMM Ewald potential needed for QM-QM Coupling using
!>      point charges
!> \param ewald_env ...
!> \param ewald_pw ...
!> \param qmmm_coupl_type ...
!> \param mm_cell ...
!> \param para_env ...
!> \param qmmm_periodic ...
!> \param print_section ...
!> \par History
!>      10.2014 created [JGH]
!> \author JGH
! *****************************************************************************
  SUBROUTINE qmmm_ewald_potential_init(ewald_env,ewald_pw,qmmm_coupl_type, mm_cell, para_env,&
                                       qmmm_periodic, print_section)
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(ewald_pw_type), POINTER             :: ewald_pw
    INTEGER, INTENT(IN)                      :: qmmm_coupl_type
    TYPE(cell_type), POINTER                 :: mm_cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: qmmm_periodic, print_section

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_ewald_potential_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ewald_type, gmax(3), iw, &
                                                o_spline, ounit
    LOGICAL                                  :: do_multipoles, failure
    REAL(KIND=dp)                            :: alpha, rcut
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: ewald_print_section, &
                                                ewald_section, poisson_section

    failure = .FALSE.
    logger => cp_get_default_logger()
    ounit = cp_logger_get_default_io_unit(logger)
    CPASSERT(.NOT.ASSOCIATED(ewald_env))
    CPASSERT(.NOT.ASSOCIATED(ewald_pw))

    ! Create Ewald environments
    poisson_section => section_vals_get_subs_vals(qmmm_periodic,"POISSON")
    CALL ewald_env_create(ewald_env,para_env)
    CALL ewald_env_set(ewald_env,poisson_section=poisson_section)
    ewald_section => section_vals_get_subs_vals(poisson_section,"EWALD")
    CALL read_ewald_section(ewald_env,ewald_section)
    ewald_print_section => section_vals_get_subs_vals(print_section,"GRID_INFORMATION")
    CALL ewald_pw_create(ewald_pw,ewald_env,mm_cell,mm_cell,print_section=ewald_print_section)

    CALL ewald_env_get(ewald_env,ewald_type=ewald_type,do_multipoles=do_multipoles,&
                       gmax=gmax,o_spline=o_spline,alpha=alpha,rcut=rcut)
    IF(do_multipoles) CALL cp_unimplemented_error(fromWhere=routineP,&
       message="No multipole force fields allowed in QM-QM Ewald long range correction",&
       error_level=cp_fatal_level)

    SELECT CASE(qmmm_coupl_type)
       CASE(do_qmmm_coulomb)
          CALL cp_unimplemented_error(fromWhere=routineP, &
            message="QM-QM long range correction not possible with COULOMB coupling", &
            error_level=cp_fatal_level)
       CASE(do_qmmm_pcharge)
          ! OK
       CASE(do_qmmm_gauss,do_qmmm_swave)
          CALL cp_unimplemented_error(fromWhere=routineP, &
            message="QM-QM long range correction not possible with GAUSS/SWAVE coupling", &
            error_level=cp_fatal_level)
       CASE DEFAULT
          ! We should never get to this point
          CPASSERT(.FALSE.)
    END SELECT

    iw=cp_print_key_unit_nr(logger,print_section,"PERIODIC_INFO",extension=".log")
    IF ( iw > 0) THEN
       WRITE (UNIT=iw,FMT="(/,T2,A)") REPEAT("-",79)
       WRITE (UNIT=iw,FMT="(T2,A,T20,A,T80,A)")"-","QMMM PERIODIC BOUNDARY CONDITION INFO","-"
       WRITE (UNIT=iw,FMT="(T2,A,T25,A,T80,A)")"-","QM-QM Long Range Correction","-"
       WRITE (UNIT=iw,FMT="(T2,A)") REPEAT("-",79)
       SELECT CASE(ewald_type)
          CASE(do_ewald_none)
             CALL cp_unimplemented_error(fromWhere=routineP, &
               message="QM-QM long range correction not compatible with Ewald=NONE", &
               error_level=cp_fatal_level)
          CASE(do_ewald_pme)
             CALL cp_unimplemented_error(fromWhere=routineP, &
               message="QM-QM long range correction not possible with Ewald=PME", &
               error_level=cp_fatal_level)
          CASE(do_ewald_ewald)
             CALL cp_unimplemented_error(fromWhere=routineP, &
               message="QM-QM long range correction not possible with Ewald method", &
               error_level=cp_fatal_level)
          CASE(do_ewald_spme)
             WRITE (UNIT=iw,FMT="(T2,A,T35,A,T75,A,T80,A)")"-","Ewald type","SPME","-"
             WRITE (UNIT=iw,FMT="(T2,A,T35,A,T61,3I6,T80,A)")"-","GMAX values",gmax,"-"
             WRITE (UNIT=iw,FMT="(T2,A,T35,A,T73,I6,T80,A)")"-","Spline order",o_spline,"-"
             WRITE (UNIT=iw,FMT="(T2,A,T35,A,T67,F12.4,T80,A)")"-","Alpha value",alpha,"-"
             WRITE (UNIT=iw,FMT="(T2,A,T35,A,T67,F12.4,T80,A)")"-","Real space cutoff value",rcut,"-"
       END SELECT
       WRITE (UNIT=iw,FMT="(T2,A)") REPEAT("-",79)
    END IF
    CALL cp_print_key_finished_output(iw,logger,print_section,"PERIODIC_INFO")

  END SUBROUTINE qmmm_ewald_potential_init

! *****************************************************************************

END MODULE qmmm_per_elpot

