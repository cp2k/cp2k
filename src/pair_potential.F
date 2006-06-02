!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/pair_potential [1.0] *
!!
!!   NAME
!!     pair_potential
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     September 2005 - Introduced the Born-Mayer-Huggins-Fumi-Tosi  Potential
!!                      BMHTF
!!     2006 - Major rewriting of the routines.. Linear scaling setup of splines
!!
!!   SOURCE
!******************************************************************************

MODULE pair_potential

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE erf_fn,                          ONLY: erfc
  USE ewald_environment_types,         ONLY: ewald_env_get,&
                                             ewald_environment_type
  USE input_constants,                 ONLY: do_ewald_none
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi
  USE pair_potential_types,            ONLY: pair_potential_pp_type,&
                                             pair_potential_single_type,&
                                             nn_type,&
                                             lj_type,&
                                             ft_type,&
                                             wl_type,&
                                             gw_type,&
                                             ip_type,&
                                             ea_type,&
                                             compare_pot,&
                                             list_pot
  USE splines,                         ONLY: init_spline,&
                                             init_splinexy,&
                                             spline_data_type,&
                                             spline_environment_type,&
                                             potential_s,&
                                             spline_data_release,&
                                             spline_data_create,&
                                             spline_env_create
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: sort
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  REAL(KIND=dp), PARAMETER :: ifourpi = 1.0_dp / ( 4.0_dp * pi )
  PUBLIC :: potential_f, potential_s
  PUBLIC :: spline_nonbond_control, get_nonbond_storage
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pair_potential'

CONTAINS

  !******************************************************************************
  !!****
  !!****s* pair_potential/ener_pot [1.0] *
  !!
  !!   NAME
  !!     ener_pot
  !!
  !!   FUNCTION
  !!     Evaluates the nonbond potential energy for the implemented FF kinds
  !!     
  !!   AUTHOR
  !!    
  !!
  !!   MODIFICATION HISTORY
  !!     
  !!
  !!   SOURCE
  !******************************************************************************
  FUNCTION ener_pot ( pot, qi, qj, r, vdw_scale14 )
    TYPE(pair_potential_single_type), &
         POINTER                                :: pot
    REAL(KIND=dp), INTENT(IN)                :: qi, qj, r
    REAL(KIND=dp)                            :: vdw_scale14
    REAL(KIND=dp)                            :: ener_pot

    INTEGER                                  :: i, index
    REAL(KIND=dp)                            :: pp, qq, scale

    IF (pot%type==lj_type) THEN
       ener_pot = vdw_scale14 * &
            4.0_dp*pot%lj%epsilon*(pot%lj%sigma12*r**(-12)-pot%lj% &
            sigma6*r**(-6)) - pot%energy_cutoff
    ELSE IF (pot%type==ip_type) THEN
       ener_pot = 0._dp
       IF ( r > pot % ipbv % rcore ) THEN
          DO i = 2, 15
             ener_pot = ener_pot + pot % ipbv % a(i)/(r**(i-1)*REAL(i-1,dp)) 
          END DO
       ELSE
          ! use a linear potential  
          ener_pot = pot % ipbv % m * r + pot % ipbv % b
       ENDIF
       ener_pot=ener_pot - pot % energy_cutoff
    ELSE IF (pot%type==wl_type) THEN
       ener_pot = pot%willis%a*EXP(-pot%willis%b*r) - pot%willis%c/r**6 &
            - pot%energy_cutoff
    ELSE IF (pot%type==gw_type) THEN
       scale = EXP(pot%goodwin%m*(-(r/pot%goodwin%dc)**pot%goodwin%mc+(pot% &
            goodwin%d/pot%goodwin%dc)**pot%goodwin%mc))
       ener_pot = scale*pot%goodwin%vr0*(pot%goodwin%d/r)**pot%goodwin%m
    ELSE IF (pot%type==ft_type) THEN
       ener_pot = pot%ft%a*EXP(-pot%ft%b*r) - pot%ft%c/r**6 - pot%ft%d/r**8 &
            - pot%energy_cutoff
    ELSE IF (pot%type==ea_type) THEN
       index = INT ( r / pot % eam % drar ) + 1
       IF ( index > 5000 ) THEN
          index = 5000
       ELSEIF ( index < 1 ) THEN
          index = 1
       ENDIF
       qq = r - pot % eam % rval (index)
       pp =  pot % eam % phi ( index ) + &
            qq * pot % eam % phip ( index )
       ener_pot = pp - pot % energy_cutoff
    ELSE 
       ener_pot = 0.0_dp
    END IF

  END FUNCTION ener_pot

  !******************************************************************************
  !!****
  !!****s* pair_potential/ener_coul [1.0] *
  !!
  !!   NAME
  !!     ener_coul
  !!
  !!   FUNCTION
  !!     Evaluates the Coulomb electrostatic energy 
  !!
  !!   AUTHOR
  !!    
  !!
  !!   MODIFICATION HISTORY
  !!     
  !!
  !!   SOURCE
  !******************************************************************************
  FUNCTION ener_coul ( pot, qi, qj, r, ewald_env, ei_scale14, my_do_14 )

    TYPE(pair_potential_single_type), &
         POINTER                                :: pot
    REAL(KIND=dp), INTENT(IN)                :: qi, qj, r
    TYPE(ewald_environment_type),POINTER     :: ewald_env
    REAL(KIND=dp)                            :: ei_scale14
    REAL(KIND=dp)                            :: ener_coul
    LOGICAL, INTENT(IN)                      :: my_do_14

    INTEGER                                  :: ewald_type
    REAL(KIND=dp)                            :: alpha, eps0

    ener_coul= 0.0_dp
    CALL ewald_env_get ( ewald_env, alpha = alpha, eps0 = eps0, &
         ewald_type = ewald_type ) 
    IF (.NOT. my_do_14 ) THEN
       IF (ewald_type /= do_ewald_none) THEN
          ener_coul = ei_scale14 * ifourpi / eps0 * qi * qj &
               * erfc ( alpha * r ) / r - pot % e_cutoff_coul
       ELSE
          CALL ewald_env_get ( ewald_env, eps0 = eps0)
          ener_coul = ei_scale14 * ifourpi / eps0 * qi * qj / r - pot % e_cutoff_coul
       END IF
    ELSE
       ! to avoid a large numberof points in the spline generation
       ! we evaluate the 1-4 interactions analytically..
       ener_coul = 0.0_dp
       !ei_scale14 * ifourpi / eps0 * qi * qj / r - pot % e_cutoff_coul
    ENDIF

  END FUNCTION ener_coul

  !******************************************************************************
  !!****
  !!****s* pair_potential/spline_nonbond_control [1.0] *
  !!
  !!   NAME
  !!     spline_nonbond_control
  !!
  !!   FUNCTION
  !!     creates the splines for the potentials
  !!
  !!   AUTHOR
  !!    
  !!
  !!   MODIFICATION HISTORY
  !!     Teo 2006.05 : Improved speed and accuracy. Linear scaling of the setup
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE spline_nonbond_control ( spline_env, potparm, atomic_kind_set, &
       ewald_env, ei_scale14, vdw_scale14, eps_spline, max_energy, rlow_nb,&
       emax_spline, iw, iw2, iw3, do_14, scale_cutoff, error)

    TYPE(spline_environment_type), POINTER   :: spline_env
    TYPE(pair_potential_pp_type), POINTER    :: potparm
    TYPE(atomic_kind_type), DIMENSION(:), &
         POINTER                             :: atomic_kind_set
    TYPE(ewald_environment_type), OPTIONAL, &
         POINTER                             :: ewald_env
    REAL(KIND=dp), OPTIONAL                  :: ei_scale14, vdw_scale14
    REAL(KIND=dp), INTENT(IN)                :: eps_spline, max_energy, rlow_nb,&
         emax_spline
    INTEGER, INTENT(IN)                      :: iw, iw2, iw3
    LOGICAL, INTENT(IN),OPTIONAL             :: do_14
    LOGICAL, INTENT(IN)                      :: scale_cutoff
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'spline_nonbond_control', &
         routineP = moduleN//':'//routineN
    INTEGER                                  :: i, j, jx, n, npoints, ntype, handle
    INTEGER                                  :: nsize, nkx, kx, ncount
    REAL(KIND=dp)                            :: dx2, e, &
         hicut, locut, qi, qj, x, x2,&
         xdum,diffmax, xdum1, xsav, ediff
    REAL(KIND=dp):: my_ei_scale14, my_vdw_scale14, eps0, myfac1, myfac2
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    LOGICAL                                  :: failure, check_bounds, my_do_14
    
    failure     = .FALSE.
    n      = 0
    ncount = 0
    my_ei_scale14 = 1.0_dp ; IF (PRESENT(ei_scale14))   my_ei_scale14 = ei_scale14
    my_vdw_scale14= 1.0_dp ; IF (PRESENT(vdw_scale14)) my_vdw_scale14 = vdw_scale14
    my_do_14= .FALSE.      ; IF (PRESENT(do_14))             my_do_14 = do_14
    ntype = SIZE ( atomic_kind_set )
    CALL timeset(routineN,'I','',handle)
    IF (iw3>0) THEN
       WRITE(iw3,'(T2,A,I0,A,I0,A)') 'SPLINE_INFO| Generating: ',&
            (ntype*(ntype+1))/2,' splines due to ',ntype,' different atomic kinds'
    ENDIF
    DO i = 1, ntype ! i:  first atom type
       atomic_kind=> atomic_kind_set(i)
       CALL get_atomic_kind(atomic_kind=atomic_kind,qeff=qi)
       DO j = 1, i  ! j:  second atom type
          IF (.NOT.potparm%pot(i,j)%pot%undef) CYCLE
          ncount = ncount + 1
          potparm%pot(i,j)%pot%undef = .FALSE.
          atomic_kind=> atomic_kind_set(j)
          CALL get_atomic_kind(atomic_kind=atomic_kind,qeff=qj)
          n = spline_env % spltab(i,j)
          locut = rlow_nb
          hicut = SQRT ( potparm %pot ( i, j ) %pot % rcutsq )
          IF (ABS(hicut) <= 1.0E-15_dp) hicut = 1.0E3_dp
          potparm%pot(i,j)%pot %energy_cutoff = 0.0_dp
          potparm%pot(i,j)%pot %e_cutoff_coul = 0.0_dp
          IF (scale_cutoff) THEN
             potparm%pot(i,j)%pot %energy_cutoff = ener_pot(potparm%pot(i,j)%pot,qi,qj,hicut,my_vdw_scale14)
             potparm%pot(i,j)%pot %e_cutoff_coul &
                  = ener_coul ( potparm%pot( i, j )%pot, qi, qj, hicut, ewald_env, my_ei_scale14, my_do_14 )
          END IF
          ! Find the real locut according emax_spline
          dx2 = (hicut-locut)/10000.0_dp
          x   = hicut
          DO jx = 1, 10000
             e = ener_pot ( potparm%pot( i, j )%pot, qi, qj, x , my_vdw_scale14 ) &
                  + ener_coul ( potparm%pot( i, j )%pot, qi, qj, x, ewald_env, my_ei_scale14, my_do_14 )
             IF (ABS(e) > emax_spline) THEN
                locut = x
                EXIT
             END IF
             x = x - dx2
          END DO
          npoints = 20
          IF (jx == 10001) npoints = 2  ! Zero spline
          DO WHILE (.TRUE.)
             CALL init_splinexy ( spline_env % spline_data_p ( n ) % spline_data, npoints+1)
             dx2 = (1.0_dp/locut**2-1.0_dp/hicut**2)/REAL(npoints,KIND=dp)
             x2  = 1.0_dp/hicut**2
             DO jx = 1, npoints + 1     ! jx: loop over 1/distance**2
                x = SQRT ( 1.0_dp/x2 )
                e = ener_pot ( potparm%pot( i, j )%pot, qi, qj, x , my_vdw_scale14 ) &
                     + ener_coul ( potparm%pot( i, j )%pot, qi, qj, x, ewald_env, my_ei_scale14, my_do_14 )
                spline_env % spline_data_p ( n ) % spline_data % x(jx) = x2
                spline_env % spline_data_p ( n ) % spline_data % y(jx) = e
                x2 = x2 + dx2
             END DO ! jx: loop over distance
             IF (my_do_14) THEN
                CALL ewald_env_get ( ewald_env, eps0 = eps0 )
                myfac1 = my_ei_scale14 * ifourpi / eps0 * qi * qj 
                myfac2 = 0.0_dp
                IF (scale_cutoff) THEN
                   myfac2 = my_ei_scale14 * ifourpi / eps0 * qi * qj / hicut
                END IF
                CALL init_spline(spline_env%spline_data_p(n)%spline_data,dx=dx2,&
                     ei_scale14=myfac1, ei_cutoff14=myfac2)
             ELSE
                CALL init_spline(spline_env%spline_data_p(n)%spline_data,dx=dx2)
             END IF
             ! This is the check for required accuracy on spline setup
             dx2 = (hicut - locut)/REAL(5*npoints+1,KIND=dp)
             x2  = locut + dx2
             diffmax = 0.0_dp
             xsav = hicut
             DO jx = 1, 5*npoints
                x = x2
                e = ener_pot ( potparm%pot( i, j )%pot, qi, qj, x , my_vdw_scale14 ) &
                     + ener_coul ( potparm%pot( i, j )%pot, qi, qj, x, ewald_env, my_ei_scale14, my_do_14 )
                IF (e<max_energy) THEN
                   xdum1   = ABS(e-potential_s(spline_env%spline_data_p(n)%spline_data,x*x,xdum))
                   diffmax = MAX(diffmax,xdum1)
                   xsav = MIN(x,xsav)
                END IF
                x2 = x2 + dx2
             END DO
             IF (npoints > 5E4) THEN
                IF (iw>0) THEN
                   WRITE(iw,'(A,I5,A,F15.9)') "SPLINE_INFO| Number of points: ",npoints,&
                        " obtained accuracy ",diffmax
                END IF
                WRITE(*,'(A)')"MM SPLINE: no convergence on required accuracy (adjust EPS_SPLINE and rerun)"
                CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
             END IF
             IF (diffmax > eps_spline ) THEN
                npoints   = CEILING(1.2_dp*REAL(npoints,KIND=dp))
             ELSE
                EXIT
             END IF
          END DO
          IF (iw>0) THEN
             WRITE(iw,*)
             WRITE(iw,'("SPLINE_INFO|",5X,A,I7)')"Spline number:",i*(i-1)/2+j
             WRITE(iw,'("SPLINE_INFO|",5X,2(A,I7,2X),A,2I7)')"Spline number",ncount,&
                  "Unique number",n," Kinds involved:",i,j
             WRITE(iw,'("SPLINE_INFO|",5X,A,I7)')"Number of points used in the Splines:",npoints
             WRITE(iw,'("SPLINE_INFO|",5X,A,F18.12)')"Accuracy requested:",eps_spline
             WRITE(iw,'("SPLINE_INFO|",5X,A,F18.12)')"Accuracy achieved: ",diffmax
             WRITE(iw,'("SPLINE_INFO|",5X,A,2F12.6)')"RMIN and RMAX  used for spline (bohr):               ",&
                  locut,hicut
             WRITE(iw,'("SPLINE_INFO|",5X,A,2F12.6)')"RMIN and RMAX used to achieve spline accuracy (bohr):",&
                  xsav,hicut
          ELSE
             IF (iw3>0) THEN
                IF (MOD(i*(i-1)/2+j,MAX(1,(ntype*(ntype+1))/(2*10)))==0) WRITE(iw3,'(T2,A3,I0)') '...',i*(i-1)/2+j
             ENDIF
          END IF
          ! print spline data on file
          IF (iw2>0) THEN
             dx2 = (hicut - locut)/REAL(5000+1,KIND=dp)
             x2  = locut + dx2
             DO jx = 1, 5000
                x = x2
                e = ener_pot ( potparm%pot( i, j )%pot, qi, qj, x , my_vdw_scale14 ) &
                     + ener_coul ( potparm%pot( i, j )%pot, qi, qj, x, ewald_env, my_ei_scale14, my_do_14 )
                WRITE(1000+n,*) x, e, potential_s(spline_env%spline_data_p(n)%spline_data,x*x,xdum)
                x2 = x2 + dx2
             END DO
             WRITE(1000+n,*)
          END IF
       END DO 
    END DO 
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
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE spline_nonbond_control

  !******************************************************************************
  !!****
  !!****s* pair_potential/get_nonbond_storage [1.0] *
  !!
  !!   NAME
  !!     get_nonbond_storage
  !!
  !!   FUNCTION
  !!     Prescreening of the effective bonds evaluations. linear scaling algorithm
  !!
  !!   AUTHOR
  !!     Teodoro Laino [tlaino] 2006.05
  !!
  !!   MODIFICATION HISTORY
  !!     
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE get_nonbond_storage ( spline_env, potparm, atomic_kind_set, &
       ei_scale14, vdw_scale14, error)

    TYPE(spline_environment_type), POINTER     :: spline_env
    TYPE(pair_potential_pp_type), POINTER      :: potparm
    TYPE(atomic_kind_type), DIMENSION(:), &
         POINTER                               :: atomic_kind_set
    REAL(KIND=dp), OPTIONAL                    :: ei_scale14, vdw_scale14
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                              :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_nonbond_storage', &
         routineP = moduleN//':'//routineN
    INTEGER                                    :: n, ntype, handle
    INTEGER                                    :: stat, nk, kk
    REAL(KIND=dp)                              :: qi, qj, my_ei_scale14, my_vdw_scale14
    TYPE(atomic_kind_type), POINTER            :: atomic_kind
    LOGICAL                                    :: failure, check
    INTEGER, DIMENSION(:,:), ALLOCATABLE       :: tmp_index

    INTEGER, DIMENSION(:), ALLOCATABLE         :: my_index, Iwork1, Iwork2
    REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: pot_par
    REAL(KIND=dp), DIMENSION(:),   ALLOCATABLE :: Rwork, Cwork, wtmp
    INTEGER                                    :: pot_target, nvar, istart, iend, i, j, k, nunique
    INTEGER                                    :: ndim, idim, tmpij(2), tmpij0(2), min_val, locij
    LOGICAL                                    :: at_least_one, first_val

    failure     = .FALSE.
    CALL timeset(routineN,'I','',handle)
    ntype = SIZE ( atomic_kind_set )
    my_ei_scale14 = 1.0_dp ; IF (PRESENT(ei_scale14))  my_ei_scale14  = ei_scale14
    my_vdw_scale14= 1.0_dp ; IF (PRESENT(vdw_scale14)) my_vdw_scale14 = vdw_scale14
    DO i = 1, ntype
       DO j = 1, i
          potparm%pot(i,j)%pot%undef = .FALSE. 
       END DO
    END DO
    ALLOCATE(tmp_index(ntype,ntype),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    !
    nunique = 0
    tmp_index = HUGE(0)    
    DO pot_target = MINVAL(list_pot), MAXVAL(list_pot)
       ndim  = 0
       first_val = .FALSE.
       n     = 0
       DO i = 1, ntype
          DO j = 1, i
             n = n + 1
             IF (potparm%pot(i,j)%pot%type==pot_target) THEN
                IF (.NOT.first_val) THEN 
                   min_val = n
                   first_val = .TRUE.
                END IF
                tmp_index(i,j) = 1
                tmp_index(j,i) = 1
                ndim = ndim + 1
             END IF
          END DO
       END DO
       IF (ndim==0) CYCLE ! No potential of this kind found
       CALL get_indexes(min_val,ntype,tmpij)
       potparm%pot(tmpij(1),tmpij(2))%pot%undef = .TRUE.
       n     = 1
       ALLOCATE(my_index(ndim),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       nvar = 0
       IF (my_ei_scale14 /=0.0_dp) nvar = nvar + 2
       SELECT CASE(pot_target)
       CASE (lj_type)
          IF (my_vdw_scale14/=0.0_dp) nvar = 3 + nvar
       CASE (wl_type)
          nvar = 3 + nvar
       CASE (gw_type)
          nvar = 5 + nvar
       CASE (ea_type)
          nvar = 3 + nvar
       CASE (ft_type)
          nvar = 4 + nvar
       CASE (ip_type)
          nvar = 3 + nvar
       CASE (nn_type)
          nvar = nvar
       CASE DEFAULT
          CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
       END SELECT
       IF (nvar /=0) THEN
          ALLOCATE(pot_par(ndim,nvar),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          n  = 0
          nk = 0
          DO i = 1, ntype
             atomic_kind=> atomic_kind_set(i)
             CALL get_atomic_kind(atomic_kind=atomic_kind,qeff=qi)
             DO j = 1, i
                atomic_kind=> atomic_kind_set(j)
                CALL get_atomic_kind(atomic_kind=atomic_kind,qeff=qj)
                n = n + 1
                IF (potparm%pot(i,j)%pot%type==pot_target) THEN
                   nk = nk + 1
                   my_index(nk) = n
                   SELECT CASE(pot_target)
                   CASE (lj_type)
                      pot_par(nk,1)  = potparm%pot(i,j)%pot%lj%epsilon
                      pot_par(nk,2)  = potparm%pot(i,j)%pot%lj%sigma6
                      pot_par(nk,3)  = potparm%pot(i,j)%pot%lj%sigma12
                   CASE (wl_type)
                      pot_par(nk,1)  = potparm%pot(i,j)%pot%willis%a
                      pot_par(nk,2)  = potparm%pot(i,j)%pot%willis%b
                      pot_par(nk,3)  = potparm%pot(i,j)%pot%willis%c
                   CASE (gw_type)
                      pot_par(nk,1)  = potparm%pot(i,j)%pot%goodwin%vr0
                      pot_par(nk,2)  = potparm%pot(i,j)%pot%goodwin%m  
                      pot_par(nk,3)  = potparm%pot(i,j)%pot%goodwin%mc 
                      pot_par(nk,4)  = potparm%pot(i,j)%pot%goodwin%d  
                      pot_par(nk,5)  = potparm%pot(i,j)%pot%goodwin%dc 
                   CASE (ea_type)
                      pot_par(nk,1)  = potparm%pot(i,j)%pot%eam%drar  
                      pot_par(nk,2)  = potparm%pot(i,j)%pot%eam%drhoar
                      pot_par(nk,3)  = potparm%pot(i,j)%pot%eam%acutal
                   CASE (ft_type)
                      pot_par(nk,1)  = potparm%pot(i,j)%pot%ft%A
                      pot_par(nk,2)  = potparm%pot(i,j)%pot%ft%B
                      pot_par(nk,3)  = potparm%pot(i,j)%pot%ft%C
                      pot_par(nk,4)  = potparm%pot(i,j)%pot%ft%D
                   CASE (ip_type)
                      pot_par(nk,1)  = potparm%pot(i,j)%pot%ipbv%rcore
                      pot_par(nk,2)  = potparm%pot(i,j)%pot%ipbv%m
                      pot_par(nk,3)  = potparm%pot(i,j)%pot%ipbv%b
                   CASE (nn_type)
                      ! no checks
                   CASE DEFAULT
                      CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
                   END SELECT
                   IF (my_ei_scale14 /=0.0_dp) THEN 
                      pot_par(nk,nvar-1)=qi
                      pot_par(nk,nvar  )=qj
                   END IF
                END IF
             END DO
          END DO
          !
          ! MAIN SORTING LOOP
          !
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
             SELECT CASE(pot_target)
             CASE(ea_type,ip_type)
                ! check the array components
                locij = my_index(j)
                CALL get_indexes(locij,ntype,tmpij)
                CALL compare_pot(potparm%pot(tmpij(1) ,tmpij(2) )%pot,&
                                 potparm%pot(tmpij0(1),tmpij0(2))%pot,&
                                 check, error)
             CASE default
                check = .TRUE.
             END SELECT
             IF (ALL(Cwork==pot_par(j,:)).AND.check) CYCLE
             Cwork = pot_par(j,:)
             ! Index for special cases
             SELECT CASE(pot_target)
             CASE(ea_type,ip_type)
                locij = my_index(j)
                CALL get_indexes(locij,ntype,tmpij0)
             END SELECT
             nunique = nunique + 1
             iend    = j - 1
             min_val = MINVAL(my_index(istart:iend))
             CALL get_indexes(min_val,ntype,tmpij)
             potparm%pot(tmpij(1),tmpij(2))%pot%undef =.TRUE.
             potparm%pot(tmpij(2),tmpij(1))%pot%undef =.TRUE.
             DO i = istart, iend
                locij = my_index(i)
                CALL get_indexes(locij,ntype,tmpij)
                tmp_index(tmpij(1),tmpij(2))=nunique
                tmp_index(tmpij(2),tmpij(1))=nunique
             END DO
             istart = j
          END DO
          nunique = nunique + 1
          iend    = ndim
          min_val = MINVAL(my_index(istart:iend)) 
          CALL get_indexes(min_val,ntype,tmpij)
          potparm%pot(tmpij(1),tmpij(2))%pot%undef =.TRUE.
          potparm%pot(tmpij(2),tmpij(1))%pot%undef =.TRUE.
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
       END IF
       DEALLOCATE(my_index,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END DO
    CALL spline_env_create ( spline_env, ntype, ntab_in=nunique, error=error )
    spline_env%spltab=tmp_index
    DEALLOCATE(tmp_index,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE get_nonbond_storage

  !******************************************************************************
  !!****
  !!****s* pair_potential/get_indexes [1.0] *
  !!
  !!   NAME
  !!     get_indexes
  !!
  !!   FUNCTION
  !!     Gives back the indices of the matrix w.r.t. the collective array index
  !!
  !!   AUTHOR
  !!     Teodoro Laino [tlaino] 2006.05
  !!
  !!   MODIFICATION HISTORY
  !!     
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE get_indexes(Inind, ndim, ij)
    IMPLICIT NONE
    INTEGER, INTENT(IN)                :: Inind, ndim
    INTEGER, DIMENSION(2), INTENT(OUT) :: ij
    
    INTEGER :: i, tmp
    
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

  !******************************************************************************
  !!****
  !!****s* pair_potential/potential_f [1.0] *
  !!
  !!   NAME
  !!     potential_f
  !!
  !!   FUNCTION
  !!    Evaluates the nonbond and electrostaitic potential analytically
  !!
  !!   AUTHOR
  !!     Teodoro Laino [tlaino] 2006.05
  !!
  !!   MODIFICATION HISTORY
  !!     
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE potential_f ( rijsq, potparm, qi, qj, iatom, jatom, energy, ewald_env )
    REAL(KIND=dp), INTENT(IN)                :: rijsq
    TYPE(pair_potential_pp_type), POINTER    :: potparm
    REAL(KIND=dp), INTENT(IN)                :: qi, qj
    INTEGER, INTENT(IN)                      :: iatom, jatom
    REAL(KIND=dp), INTENT(OUT)               :: energy
    TYPE(ewald_environment_type), OPTIONAL, &
         POINTER                                :: ewald_env

    REAL(KIND=dp)                            :: rij

    rij = SQRT ( rijsq )

    ! non-bond energy
    energy = ener_pot ( potparm%pot(iatom,jatom)%pot, qi, qj, rij, 1.0_dp)

    ! real-space ewald non-bond energy
    IF ( PRESENT ( ewald_env ) ) THEN
       energy = energy &
            + ener_coul ( potparm%pot(iatom,jatom )%pot, qi, qj, rij, ewald_env, 1.0_dp,&
            my_do_14=.FALSE.)
    END IF

  END SUBROUTINE potential_f

END MODULE pair_potential

