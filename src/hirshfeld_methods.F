!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculate Hirshfeld charges and related functions
!> \par History
!>      11.2014 created [JGH]
!> \author JGH
! *****************************************************************************
MODULE hirshfeld_methods
  USE atom_kind_orbitals,              ONLY: calculate_atomic_density
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cube_utils,                      ONLY: cube_info_type
  USE hirshfeld_types,                 ONLY: get_hirshfeld_info,&
                                             hirshfeld_type,&
                                             set_hirshfeld_info
  USE input_constants,                 ONLY: shape_function_density,&
                                             shape_function_gaussian
  USE kinds,                           ONLY: dp,&
                                             int_8
  USE mathconstants,                   ONLY: pi
  USE message_passing,                 ONLY: mp_sum
  USE particle_types,                  ONLY: particle_type
  USE periodic_table,                  ONLY: get_ptable_info
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_methods,                      ONLY: pw_integrate_function
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_type
  USE pw_types,                        ONLY: REALDATA3D,&
                                             REALSPACE,&
                                             pw_p_type,&
                                             pw_release
  USE qs_collocate_density,            ONLY: collocate_pgf_product_rspace
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_integrate_potential_low,      ONLY: integrate_pgf_product_rspace
  USE qs_kind_types,                   ONLY: get_qs_kind,&
                                             qs_kind_type
  USE qs_modify_pab_block,             ONLY: FUNC_AB
  USE qs_rho_types,                    ONLY: qs_rho_get,&
                                             qs_rho_type
  USE realspace_grid_types,            ONLY: pw2rs,&
                                             realspace_grid_desc_type,&
                                             realspace_grid_type,&
                                             rs2pw,&
                                             rs_grid_release,&
                                             rs_grid_retain,&
                                             rs_grid_zero,&
                                             rs_pw_transfer
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'hirshfeld_methods'

  PUBLIC :: create_shape_function, comp_hirshfeld_charges,&
            comp_hirshfeld_i_charges, write_hirshfeld_charges

! *****************************************************************************

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param charges ...
!> \param hirshfeld_env ...
!> \param particle_set ...
!> \param qs_kind_set ...
!> \param unit_nr ...
! *****************************************************************************
  SUBROUTINE write_hirshfeld_charges(charges,hirshfeld_env,particle_set,&
             qs_kind_set,unit_nr)
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(inout)                          :: charges
    TYPE(hirshfeld_type), POINTER            :: hirshfeld_env
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    INTEGER, INTENT(IN)                      :: unit_nr

    CHARACTER(len=*), PARAMETER :: routineN = 'write_hirshfeld_charges', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=2)                         :: element_symbol
    INTEGER                                  :: iatom, ikind, natom, nspin
    REAL(KIND=dp)                            :: refc, tc1, zeff

    natom = SIZE(charges,1)
    nspin = SIZE(charges,2)
    WRITE(unit_nr,'(/,T2,A)') '!-----------------------------------------------------------------------------!'
    WRITE (UNIT=unit_nr,FMT="(T28,A)") "Hirshfeld Charges"
    IF (nspin == 1) THEN
       WRITE (UNIT=unit_nr,FMT="(/,T3,A,A)")&
          "#Atom  Element  Kind "," Ref Charge     Population                    Net charge"
    ELSE
       WRITE (UNIT=unit_nr,FMT="(/,T3,A,A)")&
          "#Atom  Element  Kind "," Ref Charge     Population       Spin moment  Net charge"
    END IF
    tc1 = 0.0_dp
    DO iatom=1,natom
       CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
                            element_symbol=element_symbol,kind_number=ikind)
       refc = hirshfeld_env%charges(iatom)
       CALL get_qs_kind(qs_kind_set(ikind),zeff=zeff)
       IF (nspin == 1) THEN
          WRITE (UNIT=unit_nr,FMT="(i7,T15,A2,T20,i3,T27,F8.3,T42,F8.3,T72,F8.3)") &
            iatom,element_symbol,ikind,refc,charges(iatom,1),zeff-charges(iatom,1)
       ELSE
          WRITE (UNIT=unit_nr,FMT="(i7,T15,A2,T20,i3,T27,F8.3,T36,2F8.3,T61,F8.3,T72,F8.3)") &
            iatom,element_symbol,ikind,refc,charges(iatom,1),charges(iatom,2),&
            charges(iatom,1)-charges(iatom,2),zeff-SUM(charges(iatom,:))
       END IF
       tc1 = tc1 + (zeff-SUM(charges(iatom,:)))
    END DO
    WRITE (UNIT=unit_nr,FMT="(/,T3,A,T72,F8.3)") "Total Charge ",tc1
    WRITE(unit_nr,'(T2,A)') '!-----------------------------------------------------------------------------!'

  END SUBROUTINE write_hirshfeld_charges

! *****************************************************************************
!> \brief ...
!> \param hirshfeld_env ...
!> \param qs_kind_set ...
!> \param atomic_kind_set ...
! *****************************************************************************
  SUBROUTINE create_shape_function(hirshfeld_env,qs_kind_set,atomic_kind_set)
    TYPE(hirshfeld_type), POINTER            :: hirshfeld_env
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set

    CHARACTER(len=*), PARAMETER :: routineN = 'create_shape_function', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: ngto = 8

    CHARACTER(len=2)                         :: esym
    INTEGER                                  :: ikind, nkind, stat
    LOGICAL                                  :: failure, found
    REAL(KIND=dp)                            :: al, rco, zeff
    REAL(KIND=dp), DIMENSION(ngto, 2)        :: ppdens
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(qs_kind_type), POINTER              :: qs_kind

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(hirshfeld_env),cp_failure_level,routineP,failure)

    nkind = SIZE(qs_kind_set)
    ALLOCATE(hirshfeld_env%kind_shape_fn(nkind), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,failure)

    SELECT CASE(hirshfeld_env%shape_function_type)
       CASE (shape_function_gaussian)
          DO ikind=1,nkind
             hirshfeld_env%kind_shape_fn(ikind)%numexp = 1
             ALLOCATE(hirshfeld_env%kind_shape_fn(ikind)%zet(1), stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,failure)
             ALLOCATE(hirshfeld_env%kind_shape_fn(ikind)%coef(1), stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,failure)
             CALL get_qs_kind(qs_kind_set(ikind),element_symbol=esym)
             rco = 2.0_dp
             CALL get_ptable_info(symbol=esym,covalent_radius=rco,found=found)
             rco = MAX(rco,1.0_dp)
             al = 0.5_dp/rco**2
             hirshfeld_env%kind_shape_fn(ikind)%zet(1) = al
             hirshfeld_env%kind_shape_fn(ikind)%coef(1) = (al/pi)**1.5_dp
          END DO
       CASE (shape_function_density)
          ! calculate atomic density
          DO ikind=1,nkind
             atomic_kind => atomic_kind_set(ikind)
             qs_kind     => qs_kind_set(ikind)
             CALL calculate_atomic_density(ppdens(:,:),atomic_kind,qs_kind,ngto,&
                                           confine=.FALSE.)
             hirshfeld_env%kind_shape_fn(ikind)%numexp = ngto
             ALLOCATE(hirshfeld_env%kind_shape_fn(ikind)%zet(ngto), stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,failure)
             ALLOCATE(hirshfeld_env%kind_shape_fn(ikind)%coef(ngto), stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,failure)
             hirshfeld_env%kind_shape_fn(ikind)%zet(:) = ppdens(:,1)
             CALL get_qs_kind(qs_kind,zeff=zeff)
             hirshfeld_env%kind_shape_fn(ikind)%coef(:) = ppdens(:,2)/zeff
          END DO

       CASE DEFAULT
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="Unknown shape function", &
               error_level=cp_failure_level)
    END SELECT

  END SUBROUTINE create_shape_function

! *****************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param hirshfeld_env ...
!> \param charges ...
! *****************************************************************************
  SUBROUTINE comp_hirshfeld_charges(qs_env,hirshfeld_env,charges)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(hirshfeld_type), POINTER            :: hirshfeld_env
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(inout)                          :: charges

    CHARACTER(len=*), PARAMETER :: routineN = 'comp_hirshfeld_charges', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: is, stat
    LOGICAL                                  :: failure, rho_r_valid
    REAL(KIND=dp)                            :: tnfun
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: rho_r
    TYPE(pw_p_type), POINTER                 :: rhonorm
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(qs_rho_type), POINTER               :: rho

    NULLIFY(rho_r)
    failure = .FALSE.
    ! normalization function on grid
    CALL calculate_hirshfeld_normalization(qs_env,hirshfeld_env)
    ! check normalization
    tnfun = pw_integrate_function(hirshfeld_env%fnorm%pw)
    tnfun = ABS(tnfun - SUM(hirshfeld_env%charges))
    ! 
    ALLOCATE(rhonorm, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,failure)
    !
    CALL get_qs_env(qs_env=qs_env,pw_env=pw_env,rho=rho)
    CALL qs_rho_get(rho, rho_r=rho_r, rho_r_valid=rho_r_valid)
    CALL pw_env_get(pw_env=pw_env,auxbas_pw_pool=auxbas_pw_pool)
    CALL pw_pool_create_pw(auxbas_pw_pool,rhonorm%pw,use_data=REALDATA3D)
    ! loop over spins
    DO is=1,SIZE(rho_r)
      IF(rho_r_valid) THEN
        CALL hfun_scale(rhonorm%pw%cr3d,rho_r(is)%pw%cr3d,&
                        hirshfeld_env%fnorm%pw%cr3d)
      ELSE
        CALL cp_unimplemented_error(fromWhere=routineP, &
             message="We need rho in real space",error_level=cp_failure_level)
      END IF
      CALL hirshfeld_integration(qs_env,hirshfeld_env,rhonorm,charges(:,is))
      charges(:,is) = charges(:,is)*hirshfeld_env%charges(:)
    END DO
    CALL pw_pool_give_back_pw(auxbas_pw_pool,rhonorm%pw)
    !
    DEALLOCATE(rhonorm, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,failure)

  END SUBROUTINE comp_hirshfeld_charges
! *****************************************************************************
!> \brief Calculate fout = fun1/fun2
!> \param fout ...
!> \param fun1 ...
!> \param fun2 ...
! *****************************************************************************
  SUBROUTINE hfun_scale(fout,fun1,fun2)
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(OUT)                            :: fout
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(IN)                             :: fun1, fun2

    CHARACTER(len=*), PARAMETER :: routineN = 'hfun_scale', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: small = 1.0e-12_dp

    INTEGER                                  :: i1, i2, i3, n1, n2, n3
    LOGICAL                                  :: failure

    failure = .FALSE.
    n1 = SIZE(fout,1)
    n2 = SIZE(fout,2)
    n3 = SIZE(fout,3)
    CPPrecondition(n1==SIZE(fun1,1),cp_failure_level,routineP,failure)
    CPPrecondition(n2==SIZE(fun1,2),cp_failure_level,routineP,failure)
    CPPrecondition(n3==SIZE(fun1,3),cp_failure_level,routineP,failure)
    CPPrecondition(n1==SIZE(fun2,1),cp_failure_level,routineP,failure)
    CPPrecondition(n2==SIZE(fun2,2),cp_failure_level,routineP,failure)
    CPPrecondition(n3==SIZE(fun2,3),cp_failure_level,routineP,failure)

    DO i3=1,n3
      DO i2=1,n2
        DO i1=1,n1
          IF(fun2(i1,i2,i3) > small) THEN
            fout(i1,i2,i3) = fun1(i1,i2,i3)/fun2(i1,i2,i3)
          ELSE
            fout(i1,i2,i3) = 0.0_dp
          END IF
        END DO
      END DO
    END DO

  END SUBROUTINE hfun_scale

! *****************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param hirshfeld_env ...
!> \param charges ...
!> \param ounit ...
! *****************************************************************************
  SUBROUTINE comp_hirshfeld_i_charges(qs_env,hirshfeld_env,charges,ounit)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(hirshfeld_type), POINTER            :: hirshfeld_env
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(inout)                          :: charges
    INTEGER, INTENT(IN)                      :: ounit

    CHARACTER(len=*), PARAMETER :: routineN = 'comp_hirshfeld_i_charges', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: maxloop = 100
    REAL(KIND=dp), PARAMETER                 :: maxres = 1.0e-2_dp

    CHARACTER(len=3)                         :: yesno
    INTEGER                                  :: iat, iloop, is, natom, stat
    LOGICAL                                  :: failure, rho_r_valid
    REAL(KIND=dp)                            :: res, tnfun
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: rho_r
    TYPE(pw_p_type), POINTER                 :: rhonorm
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(qs_rho_type), POINTER               :: rho

    NULLIFY(rho_r)
    failure = .FALSE.

    natom = SIZE(charges,1)

    IF(ounit>0) WRITE(ounit,"(/,T2,A)") "Hirshfeld charge iterations: Residuals ..."
    !
    ALLOCATE(rhonorm, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,failure)
    !
    CALL get_qs_env(qs_env=qs_env,pw_env=pw_env,rho=rho)
    CALL qs_rho_get(rho, rho_r=rho_r, rho_r_valid=rho_r_valid)
    CALL pw_env_get(pw_env=pw_env,auxbas_pw_pool=auxbas_pw_pool)
    CALL pw_pool_create_pw(auxbas_pw_pool,rhonorm%pw,use_data=REALDATA3D)
    !
    DO iloop=1,maxloop

       ! normalization function on grid
       CALL calculate_hirshfeld_normalization(qs_env,hirshfeld_env)
       ! check normalization
       tnfun = pw_integrate_function(hirshfeld_env%fnorm%pw)
       tnfun = ABS(tnfun - SUM(hirshfeld_env%charges))
       ! loop over spins
       DO is=1,SIZE(rho_r)
         IF(rho_r_valid) THEN
           CALL hfun_scale(rhonorm%pw%cr3d,rho_r(is)%pw%cr3d,&
                           hirshfeld_env%fnorm%pw%cr3d)
         ELSE
           CALL cp_unimplemented_error(fromWhere=routineP, &
                message="We need rho in real space",error_level=cp_failure_level)
         END IF
         CALL hirshfeld_integration(qs_env,hirshfeld_env,rhonorm,charges(:,is))
         charges(:,is) = charges(:,is)*hirshfeld_env%charges(:)
       END DO
       ! residual
       res=0.0_dp
       DO iat=1,natom
          res = res + (SUM(charges(iat,:)) - hirshfeld_env%charges(iat))**2
       END DO
       res = SQRT(res/REAL(natom,KIND=dp))
       IF(ounit>0) THEN
          yesno = "NO "
          IF(MOD(iloop,10)==0) yesno = "YES"
          WRITE(ounit,FMT="(F8.3)",ADVANCE=yesno) res
       END IF
       ! update
       DO iat=1,natom
          hirshfeld_env%charges(iat) = SUM(charges(iat,:))
       END DO
       IF(res<maxres) EXIT

    END DO
    !
    CALL pw_pool_give_back_pw(auxbas_pw_pool,rhonorm%pw)
    !
    DEALLOCATE(rhonorm, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,failure)

  END SUBROUTINE comp_hirshfeld_i_charges

! *****************************************************************************
!> \brief 
!> \param qs_env ...
!> \param hirshfeld_env ...
! *****************************************************************************
  SUBROUTINE calculate_hirshfeld_normalization(qs_env,hirshfeld_env)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(hirshfeld_type), POINTER            :: hirshfeld_env

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'calculate_hirshfeld_normalization', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: atom_a, handle, iatom, ierr, &
                                                iex, ikind, ithread, j, &
                                                natom, npme, nthread, numexp
    INTEGER(KIND=int_8)                      :: subpatch_pattern
    INTEGER, DIMENSION(:), POINTER           :: atom_list, cores
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: alpha, coef, eps_rho_rspace
    REAL(KIND=dp), DIMENSION(3)              :: ra
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: pab
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cube_info_type)                     :: cube_info
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type), POINTER                 :: fnorm
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(realspace_grid_desc_type), POINTER  :: auxbas_rs_desc
    TYPE(realspace_grid_type), POINTER       :: rs_rho

    CALL timeset(routineN,handle)
    failure = .FALSE.

    CALL get_qs_env(qs_env=qs_env,atomic_kind_set=atomic_kind_set,cell=cell,&
         dft_control=dft_control,particle_set=particle_set,pw_env=pw_env)
    CALL pw_env_get(pw_env,auxbas_rs_desc=auxbas_rs_desc,auxbas_rs_grid=rs_rho,&
         auxbas_pw_pool=auxbas_pw_pool)
    cube_info=pw_env%cube_info(1)
    ! be careful in parallel nsmax is choosen with multigrid in mind!
    CALL rs_grid_retain(rs_rho)
    CALL rs_grid_zero(rs_rho)

    eps_rho_rspace = dft_control%qs_control%eps_rho_rspace
    ALLOCATE(pab(1,1),STAT=ierr)
    CPPostcondition(ierr==0,cp_failure_level,routineP,failure)
    nthread = 1
    ithread = 0

    DO ikind=1,SIZE(atomic_kind_set)
      numexp = hirshfeld_env%kind_shape_fn(ikind)%numexp
      IF ( numexp <= 0 ) CYCLE
      CALL get_atomic_kind(atomic_kind_set(ikind), natom=natom, atom_list=atom_list)
      ALLOCATE(cores(natom),STAT=ierr)
      CPPostcondition(ierr==0,cp_failure_level,routineP,failure)

      DO iex=1,numexp
         alpha=hirshfeld_env%kind_shape_fn(ikind)%zet(iex)
         coef =hirshfeld_env%kind_shape_fn(ikind)%coef(iex)
         npme = 0
         cores = 0
         DO iatom = 1, natom
            atom_a = atom_list(iatom)
            ra(:) = pbc(particle_set(atom_a)%r,cell)
            IF(rs_rho%desc%parallel .AND. .NOT. rs_rho%desc%distributed) THEN
               ! replicated realspace grid, split the atoms up between procs
               IF (MODULO(iatom,rs_rho%desc%group_size) == rs_rho % desc % my_pos ) THEN
                  npme = npme + 1
                  cores (npme) = iatom
               ENDIF
            ELSE
               npme = npme + 1
               cores (npme) = iatom
            ENDIF
         END DO
         DO j=1,npme
            iatom = cores(j)
            atom_a = atom_list(iatom)
            pab(1,1) = hirshfeld_env%charges(atom_a) * coef
            ra(:) = pbc(particle_set(atom_a)%r,cell)
            subpatch_pattern=0
            CALL collocate_pgf_product_rspace(0,alpha,0,0,0.0_dp,0,ra,&
               (/0.0_dp,0.0_dp,0.0_dp/),0.0_dp,1.0_dp,pab,0,0,rs_rho,&
               cell,cube_info,eps_rho_rspace,ga_gb_function=FUNC_AB,&
               ithread=ithread,use_subpatch=.TRUE.,subpatch_pattern=subpatch_pattern)
         END DO
      END DO

      DEALLOCATE(cores,STAT=ierr)
      CPPostcondition(ierr==0,cp_failure_level,routineP,failure)
    END DO
    DEALLOCATE (pab,STAT=ierr)
    CPPostcondition(ierr==0,cp_failure_level,routineP,failure)

    NULLIFY(fnorm)
    CALL get_hirshfeld_info(hirshfeld_env,fnorm=fnorm)
    IF (ASSOCIATED(fnorm)) THEN
       CALL pw_release(fnorm%pw)
       DEALLOCATE(fnorm,stat=ierr)
       CPPostconditionNoFail(ierr==0,cp_failure_level,routineP)
    ENDIF
    ALLOCATE(fnorm,stat=ierr)
    CPPostcondition(ierr==0,cp_failure_level,routineP,failure)
    CALL pw_pool_create_pw(auxbas_pw_pool,fnorm%pw,use_data=REALDATA3D)
    fnorm%pw%in_space=REALSPACE
    CALL set_hirshfeld_info(hirshfeld_env,fnorm=fnorm)

    CALL rs_pw_transfer(rs_rho,fnorm%pw,rs2pw)
    CALL rs_grid_release(rs_rho)

    CALL timestop(handle)

  END SUBROUTINE calculate_hirshfeld_normalization

! *****************************************************************************
!> \brief 
!> \param qs_env ...
!> \param hirshfeld_env ...
!> \param rfun ...
!> \param fval ...
!> \param fderiv ...
! *****************************************************************************
  SUBROUTINE hirshfeld_integration(qs_env,hirshfeld_env,rfun,fval,fderiv)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(hirshfeld_type), POINTER            :: hirshfeld_env
    TYPE(pw_p_type), POINTER                 :: rfun
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(inout)                          :: fval
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(inout), OPTIONAL                :: fderiv

    CHARACTER(len=*), PARAMETER :: routineN = 'hirshfeld_integration', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: atom_a, handle, iatom, iex, &
                                                ikind, ithread, j, natom, &
                                                npme, nthread, numexp, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: cores
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: do_force, failure
    REAL(KIND=dp)                            :: alpha, coef, dvol, &
                                                eps_rho_rspace
    REAL(KIND=dp), DIMENSION(3)              :: force_a, force_b, ra
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: hab, pab
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(realspace_grid_desc_type), POINTER  :: auxbas_rs_desc
    TYPE(realspace_grid_type), POINTER       :: rs_v

    CALL timeset(routineN,handle)

    failure=.FALSE.
    do_force = PRESENT(fderiv)
    fval = 0.0_dp
    dvol = rfun%pw%pw_grid%dvol

    NULLIFY(pw_env,auxbas_rs_desc)
    CALL get_qs_env(qs_env=qs_env,pw_env=pw_env)
    CALL pw_env_get(pw_env=pw_env,auxbas_rs_desc=auxbas_rs_desc, &
                    auxbas_rs_grid=rs_v)
    CALL rs_grid_retain(rs_v)
    CALL rs_pw_transfer(rs_v,rfun%pw,pw2rs)

    CALL get_qs_env(qs_env=qs_env,atomic_kind_set=atomic_kind_set,cell=cell,&
         dft_control=dft_control,particle_set=particle_set)
    eps_rho_rspace = dft_control%qs_control%eps_rho_rspace

    nthread = 1
    ithread = 0
    ALLOCATE(hab(1,1),pab(1,1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,failure)

    DO ikind=1,SIZE(atomic_kind_set)
      numexp = hirshfeld_env%kind_shape_fn(ikind)%numexp
      IF ( numexp <= 0 ) CYCLE
      CALL get_atomic_kind(atomic_kind_set(ikind),natom=natom,atom_list=atom_list)
      ALLOCATE(cores(natom),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,failure)

      DO iex=1,numexp
         alpha=hirshfeld_env%kind_shape_fn(ikind)%zet(iex)
         coef =hirshfeld_env%kind_shape_fn(ikind)%coef(iex)
         npme = 0
         cores = 0
         DO iatom = 1, natom
            atom_a = atom_list(iatom)
            ra(:) = pbc(particle_set(atom_a)%r,cell)
            IF(rs_v%desc%parallel .AND. .NOT. rs_v%desc%distributed) THEN
                ! replicated realspace grid, split the atoms up between procs
                IF (MODULO(iatom,rs_v%desc%group_size) == rs_v % desc % my_pos ) THEN
                   npme = npme + 1
                   cores (npme) = iatom
                ENDIF
             ELSE
                npme = npme + 1
                cores (npme) = iatom
             ENDIF
         END DO

         DO j=1,npme
           iatom = cores(j)
           atom_a = atom_list(iatom)
           ra(:) = pbc(particle_set(atom_a)%r,cell)
           pab(1,1) = coef
           hab(1,1) = 0.0_dp
           force_a(:) = 0.0_dp
           force_b(:) = 0.0_dp
           !
           CALL integrate_pgf_product_rspace(0,alpha,0,&
                0,0.0_dp,0,ra,(/0.0_dp,0.0_dp,0.0_dp/),0.0_dp,&
                rs_v,cell,pw_env%cube_info(1),hab,pab=pab,o1=0,o2=0,&
                eps_gvg_rspace=eps_rho_rspace,calculate_forces=do_force,&
                force_a=force_a,force_b=force_b,use_virial=.FALSE.,&
                use_subpatch=.TRUE.,subpatch_pattern=0_int_8)
           fval(atom_a) = fval(atom_a) + hab(1,1)*dvol*coef
           IF(do_force) THEN
              fderiv(:,atom_a) = fderiv(:,atom_a) + force_a(:)*dvol
           END IF
         END DO

      END DO
      DEALLOCATE(cores,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,failure)

    END DO

    CALL rs_grid_release(rs_v)
    DEALLOCATE(hab,pab,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,failure)

    CALL get_qs_env(qs_env=qs_env,para_env=para_env)
    CALL mp_sum(fval,para_env%group)

    CALL timestop(handle)

  END SUBROUTINE hirshfeld_integration

END MODULE hirshfeld_methods
