!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2003 CP2K developers group                                  !
!-----------------------------------------------------------------------------!
!!****** cp2k/level_shift [1.0] *
!!
!!   NAME
!!     level_shift
!!
!!   FUNCTION
!!     eigenvalue shifting
!!
!!   AUTHOR
!!     Joost VandeVondele (10.2003)
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

MODULE level_shift
! *****************************************************************************
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE coefficient_types,               ONLY: coeff_transform_space,&
                                             coeff_type, &
                                             coeff_sumup, &
                                             coeff_add
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_debug,&
                                             cp_error_dealloc_ref,&
                                             cp_error_get_logger,&
                                             cp_error_init,&
                                             cp_error_message,&
                                             cp_error_type,&
                                             cp_internal_error,&
                                             cp_unimplemented_error
  USE cp_rs_pool_types,                ONLY: cp_rs_pool_type,&
                                             rs_pool_create_rs,&
                                             rs_pool_give_back_rs
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: twopi
  USE particle_types,                  ONLY: particle_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             COMPLEXDATA3D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_type
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_pool_types,                   ONLY: pw_pool_give_back_coeff,&
                                             pw_pool_init_coeff,&
                                             pw_pool_type
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE realspace_grid_types,            ONLY: realspace_grid_type,&
                                             rs_pw_transfer
  USE simulation_cell,                 ONLY: cell_type,&
                                             pbc
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  IMPLICIT NONE

  PRIVATE

  PUBLIC level_shift_testing

!!***
! *****************************************************************************

CONTAINS

  SUBROUTINE level_shift_testing(atomic_kind_set, particle_set, cell, centers, pw_env, rho, rho_core, error )
    TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), POINTER    :: particle_set
    TYPE(cell_type), POINTER                      :: cell
    TYPE(pw_env_type), POINTER                    :: pw_env
    TYPE(coeff_type), POINTER                     :: rho_core
    TYPE(qs_rho_type)                             :: rho
    REAL(KIND = dp), DIMENSION(:,:), POINTER             :: centers
    TYPE(cp_error_type), OPTIONAL, INTENT(INOUT)  :: error

    TYPE(coeff_type)            :: rho_tot_rspace, rho_tot_gspace
    TYPE(pw_pool_type), POINTER :: auxbas_pw_pool
    INTEGER :: ispin, handle, power
    REAL(KIND = dp) :: origin(3),electronic_eshift,core_eshift

    CALL timeset("level_shift_testing","I","",handle)

    CALL wannier_part(centers,electronic_eshift)
    CALL core_part(atomic_kind_set, particle_set, cell, core_eshift)

    write(6,*) "electronic_eshift       ",electronic_eshift
    write(6,*) "core_eshift             ",core_eshift
    write(6,*) "total eshift            ",core_eshift+electronic_eshift
    write(6,*) "eval shift in a.u. ", &
                   (core_eshift+electronic_eshift)*twopi / ( 3.0_dp * cell%deth )

    ! realspace based approach 
    ! generate total density in real space
    NULLIFY(auxbas_pw_pool)
    CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool) 
    CALL pw_pool_init_coeff(auxbas_pw_pool,rho_tot_rspace,&
            use_data = REALDATA3D, in_space = REALSPACE, error=error )
    CALL pw_pool_init_coeff(auxbas_pw_pool,rho_tot_gspace,&
            use_data = COMPLEXDATA1D, in_space = RECIPROCALSPACE, error=error )
    CALL coeff_transform_space(rho_core,rho_tot_rspace)
    CALL coeff_sumup(rho%rho_r(1),rho_tot_rspace)
    DO ispin=2,SIZE(rho%rho_r)
       CALL coeff_sumup(rho%rho_r(ispin),rho_tot_rspace)
    END DO
    CALL coeff_transform_space(rho_tot_rspace,rho_tot_gspace)

    ! start realspace analysis
    origin=0
    power=1
    CALL  realspace_approach(rho_tot_rspace,pw_env,cell,origin,power,.TRUE.)

    ! start g-space analysis
    CALL  gspace_approach(rho_tot_gspace,pw_env,cell)

    ! give back coeffs and pools
    CALL pw_pool_give_back_coeff(auxbas_pw_pool,rho_tot_rspace)
    CALL pw_pool_give_back_coeff(auxbas_pw_pool,rho_tot_gspace)

    CALL timestop(0.0_dp,handle)
 
  END SUBROUTINE level_shift_testing

  SUBROUTINE wannier_part(centers,electronic_eshift)
    REAL(KIND = dp), DIMENSION(:,:), POINTER :: centers
    REAL(KIND = dp) ::electronic_eshift

    INTEGER :: i
    electronic_eshift=0.0_dp
    DO i=1,SIZE(centers,1)
       electronic_eshift=electronic_eshift&
                         -2.0_dp*SUM(centers(i,1:3)**2+centers(i,4:6)**2) ! NOT FOR LSD
    ENDDO
  END SUBROUTINE wannier_part

  SUBROUTINE core_part(atomic_kind_set, particle_set, cell, core_eshift)
    REAL(KIND = dp) :: core_eshift
    TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), POINTER    :: particle_set
    TYPE(cell_type), POINTER        :: cell

    TYPE(atomic_kind_type), POINTER :: atomic_kind
    INTEGER, DIMENSION(:), POINTER  :: atom_list
    INTEGER :: ikind,natom,iatom,atom_a,zeff
    REAL(KIND = dp)    :: alpha
    REAL(KIND = dp), DIMENSION(3) :: ra, center
    REAL(KIND = dp) :: pos_part, spread_part

    NULLIFY(atomic_kind,atom_list)
    core_eshift=0.0_dp
    pos_part=0.0_dp
    spread_part=0.0_dp
    center=0.0_dp
    DO ikind=1,SIZE(atomic_kind_set)
      atomic_kind => atomic_kind_set(ikind)
      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           natom=natom,&
                           atom_list=atom_list,&
                           alpha_core_charge=alpha,&
                           zeff=zeff)
      DO iatom=1,natom
          atom_a = atom_list(iatom)
          ra(:) = pbc(particle_set(atom_a)%r,cell)
          pos_part=pos_part+Zeff*DOT_PRODUCT(ra,ra)
          spread_part=spread_part+3*Zeff/(2*alpha)
          center=center+Zeff*ra
      ENDDO
    ENDDO
    !write(6,*) "spread part ",spread_part
    !write(6,*) "pos part ",pos_part
    !write(6,*) "center ",center
    core_eshift=pos_part+spread_part
  END SUBROUTINE core_part

  FUNCTION get_value(rho_g,g) result(RES)
    TYPE(coeff_type)                      :: rho_g
    INTEGER :: g(3)
    COMPLEX(dp) :: RES
   
    INTEGER :: i 
    RES=(0.0_dp,0.0_dp)
    DO i=1,rho_g%pw%pw_grid%ngpts_local
       IF (ALL(rho_g%pw%pw_grid%g_hat(:,i).EQ.g)) THEN
           RES=rho_g%pw%cc(i)
           EXIT
       ENDIF
    ENDDO
  END FUNCTION get_value
 
  SUBROUTINE gspace_approach(rho_g,pw_env,cell)
    TYPE(coeff_type)                      :: rho_g
    TYPE(pw_env_type), POINTER            :: pw_env
    TYPE(cell_type), POINTER              :: cell

    INTEGER :: i
    INTEGER :: g(3)
    REAL(KIND=dp) :: resp,resm,rest
    rest=(0.0_8,0.0_8)
    DO i=1,3
       g=0
       g(i)=1
       write(6,*) "g_",i,get_value(rho_g,g)
       resp=ABS(get_value(rho_g,g))
       g(i)=-1
       write(6,*) "g_",i,get_value(rho_g,g)
       resm=ABS(get_value(rho_g,g))
       rest=rest+cell%deth*(resp+resm)/(twopi/cell%hmat(i,i))**2
    ENDDO
    write(6,*) "gspace number",rest
    
  END SUBROUTINE

  SUBROUTINE realspace_approach(rho_r,pw_env,cell,center,power,ionode)
    TYPE(coeff_type)                      :: rho_r
    TYPE(pw_env_type), POINTER            :: pw_env
    logical :: ionode
    REAL(KIND = dp)                         :: center(3)
    TYPE(cell_type), POINTER              :: cell
    INTEGER :: power

    TYPE(realspace_grid_type), POINTER    :: rs_rho
    TYPE(cp_rs_pool_type), POINTER        :: auxbas_rs_pool

    REAL(KIND = dp)   :: monopole,dipole(3),quadrupole(3,3),r(3),boxl(3)

    INTEGER :: i,j,k,m,n,shift
    INTEGER, DIMENSION(3) :: gvec

    NULLIFY(rs_rho,auxbas_rs_pool)
    CALL pw_env_get(pw_env, auxbas_rs_pool=auxbas_rs_pool)
    CALL rs_pool_create_rs(auxbas_rs_pool,rs_rho)
    CALL rs_pw_transfer(rs_rho,rho_r%pw,"BACKWARD")

    ! could all be done in parallel using the pw_grids I guess
    monopole    =0.0_dp
    dipole    =0.0_dp
    quadrupole=0.0_dp
    boxl(1)=(UBOUND(rs_rho%r,1)-LBOUND(rs_rho%r,1)+1)*rs_rho%dr(1)
    boxl(2)=(UBOUND(rs_rho%r,2)-LBOUND(rs_rho%r,2)+1)*rs_rho%dr(2)
    boxl(3)=(UBOUND(rs_rho%r,3)-LBOUND(rs_rho%r,3)+1)*rs_rho%dr(3)
    DO k=LBOUND(rs_rho%r,3),UBOUND(rs_rho%r,3)
       r(3)=(k-LBOUND(rs_rho%r,3))*rs_rho%dr(3)
    DO j=LBOUND(rs_rho%r,2),UBOUND(rs_rho%r,2)
       r(2)=(j-LBOUND(rs_rho%r,2))*rs_rho%dr(2)
    DO i=LBOUND(rs_rho%r,1),UBOUND(rs_rho%r,1)
       r(1)=(i-LBOUND(rs_rho%r,1))*rs_rho%dr(1)
       r=r-center
       r=pbc(r,cell)
       r=r+center
       monopole    =monopole+rs_rho%r(i,j,k)**power
       DO m=1,3
          dipole(m)    =dipole(m)+rs_rho%r(i,j,k)**power*r(m)
          DO n=1,3
             quadrupole(n,m)=quadrupole(n,m)+rs_rho%r(i,j,k)**power*r(m)*r(n)
          ENDDO
       ENDDO
    ENDDO
    ENDDO
    ENDDO
    monopole=monopole*PRODUCT(rs_rho%dr)
    dipole=dipole*PRODUCT(rs_rho%dr)
    quadrupole=quadrupole*PRODUCT(rs_rho%dr)
    IF (IONODE) THEN
    write(6,'(A)') "origin centered moments"
    write(6,'(A,F16.8)') "zero", monopole
    write(6,'(A,3F16.8)') "one", dipole
    write(6,'(A,3F16.8)') "two", quadrupole(1,1:3)
    write(6,'(A,3F16.8)') "   ", quadrupole(2,1:3)
    write(6,'(A,3F16.8)') "   ", quadrupole(3,1:3)
    write(6,'(A,F16.8)') "two trace", quadrupole(1,1)+quadrupole(2,2)+quadrupole(3,3)
    ENDIF
  END SUBROUTINE realspace_approach
END MODULE level_shift
