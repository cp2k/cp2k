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
  USE kinds,                           ONLY: wp => dp
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
    REAL(wp), DIMENSION(:,:), POINTER             :: centers
    TYPE(cp_error_type), OPTIONAL, INTENT(INOUT)  :: error

    TYPE(coeff_type)            :: rho_tot_rspace
    TYPE(pw_pool_type), POINTER :: auxbas_pw_pool
    INTEGER :: ispin, handle, power
    REAL(wp) :: origin(3),electronic_eshift,core_eshift

    CALL timeset("level_shift_testing","I","",handle)

    CALL wannier_part(centers,electronic_eshift)
    CALL core_part(atomic_kind_set, particle_set, cell, core_eshift)

    write(6,*) "electronic_eshift       ",electronic_eshift
    write(6,*) "core_eshift             ",core_eshift
    write(6,*) "total eshift            ",core_eshift+electronic_eshift
    write(6,*) "eval shift in a.u. ", &
                   (core_eshift+electronic_eshift)*twopi / ( 3.0_wp * cell%deth )

    ! realspace based approach 
    ! generate total density in real space
    NULLIFY(auxbas_pw_pool)
    CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool) 
    CALL pw_pool_init_coeff(auxbas_pw_pool,rho_tot_rspace,&
            use_data = REALDATA3D, in_space = REALSPACE, error=error )
    CALL coeff_transform_space(rho_core,rho_tot_rspace)
    CALL coeff_sumup(rho%rho_r(1),rho_tot_rspace)
    DO ispin=2,SIZE(rho%rho_r)
       CALL coeff_sumup(rho%rho_r(ispin),rho_tot_rspace)
    END DO
    ! start analysis
    origin=0
    power=1
    CALL  realspace_approach(rho_tot_rspace,pw_env,cell,origin,power,.TRUE.)
    ! give back coeffs and pools
    CALL pw_pool_give_back_coeff(auxbas_pw_pool,rho_tot_rspace)

    CALL timestop(0.0_wp,handle)
 
  END SUBROUTINE level_shift_testing

  SUBROUTINE wannier_part(centers,electronic_eshift)
    REAL(wp), DIMENSION(:,:), POINTER :: centers
    REAL(wp) ::electronic_eshift

    INTEGER :: i
    electronic_eshift=0.0_wp
    DO i=1,SIZE(centers,1)
       electronic_eshift=electronic_eshift&
                         -2.0_wp*SUM(centers(i,4:6)**2) ! NOT FOR LSD
    ENDDO
    write(6,*) "spread part ",electronic_eshift
    electronic_eshift=0.0_wp
    DO i=1,SIZE(centers,1)
       electronic_eshift=electronic_eshift&
                         -2.0_wp*SUM(centers(i,1:3)**2) ! NOT FOR LSD
    ENDDO
    write(6,*) "center ",(-2.0_wp*SUM(centers(:,i)),i=1,3)

    electronic_eshift=0.0_wp
    DO i=1,SIZE(centers,1)
       electronic_eshift=electronic_eshift&
                         -2.0_wp*SUM(centers(i,1:3)**2+centers(i,4:6)**2) ! NOT FOR LSD
    ENDDO
  END SUBROUTINE wannier_part

  SUBROUTINE core_part(atomic_kind_set, particle_set, cell, core_eshift)
    REAL(wp) :: core_eshift
    TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), POINTER    :: particle_set
    TYPE(cell_type), POINTER        :: cell

    TYPE(atomic_kind_type), POINTER :: atomic_kind
    INTEGER, DIMENSION(:), POINTER  :: atom_list
    INTEGER :: ikind,natom,iatom,atom_a,zeff
    REAL(wp)    :: alpha
    REAL(wp), DIMENSION(3) :: ra, center
    REAL(wp) :: pos_part, spread_part

    NULLIFY(atomic_kind,atom_list)
    core_eshift=0.0_wp
    pos_part=0.0_wp
    spread_part=0.0_wp
    center=0.0_wp
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
    write(6,*) "spread part ",spread_part
    write(6,*) "pos part ",pos_part
    write(6,*) "center ",center
    core_eshift=pos_part+spread_part
  END SUBROUTINE core_part

  SUBROUTINE realspace_approach(rho_r,pw_env,cell,center,power,ionode)
    TYPE(coeff_type)                      :: rho_r
    TYPE(pw_env_type), POINTER            :: pw_env
    logical :: ionode
    REAL(kind=wp)                         :: center(3)
    TYPE(cell_type), POINTER              :: cell
    INTEGER :: power

    TYPE(realspace_grid_type), POINTER    :: rs_rho
    TYPE(cp_rs_pool_type), POINTER        :: auxbas_rs_pool

    REAL(kind=wp)   :: monopole,dipole(3),quadrupole(3,3),r(3),boxl(3)

    INTEGER :: i,j,k,m,n,shift
    INTEGER, DIMENSION(3) :: gvec

    NULLIFY(rs_rho,auxbas_rs_pool)
    CALL pw_env_get(pw_env, auxbas_rs_pool=auxbas_rs_pool)
    CALL rs_pool_create_rs(auxbas_rs_pool,rs_rho)
    CALL rs_pw_transfer(rs_rho,rho_r%pw,"BACKWARD")

    ! could all be done in parallel using the pw_grids I guess
    monopole    =0.0_wp
    dipole    =0.0_wp
    quadrupole=0.0_wp
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
