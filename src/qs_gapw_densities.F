!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

!!****h* cp2k/qs_gapw_densities [1.0] *
!!
!!   NAME
!!     qs_gapw_densities
!!
!!   FUNCTION
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     
!!
!!
!!   SOURCE
!****************************************************************************
MODULE qs_gapw_densities
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_error_handling,               ONLY: cp_error_type
  USE kinds,                           ONLY: dp
  USE memory_utilities,                ONLY: reallocate
  USE qs_charges_types,                ONLY: qs_charges_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_grid_atom,                    ONLY: grid_atom_type
  USE qs_harmonics_atom,               ONLY: harmonics_atom_type
  USE qs_p_env_types,                  ONLY: qs_p_env_type
  USE qs_rho0_ggrid,                   ONLY: put_rho0_s_on_grid
  USE qs_rho0_methods,                 ONLY: calculate_rho0_atom
  USE qs_rho0_types,                   ONLY: rho0_mpole_type, &
                                             rho0_atom_type
  USE qs_rho_atom_methods,             ONLY: calculate_rho_atom
  USE qs_rho_atom_types,               ONLY: rho_atom_type
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PRIVATE, PARAMETER :: module_name='qs_gapw_densities'

  PUBLIC :: prepare_gapw_den 


!******************************************************************************


CONTAINS

!******************************************************************************

 SUBROUTINE prepare_gapw_den(qs_env,tddft,p_env,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(IN), OPTIONAL            :: tddft
    TYPE(qs_p_env_type), POINTER, OPTIONAL   :: p_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER         :: routine_name = "prepare_gapw_den"

    TYPE(atomic_kind_type), DIMENSION(:), &
                                 POINTER     :: atomic_kind_set
    TYPE(atomic_kind_type),      POINTER     :: atom_kind
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(qs_charges_type), POINTER           :: qs_charges
    TYPE(grid_atom_type),      POINTER       :: grid_atom
    TYPE(rho0_mpole_type),      POINTER      :: rho0_mpole
    INTEGER                                  :: handle,ispin, ikind, natom, &
                                                nkind, nspins
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: paw_atom, my_tddft
    TYPE(harmonics_atom_type), POINTER       :: harmonics
    TYPE(rho_atom_type), DIMENSION(:), &
      POINTER                                :: rho_atom_set
    TYPE(rho0_atom_type), DIMENSION(:),&
                           POINTER           :: rho0_atom_set
    REAL(dp)                                 :: rho0_h_tot, rho0_s_tot,&
                                                tot_rs_int
    REAL(dp), DIMENSION(:), POINTER          :: rho1_h_tot,rho1_s_tot
!******************************************************************************

    CALL timeset("prepare_gapw_den","I"," ",handle)

    NULLIFY(atomic_kind_set)
    NULLIFY(atom_kind)
    NULLIFY(dft_control)
    NULLIFY(atom_list)
    NULLIFY(grid_atom)
    NULLIFY(rho0_mpole)
    NULLIFY(qs_charges)
    NULLIFY(rho1_h_tot,rho1_s_tot)
    NULLIFY(rho_atom_set)
    NULLIFY(rho0_atom_set)

    my_tddft = .FALSE.
    IF (PRESENT(tddft)) my_tddft = tddft

    CALL get_qs_env(qs_env=qs_env,dft_control=dft_control,&
                    qs_charges=qs_charges,&
                    atomic_kind_set=atomic_kind_set,&
                    rho0_mpole=rho0_mpole,&
                    rho_atom_set=rho_atom_set,&
                    rho0_atom_set=rho0_atom_set)

    IF (my_tddft) THEN
       rho0_mpole   => p_env%local_rho_set%rho0_mpole
       rho_atom_set => p_env%local_rho_set%rho_atom_set
       rho0_atom_set => p_env%local_rho_set%rho0_atom_set
    END IF

    nkind = SIZE(atomic_kind_set)
    nspins = dft_control%nspins

    rho0_h_tot = 0.0_dp
    rho0_s_tot = 0.0_dp
    CALL reallocate(rho1_h_tot,1,nspins)
    CALL reallocate(rho1_s_tot,1,nspins)

    DO ikind = 1,nkind
      atom_kind => atomic_kind_set(ikind)
      CALL get_atomic_kind(atomic_kind=atom_kind, atom_list=atom_list,&
                           grid_atom=grid_atom,natom=natom,&
                           paw_atom=paw_atom, harmonics=harmonics)

!     Calculate rho1_h and rho1_s on the radial gris centered on the atomic position
      IF(paw_atom) &
        CALL calculate_rho_atom(rho_atom_set,atom_kind,atom_list,grid_atom,natom,&
                                nspins,rho1_h_tot,rho1_s_tot)

!     Calculate rho0_h and rho0_s on the radial gris centered on the atomic position
      CALL calculate_rho0_atom(rho_atom_set,rho0_atom_set,rho0_mpole,atom_list,grid_atom,&
                               paw_atom,natom,ikind,harmonics,rho0_h_tot,rho0_s_tot)

    ENDDO

!dbg
!  stop 'qlm'
    DO ispin = 1,nspins
      qs_charges%total_rho1_hard(ispin) = - rho1_h_tot(ispin)
      qs_charges%total_rho1_soft(ispin) = - rho1_s_tot(ispin)
    END DO 

    rho0_mpole%total_rho0_h = -rho0_h_tot 
    rho0_mpole%total_rho0_s = -rho0_s_tot

!   Put the rho0_soft on the global grid
    CALL put_rho0_s_on_grid(qs_env,atomic_kind_set,rho0_mpole,tot_rs_int,error)

    WRITE(*,*)  -rho0_s_tot,  tot_rs_int,tot_rs_int/rho0_s_tot 

    IF(ABS(1.0_dp-ABS(tot_rs_int/rho0_s_tot)).GT.1.0E-6_dp) THEN
        WRITE(6,'(/,72("*"))')
        WRITE(6,'(T2,A,T66,1E20.8)') & 
          "WARNING: rho0_soft calculated on the local grid is  :", -rho0_s_tot,&
          "         rho0_soft calculated on the global grid is :", tot_rs_int
        WRITE(6,'(T2,A)') &
          "         bad integration"
        WRITE(6,'(72("*"),/)')    
        CALL stop_program(routine_name,"check the Qlm coefficients of rho collocation")
   END IF

   qs_charges%total_rho0_soft_rspace = tot_rs_int
   qs_charges%total_rho0_hard_lebedev = rho0_h_tot

   CALL timestop(0.0_dp,handle)

 END SUBROUTINE  prepare_gapw_den
!******************************************************************************

!******************************************************************************
 END MODULE  qs_gapw_densities 
!******************************************************************************
