!-----------------------------------------------------------------------------
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/kg_force_numer [1.0] *
!!
!!   NAME
!!    kg_force_numer
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     gt 
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

 MODULE kg_force_numer
!------------------------------------------------------------------------------!
!
  USE atomic_kinds, ONLY: kind_info_type
  USE coefficient_types, ONLY : coeff_type, coeff_allocate, &
      coeff_deallocate, coeff_transform_space, coeff_sumup, &
      coeff_zero
  USE cube_utils, ONLY : cube_info_type
  USE gaussian_gridlevels, ONLY: gridlevel_info_type
  USE global_types, ONLY : global_environment_type
  USE hartree,                  ONLY: calculate_hartree
  USE l_utils, ONLY : l_info_type
  USE kg_core, ONLY : calculate_ecore_rspace
  USE kg_density, ONLY : calculate_density, calculate_epc_density,&
                     calculate_p_density 
  USE kg_environment_types, ONLY : get_kg_env, kg_environment_type
  USE kg_energy_types, ONLY : kg_energy_type 
!gt
! kg intra temporarily removed from OBJECT file
!    USE kg_intra, ONLY : calculate_ebond_corr
!gt
  USE kg_kxc, ONLY : calculate_kxc_potential
  USE kg_neighbor_lists, ONLY : build_kg_neighbor_lists
  USE kg_ppl, ONLY : calculate_rho0_ppl, calculate_drho_ppl
  USE kinds, ONLY : dbl
  USE pw_grid_types, ONLY : pw_grid_type
  USE particle_types, ONLY : particle_type
  USE pw_pool_types, ONLY : pw_pool_type
  USE pw_types,  ONLY: COMPLEXDATA1D, COMPLEXDATA3D,&
       REALDATA3D, REALSPACE, RECIPROCALSPACE
  USE simulation_cell, ONLY : cell_type
  USE termination, ONLY : stop_memory, stop_program

  PRIVATE
  PUBLIC :: fpart_gspace_numer, fpart_bc_numer, fpart_xc_numer, &
            fpart_core_rspace_numer, fpart_ppl_numer, fpart_p_xc_numer, &
            fpart_p_gspace_numer, fcoef_p_xc_numer, fcoef_p_gspace_numer,&
            fpart_p_ppl_numer, fcoef_p_ppl_numer
!!***
!------------------------------------------------------------------------------!
!
   CONTAINS
!
!------------------------------------------------------------------------------!

 SUBROUTINE  fpart_gspace_numer  ( delta, kg_env, pw_grid, iref, &
                                   gridlevel_info, l_info, cube_info, &
                                   f_part, energy_numer, globenv )
   IMPLICIT NONE
!
! begin passed variable declaration:
!
   TYPE ( cube_info_type ), DIMENSION(:), POINTER :: cube_info
   TYPE ( global_environment_type ), INTENT ( INOUT ) :: globenv
   TYPE ( gridlevel_info_type ), INTENT(IN)       :: gridlevel_info
   TYPE ( kg_environment_type ), POINTER :: kg_env
   TYPE ( l_info_type ), INTENT(IN)               :: l_info
   TYPE ( pw_grid_type ), INTENT ( IN ), DIMENSION ( : ) :: pw_grid
   INTEGER, INTENT ( IN ) :: iref
   REAL ( dbl ), INTENT ( IN ) :: delta
   REAL ( dbl ), INTENT ( OUT ) :: energy_numer
   REAL ( dbl ), INTENT ( OUT ), DIMENSION ( :, : ) :: f_part

! begin local variable declaration:
!
   CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE fpart_gspace_numer (MODULE kg_force_numer)"
   INTEGER :: i, j, ngrid_level, igrid_level  
   INTEGER :: id,  natoms, istat
   REAL (dbl) :: energy
   REAL (dbl) :: energy_plus, energy_minus
   REAL (dbl) :: total_rho0_rspace, total_rho_core_rspace
   TYPE (particle_type), DIMENSION(:), POINTER :: particle_set
   TYPE ( coeff_type ) :: rho_gspace, rho_elec_rspace, rho_core, v_gspace
   TYPE ( coeff_type ), POINTER :: mgrid_gspace (:)
   TYPE ( coeff_type ), POINTER :: mgrid_rspace (:)

   energy_numer = 0.0_dbl
   f_part = 0.0_dbl
   energy_plus = 0.0_dbl
   energy_minus = 0.0_dbl
   NULLIFY(mgrid_gspace, mgrid_rspace)
   NULLIFY (particle_set)
   ngrid_level = SIZE ( pw_grid )

   CALL get_kg_env(kg_env=kg_env,particle_set=particle_set)
! allocate densities on the grids

   CALL coeff_allocate(rho_elec_rspace,pw_grid(iref),REALDATA3D)
   rho_elec_rspace%pw%in_space = REALSPACE
   CALL coeff_allocate(rho_gspace,pw_grid(iref),COMPLEXDATA1D)
   rho_gspace%pw%in_space = RECIPROCALSPACE
   CALL coeff_allocate(rho_core,pw_grid(iref),COMPLEXDATA1D)
   rho_core%pw%in_space = RECIPROCALSPACE
   CALL coeff_allocate(v_gspace,pw_grid(iref),COMPLEXDATA1D)
   v_gspace%pw%in_space = RECIPROCALSPACE

! allocate the multi-grid rho's
   ALLOCATE (mgrid_rspace(ngrid_level),STAT=istat)
   IF (istat /= 0) CALL stop_memory(routine,"mgrid_rspace",0)
   ALLOCATE (mgrid_gspace(ngrid_level),STAT=istat)
   IF (istat /= 0) CALL stop_memory(routine,"mgrid_gspace",0)
   DO igrid_level=1,ngrid_level
     CALL coeff_allocate(mgrid_rspace(igrid_level),pw_grid(igrid_level),&
          REALDATA3D)
     mgrid_rspace(igrid_level)%pw%in_space = REALSPACE
     CALL coeff_allocate(mgrid_gspace(igrid_level),pw_grid(igrid_level),&
          COMPLEXDATA1D)
     mgrid_gspace(igrid_level)%pw%in_space = RECIPROCALSPACE
   END DO

!
! starting the force loop
!
   natoms = SIZE ( particle_set )

   DO i = 1, natoms

     DO id = 1,3

        particle_set ( i ) % r ( id ) = particle_set ( i ) % r ( id ) + delta
        CALL coeff_zero(rho_core)
        CALL coeff_zero(rho_elec_rspace)
        CALL coeff_zero(rho_gspace)
        CALL coeff_zero (v_gspace)
        DO igrid_level=1,ngrid_level
          CALL coeff_zero (mgrid_rspace(igrid_level))
          CALL coeff_zero (mgrid_gspace(igrid_level))
        END DO
! Calculate the density generated by the effective point core charges qeff
        CALL calculate_epc_density(rho_core, cube_info(1),&
                         l_info,&
                         total_rho_core_rspace, kg_env, globenv)
! calculate the frozen electronic pseudo density
        CALL calculate_density ( rho_elec_rspace, rho_gspace, &
            mgrid_rspace, mgrid_gspace,&
            gridlevel_info, cube_info,&
            l_info, total_rho0_rspace, kg_env, globenv )
! get the total gspace density
        CALL coeff_sumup(rho_core,rho_gspace)
        CALL calculate_hartree(rho_gspace,energy_plus,v_gspace)

        particle_set ( i ) % r ( id ) = particle_set ( i ) % r ( id ) - 2._dbl * delta
        CALL coeff_zero(rho_core)
        CALL coeff_zero(rho_elec_rspace)
        CALL coeff_zero(rho_gspace)
        CALL coeff_zero (v_gspace)
        DO igrid_level=1,ngrid_level
          CALL coeff_zero (mgrid_rspace(igrid_level))
          CALL coeff_zero (mgrid_gspace(igrid_level))
        END DO
! Calculate the density generated by the effective point core charges qeff
       CALL calculate_epc_density(rho_core, cube_info(1),&
                         l_info,&
                         total_rho_core_rspace, kg_env, globenv)
! calculate the frozen electronic pseudo density
       CALL calculate_density ( rho_elec_rspace, rho_gspace, &
            mgrid_rspace, mgrid_gspace,&
            gridlevel_info, cube_info,&
            l_info, total_rho0_rspace, kg_env, globenv )
! get the total gspace density
        CALL coeff_sumup(rho_core,rho_gspace)
        CALL calculate_hartree(rho_gspace,energy_minus,v_gspace)

        f_part ( id, i ) = energy_minus - energy_plus

        particle_set ( i ) % r ( id ) = particle_set ( i ) % r ( id ) + delta

     END DO

  END DO 

  f_part = f_part / 2._dbl / delta
! get ener_numer
  CALL coeff_zero(rho_core)
  CALL coeff_zero(rho_elec_rspace)
  CALL coeff_zero(rho_gspace)
  CALL coeff_zero (v_gspace)
  DO igrid_level=1,ngrid_level
    CALL coeff_zero (mgrid_rspace(igrid_level))
    CALL coeff_zero (mgrid_gspace(igrid_level))
  END DO
! Calculate the density generated by the effective point core charges qeff
  CALL calculate_epc_density(rho_core, cube_info(1),&
                   l_info,&
                   total_rho_core_rspace, kg_env, globenv)
! calculate the frozen electronic pseudo density
  CALL calculate_density ( rho_elec_rspace, rho_gspace, &
     mgrid_rspace, mgrid_gspace,&
     gridlevel_info, cube_info,&
     l_info, total_rho0_rspace, kg_env, globenv )
! get the total gspace density
  CALL coeff_sumup(rho_core,rho_gspace)
  CALL calculate_hartree(rho_gspace,energy_numer,v_gspace)

 
! deallocate densities

  CALL coeff_deallocate(rho_elec_rspace)
  CALL coeff_deallocate(rho_gspace)
  CALL coeff_deallocate(rho_core)
  CALL coeff_deallocate(v_gspace)
  DO igrid_level=1,ngrid_level
     CALL coeff_deallocate(mgrid_rspace(igrid_level))
     CALL coeff_deallocate(mgrid_gspace(igrid_level))
  END DO
  DEALLOCATE (mgrid_rspace, STAT=istat)
  IF (istat /= 0) CALL stop_memory(routine,"mgrid_rspace",0)
  DEALLOCATE (mgrid_gspace,STAT=istat)
  IF (istat /= 0) CALL stop_memory(routine,"mgrid_gspace",0)

  RETURN
 
END SUBROUTINE  fpart_gspace_numer 

!------------------------------------------------------------------------------!

 SUBROUTINE  fpart_xc_numer  ( delta, kg_env, pw_grid, pw_pool, iref, &
                                   gridlevel_info, l_info, cube_info, &
                                   f_part, energy_numer, globenv )
   IMPLICIT NONE
!
! begin passed variable declaration:
!
   TYPE ( cube_info_type ), DIMENSION(:), POINTER :: cube_info
   TYPE ( global_environment_type ), INTENT ( INOUT ) :: globenv
   TYPE ( gridlevel_info_type ), INTENT(IN)       :: gridlevel_info
   TYPE ( kg_environment_type ), INTENT ( INOUT ) :: kg_env
   TYPE ( l_info_type ), INTENT(IN)               :: l_info
   TYPE ( pw_grid_type ), INTENT ( IN ), DIMENSION ( : ) :: pw_grid
   TYPE ( pw_pool_type ), pointer :: pw_pool
   INTEGER, INTENT ( IN ) :: iref
   REAL ( dbl ), INTENT ( IN ) :: delta
   REAL ( dbl ), INTENT ( OUT ) :: energy_numer
   REAL ( dbl ), INTENT ( OUT ), DIMENSION ( :, : ) :: f_part

! begin local variable declaration:
!
   CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE fpart_xc_numer (MODULE kg_force_numer)"
   INTEGER :: i, j, ngrid_level, igrid_level  
   INTEGER :: id,  natoms, istat
   REAL (dbl) :: energy
   REAL (dbl) :: energy_plus, energy_minus
   REAL (dbl) :: total_rho0_rspace, total_rho_core_rspace
   TYPE (particle_type), DIMENSION(:), POINTER :: particle_set
   TYPE ( coeff_type ) :: rho_gspace, rho_elec_rspace, v_xc_rspace
   TYPE ( coeff_type ), POINTER :: mgrid_gspace (:)
   TYPE ( coeff_type ), POINTER :: mgrid_rspace (:)

   energy_numer = 0.0_dbl
   energy_plus = 0.0_dbl
   energy_minus = 0.0_dbl
   f_part = 0.0_dbl
   NULLIFY(mgrid_gspace, mgrid_rspace)
   NULLIFY (particle_set)
   ngrid_level = SIZE ( pw_grid )
   CALL get_kg_env(kg_env=kg_env,particle_set=particle_set)

! allocate densities on the grids

   CALL coeff_allocate(rho_elec_rspace,pw_grid(iref),REALDATA3D)
   rho_elec_rspace%pw%in_space = REALSPACE
   CALL coeff_allocate(rho_gspace,pw_grid(iref),COMPLEXDATA1D)
   rho_gspace%pw%in_space = RECIPROCALSPACE
   CALL coeff_allocate(v_xc_rspace,pw_grid(iref),REALDATA3D)
   v_xc_rspace%pw%in_space = REALSPACE

! allocate the multi-grid rho's
   ALLOCATE (mgrid_rspace(ngrid_level),STAT=istat)
   IF (istat /= 0) CALL stop_memory(routine,"mgrid_rspace",0)
   ALLOCATE (mgrid_gspace(ngrid_level),STAT=istat)
   IF (istat /= 0) CALL stop_memory(routine,"mgrid_gspace",0)
   DO igrid_level=1,ngrid_level
     CALL coeff_allocate(mgrid_rspace(igrid_level),pw_grid(igrid_level),&
          REALDATA3D)
     mgrid_rspace(igrid_level)%pw%in_space = REALSPACE
     CALL coeff_allocate(mgrid_gspace(igrid_level),pw_grid(igrid_level),&
          COMPLEXDATA1D)
     mgrid_gspace(igrid_level)%pw%in_space = RECIPROCALSPACE
   END DO

!
! starting the force loop
!
   natoms = SIZE ( particle_set )

   DO i = 1, natoms

     DO id = 1,3

        particle_set ( i ) % r ( id ) = particle_set ( i ) % r ( id ) + delta
        CALL coeff_zero(rho_elec_rspace)
        CALL coeff_zero(rho_gspace)
        CALL coeff_zero (v_xc_rspace)
        DO igrid_level=1,ngrid_level
          CALL coeff_zero (mgrid_rspace(igrid_level))
          CALL coeff_zero (mgrid_gspace(igrid_level))
        END DO
! calculate the frozen electronic pseudo density
        CALL calculate_density ( rho_elec_rspace, rho_gspace, &
            mgrid_rspace, mgrid_gspace,&
            gridlevel_info, cube_info,&
            l_info, total_rho0_rspace, kg_env, globenv )
        CALL calculate_kxc_potential(pw_pool, rho_elec_rspace, rho_gspace, v_xc_rspace, &
                            energy_plus, kg_env%dft_control)

        particle_set ( i ) % r ( id ) = particle_set ( i ) % r ( id ) - 2._dbl * delta
        CALL coeff_zero(rho_elec_rspace)
        CALL coeff_zero(rho_gspace)
        CALL coeff_zero (v_xc_rspace)
        DO igrid_level=1,ngrid_level
          CALL coeff_zero (mgrid_rspace(igrid_level))
          CALL coeff_zero (mgrid_gspace(igrid_level))
        END DO
! calculate the frozen electronic pseudo density
        CALL calculate_density ( rho_elec_rspace, rho_gspace, &
            mgrid_rspace, mgrid_gspace,&
            gridlevel_info, cube_info,&
            l_info, total_rho0_rspace, kg_env, globenv )
        CALL calculate_kxc_potential(pw_pool,rho_elec_rspace, rho_gspace, v_xc_rspace, &
                            energy_minus, kg_env%dft_control)
        f_part ( id, i ) = energy_minus - energy_plus

        particle_set ( i ) % r ( id ) = particle_set ( i ) % r ( id ) + delta

     END DO

  END DO 

  f_part = f_part / 2._dbl / delta
! get ener_numer
  CALL coeff_zero(rho_elec_rspace)
  CALL coeff_zero(rho_gspace)
  CALL coeff_zero (v_xc_rspace)
  DO igrid_level=1,ngrid_level
    CALL coeff_zero (mgrid_rspace(igrid_level))
    CALL coeff_zero (mgrid_gspace(igrid_level))
  END DO
! Calculate the density generated by the frozen electronic denisty         
  CALL calculate_density ( rho_elec_rspace, rho_gspace, &
                           mgrid_rspace, mgrid_gspace,&
                           gridlevel_info, cube_info,&
                           l_info, total_rho0_rspace, kg_env, globenv )
  CALL calculate_kxc_potential(pw_pool,rho_elec_rspace, rho_gspace, v_xc_rspace, &
                           energy_numer, kg_env%dft_control)


! deallocate densities

  CALL coeff_deallocate(rho_elec_rspace)
  CALL coeff_deallocate(rho_gspace)
  CALL coeff_deallocate(v_xc_rspace)
  DO igrid_level=1,ngrid_level
     CALL coeff_deallocate(mgrid_rspace(igrid_level))
     CALL coeff_deallocate(mgrid_gspace(igrid_level))
  END DO
  DEALLOCATE (mgrid_rspace, STAT=istat)
  IF (istat /= 0) CALL stop_memory(routine,"mgrid_rspace",0)
  DEALLOCATE (mgrid_gspace,STAT=istat)
  IF (istat /= 0) CALL stop_memory(routine,"mgrid_gspace",0)

  RETURN
 
END SUBROUTINE  fpart_xc_numer 
!------------------------------------------------------------------------------!

 SUBROUTINE  fpart_core_rspace_numer  ( delta, kg_env,  &
                                   f_part, energy_numer, globenv )
   IMPLICIT NONE
!
! begin passed variable declaration:
!
   TYPE ( global_environment_type ), INTENT ( INOUT ) :: globenv
   TYPE ( kg_environment_type ), INTENT ( INOUT ) :: kg_env
   REAL ( dbl ), INTENT ( IN ) :: delta
   REAL ( dbl ), INTENT ( OUT ) :: energy_numer
   REAL ( dbl ), INTENT ( OUT ), DIMENSION ( :, : ) :: f_part

! begin local variable declaration:
!
   CHARACTER(LEN=*), PARAMETER :: routine =&
     "SUBROUTINE fpart_core_rspace_numer (MODULE kg_force_numer)"
   INTEGER :: i, j
   INTEGER :: id,  natoms, istat
   TYPE (particle_type), DIMENSION(:), POINTER :: particle_set
   REAL (dbl) :: energy
   REAL (dbl) :: energy_plus, energy_minus

   energy_numer = 0.0_dbl
   energy_plus = 0.0_dbl
   energy_minus = 0.0_dbl
   f_part = 0.0_dbl
   NULLIFY (particle_set)
   CALL get_kg_env(kg_env=kg_env,particle_set=particle_set)
!
! starting the force loop
!
   natoms = SIZE ( particle_set )

   DO i = 1, natoms

     DO id = 1,3

       particle_set ( i ) % r ( id ) = particle_set ( i ) % r ( id ) + delta
       CALL calculate_ecore_rspace(kg_env,globenv,energy_plus)

       particle_set ( i ) % r ( id ) = particle_set ( i ) % r ( id ) - 2._dbl * delta
       CALL calculate_ecore_rspace(kg_env,globenv,energy_minus)

       f_part ( id, i ) = energy_minus - energy_plus

       particle_set ( i ) % r ( id ) = particle_set ( i ) % r ( id ) + delta

     END DO

  END DO

  f_part = f_part / 2._dbl / delta
  CALL calculate_ecore_rspace(kg_env,globenv,energy_numer)

  RETURN

END SUBROUTINE  fpart_core_rspace_numer


!------------------------------------------------------------------------------!

 SUBROUTINE  fpart_ppl_numer  ( delta, kg_env,  &
                                   f_part, energy_numer, globenv )
   IMPLICIT NONE
!
! begin passed variable declaration:
!
   TYPE ( global_environment_type ), INTENT ( INOUT ) :: globenv
   TYPE ( kg_environment_type ), POINTER          :: kg_env
   REAL ( dbl ), INTENT ( IN ) :: delta
   REAL ( dbl ), INTENT ( OUT ) :: energy_numer
   REAL ( dbl ), INTENT ( OUT ), DIMENSION ( :, : ) :: f_part

! begin local variable declaration:
!
   CHARACTER(LEN=*), PARAMETER :: routine =&
     "SUBROUTINE fpart_ppl_numer (MODULE kg_force_numer)"
   INTEGER :: i, j
   INTEGER :: id,  natoms, istat
   REAL (dbl) :: energy_plus, energy_minus
   TYPE (particle_type), DIMENSION(:), POINTER :: particle_set
   TYPE(kg_energy_type), POINTER :: energy
   LOGICAL :: calculate_forces

   energy_numer = 0.0_dbl
   energy_plus = 0.0_dbl
   energy_minus = 0.0_dbl
   f_part = 0.0_dbl
   calculate_forces = .FALSE.
   NULLIFY (particle_set,energy)
   CALL get_kg_env(kg_env=kg_env,particle_set=particle_set)
!
! starting the force loop
!
   natoms = SIZE ( particle_set )

   DO i = 1, natoms

     DO id = 1,3

       particle_set ( i ) % r ( id ) = particle_set ( i ) % r ( id ) + delta
       CALL build_kg_neighbor_lists(kg_env,globenv)
       CALL calculate_rho0_ppl(kg_env,globenv,calculate_forces)
       CALL get_kg_env(kg_env=kg_env,energy=energy)
       energy_plus=energy%pseudo

       particle_set ( i ) % r ( id ) = particle_set ( i ) % r ( id ) - 2._dbl * delta
       CALL build_kg_neighbor_lists(kg_env,globenv)
       CALL calculate_rho0_ppl(kg_env,globenv,calculate_forces)
       CALL get_kg_env(kg_env=kg_env,energy=energy)
       energy_minus=energy%pseudo

       f_part ( id, i ) = energy_minus - energy_plus

       particle_set ( i ) % r ( id ) = particle_set ( i ) % r ( id ) + delta

     END DO

  END DO

  f_part = f_part / 2._dbl / delta
  CALL build_kg_neighbor_lists(kg_env,globenv)
  CALL calculate_rho0_ppl(kg_env,globenv,calculate_forces)
  CALL get_kg_env(kg_env=kg_env,energy=energy)
  energy_numer=energy%pseudo

  RETURN

END SUBROUTINE  fpart_ppl_numer

!------------------------------------------------------------------------------!

 SUBROUTINE  fpart_p_ppl_numer  ( delta, kg_env,  &
                                   f_part, energy_numer, globenv )
   IMPLICIT NONE
!
! begin passed variable declaration:
!
   TYPE ( global_environment_type ), INTENT ( INOUT ) :: globenv
   TYPE ( kg_environment_type ), POINTER          :: kg_env
   REAL ( dbl ), INTENT ( IN ) :: delta
   REAL ( dbl ), INTENT ( OUT ) :: energy_numer
   REAL ( dbl ), INTENT ( OUT ), DIMENSION ( :, : ) :: f_part

! begin local variable declaration:
!
   CHARACTER(LEN=*), PARAMETER :: routine =&
     "SUBROUTINE fpart_ppl_numer (MODULE kg_force_numer)"
   INTEGER :: i, j
   INTEGER :: id,  natoms, istat
   REAL (dbl) :: energy2
   REAL (dbl) :: energy_plus, energy_minus
   LOGICAL :: calculate_forces
   TYPE (kg_energy_type), POINTER :: energy
   TYPE (particle_type), DIMENSION(:), POINTER :: particle_set
   TYPE ( coeff_type ), POINTER :: coef_pos
   TYPE ( coeff_type ), POINTER :: coef_vel
   TYPE ( coeff_type ), POINTER :: coef_force

!   coef_pos => kg_env % coef_pos (1)
!   coef_vel => kg_env % coef_vel (1)
!   coef_force => kg_env % coef_force (1)

   NULLIFY (particle_set,energy)
   energy_numer = 0.0_dbl
   energy_plus = 0.0_dbl
   energy_minus = 0.0_dbl
   f_part = 0.0_dbl
   calculate_forces = .FALSE.
   CALL get_kg_env(kg_env=kg_env,particle_set=particle_set)
!
! starting the force loop
!
   natoms = SIZE ( particle_set )

   DO i = 1, natoms

     DO id = 1,3

       particle_set ( i ) % r ( id ) = particle_set ( i ) % r ( id ) + delta
       CALL build_kg_neighbor_lists(kg_env,globenv)
       CALL calculate_rho0_ppl(kg_env,globenv,calculate_forces)
       CALL get_kg_env(kg_env=kg_env,energy=energy)
       CALL calculate_drho_ppl(coef_pos%ao%cr,coef_force%ao%cr,&
            kg_env,globenv,energy2,calculate_forces)
       energy_plus=energy%pseudo+energy2
       

       particle_set ( i ) % r ( id ) = particle_set ( i ) % r ( id ) - 2._dbl * delta
       CALL build_kg_neighbor_lists(kg_env,globenv)
       CALL calculate_rho0_ppl(kg_env,globenv,calculate_forces)
       CALL get_kg_env(kg_env=kg_env,energy=energy)
       CALL calculate_drho_ppl(coef_pos%ao%cr,coef_force%ao%cr,&
                                kg_env,globenv,energy2,calculate_forces)
       energy_minus=energy%pseudo+energy2

       f_part ( id, i ) = energy_minus - energy_plus

       particle_set ( i ) % r ( id ) = particle_set ( i ) % r ( id ) + delta

     END DO

  END DO

  f_part = f_part / 2._dbl / delta
  CALL build_kg_neighbor_lists(kg_env,globenv)
  CALL calculate_rho0_ppl(kg_env,globenv,calculate_forces)
  CALL get_kg_env(kg_env=kg_env,energy=energy)
  CALL calculate_drho_ppl(coef_pos%ao%cr,coef_force%ao%cr,&
                        kg_env,globenv,energy2,calculate_forces)
  energy_numer=energy%pseudo +energy2

  RETURN

END SUBROUTINE  fpart_p_ppl_numer


!------------------------------------------------------------------------------!

 SUBROUTINE  fpart_bc_numer  ( delta, kg_env,  &
                                   f_part, energy_numer, globenv )
   IMPLICIT NONE
!
! begin passed variable declaration:
!
   TYPE ( global_environment_type ), INTENT ( INOUT ) :: globenv
   TYPE ( kg_environment_type ), INTENT ( INOUT ) :: kg_env
   REAL ( dbl ), INTENT ( IN ) :: delta
   REAL ( dbl ), INTENT ( OUT ) :: energy_numer
   REAL ( dbl ), INTENT ( OUT ), DIMENSION ( :, : ) :: f_part

! begin local variable declaration:
!
   CHARACTER(LEN=*), PARAMETER :: routine =&
     "SUBROUTINE fpart_bc_numer (MODULE kg_force_numer)"
   INTEGER :: i, j
   INTEGER :: id,  natoms, istat
   TYPE (particle_type), DIMENSION(:), POINTER :: particle_set
   REAL (dbl) :: energy
   REAL (dbl) :: energy_plus, energy_minus


   NULLIFY (particle_set)
   CALL get_kg_env(kg_env=kg_env,particle_set=particle_set)
   energy_numer = 0.0_dbl
   energy_plus = 0.0_dbl
   energy_minus = 0.0_dbl
   f_part = 0.0_dbl
!
! starting the force loop
!
   natoms = SIZE ( particle_set )

   DO i = 1, natoms

     DO id = 1,3

       particle_set ( i ) % r ( id ) = particle_set ( i ) % r ( id ) + delta
!       CALL calculate_ebond_corr(kg_env,globenv,energy_plus)

       particle_set ( i ) % r ( id ) = particle_set ( i ) % r ( id ) - 2._dbl * delta
!       CALL calculate_ebond_corr(kg_env,globenv,energy_minus)

       f_part ( id, i ) = energy_minus - energy_plus

       particle_set ( i ) % r ( id ) = particle_set ( i ) % r ( id ) + delta

     END DO

  END DO

  f_part = f_part / 2._dbl / delta
!  CALL calculate_ebond_corr(kg_env,globenv,energy_numer)

  RETURN

END SUBROUTINE  fpart_bc_numer
!------------------------------------------------------------------------------!

 SUBROUTINE  fpart_p_xc_numer  ( delta, kg_env, pw_grid, pw_pool, iref, &
                                   gridlevel_info, l_info, cube_info, &
                                   f_part, energy_numer, globenv )
   IMPLICIT NONE
!
! begin passed variable declaration:
!
   TYPE ( cube_info_type ), DIMENSION(:), POINTER :: cube_info
   TYPE ( global_environment_type ), INTENT ( INOUT ) :: globenv
   TYPE ( gridlevel_info_type ), INTENT(IN)       :: gridlevel_info
   TYPE ( kg_environment_type ), INTENT ( INOUT ) :: kg_env
   TYPE ( l_info_type ), INTENT(IN)               :: l_info
   TYPE ( pw_grid_type ), INTENT ( IN ), DIMENSION ( : ) :: pw_grid
   TYPE ( pw_pool_type ), pointer :: pw_pool
   INTEGER, INTENT ( IN ) :: iref
   REAL ( dbl ), INTENT ( IN ) :: delta
   REAL ( dbl ), INTENT ( OUT ) :: energy_numer
   REAL ( dbl ), INTENT ( OUT ), DIMENSION ( :, : ) :: f_part

! begin local variable declaration:
!
   CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE fpart_xc_numer (MODULE kg_force_numer)"
   INTEGER :: i, j, ngrid_level, igrid_level  
   INTEGER :: id,  natoms, istat
   REAL (dbl) :: energy, exc, e1_xc, e2_xc
   REAL (dbl) :: energy_plus, energy_minus
   REAL (dbl) :: total_rho0_rspace, total_rho_core_rspace
   REAL (dbl) :: total_rhop_rspace
   TYPE ( coeff_type ) :: rho_gspace, rho_elec_rspace, v_xc_rspace
   TYPE ( coeff_type ) :: rho_p_gspace, rho_p_rspace
   TYPE ( coeff_type ) :: v1_xc_rspace, v2_xc_rspace
   TYPE ( coeff_type ), POINTER :: mgrid_gspace (:)
   TYPE ( coeff_type ), POINTER :: mgrid_rspace (:)
   TYPE (particle_type), DIMENSION(:), POINTER :: particle_set
   TYPE ( coeff_type ), POINTER :: coef_pos
   TYPE ( coeff_type ), POINTER :: coef_vel
   TYPE ( coeff_type ), POINTER :: coef_force

!   coef_pos => kg_env % coef_pos (1)
!   coef_vel => kg_env % coef_vel (1)
!   coef_force => kg_env % coef_force (1)
   energy_numer = 0.0_dbl
   energy_plus = 0.0_dbl
   energy_minus = 0.0_dbl
   f_part = 0.0_dbl
   NULLIFY(mgrid_gspace, mgrid_rspace)
   NULLIFY(particle_set)
   ngrid_level = SIZE ( pw_grid )
   CALL get_kg_env(kg_env=kg_env,particle_set=particle_set)

! allocate densities on the grids

   CALL coeff_allocate(rho_elec_rspace,pw_grid(iref),REALDATA3D)
   rho_elec_rspace%pw%in_space = REALSPACE
   CALL coeff_allocate(rho_p_rspace,pw_grid(iref),REALDATA3D)
   rho_p_rspace%pw%in_space = REALSPACE
   CALL coeff_allocate(rho_gspace,pw_grid(iref),COMPLEXDATA1D)
   rho_gspace%pw%in_space = RECIPROCALSPACE
   CALL coeff_allocate(rho_p_gspace,pw_grid(iref),COMPLEXDATA1D)
   rho_p_gspace%pw%in_space = RECIPROCALSPACE
   CALL coeff_allocate(v_xc_rspace,pw_grid(iref),REALDATA3D)
   v_xc_rspace%pw%in_space = REALSPACE
   CALL coeff_allocate(v1_xc_rspace,pw_grid(iref),REALDATA3D)
   v1_xc_rspace%pw%in_space = REALSPACE
   CALL coeff_allocate(v2_xc_rspace,pw_grid(iref),REALDATA3D)
   v2_xc_rspace%pw%in_space = REALSPACE


! allocate the multi-grid rho's
   ALLOCATE (mgrid_rspace(ngrid_level),STAT=istat)
   IF (istat /= 0) CALL stop_memory(routine,"mgrid_rspace",0)
   ALLOCATE (mgrid_gspace(ngrid_level),STAT=istat)
   IF (istat /= 0) CALL stop_memory(routine,"mgrid_gspace",0)
   DO igrid_level=1,ngrid_level
     CALL coeff_allocate(mgrid_rspace(igrid_level),pw_grid(igrid_level),&
          REALDATA3D)
     mgrid_rspace(igrid_level)%pw%in_space = REALSPACE
     CALL coeff_allocate(mgrid_gspace(igrid_level),pw_grid(igrid_level),&
          COMPLEXDATA1D)
     mgrid_gspace(igrid_level)%pw%in_space = RECIPROCALSPACE
   END DO

!
! starting the force loop
!
   natoms = SIZE ( particle_set )

   DO i = 1, natoms

     DO id = 1,3

        particle_set ( i ) % r ( id ) = particle_set ( i ) % r ( id ) + delta
        CALL coeff_zero(rho_elec_rspace)
        CALL coeff_zero(rho_gspace)
        CALL coeff_zero (v_xc_rspace)
        CALL coeff_zero(rho_p_rspace)
        CALL coeff_zero(rho_p_gspace)
        CALL coeff_zero (v1_xc_rspace)
        CALL coeff_zero (v2_xc_rspace)
        DO igrid_level=1,ngrid_level
          CALL coeff_zero (mgrid_rspace(igrid_level))
          CALL coeff_zero (mgrid_gspace(igrid_level))
        END DO
! calculate the frozen electronic pseudo density
        CALL calculate_density ( rho_elec_rspace, rho_gspace, &
            mgrid_rspace, mgrid_gspace,&
            gridlevel_info, cube_info,&
            l_info, total_rho0_rspace, kg_env, globenv )
! calculate the polarization density
        CALL calculate_p_density ( coef_pos % ao % cr, &
            rho_p_rspace, rho_p_gspace, &
            mgrid_rspace, mgrid_gspace,&
            gridlevel_info, cube_info,&
            l_info, total_rhop_rspace, kg_env, globenv )
        CALL calculate_kxc_potential (pw_pool,rho_elec_rspace, rho_gspace, &
                                rho_p_rspace, v_xc_rspace, &
                                v1_xc_rspace, v2_xc_rspace, &
                                exc, e1_xc, e2_xc,kg_env% dft_control)
         energy_plus = exc + e1_xc + e2_xc
!         energy_plus =e2_xc 

        particle_set ( i ) % r ( id ) = particle_set ( i ) % r ( id ) - 2._dbl * delta
        CALL coeff_zero(rho_elec_rspace)
        CALL coeff_zero(rho_gspace)
        CALL coeff_zero (v_xc_rspace)
        CALL coeff_zero(rho_p_rspace)
        CALL coeff_zero(rho_p_gspace)
        CALL coeff_zero (v1_xc_rspace)
        CALL coeff_zero (v2_xc_rspace)
        DO igrid_level=1,ngrid_level
          CALL coeff_zero (mgrid_rspace(igrid_level))
          CALL coeff_zero (mgrid_gspace(igrid_level))
        END DO
! calculate the frozen electronic pseudo density
        CALL calculate_density ( rho_elec_rspace, rho_gspace, &
            mgrid_rspace, mgrid_gspace,&
            gridlevel_info, cube_info,&
            l_info, total_rho0_rspace, kg_env, globenv )
        CALL calculate_p_density ( coef_pos % ao % cr, &
            rho_p_rspace, rho_p_gspace, &
            mgrid_rspace, mgrid_gspace,&
            gridlevel_info, cube_info,&
            l_info, total_rhop_rspace, kg_env, globenv )
        CALL calculate_kxc_potential (pw_pool,rho_elec_rspace, rho_gspace, &
                                rho_p_rspace, v_xc_rspace, &
                                v1_xc_rspace, v2_xc_rspace, &
                                exc, e1_xc, e2_xc,kg_env% dft_control)
         energy_minus = exc + e1_xc + e2_xc
!        energy_minus = e2_xc 

        f_part ( id, i ) = energy_minus - energy_plus

        particle_set ( i ) % r ( id ) = particle_set ( i ) % r ( id ) + delta

     END DO

  END DO 

  f_part = f_part / 2._dbl / delta
! get ener_numer
  CALL coeff_zero(rho_elec_rspace)
  CALL coeff_zero(rho_gspace)
  CALL coeff_zero (v_xc_rspace)
        CALL coeff_zero(rho_p_rspace)
        CALL coeff_zero(rho_p_gspace)
        CALL coeff_zero (v1_xc_rspace)
        CALL coeff_zero (v2_xc_rspace)
  DO igrid_level=1,ngrid_level
    CALL coeff_zero (mgrid_rspace(igrid_level))
    CALL coeff_zero (mgrid_gspace(igrid_level))
  END DO
! Calculate the density generated by the frozen electronic denisty         
  CALL calculate_density ( rho_elec_rspace, rho_gspace, &
                           mgrid_rspace, mgrid_gspace,&
                           gridlevel_info, cube_info,&
                           l_info, total_rho0_rspace, kg_env, globenv )
  CALL calculate_p_density ( coef_pos % ao % cr, &
      rho_p_rspace, rho_p_gspace, &
      mgrid_rspace, mgrid_gspace,&
      gridlevel_info, cube_info,&
      l_info, total_rhop_rspace, kg_env, globenv )
  CALL calculate_kxc_potential (pw_pool,rho_elec_rspace, rho_gspace, &
                          rho_p_rspace, v_xc_rspace, &
                          v1_xc_rspace, v2_xc_rspace, &
                          exc, e1_xc, e2_xc,kg_env% dft_control)
   energy_numer = exc + e1_xc + e2_xc
!  energy_numer = e2_xc 

! deallocate densities

  CALL coeff_deallocate(rho_elec_rspace)
  CALL coeff_deallocate(rho_p_rspace)
  CALL coeff_deallocate(rho_gspace)
  CALL coeff_deallocate(rho_p_gspace)
  CALL coeff_deallocate(v_xc_rspace)
  CALL coeff_deallocate(v1_xc_rspace)
  CALL coeff_deallocate(v2_xc_rspace)
  DO igrid_level=1,ngrid_level
     CALL coeff_deallocate(mgrid_rspace(igrid_level))
     CALL coeff_deallocate(mgrid_gspace(igrid_level))
  END DO
  DEALLOCATE (mgrid_rspace, STAT=istat)
  IF (istat /= 0) CALL stop_memory(routine,"mgrid_rspace",0)
  DEALLOCATE (mgrid_gspace,STAT=istat)
  IF (istat /= 0) CALL stop_memory(routine,"mgrid_gspace",0)

  RETURN
 
END SUBROUTINE  fpart_p_xc_numer 
!------------------------------------------------------------------------------!

 SUBROUTINE  fpart_p_gspace_numer  ( delta, kg_env, pw_grid, iref, &
                                   gridlevel_info, l_info, cube_info, &
                                   f_part, energy_numer, globenv )
   IMPLICIT NONE
!
! begin passed variable declaration:
!
   TYPE ( cube_info_type ), DIMENSION(:), POINTER :: cube_info
   TYPE ( global_environment_type ), INTENT ( INOUT ) :: globenv
   TYPE ( gridlevel_info_type ), INTENT(IN)       :: gridlevel_info
   TYPE ( kg_environment_type ), INTENT ( INOUT ) :: kg_env
   TYPE ( l_info_type ), INTENT(IN)               :: l_info
   TYPE ( pw_grid_type ), INTENT ( IN ), DIMENSION ( : ) :: pw_grid
   INTEGER, INTENT ( IN ) :: iref
   REAL ( dbl ), INTENT ( IN ) :: delta
   REAL ( dbl ), INTENT ( OUT ) :: energy_numer
   REAL ( dbl ), INTENT ( OUT ), DIMENSION ( :, : ) :: f_part

! begin local variable declaration:
!
   CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE fpart_gspace_numer (MODULE kg_force_numer)"
   INTEGER :: i, j, ngrid_level, igrid_level  
   INTEGER :: id,  natoms, istat
   REAL (dbl) :: energy
   REAL (dbl) :: energy_plus, energy_minus
   REAL (dbl) :: total_rho0_rspace, total_rho_core_rspace
   REAL (dbl) :: total_rhop_rspace
   type (particle_type), DIMENSION(:), POINTER :: particle_set
   TYPE ( coeff_type ) :: rho_gspace, rho_elec_rspace, rho_core, v_gspace
   TYPE ( coeff_type ) :: rho_p_gspace, rho_p_rspace
   TYPE ( coeff_type ), POINTER :: mgrid_gspace (:)
   TYPE ( coeff_type ), POINTER :: mgrid_rspace (:)
   TYPE ( coeff_type ), POINTER :: coef_pos
   TYPE ( coeff_type ), POINTER :: coef_vel
   TYPE ( coeff_type ), POINTER :: coef_force

!   coef_pos => kg_env % coef_pos (1)
!   coef_vel => kg_env % coef_vel (1)
!   coef_force => kg_env % coef_force (1)


   energy_numer = 0.0_dbl
   f_part = 0.0_dbl
   energy_plus = 0.0_dbl
   energy_minus = 0.0_dbl
   NULLIFY(mgrid_gspace, mgrid_rspace)
   NULLIFY (particle_set)
   ngrid_level = SIZE ( pw_grid )

   CALL get_kg_env(kg_env=kg_env,particle_set=particle_set)

! allocate densities on the grids

   CALL coeff_allocate(rho_elec_rspace,pw_grid(iref),REALDATA3D)
   rho_elec_rspace%pw%in_space = REALSPACE
   CALL coeff_allocate(rho_p_rspace,pw_grid(iref),REALDATA3D)
   rho_p_rspace%pw%in_space = REALSPACE
   CALL coeff_allocate(rho_gspace,pw_grid(iref),COMPLEXDATA1D)
   rho_gspace%pw%in_space = RECIPROCALSPACE
   CALL coeff_allocate(rho_p_gspace,pw_grid(iref),COMPLEXDATA1D)
   rho_p_gspace%pw%in_space = RECIPROCALSPACE
   CALL coeff_allocate(rho_core,pw_grid(iref),COMPLEXDATA1D)
   rho_core%pw%in_space = RECIPROCALSPACE
   CALL coeff_allocate(v_gspace,pw_grid(iref),COMPLEXDATA1D)
   v_gspace%pw%in_space = RECIPROCALSPACE

! allocate the multi-grid rho's
   ALLOCATE (mgrid_rspace(ngrid_level),STAT=istat)
   IF (istat /= 0) CALL stop_memory(routine,"mgrid_rspace",0)
   ALLOCATE (mgrid_gspace(ngrid_level),STAT=istat)
   IF (istat /= 0) CALL stop_memory(routine,"mgrid_gspace",0)
   DO igrid_level=1,ngrid_level
     CALL coeff_allocate(mgrid_rspace(igrid_level),pw_grid(igrid_level),&
          REALDATA3D)
     mgrid_rspace(igrid_level)%pw%in_space = REALSPACE
     CALL coeff_allocate(mgrid_gspace(igrid_level),pw_grid(igrid_level),&
          COMPLEXDATA1D)
     mgrid_gspace(igrid_level)%pw%in_space = RECIPROCALSPACE
   END DO

!
! starting the force loop
!
   natoms = SIZE ( particle_set )

   DO i = 1, natoms

     DO id = 1,3

        particle_set ( i ) % r ( id ) = particle_set ( i ) % r ( id ) + delta
        CALL coeff_zero(rho_core)
        CALL coeff_zero(rho_elec_rspace)
        CALL coeff_zero(rho_gspace)
        CALL coeff_zero(rho_p_gspace)
        CALL coeff_zero(rho_p_rspace)
        CALL coeff_zero (v_gspace)
        DO igrid_level=1,ngrid_level
          CALL coeff_zero (mgrid_rspace(igrid_level))
          CALL coeff_zero (mgrid_gspace(igrid_level))
        END DO
! Calculate the density generated by the effective point core charges qeff
        CALL calculate_epc_density(rho_core, cube_info(1),&
                         l_info,&
                         total_rho_core_rspace, kg_env, globenv)
! calculate the frozen electronic pseudo density
        CALL calculate_density ( rho_elec_rspace, rho_gspace, &
            mgrid_rspace, mgrid_gspace,&
            gridlevel_info, cube_info,&
            l_info, total_rho0_rspace, kg_env, globenv )
! calculate the polarization density
        CALL calculate_p_density ( coef_pos % ao % cr, &
            rho_p_rspace, rho_p_gspace, &
            mgrid_rspace, mgrid_gspace,&
            gridlevel_info, cube_info,&
            l_info, total_rhop_rspace, kg_env, globenv )
! get the total gspace density
        CALL coeff_sumup(rho_core,rho_gspace)
!dbg gt
!      rho_p_gspace%pw%cc(:)=-1.0_dbl*rho_p_gspace%pw%cc(:)
        CALL coeff_sumup(rho_p_gspace,rho_gspace)
        CALL calculate_hartree(rho_gspace,energy_plus,v_gspace)

        particle_set ( i ) % r ( id ) = particle_set ( i ) % r ( id ) - 2._dbl * delta
        CALL coeff_zero(rho_core)
        CALL coeff_zero(rho_elec_rspace)
        CALL coeff_zero(rho_gspace)
        CALL coeff_zero(rho_p_gspace)
        CALL coeff_zero(rho_p_rspace)
        CALL coeff_zero (v_gspace)
        DO igrid_level=1,ngrid_level
          CALL coeff_zero (mgrid_rspace(igrid_level))
          CALL coeff_zero (mgrid_gspace(igrid_level))
        END DO
! Calculate the density generated by the effective point core charges qeff
       CALL calculate_epc_density(rho_core, cube_info(1),&
                         l_info,&
                         total_rho_core_rspace, kg_env, globenv)
! calculate the frozen electronic pseudo density
       CALL calculate_density ( rho_elec_rspace, rho_gspace, &
            mgrid_rspace, mgrid_gspace,&
            gridlevel_info, cube_info,&
            l_info, total_rho0_rspace, kg_env, globenv )
! calculate the polarization density
        CALL calculate_p_density ( coef_pos % ao % cr, &
            rho_p_rspace, rho_p_gspace, &
            mgrid_rspace, mgrid_gspace,&
            gridlevel_info, cube_info,&
            l_info, total_rhop_rspace, kg_env, globenv )
! get the total gspace density
        CALL coeff_sumup(rho_core,rho_gspace)
!      rho_p_gspace%pw%cc(:)=-1.0_dbl*rho_p_gspace%pw%cc(:)
        CALL coeff_sumup(rho_p_gspace,rho_gspace)
        CALL calculate_hartree(rho_gspace,energy_minus,v_gspace)

        f_part ( id, i ) = energy_minus - energy_plus

        particle_set ( i ) % r ( id ) = particle_set ( i ) % r ( id ) + delta

     END DO

  END DO 

  f_part = f_part / 2._dbl / delta
! get ener_numer
  CALL coeff_zero(rho_core)
  CALL coeff_zero(rho_elec_rspace)
  CALL coeff_zero(rho_gspace)
  CALL coeff_zero(rho_p_gspace)
  CALL coeff_zero(rho_p_rspace)
  CALL coeff_zero (v_gspace)
  DO igrid_level=1,ngrid_level
    CALL coeff_zero (mgrid_rspace(igrid_level))
    CALL coeff_zero (mgrid_gspace(igrid_level))
  END DO
! Calculate the density generated by the effective point core charges qeff
  CALL calculate_epc_density(rho_core, cube_info(1),&
                   l_info,&
                   total_rho_core_rspace, kg_env, globenv)
! calculate the frozen electronic pseudo density
  CALL calculate_density ( rho_elec_rspace, rho_gspace, &
     mgrid_rspace, mgrid_gspace,&
     gridlevel_info, cube_info,&
     l_info, total_rho0_rspace, kg_env, globenv )
! calculate the polarization density
  CALL calculate_p_density ( coef_pos % ao % cr, &
     rho_p_rspace, rho_p_gspace, &
     mgrid_rspace, mgrid_gspace,&
     gridlevel_info, cube_info,&
    l_info, total_rhop_rspace, kg_env, globenv )
! get the total gspace density
  CALL coeff_sumup(rho_core,rho_gspace)
!      rho_p_gspace%pw%cc(:)=-1.0_dbl*rho_p_gspace%pw%cc(:)
  CALL coeff_sumup(rho_p_gspace,rho_gspace)
  CALL calculate_hartree(rho_gspace,energy_numer,v_gspace)

 
! deallocate densities

  CALL coeff_deallocate(rho_elec_rspace)
  CALL coeff_deallocate(rho_gspace)
  CALL coeff_deallocate(rho_p_gspace)
  CALL coeff_deallocate(rho_p_rspace)
  CALL coeff_deallocate(rho_core)
  CALL coeff_deallocate(v_gspace)
  DO igrid_level=1,ngrid_level
     CALL coeff_deallocate(mgrid_rspace(igrid_level))
     CALL coeff_deallocate(mgrid_gspace(igrid_level))
  END DO
  DEALLOCATE (mgrid_rspace, STAT=istat)
  IF (istat /= 0) CALL stop_memory(routine,"mgrid_rspace",0)
  DEALLOCATE (mgrid_gspace,STAT=istat)
  IF (istat /= 0) CALL stop_memory(routine,"mgrid_gspace",0)

  RETURN
 
END SUBROUTINE  fpart_p_gspace_numer 

!------------------------------------------------------------------------------!

 SUBROUTINE  fcoef_p_xc_numer  ( delta, kg_env, pw_grid, pw_pool, iref, &
                                   gridlevel_info, l_info, cube_info, &
                                   f_coef, energy_numer, globenv )
   IMPLICIT NONE
!
! begin passed variable declaration:
!
   TYPE ( cube_info_type ), DIMENSION(:), POINTER :: cube_info
   TYPE ( global_environment_type ), INTENT ( INOUT ) :: globenv
   TYPE ( gridlevel_info_type ), INTENT(IN)       :: gridlevel_info
   TYPE ( kg_environment_type ), INTENT ( INOUT ) :: kg_env
   TYPE ( l_info_type ), INTENT(IN)               :: l_info
   TYPE ( pw_grid_type ), INTENT ( IN ), DIMENSION ( : ) :: pw_grid
   TYPE ( pw_pool_type ), pointer  :: pw_pool
   INTEGER, INTENT ( IN ) :: iref
   REAL ( dbl ), INTENT ( IN ) :: delta
   REAL ( dbl ), INTENT ( OUT ) :: energy_numer
   REAL ( dbl ), INTENT ( OUT ), DIMENSION ( : ) :: f_coef

! begin local variable declaration:
!
   CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE fcoef_p_xc_numer (MODULE kg_force_numer)"
   INTEGER :: i, j, ngrid_level, igrid_level  
   INTEGER :: id,  ncoef, istat
   REAL (dbl) :: energy, exc, e1_xc, e2_xc
   REAL (dbl) :: energy_plus, energy_minus
   REAL (dbl) :: total_rho0_rspace, total_rho_core_rspace
   REAL (dbl) :: total_rhop_rspace
   TYPE ( coeff_type ) :: rho_gspace, rho_elec_rspace, v_xc_rspace
   TYPE ( coeff_type ) :: rho_p_gspace, rho_p_rspace
   TYPE ( coeff_type ) :: v1_xc_rspace, v2_xc_rspace
   TYPE ( coeff_type ), POINTER :: mgrid_gspace (:)
   TYPE ( coeff_type ), POINTER :: mgrid_rspace (:)
   TYPE ( coeff_type ), POINTER :: coef_pos
   TYPE ( coeff_type ), POINTER :: coef_vel
   TYPE ( coeff_type ), POINTER :: coef_force

!   coef_pos => kg_env % coef_pos (1)
!   coef_vel => kg_env % coef_vel (1)
!   coef_force => kg_env % coef_force (1)
   energy_numer = 0.0_dbl
   energy_plus = 0.0_dbl
   energy_minus = 0.0_dbl
   f_coef = 0.0_dbl
   NULLIFY(mgrid_gspace, mgrid_rspace)
   ngrid_level = SIZE ( pw_grid )

! allocate densities on the grids

   CALL coeff_allocate(rho_elec_rspace,pw_grid(iref),REALDATA3D)
   rho_elec_rspace%pw%in_space = REALSPACE
   CALL coeff_allocate(rho_p_rspace,pw_grid(iref),REALDATA3D)
   rho_p_rspace%pw%in_space = REALSPACE
   CALL coeff_allocate(rho_gspace,pw_grid(iref),COMPLEXDATA1D)
   rho_gspace%pw%in_space = RECIPROCALSPACE
   CALL coeff_allocate(rho_p_gspace,pw_grid(iref),COMPLEXDATA1D)
   rho_p_gspace%pw%in_space = RECIPROCALSPACE
   CALL coeff_allocate(v_xc_rspace,pw_grid(iref),REALDATA3D)
   v_xc_rspace%pw%in_space = REALSPACE
   CALL coeff_allocate(v1_xc_rspace,pw_grid(iref),REALDATA3D)
   v1_xc_rspace%pw%in_space = REALSPACE
   CALL coeff_allocate(v2_xc_rspace,pw_grid(iref),REALDATA3D)
   v2_xc_rspace%pw%in_space = REALSPACE


! allocate the multi-grid rho's
   ALLOCATE (mgrid_rspace(ngrid_level),STAT=istat)
   IF (istat /= 0) CALL stop_memory(routine,"mgrid_rspace",0)
   ALLOCATE (mgrid_gspace(ngrid_level),STAT=istat)
   IF (istat /= 0) CALL stop_memory(routine,"mgrid_gspace",0)
   DO igrid_level=1,ngrid_level
     CALL coeff_allocate(mgrid_rspace(igrid_level),pw_grid(igrid_level),&
          REALDATA3D)
     mgrid_rspace(igrid_level)%pw%in_space = REALSPACE
     CALL coeff_allocate(mgrid_gspace(igrid_level),pw_grid(igrid_level),&
          COMPLEXDATA1D)
     mgrid_gspace(igrid_level)%pw%in_space = RECIPROCALSPACE
   END DO

!
! starting the force loop
!
   ncoef = SIZE ( coef_pos%ao%cr )

   DO i = 1, ncoef

      coef_pos % ao % cr ( i ) = coef_pos % ao % cr ( i ) + delta
      CALL coeff_zero(rho_elec_rspace)
      CALL coeff_zero(rho_gspace)
      CALL coeff_zero (v_xc_rspace)
      CALL coeff_zero(rho_p_rspace)
      CALL coeff_zero(rho_p_gspace)
      CALL coeff_zero (v1_xc_rspace)
      CALL coeff_zero (v2_xc_rspace)
      DO igrid_level=1,ngrid_level
        CALL coeff_zero (mgrid_rspace(igrid_level))
        CALL coeff_zero (mgrid_gspace(igrid_level))
      END DO
! calculate the frozen electronic pseudo density
      CALL calculate_density ( rho_elec_rspace, rho_gspace, &
          mgrid_rspace, mgrid_gspace,&
          gridlevel_info, cube_info,&
          l_info, total_rho0_rspace, kg_env, globenv )
! calculate the polarization density
      CALL calculate_p_density ( coef_pos % ao % cr, &
          rho_p_rspace, rho_p_gspace, &
          mgrid_rspace, mgrid_gspace,&
          gridlevel_info, cube_info,&
          l_info, total_rhop_rspace, kg_env, globenv )
      CALL calculate_kxc_potential (pw_pool,rho_elec_rspace, rho_gspace, &
                              rho_p_rspace, v_xc_rspace, &
                              v1_xc_rspace, v2_xc_rspace, &
                              exc, e1_xc, e2_xc,kg_env% dft_control)
       energy_plus = exc + e1_xc + e2_xc
!      energy_plus = e2_xc 

      coef_pos % ao % cr ( i ) = coef_pos % ao % cr ( i ) - 2._dbl * delta
      CALL coeff_zero(rho_elec_rspace)
      CALL coeff_zero(rho_gspace)
      CALL coeff_zero (v_xc_rspace)
      CALL coeff_zero(rho_p_rspace)
      CALL coeff_zero(rho_p_gspace)
      CALL coeff_zero (v1_xc_rspace)
      CALL coeff_zero (v2_xc_rspace)
      DO igrid_level=1,ngrid_level
        CALL coeff_zero (mgrid_rspace(igrid_level))
        CALL coeff_zero (mgrid_gspace(igrid_level))
      END DO
! cculate the frozen electronic pseudo density
      CALL calculate_density ( rho_elec_rspace, rho_gspace, &
          mgrid_rspace, mgrid_gspace,&
          gridlevel_info, cube_info,&
          l_info, total_rho0_rspace, kg_env, globenv )
      CALL calculate_p_density ( coef_pos % ao % cr, &
          rho_p_rspace, rho_p_gspace, &
          mgrid_rspace, mgrid_gspace,&
          gridlevel_info, cube_info,&
          l_info, total_rhop_rspace, kg_env, globenv )
      CALL calculate_kxc_potential (pw_pool,rho_elec_rspace, rho_gspace, &
                              rho_p_rspace, v_xc_rspace, &
                              v1_xc_rspace, v2_xc_rspace, &
                              exc, e1_xc, e2_xc,kg_env% dft_control)
       energy_minus = exc + e1_xc + e2_xc
!      energy_minus = e2_xc 

      f_coef ( i ) = energy_minus - energy_plus

      coef_pos % ao % cr (i) = coef_pos % ao % cr (i) + delta

  END DO 

  f_coef = f_coef / 2._dbl / delta
! get ener_numer
  CALL coeff_zero(rho_elec_rspace)
  CALL coeff_zero(rho_gspace)
  CALL coeff_zero (v_xc_rspace)
        CALL coeff_zero(rho_p_rspace)
        CALL coeff_zero(rho_p_gspace)
        CALL coeff_zero (v1_xc_rspace)
        CALL coeff_zero (v2_xc_rspace)
  DO igrid_level=1,ngrid_level
    CALL coeff_zero (mgrid_rspace(igrid_level))
    CALL coeff_zero (mgrid_gspace(igrid_level))
  END DO
! Calculate the density generated by the frozen electronic denisty         
  CALL calculate_density ( rho_elec_rspace, rho_gspace, &
                           mgrid_rspace, mgrid_gspace,&
                           gridlevel_info, cube_info,&
                           l_info, total_rho0_rspace, kg_env, globenv )
  CALL calculate_p_density ( coef_pos % ao % cr, &
      rho_p_rspace, rho_p_gspace, &
      mgrid_rspace, mgrid_gspace,&
      gridlevel_info, cube_info,&
      l_info, total_rhop_rspace, kg_env, globenv )
  CALL calculate_kxc_potential (pw_pool,rho_elec_rspace, rho_gspace, &
                          rho_p_rspace, v_xc_rspace, &
                          v1_xc_rspace, v2_xc_rspace, &
                          exc, e1_xc, e2_xc,kg_env% dft_control)
  energy_numer = exc + e1_xc + e2_xc         
!   energy_numer = e2_xc 

! deallocate densities

  CALL coeff_deallocate(rho_elec_rspace)
  CALL coeff_deallocate(rho_p_rspace)
  CALL coeff_deallocate(rho_gspace)
  CALL coeff_deallocate(rho_p_gspace)
  CALL coeff_deallocate(v_xc_rspace)
  CALL coeff_deallocate(v1_xc_rspace)
  CALL coeff_deallocate(v2_xc_rspace)
  DO igrid_level=1,ngrid_level
     CALL coeff_deallocate(mgrid_rspace(igrid_level))
     CALL coeff_deallocate(mgrid_gspace(igrid_level))
  END DO
  DEALLOCATE (mgrid_rspace, STAT=istat)
  IF (istat /= 0) CALL stop_memory(routine,"mgrid_rspace",0)
  DEALLOCATE (mgrid_gspace,STAT=istat)
  IF (istat /= 0) CALL stop_memory(routine,"mgrid_gspace",0)

  RETURN
 
END SUBROUTINE  fcoef_p_xc_numer 
!------------------------------------------------------------------------------!

 SUBROUTINE  fcoef_p_gspace_numer  ( delta, kg_env, pw_grid, iref, &
                                   gridlevel_info, l_info, cube_info, &
                                   f_coef, energy_numer, globenv )
   IMPLICIT NONE
!
! begin passed variable declaration:
!
   TYPE ( cube_info_type ), DIMENSION(:), POINTER :: cube_info
   TYPE ( global_environment_type ), INTENT ( INOUT ) :: globenv
   TYPE ( gridlevel_info_type ), INTENT(IN)       :: gridlevel_info
   TYPE ( kg_environment_type ), INTENT ( INOUT ) :: kg_env
   TYPE ( l_info_type ), INTENT(IN)               :: l_info
   TYPE ( pw_grid_type ), INTENT ( IN ), DIMENSION ( : ) :: pw_grid
   INTEGER, INTENT ( IN ) :: iref
   REAL ( dbl ), INTENT ( IN ) :: delta
   REAL ( dbl ), INTENT ( OUT ) :: energy_numer
   REAL ( dbl ), INTENT ( OUT ), DIMENSION ( : ) :: f_coef

! begin local variable declaration:
!
   CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE fcoef_gspace_numer (MODULE kg_force_numer)"
   INTEGER :: i, j, ngrid_level, igrid_level  
   INTEGER :: id,  ncoef , istat
   REAL (dbl) :: energy
   REAL (dbl) :: energy_plus, energy_minus
   REAL (dbl) :: total_rho0_rspace, total_rho_core_rspace
   REAL (dbl) :: total_rhop_rspace
   TYPE ( coeff_type ) :: rho_gspace, rho_elec_rspace, rho_core, v_gspace
   TYPE ( coeff_type ) :: rho_p_gspace, rho_p_rspace
   TYPE ( coeff_type ), POINTER :: mgrid_gspace (:)
   TYPE ( coeff_type ), POINTER :: mgrid_rspace (:)
   TYPE ( coeff_type ), POINTER :: coef_pos
   TYPE ( coeff_type ), POINTER :: coef_vel
   TYPE ( coeff_type ), POINTER :: coef_force

!   coef_pos => kg_env % coef_pos (1)
!   coef_vel => kg_env % coef_vel (1)
!   coef_force => kg_env % coef_force (1)


   energy_numer = 0.0_dbl
   f_coef = 0.0_dbl
   energy_plus = 0.0_dbl
   energy_minus = 0.0_dbl
   NULLIFY(mgrid_gspace, mgrid_rspace)
   ngrid_level = SIZE ( pw_grid )

! allocate densities on the grids

   CALL coeff_allocate(rho_elec_rspace,pw_grid(iref),REALDATA3D)
   rho_elec_rspace%pw%in_space = REALSPACE
   CALL coeff_allocate(rho_p_rspace,pw_grid(iref),REALDATA3D)
   rho_p_rspace%pw%in_space = REALSPACE
   CALL coeff_allocate(rho_gspace,pw_grid(iref),COMPLEXDATA1D)
   rho_gspace%pw%in_space = RECIPROCALSPACE
   CALL coeff_allocate(rho_p_gspace,pw_grid(iref),COMPLEXDATA1D)
   rho_p_gspace%pw%in_space = RECIPROCALSPACE
   CALL coeff_allocate(rho_core,pw_grid(iref),COMPLEXDATA1D)
   rho_core%pw%in_space = RECIPROCALSPACE
   CALL coeff_allocate(v_gspace,pw_grid(iref),COMPLEXDATA1D)
   v_gspace%pw%in_space = RECIPROCALSPACE

! allocate the multi-grid rho's
   ALLOCATE (mgrid_rspace(ngrid_level),STAT=istat)
   IF (istat /= 0) CALL stop_memory(routine,"mgrid_rspace",0)
   ALLOCATE (mgrid_gspace(ngrid_level),STAT=istat)
   IF (istat /= 0) CALL stop_memory(routine,"mgrid_gspace",0)
   DO igrid_level=1,ngrid_level
     CALL coeff_allocate(mgrid_rspace(igrid_level),pw_grid(igrid_level),&
          REALDATA3D)
     mgrid_rspace(igrid_level)%pw%in_space = REALSPACE
     CALL coeff_allocate(mgrid_gspace(igrid_level),pw_grid(igrid_level),&
          COMPLEXDATA1D)
     mgrid_gspace(igrid_level)%pw%in_space = RECIPROCALSPACE
   END DO

!
! starting the force loop
!
   ncoef = SIZE (coef_pos%ao%cr)

   DO i = 1, ncoef 

      coef_pos % ao % cr ( i )  = coef_pos % ao % cr (i )  + delta            
      CALL coeff_zero(rho_core)
      CALL coeff_zero(rho_elec_rspace)
      CALL coeff_zero(rho_gspace)
      CALL coeff_zero(rho_p_gspace)
      CALL coeff_zero(rho_p_rspace)
      CALL coeff_zero (v_gspace)
      DO igrid_level=1,ngrid_level
        CALL coeff_zero (mgrid_rspace(igrid_level))
        CALL coeff_zero (mgrid_gspace(igrid_level))
      END DO
! Calculate the density generated by the effective point core charges qeff
      CALL calculate_epc_density(rho_core, cube_info(1),&
                       l_info,&
                       total_rho_core_rspace, kg_env, globenv)
! calculate the frozen electronic pseudo density
      CALL calculate_density ( rho_elec_rspace, rho_gspace, &
          mgrid_rspace, mgrid_gspace,&
          gridlevel_info, cube_info,&
          l_info, total_rho0_rspace, kg_env, globenv )
! calculate the polarization density
      CALL calculate_p_density ( coef_pos % ao % cr, &
          rho_p_rspace, rho_p_gspace, &
          mgrid_rspace, mgrid_gspace,&
          gridlevel_info, cube_info,&
          l_info, total_rhop_rspace, kg_env, globenv )
! get the total gspace density
      CALL coeff_sumup(rho_core,rho_gspace)
!     rho_p_gspace%pw%cc(:)=-1.0_dbl*rho_p_gspace%pw%cc(:)
      CALL coeff_sumup(rho_p_gspace,rho_gspace)
      CALL calculate_hartree(rho_gspace,energy_plus,v_gspace)

      coef_pos % ao % cr ( i )  = coef_pos % ao % cr (i ) - 2._dbl * delta
      CALL coeff_zero(rho_core)
      CALL coeff_zero(rho_elec_rspace)
      CALL coeff_zero(rho_gspace)
      CALL coeff_zero(rho_p_gspace)
      CALL coeff_zero(rho_p_rspace)
      CALL coeff_zero (v_gspace)
      DO igrid_level=1,ngrid_level
        CALL coeff_zero (mgrid_rspace(igrid_level))
        CALL coeff_zero (mgrid_gspace(igrid_level))
      END DO
! Calculate the density generated by the effective point core charges qeff
     CALL calculate_epc_density(rho_core, cube_info(1),&
                       l_info,&
                       total_rho_core_rspace, kg_env, globenv)
! calculate the frozen electronic pseudo density
     CALL calculate_density ( rho_elec_rspace, rho_gspace, &
          mgrid_rspace, mgrid_gspace,&
          gridlevel_info, cube_info,&
          l_info, total_rho0_rspace, kg_env, globenv )
! calculate the polarization density
      CALL calculate_p_density ( coef_pos % ao % cr, &
          rho_p_rspace, rho_p_gspace, &
          mgrid_rspace, mgrid_gspace,&
          gridlevel_info, cube_info,&
          l_info, total_rhop_rspace, kg_env, globenv )
! get the total gspace density

      CALL coeff_sumup(rho_core,rho_gspace)
!      rho_p_gspace%pw%cc(:)=-1.0_dbl*rho_p_gspace%pw%cc(:)
      CALL coeff_sumup(rho_p_gspace,rho_gspace)
      CALL calculate_hartree(rho_gspace,energy_minus,v_gspace)

      f_coef ( i ) = energy_minus - energy_plus

      coef_pos % ao % cr ( i )  = coef_pos % ao % cr (i)   + delta

  END DO 

  f_coef = f_coef / 2._dbl / delta
! get ener_numer
  CALL coeff_zero(rho_core)
  CALL coeff_zero(rho_elec_rspace)
  CALL coeff_zero(rho_gspace)
  CALL coeff_zero(rho_p_gspace)
  CALL coeff_zero(rho_p_rspace)
  CALL coeff_zero (v_gspace)
  DO igrid_level=1,ngrid_level
    CALL coeff_zero (mgrid_rspace(igrid_level))
    CALL coeff_zero (mgrid_gspace(igrid_level))
  END DO
! Calculate the density generated by the effective point core charges qeff
  CALL calculate_epc_density(rho_core, cube_info(1),&
                   l_info,&
                   total_rho_core_rspace, kg_env, globenv)
! calculate the frozen electronic pseudo density
  CALL calculate_density ( rho_elec_rspace, rho_gspace, &
     mgrid_rspace, mgrid_gspace,&
     gridlevel_info, cube_info,&
     l_info, total_rho0_rspace, kg_env, globenv )
! calculate the polarization density
  CALL calculate_p_density ( coef_pos % ao % cr, &
     rho_p_rspace, rho_p_gspace, &
     mgrid_rspace, mgrid_gspace,&
     gridlevel_info, cube_info,&
    l_info, total_rhop_rspace, kg_env, globenv )
! get the total gspace density
  CALL coeff_sumup(rho_core,rho_gspace)
!      rho_p_gspace%pw%cc(:)=-1.0_dbl*rho_p_gspace%pw%cc(:)
  CALL coeff_sumup(rho_p_gspace,rho_gspace)
  CALL calculate_hartree(rho_gspace,energy_numer,v_gspace)

 
! deallocate densities

  CALL coeff_deallocate(rho_elec_rspace)
  CALL coeff_deallocate(rho_gspace)
  CALL coeff_deallocate(rho_p_gspace)
  CALL coeff_deallocate(rho_p_rspace)
  CALL coeff_deallocate(rho_core)
  CALL coeff_deallocate(v_gspace)
  DO igrid_level=1,ngrid_level
     CALL coeff_deallocate(mgrid_rspace(igrid_level))
     CALL coeff_deallocate(mgrid_gspace(igrid_level))
  END DO
  DEALLOCATE (mgrid_rspace, STAT=istat)
  IF (istat /= 0) CALL stop_memory(routine,"mgrid_rspace",0)
  DEALLOCATE (mgrid_gspace,STAT=istat)
  IF (istat /= 0) CALL stop_memory(routine,"mgrid_gspace",0)

  RETURN
 
END SUBROUTINE  fcoef_p_gspace_numer 
!------------------------------------------------------------------------------!

 SUBROUTINE  fcoef_p_ppl_numer  ( delta, kg_env,  &
                                   f_coef, energy_numer, globenv )
   IMPLICIT NONE
!
! begin passed variable declaration:
!
   TYPE ( global_environment_type ), INTENT ( INOUT ) :: globenv
   TYPE ( kg_environment_type ), POINTER          :: kg_env
   REAL ( dbl ), INTENT ( IN ) :: delta
   REAL ( dbl ), INTENT ( OUT ) :: energy_numer
   REAL ( dbl ), INTENT ( OUT ), DIMENSION ( : ) :: f_coef

! begin local variable declaration:
!
   CHARACTER(LEN=*), PARAMETER :: routine =&
     "SUBROUTINE fcoef_ppl_numer (MODULE kg_force_numer)"
   INTEGER :: i, j
   INTEGER :: id,  ncoef, istat
   REAL (dbl) :: energy,energy2
   REAL (dbl) :: energy_plus, energy_minus
   LOGICAL :: calculate_forces
   TYPE ( coeff_type ), POINTER :: coef_pos
   TYPE ( coeff_type ), POINTER :: coef_vel
   TYPE ( coeff_type ), POINTER :: coef_force

!   coef_pos => kg_env % coef_pos (1)
!   coef_vel => kg_env % coef_vel (1)
!    coef_force => kg_env % coef_force (1)

   energy_numer = 0.0_dbl
   energy_plus = 0.0_dbl
   energy_minus = 0.0_dbl
   f_coef = 0.0_dbl
   calculate_forces = .FALSE.
!
! starting the force loop
!
   ncoef = SIZE ( coef_pos%ao%cr )

   DO i = 1, ncoef

      coef_pos % ao % cr ( i ) = coef_pos % ao % cr ( i ) + delta

      CALL build_kg_neighbor_lists(kg_env,globenv)
      CALL calculate_drho_ppl(coef_pos%ao%cr,coef_force%ao%cr,kg_env, &
                globenv,energy_plus,calculate_forces)
       
      coef_pos % ao % cr ( i ) = coef_pos % ao % cr ( i ) - 2._dbl*delta

      CALL build_kg_neighbor_lists(kg_env,globenv)
      CALL calculate_drho_ppl(coef_pos%ao%cr,coef_force%ao%cr,kg_env,&
                 globenv,energy_minus,calculate_forces)

      f_coef ( i ) = energy_minus - energy_plus

      coef_pos % ao % cr ( i ) = coef_pos % ao % cr ( i ) + delta

   END DO

   f_coef = f_coef / 2._dbl / delta
   CALL build_kg_neighbor_lists(kg_env,globenv)
   CALL calculate_drho_ppl(coef_pos%ao%cr,coef_force%ao%cr,&
         kg_env,globenv,energy_numer,calculate_forces)

  RETURN

END SUBROUTINE  fcoef_p_ppl_numer


!------------------------------------------------------------------------------!
END MODULE kg_force_numer
!------------------------------------------------------------------------------!
