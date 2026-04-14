!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
! **************************************************************************************************
!> \author  Ziwei Chai
!> \date    09.04.2026
!> \version 1.0
! **************************************************************************************************
MODULE gce_methods

   USE cp_control_types,                ONLY: pcc_control_type, &
                                              dft_control_type
   USE cp_dbcsr_output,                 ONLY: cp_dbcsr_write_sparse_matrix
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_realspace_grid_cube,          ONLY: cp_pw_to_cube
   USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                              cp_subsys_type
   USE cp_output_handling,              ONLY: cp_print_key_finished_output, &
                                              cp_print_key_unit_nr
   USE cp_log_handling,                 ONLY: cp_get_default_logger, &
                                              cp_logger_get_default_io_unit
   USE input_section_types,             ONLY: section_vals_type,&
                                              section_get_ivals
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE particle_list_types,             ONLY: particle_list_type
   USE physcon,                         ONLY: angstrom,&
                                              evolt
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_integrate_function,&
                                              pw_scale,&
                                              pw_transfer,&
                                              pw_zero
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type


#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gce_methods'

   PUBLIC :: planar_averaged_v_hartree_3d, planar_counter_charge

CONTAINS

! **************************************************************************************************
!> \brief calculate the planar averaged real space potential (e.g. Hartree potential) along the 
!>        surface normal of a slab model and the corresponding reference (electrostatic) energy
!>        near the simulation cell edge. 
!> \author Ziwei Chai
! **************************************************************************************************
   SUBROUTINE planar_averaged_v_hartree_3d(v_rspace, dft_control, do_gce, ref_esp, para_env)

      TYPE(pw_r3d_rs_type), POINTER                      :: v_rspace
      TYPE(dft_control_type), POINTER                    :: dft_control
      LOGICAL, INTENT(IN)                                :: do_gce
      REAL(KIND=dp), INTENT(INOUT)                       :: ref_esp
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'planar_averaged_v_hartree_3d'
      INTEGER                                            :: handle, a, b, c, x, y, z, my_io_unit
      REAL(KIND=dp)                                      :: ngrids_in_plane, tmp, my_ref_esp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: profile

      CALL timeset(routineN, handle)
      CPASSERT(ASSOCIATED(v_rspace))

      IF (do_gce) THEN

         ! a and b cell vectors are paralle to the surface, c vector is not.
         IF (dft_control%pcc_control%surf_normal == 3) THEN
            a = 1
            b = 2
            c = 3
         ELSE IF (dft_control%pcc_control%surf_normal == 2) THEN
            c = 1
            a = 2
            b = 3
         ELSE
            b = 1
            c = 2
            a = 3
         END IF
   
         my_ref_esp = 0.0_dp
         z =v_rspace%pw_grid%npts(a)*v_rspace%pw_grid%npts(b)
         ngrids_in_plane = REAL(z, dp)
   
!$OMP    PARALLEL DO DEFAULT(NONE) &
!$OMP                PRIVATE(x,y,z) &
!$OMP                SHARED(v_rspace,a,b,c)&
!$OMP                REDUCTION(+:my_ref_esp)
         DO y = v_rspace%pw_grid%bounds_local(1, b), v_rspace%pw_grid%bounds_local(2, b)
            DO x = v_rspace%pw_grid%bounds_local(1, a), v_rspace%pw_grid%bounds_local(2, a)
               DO z = v_rspace%pw_grid%bounds_local(1, c), v_rspace%pw_grid%bounds_local(2, c)
                  IF (z == v_rspace%pw_grid%bounds(2, c)) THEN
                     IF (c == 3) THEN
                        my_ref_esp = v_rspace%array(x, y, z) + my_ref_esp
                     ELSE IF (c == 2) THEN
                        my_ref_esp = v_rspace%array(y, z, x) + my_ref_esp
                     ELSE
                        my_ref_esp = v_rspace%array(z, x, y) + my_ref_esp
                     END IF
                  END IF
               END DO
            END DO
         END DO
!$OMP    END PARALLEL DO
   
         CALL para_env%sum(my_ref_esp)
         my_ref_esp = my_ref_esp / ngrids_in_plane
         ref_esp = my_ref_esp

      END IF

      IF (dft_control%do_paep) THEN

         IF (dft_control%paep_control%surf_normal == 3) THEN
            a = 1
            b = 2
            c = 3
         ELSE IF (dft_control%paep_control%surf_normal == 2) THEN
            c = 1
            a = 2
            b = 3
         ELSE
            b = 1
            c = 2
            a = 3
         END IF

         z =v_rspace%pw_grid%npts(a)*v_rspace%pw_grid%npts(b)
         ngrids_in_plane = REAL(z, dp)

         ALLOCATE(profile(v_rspace%pw_grid%bounds(1, c):v_rspace%pw_grid%bounds(2, c))) 
         profile(:) = 0.0_dp

!$OMP    PARALLEL DO DEFAULT(NONE) &
!$OMP                PRIVATE(x,y,z) &
!$OMP                SHARED(v_rspace,a,b,c) &
!$OMP                REDUCTION(+:profile)
         DO z = v_rspace%pw_grid%bounds_local(1, c), v_rspace%pw_grid%bounds_local(2, c)
            DO y = v_rspace%pw_grid%bounds_local(1, b), v_rspace%pw_grid%bounds_local(2, b)
               DO x = v_rspace%pw_grid%bounds_local(1, a), v_rspace%pw_grid%bounds_local(2, a)
                  profile(z) = profile(z) + v_rspace%array(x, y, z)
               END DO
            END DO
         END DO
!$OMP    END PARALLEL DO

         CALL para_env%sum(profile(:))
         profile(:) = profile(:) / ngrids_in_plane

         my_io_unit = cp_logger_get_default_io_unit()

! The following implementation is not pretty enough.
! A cleaner way of printing the curve is still not implemented yet.
! Ziwei Chai

!         IF(my_io_unit>0)THEN
!           WRITE(UNIT=my_io_unit, FMT=*) "PLANAR AVERAGED HARTREE POTENTIAL: POSITION(ANGSTROM), POTENTIAL(eV)"
!           DO x = v_rspace%pw_grid%bounds(1, c), v_rspace%pw_grid%bounds(2, c)
!              WRITE(UNIT=my_io_unit, FMT=*) real(x, dp)*v_rspace%pw_grid%dr(c)*angstrom, evolt*profile(x)
!           END DO
!         END IF

         DEALLOCATE(profile)

      END IF

      CALL timestop(handle)

   END SUBROUTINE planar_averaged_v_hartree_3d


! **************************************************************************************************
!> \brief add the planar counter charge density to the total charge density
!> \author Ziwei Chai
! **************************************************************************************************
   SUBROUTINE planar_counter_charge(qs_env, rho_tot_gspace, pcc_env, auxbas_pw_pool)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(pw_c1d_gs_type), INTENT(INOUT)                :: rho_tot_gspace
      TYPE(pcc_control_type), POINTER                    :: pcc_env
      TYPE(pw_pool_type), POINTER, INTENT(IN)            :: auxbas_pw_pool
      REAL(KIND=dp)                                      :: ref_esp

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'planar_counter_charge'
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_subsys_type), POINTER                      :: cp_subsys
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(pw_r3d_rs_type)                               :: pcc_density_r
      TYPE(pw_c1d_gs_type)                               :: pcc_density_g
      TYPE(section_vals_type), POINTER                   :: input
      INTEGER                                            :: handle, a, b, c, z, surf_normal, ucenter, lcenter,&
                                                            cube_unit, my_io_unit, x, y
      LOGICAL                                            :: my_mpi_io
      REAL(KIND=dp)                                      :: gau_a, tot_charge, dz, val
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: profile

      CALL timeset(routineN, handle)

      my_io_unit = cp_logger_get_default_io_unit()

      ! a and b cell vectors are paralle to the surface, c vector is not.
      IF (pcc_env%surf_normal == 3) THEN
         a = 1
         b = 2
         c = 3
      ELSE IF (pcc_env%surf_normal == 1) THEN
         c = 1
         a = 2
         b = 3
      ELSE IF (pcc_env%surf_normal == 2) THEN
         b = 1
         c = 2
         a = 3
      END IF

      CALL auxbas_pw_pool%create_pw(pcc_density_r)
      CALL auxbas_pw_pool%create_pw(pcc_density_g)
      CALL pw_zero(pcc_density_r)
      CALL pw_zero(pcc_density_g)

      tot_charge = pw_integrate_function(rho_tot_gspace, isign=-1)

      dz = pcc_density_r%pw_grid%dr(c)

      ucenter = pcc_density_r%pw_grid%bounds(2, c)- INT(pcc_env%dist_edge / dz)
      lcenter = pcc_density_r%pw_grid%bounds(1, c)+ INT(pcc_env%dist_edge / dz) 
      gau_a = 1.0_dp
!      pcc_env%gau_c = 0.1_dp
!$OMP    PARALLEL DO DEFAULT(NONE) &
!$OMP                PRIVATE(x,y,z,val) &
!$OMP                SHARED(pcc_density_r,gau_a,ucenter,lcenter,a,b,c,dz,pcc_env)
      DO z = pcc_density_r%pw_grid%bounds_local(1, c), pcc_density_r%pw_grid%bounds_local(2, c)
         val = gau_a * EXP(-(dz * REAL(z-ucenter, dp))**2.0_dp / (2.0_dp * pcc_env%gau_c**2.0_dp))
         val = val + gau_a * EXP(-(dz * REAL(z-lcenter, dp))**2.0_dp / (2.0_dp * pcc_env%gau_c**2.0_dp))
         IF (val<1.0E-15_dp) val = 0.0_dp
         DO y = pcc_density_r%pw_grid%bounds_local(1, b), pcc_density_r%pw_grid%bounds_local(2, b)
            DO x = pcc_density_r%pw_grid%bounds_local(1, a), pcc_density_r%pw_grid%bounds_local(2, a)
               IF(c==3) pcc_density_r%array(x, y, z) = val
               IF(c==2) pcc_density_r%array(y, z, x) = val
               IF(c==1) pcc_density_r%array(z, x, y) = val
            END DO
         END DO
      END DO
!$OMP    END PARALLEL DO

      tot_charge = tot_charge / pw_integrate_function(pcc_density_r, isign=1)
      CALL pw_scale(pcc_density_r, tot_charge)
      CALL pw_transfer(pcc_density_r, pcc_density_g)
      CALL pw_axpy(pcc_density_g, rho_tot_gspace)
      CALL pw_transfer(rho_tot_gspace, pcc_density_r)

      tot_charge = pw_integrate_function(rho_tot_gspace, isign=-1)

      CALL auxbas_pw_pool%give_back_pw(pcc_density_r)
      CALL auxbas_pw_pool%give_back_pw(pcc_density_g)


      CALL timestop(handle)

   END SUBROUTINE planar_counter_charge

END MODULE gce_methods
