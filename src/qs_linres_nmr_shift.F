!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief from the response current density calculates the shift tensor
!>      and the susceptibility
!> \par History
!>      created 02-2006 [MI]
!> \author MI
! *****************************************************************************
MODULE qs_linres_nmr_shift

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE cp_array_i_utils,                ONLY: cp_2d_i_p_type
  USE cp_array_r_utils,                ONLY: cp_2d_r_p_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_fm_basic_linalg,              ONLY: cp_fm_trace
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_submatrix,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_set_submatrix,&
                                             cp_fm_type
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_rs_pool_types,                ONLY: cp_rs_pool_type
  USE cp_sm_fm_interactions,           ONLY: cp_sm_fm_multiply
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE mathconstants,                   ONLY: gaussi,&
                                             twopi
  USE mathlib,                         ONLY: diamat_all
  USE message_passing,                 ONLY: mp_sum
  USE particle_types,                  ONLY: particle_type
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_grid_types,                   ONLY: pw_grid_type
  USE pw_methods,                      ONLY: pw_axpy,&
                                             pw_scale,&
                                             pw_zero
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_p_type,&
                                             pw_pool_type
  USE pw_spline_utils,                 ONLY: Eval_Interp_Spl3_pbc,&
                                             find_coeffs,&
                                             pw_spline_do_precond,&
                                             pw_spline_precond_create,&
                                             pw_spline_precond_release,&
                                             pw_spline_precond_set_kind,&
                                             pw_spline_precond_type,&
                                             spl3_pbc
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_p_type,&
                                             pw_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_grid_atom,                    ONLY: grid_atom_type
  USE qs_harmonics_atom,               ONLY: harmonics_atom_type
  USE qs_linres_nmr_epr_common_utils,  ONLY: mult_G_ov_G2_grid
  USE qs_linres_op,                    ONLY: fac_vecp,&
                                             ind_m2,&
                                             set_vecp,&
                                             set_vecp_rev
  USE qs_linres_types,                 ONLY: current_env_type,&
                                             get_current_env,&
                                             get_nmr_env,&
                                             jrho_atom_type,&
                                             nmr_env_type
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
  USE qs_operators_ao,                 ONLY: rRc_xyz_der_ao,&
                                             set_up_op_sm
  USE sparse_matrix_types,             ONLY: allocate_matrix_set,&
                                             deallocate_matrix_set,&
                                             real_matrix_p_type,&
                                             replicate_matrix_structure,&
                                             set_matrix
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: get_limit
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  ! *** Public subroutines ***
  PUBLIC :: interpolate_shift_pwgrid, gsum_shift_pwgrid, nmr_shift_print,&
       &    nmr_shift

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_linres_nmr_shift'

! **************

CONTAINS

  SUBROUTINE nmr_shift(nmr_env,current_env,qs_env,error)

    TYPE(nmr_env_type)                       :: nmr_env
    TYPE(current_env_type)                   :: current_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'nmr_shift', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iB, idir, idir2, &
                                                idir3, iiB, iiiB, ispin, &
                                                istat, natom, nspins
    LOGICAL                                  :: failure, gapw
    REAL(dp)                                 :: scale_fac
    REAL(dp), DIMENSION(3, 3)                :: chi_tensor, chi_tensor_loc
    REAL(dp), DIMENSION(:, :, :), POINTER    :: chemical_shift, &
                                                chemical_shift_loc, &
                                                chemical_shift_nics, &
                                                chemical_shift_nics_loc
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_rs_pool_type), POINTER           :: auxbas_rs_pool
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type)                          :: pw_gspace_work, &
                                                shift_pw_rspace
    TYPE(pw_p_type), DIMENSION(:, :), &
      POINTER                                :: shift_pw_gspace
    TYPE(pw_p_type), POINTER                 :: rho_gspace, rho_rspace
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool

    failure = .FALSE.
    CALL timeset(routineN,"I"," ",handle)

    NULLIFY(chemical_shift,chemical_shift_loc,chemical_shift_nics,chemical_shift_nics_loc,&
         &  cell,dft_control,pw_env,auxbas_rs_pool,auxbas_pw_pool,rho_gspace,rho_rspace,&
         &  pw_pools,particle_set)

    CALL get_qs_env(qs_env=qs_env,cell=cell,dft_control=dft_control,&
         &          particle_set=particle_set,error=error)

    gapw = dft_control%qs_control%gapw
    natom = SIZE(particle_set,1)
    nspins = dft_control%nspins

    CALL get_nmr_env(nmr_env=nmr_env,chemical_shift=chemical_shift,&
                     chemical_shift_loc=chemical_shift_loc,&
                     chemical_shift_nics=chemical_shift_nics,&
                     chemical_shift_nics_loc=chemical_shift_nics_loc,&
                     chi_tensor=chi_tensor,chi_tensor_loc=chi_tensor_loc,&
                     error=error)

    CALL get_current_env(current_env=current_env,error=error)


    CALL get_qs_env(qs_env=qs_env,pw_env=pw_env,error=error)
    CALL pw_env_get(pw_env,auxbas_rs_pool=auxbas_rs_pool,&
         &          auxbas_pw_pool=auxbas_pw_pool,pw_pools=pw_pools,&
         &          error=error)
    !
    ! Initialize
    ! Allocate grids for the calculation of jrho and the shift
    ALLOCATE(shift_pw_gspace(3,nspins),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DO ispin = 1,nspins
       DO idir = 1,3
          CALL pw_pool_create_pw(auxbas_pw_pool,shift_pw_gspace(idir,ispin)%pw,&
               &                  use_data=COMPLEXDATA1D,&
               &                  in_space=RECIPROCALSPACE,error=error)
          CALL pw_zero(shift_pw_gspace(idir,ispin)%pw)
       ENDDO
    ENDDO
    CALL pw_pool_create_pw(auxbas_pw_pool,shift_pw_rspace%pw,&
         &                 use_data=REALDATA3D,in_space=REALSPACE,error=error)
    CALL pw_zero(shift_pw_rspace%pw)
    CALL pw_pool_create_pw(auxbas_pw_pool,pw_gspace_work%pw,&
            &              use_data=COMPLEXDATA1D,&
            &              in_space=RECIPROCALSPACE,error=error)

    CALL pw_zero(pw_gspace_work%pw)
    !
    nmr_env%chi_tensor(:,:)     = 0.0_dp
    nmr_env%chi_tensor_loc(:,:) = 0.0_dp
    nmr_env%chemical_shift(:,:,:)     = 0.0_dp
    nmr_env%chemical_shift_loc(:,:,:) = 0.0_dp
    IF(nmr_env%do_nics) THEN
       nmr_env%chemical_shift_nics(:,:,:)     = 0.0_dp
       nmr_env%chemical_shift_nics_loc(:,:,:) = 0.0_dp
    ENDIF
    !
    DO iB = 1,3
       !
       CALL set_vecp(iB,iiB,iiiB)
       !
       DO ispin = 1,nspins
          !
          CALL chi_analytic_new(nmr_env,current_env,qs_env,ispin,iB,error=error)
          DO idir3=1,3
             ! set to zero for the calculation of the shift
             CALL pw_zero(shift_pw_gspace(idir3,ispin)%pw)
          ENDDO
          !
          DO idir = 1,3
             rho_rspace => current_env%jrho1_set(idir,iB)%rho%rho_r(ispin)
             rho_gspace => current_env%jrho1_set(idir,iB)%rho%rho_g(ispin)
             ! Field gradient
             ! loop over the Gvec  components: x,y,z
             DO idir2 = 1,3
                IF(idir /= idir2) THEN
                   ! in reciprocal space multiply (G_idir2(i)/G(i)^2)J_(idir)(G(i))
                   CALL mult_G_ov_G2_grid(cell,auxbas_pw_pool,rho_gspace,&
                        &                 pw_gspace_work,idir2,0.0_dp,error=error)
                   !
                   ! scale and add to the correct component of the shift column
                   CALL set_vecp_rev(idir,idir2,idir3)
                   scale_fac=fac_vecp(idir3,idir2,idir)
                   CALL pw_scale(pw_gspace_work%pw,scale_fac)
                   CALL pw_axpy(pw_gspace_work%pw,shift_pw_gspace(idir3,ispin)%pw)
                ENDIF
             ENDDO
             !
          ENDDO ! idir
       ENDDO ! ispin
       !
       ! compute shildings
       DO ispin = 1,nspins
          DO idir = 1,3
             ! Here the shifts are computed from summation of the coeff on the G-grip .
             CALL gsum_shift_pwgrid(nmr_env,pw_env,particle_set,cell,&
                  &                 shift_pw_gspace(idir,ispin),iB,idir,&
                  &                 error=error)
             ! Here first G->R and then interpolation to get the shifts.
             ! The interpolation doesnt work in parallel yet.
             ! The choice between both methods should be left to the user.
             !CALL pw_transfer(shift_pw_gspace(idir)%pw,shift_pw_rspace%pw)
             !CALL interpolate_shift_pwgrid(nmr_env,pw_env,particle_set,cell,shift_pw_rspace,&
             !             i_B,idir,nmr_section,error=error)
          ENDDO
       ENDDO
       !
       IF(gapw) THEN
          DO idir = 1,3
             ! Finally the radial functions are multiplied by the YLM and properly summed
             ! The resulting array is J on the local grid. One array per atom.
             ! Local contributions by numerical integration over the spherical grids
             CALL nmr_shift_gapw(nmr_env,current_env,qs_env,iB,idir,error=error)
          ENDDO ! idir
       ENDIF
    ENDDO
    !
    ! Dellocate grids for the calculation of jrho and the shift
    CALL pw_pool_give_back_pw(auxbas_pw_pool,pw_gspace_work%pw,error=error)
    DO ispin = 1,dft_control%nspins
       DO idir = 1,3
          CALL pw_pool_give_back_pw(auxbas_pw_pool,shift_pw_gspace(idir,ispin)%pw,error=error)
       ENDDO
    ENDDO
    DEALLOCATE(shift_pw_gspace,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    CALL pw_pool_give_back_pw(auxbas_pw_pool,shift_pw_rspace%pw,error=error)
    !
    ! Finalize
    CALL timestop(0.0_dp,handle)
    !
  END SUBROUTINE nmr_shift

  SUBROUTINE nmr_shift_gapw(nmr_env,current_env,qs_env,iB,idir,error)

    TYPE(nmr_env_type)                       :: nmr_env
    TYPE(current_env_type)                   :: current_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    INTEGER, INTENT(IN)                      :: IB, idir
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'nmr_shift_gapw', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, ia, iat, iatom, idir2, idir3, ikind, ir, ispin, istat, &
      jatom, mepos, n_nics, na, natom, natom_tot, nkind, nr, nspins, num_pe, &
      output_unit
    INTEGER, DIMENSION(2)                    :: bo
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: failure, jatomisneighbor, &
                                                paw_atom
    REAL(dp) :: ddiff, distance, itegrated_jrho, jrho1_tot, r_iatom(3), &
      r_jatom(3), rdiff(3), rij(3), rvec(3), scale_fac, threshold
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: chemical_shift_loc_tmp
    REAL(dp), DIMENSION(:, :), POINTER       :: jrho_h_grid, jrho_s_grid, &
                                                r_nics
    REAL(dp), DIMENSION(:, :, :), POINTER    :: chemical_shift_loc, &
                                                chemical_shift_nics_loc
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atom_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(grid_atom_type), POINTER            :: grid_atom
    TYPE(harmonics_atom_type), POINTER       :: harmonics
    TYPE(jrho_atom_type), DIMENSION(:), &
      POINTER                                :: jrho1_atom_set
    TYPE(jrho_atom_type), POINTER            :: jrho1_atom
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

!
!

    CALL timeset(routineN,"I"," ",handle)
    !
    failure =.FALSE.
    NULLIFY(atomic_kind_set,cell,dft_control,para_env,particle_set,&
         &  chemical_shift_loc,chemical_shift_nics_loc,jrho1_atom_set,&
         &  jrho1_atom,r_nics,jrho_h_grid,jrho_s_grid,atom_kind,&
         &  atom_list,grid_atom,harmonics,logger)
    !
    logger => cp_error_get_logger(error)
    output_unit= cp_logger_get_default_unit_nr(logger)
    !
    CALL get_qs_env(qs_env=qs_env,&
         &          atomic_kind_set=atomic_kind_set,&
         &          cell=cell,&
         &          dft_control=dft_control,&
         &          para_env=para_env,&
         &          particle_set=particle_set,&
         &          error=error)

    CALL get_nmr_env(nmr_env=nmr_env,&
         &           chemical_shift_loc=chemical_shift_loc,&
         &           chemical_shift_nics_loc=chemical_shift_nics_loc,&
         &           n_nics=n_nics,&
         &           r_nics=r_nics,&
         &           error=error)

    CALL get_current_env(current_env=current_env,&
         &               jrho1_atom_set=jrho1_atom_set,&
         &               error=error)
    !
    nkind = SIZE(atomic_kind_set,1)
    natom_tot = SIZE(particle_set,1)
    nspins = dft_control%nspins
    itegrated_jrho = 0.0_dp
    !
    ALLOCATE(chemical_shift_loc_tmp(3,natom_tot),STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    chemical_shift_loc_tmp = 0.0_dp
    !
    ! Loop over atoms to collocate the current on each atomic grid, JA
    ! Per each JA, loop over the points where the shift needs to be computed
    DO ikind = 1,nkind

       NULLIFY(atom_kind,atom_list,grid_atom,harmonics)
       atom_kind => atomic_kind_set(ikind)
       CALL get_atomic_kind(atomic_kind=atom_kind,&
            &               atom_list=atom_list,&
            &               natom=natom,&
            &               paw_atom=paw_atom,&
            &               harmonics=harmonics,&
            &               grid_atom=grid_atom)
       !
       na = grid_atom%ng_sphere
       nr = grid_atom%nr
       !
       ! Quick cycle if needed
       IF(.NOT.paw_atom) CYCLE
       !
       ! Distribute the atoms of this kind
       num_pe = para_env%num_pe
       mepos  = para_env%mepos
       bo = get_limit( natom, num_pe, mepos )
       !
       DO iat = bo(1),bo(2)
          iatom = atom_list(iat)
          r_iatom(1:3) = pbc(particle_set(iatom)%r,cell)
          !
          NULLIFY(jrho1_atom,jrho_h_grid,jrho_s_grid)
          jrho1_atom => jrho1_atom_set(iatom)
          !
          DO ispin = 1,nspins
             jrho_h_grid => jrho1_atom%jrho_vec_rad_h(idir,iB,ispin)%r_coef
             jrho_s_grid => jrho1_atom%jrho_vec_rad_s(idir,iB,ispin)%r_coef
             !
             ! loop over the atoms neighbors of iatom in terms of the current density
             ! for each compute the contribution to the shift coming from the
             ! local current density at iatom
             DO ia = 1,na
                DO ir = 1,nr
                   !
                   rvec(1) = grid_atom%rad(ir)*harmonics%a(1,ia)
                   rvec(2) = grid_atom%rad(ir)*harmonics%a(2,ia)
                   rvec(3) = grid_atom%rad(ir)*harmonics%a(3,ia)
                   !
                   !vw Can we compute that in calculate_jrho_atom_ang? wv!
                   jrho1_tot = (jrho_h_grid(ir,ia)-jrho_s_grid(ir,ia)) * grid_atom%weight(ia,ir)
                   itegrated_jrho = itegrated_jrho+jrho1_tot
                   !
                   ! local contribution to chi (NEEDS SOME MODIFICATIONS FOR PBC)
                   DO idir2 = 1,3
                      IF(idir2/=idir) THEN
                         CALL set_vecp_rev(idir,idir2,idir3)
                         scale_fac=fac_vecp(idir3,idir2,idir)
                         nmr_env%chi_tensor_loc(idir3,iB)=nmr_env%chi_tensor_loc(idir3,iB)+&
                              & scale_fac*(jrho1_tot)*nmr_env%chi_factor*&
                              & (rvec(idir2)+particle_set(iatom)%r(idir2))
                      ENDIF
                   ENDDO
                   !
                   !
                   DO jatom = 1,natom_tot
                      ! be sure you take the minimum image
                      rij(1:3) = pbc(particle_set(iatom)%r,particle_set(jatom)%r,cell)
                      distance = SQRT(rij(1)*rij(1)+rij(2)*rij(2)+rij(3)*rij(3))
                      ! the right thing to do is to build a list of neighbors, but for the time being
                      ! let's simply use a cutoff here

! my_minimum_image = .FALSE.
! IF(PRESENT(minimum_image)) THEN
!    my_minimum_image=minimum_image
!    Lxo2 = SQRT ( SUM ( cell % hmat ( :, 1 ) ** 2 ) )/2.0_dp
!    Lyo2 = SQRT ( SUM ( cell % hmat ( :, 2 ) ** 2 ) )/2.0_dp
!    Lzo2 = SQRT ( SUM ( cell % hmat ( :, 3 ) ** 2 ) )/2.0_dp
! END IF
!         IF(my_minimum_image) THEN
!           IF(ABS(rab(1)) > Lxo2 .OR. ABS(rab(2)) > Lyo2 .OR. ABS(rab(3)) > Lzo2) THEN
!             sab_orb_neighbor_node => next(sab_orb_neighbor_node)
!             CYCLE
!           END IF
!         END IF

                      !threshold = 60.0_dp
                      jatomisneighbor = .TRUE. !(distance .LE. threshold)
                      IF(jatomisneighbor) THEN

                         rdiff(1) = rij(1) - rvec(1)
                         rdiff(2) = rij(2) - rvec(2)
                         rdiff(3) = rij(3) - rvec(3)
                         ddiff = SQRT( rdiff(1)*rdiff(1)+ rdiff(2)*rdiff(2) + rdiff(3)*rdiff(3))
                         ! the component idir of J contributes to the components idir2/=idir of the shift
                         DO idir2 = 1,3

                            IF(idir2/=idir) THEN
                               CALL set_vecp_rev(idir,idir2,idir3)
                               scale_fac=fac_vecp(idir3,idir2,idir)

                               IF(ddiff.GT.1.0E-12_dp) THEN
                                  chemical_shift_loc_tmp(idir3,jatom) = &
                                       chemical_shift_loc_tmp(idir3,jatom) + &
                                       scale_fac*rdiff(idir2) * jrho1_tot/(ddiff*ddiff*ddiff)
                               ENDIF

                            ENDIF
                         ENDDO ! idir2
                         
                      ENDIF

                   ENDDO ! jatom

                   IF(nmr_env%do_nics) THEN

                      DO jatom = 1,n_nics
                         r_jatom(1:3) = r_nics(1:3,jatom)
                         rij(1:3) = pbc(particle_set(iatom)%r,r_jatom,cell)

                         distance = SQRT(rij(1)*rij(1)+rij(2)*rij(2)+rij(3)*rij(3))
                         jatomisneighbor = (distance .LE. threshold)

                         IF(jatomisneighbor) THEN

                            rdiff(1) = rij(1) - rvec(1)
                            rdiff(2) = rij(2) - rvec(2)
                            rdiff(3) = rij(3) - rvec(3)
                            ddiff = SQRT( rdiff(1)*rdiff(1)+ rdiff(2)*rdiff(2) + rdiff(3)*rdiff(3))
                            ! the component idir of J contributes to the components idir2/=idir of the shift
                            DO idir2 = 1,3

                               IF(idir2/=idir) THEN
                                  CALL set_vecp_rev(idir,idir2,idir3)
                                  scale_fac=fac_vecp(idir3,idir2,idir)

                                  IF(ddiff .GT. EPSILON(0.0_dp)) THEN
                                     chemical_shift_nics_loc(idir3,iB,jatom) = &
                                          chemical_shift_nics_loc(idir3,iB,jatom) + &
                                          scale_fac*rdiff(idir2)*jrho1_tot/(ddiff*ddiff*ddiff)
                                  ENDIF

                               ENDIF
                            ENDDO ! idir2
                            
                         ENDIF

                      ENDDO ! jatom

                   ENDIF ! do_nics
                ENDDO ! ir
             ENDDO ! ia
          ENDDO ! ispin
       ENDDO ! iat
    ENDDO ! ikind
    !
    IF(output_unit>0) THEN
       WRITE(output_unit,'(T2,A,E24.16)') 'Integrated local j_'&
            &//ACHAR(idir+119)//ACHAR(iB+119)//'(r)=',itegrated_jrho
    ENDIF
    !
    CALL mp_sum(chemical_shift_loc_tmp,para_env%group)
    chemical_shift_loc(:,iB,:) = chemical_shift_loc(:,iB,:) &
         & - nmr_env%shift_factor_gapw * chemical_shift_loc_tmp(:,:) / 2.0_dp
    IF(nmr_env%do_nics) THEN
       !here need a chemical_shift_nics_loc_tmp... also for chi_loc
       CALL stop_program(routineP,"local nics needs some more work")
       CALL mp_sum(chemical_shift_nics_loc,para_env%group)
    ENDIF
    !
    DEALLOCATE(chemical_shift_loc_tmp,STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    !
    CALL timestop(0.0_dp,handle)
    !
  END SUBROUTINE nmr_shift_gapw

! *****************************************************************************
  SUBROUTINE chi_analytic_new(nmr_env,current_env,qs_env,ispin,iB,error)
    !
    TYPE(nmr_env_type)                       :: nmr_env
    TYPE(current_env_type)                   :: current_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    INTEGER, INTENT(IN)                      :: ispin, iB
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'chi_analytic_new', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, icenter, idir, &
                                                idir2, ii, iiB, iii, iiiB, &
                                                istat, istate, j, nao, &
                                                nbr_center(2), output_unit
    LOGICAL                                  :: failure, gapw, ionode
    REAL(dp)                                 :: chi(3), contrib, contrib2, &
                                                dk(3)
    REAL(dp), DIMENSION(:), POINTER          :: occupation
    REAL(dp), DIMENSION(:, :), POINTER       :: vecbuf_0, vecbuf_1_d, &
                                                vecbuf_1_rxp
    TYPE(cp_2d_i_p_type), DIMENSION(:), &
      POINTER                                :: center_list
    TYPE(cp_2d_r_p_type)                     :: vecbuf_1_p(3)
    TYPE(cp_2d_r_p_type), DIMENSION(:), &
      POINTER                                :: centers_set
    TYPE(cp_fm_p_type), DIMENSION(:, :), &
      POINTER                                :: psi1_D, psi1_p, psi1_rxp
    TYPE(cp_fm_struct_type), POINTER         :: tmp_fm_struct
    TYPE(cp_fm_type), POINTER                :: mo_coeff, momxpsi_istate, &
                                                psi0_istate, psi_p1_istate, &
                                                psi_p2_istate, psi_rxp_istate
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s, op_mom_ao, op_p_ao
    TYPE(real_matrix_p_type), &
      DIMENSION(:, :), POINTER               :: op_mom_der_ao

!
!

    failure = .FALSE.
    !
    CALL timeset(routineN,"I"," ",handle)
    !
    NULLIFY(dft_control,matrix_s,mos,para_env,mo_coeff,occupation,op_mom_ao,&
         &  op_mom_der_ao,center_list,centers_set,psi0_istate,psi_rxp_istate,&
         &  psi_p1_istate,psi_p2_istate,momxpsi_istate,op_p_ao,psi1_p,psi1_rxp,&
         &  psi1_D)

    logger => cp_error_get_logger(error)
    ionode = logger%para_env%mepos==logger%para_env%source
    output_unit = -1
    IF (ionode) THEN
      output_unit= cp_logger_get_default_unit_nr(logger)
    END IF

    CALL get_qs_env(qs_env=qs_env,&
         &          dft_control=dft_control,&
         &          matrix_s=matrix_s,&
         &          mos=mos,&
         &          para_env=para_env,&
         &          error=error)

    gapw = dft_control%qs_control%gapw
    IF (output_unit>0) THEN
       WRITE(output_unit,*) 'WARNING WARNING WARNING WARNING WARNING WARNING WARNING'
       WRITE(output_unit,*) 'chi_analytic: gapw=',gapw!,' reset to .false.'
       WRITE(output_unit,*) 'WARNING WARNING WARNING WARNING WARNING WARNING WARNING'
    END IF
    !gapw=.FALSE.

    CALL get_current_env(current_env=current_env,&
         &               nao=nao,&
         &               nbr_center=nbr_center,&
         &               center_list=center_list,&
         &               centers_set=centers_set,&
         &               op_p_ao=op_p_ao,&
         &               psi1_p=psi1_p,&
         &               psi1_rxp=psi1_rxp,&
         &               psi1_D=psi1_D,&
         &               error=error)
    !
    ! Allocate sparse matrices for dipole, quadrupole and their derivatives => 9x3
    ! Remember the derivatives are antisymmetric
    CALL allocate_matrix_set(op_mom_ao,9,error=error)
    CALL allocate_matrix_set(op_mom_der_ao,9,3,error=error)
    CALL set_up_op_sm(op_mom_ao(1)%matrix,qs_env,symmetry="none",&
         &            name="op_mom",error=error)
    DO idir2=1,3
       CALL replicate_matrix_structure(op_mom_ao(1)%matrix, &
            &                          op_mom_der_ao(1,idir2)%matrix,&
            &                          "op_mom_der_ao"//"-"//TRIM(&
            &                          ADJUSTL(cp_to_string(idir2))),&
            &                          target_symmetry="none",&
            &                          error=error)
       CALL set_matrix(op_mom_der_ao(1,idir2)%matrix,0.0_dp)
    ENDDO

    DO idir = 2,SIZE(op_mom_ao,1)
       CALL replicate_matrix_structure(op_mom_ao(1)%matrix, &
            &                          op_mom_ao(idir)%matrix,&
            &                          "op_mom_ao"//"-"//TRIM(&
            &                          ADJUSTL(cp_to_string(idir))),&
            &                          target_symmetry="none",&
            &                          error=error)
       CALL set_matrix(op_mom_ao(idir)%matrix,0.0_dp)
       DO idir2=1,3
          CALL replicate_matrix_structure(op_mom_ao(1)%matrix,&
               &                          op_mom_der_ao(idir,idir2)%matrix,&
               &                          "op_mom_der_ao"//"-"//TRIM(&
               &                          ADJUSTL(cp_to_string(idir*idir2))),&
               &                          target_symmetry="none",&
               &                          error=error)
          CALL set_matrix(op_mom_der_ao(idir,idir2)%matrix,0.0_dp)
       ENDDO
    ENDDO
    !
    ! get iiB and iiiB
    CALL set_vecp(iB,iiB,iiiB)
    !
    ! get ground state MOS
    CALL get_mo_set(mo_set=mos(ispin)%mo_set,mo_coeff=mo_coeff,&
         &          occupation_numbers=occupation)
    !
    ! Allocate full matrices for only one vector
    NULLIFY(tmp_fm_struct)
    CALL cp_fm_struct_create(tmp_fm_struct,nrow_global=nao,&
         &                   ncol_global=1,para_env=para_env,&
         &                   context=mo_coeff%matrix_struct%context,&
         &                   error=error)
    CALL cp_fm_create(psi0_istate, tmp_fm_struct ,error=error)
    CALL cp_fm_create(psi_rxp_istate, tmp_fm_struct ,error=error)
    CALL cp_fm_create(psi_p1_istate, tmp_fm_struct ,error=error)
    CALL cp_fm_create(psi_p2_istate, tmp_fm_struct ,error=error)
    CALL cp_fm_create(momxpsi_istate, tmp_fm_struct ,error=error)
    CALL cp_fm_struct_release( tmp_fm_struct ,error=error)
    !
    ! Allocate buffer vectors
    NULLIFY(vecbuf_0,vecbuf_1_rxp,vecbuf_1_d)
    ALLOCATE(vecbuf_0(1,nao),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(vecbuf_1_rxp(1,nao),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(vecbuf_1_d(1,nao),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DO idir=1,3
       NULLIFY(vecbuf_1_p(idir)%array)
       ALLOCATE(vecbuf_1_p(idir)%array(1,nao),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ENDDO
    !
    ! Initialize the temporary vector chi
    chi = 0.0_dp
    !
    ! Start loop over the occupied  states
    DO icenter = 1,nbr_center(ispin)
       !
       ! Get the Wannier center of the istate-th ground state orbital
       dk(1:3) = centers_set(ispin)%array(1:3,icenter)
       !
       ! Compute the multipole integrals for the state istate,
       ! using as reference center the corresponding Wannier center
       DO idir = 1,9
          CALL set_matrix(op_mom_ao(idir)%matrix, 0.0_dp)
          DO idir2 = 1,3
             CALL set_matrix(op_mom_der_ao(idir,idir2)%matrix, 0.0_dp)
          ENDDO
       ENDDO
       !
       !WRITE(*,*) ' state ', istate, dk(1:3)
       CALL rRc_xyz_der_ao(op_mom_ao,op_mom_der_ao,qs_env,dk,order=2,&
            &              minimum_image=.FALSE.,soft=gapw,error=error)


       DO j = center_list(ispin)%array(1,icenter),center_list(ispin)%array(1,icenter+1)-1
          istate = center_list(ispin)%array(2,j)

          !
          ! EXTRACT from the ground state full MOS matrix the coefficients of istate
          CALL cp_fm_get_submatrix(mo_coeff,vecbuf_0,1,istate,nao,1,&
               &                   transpose=.TRUE.,error=error)
          vecbuf_0(1,1:nao) = vecbuf_0(1,1:nao) * occupation(istate)
          CALL cp_fm_set_submatrix(psi0_istate,vecbuf_0,1,1,nao,1,&
               &                   transpose=.TRUE.,error=error)
          !
          ! Extract from the full MOS matrices of the response functions the coefficients
          ! referring to istate
          ! psi1_rxp_iB_istate
          CALL cp_fm_get_submatrix(psi1_rxp(ispin,iB)%matrix,vecbuf_1_rxp,1,&
               &                   istate,nao,1,transpose=.TRUE.,error=error)
       
          IF(current_env%full) THEN
             ! psi1_d_iB_istate
             CALL cp_fm_get_submatrix(psi1_D(ispin,iB)%matrix,vecbuf_1_d,1,istate,&
                  &                   nao,1,transpose=.TRUE.,error=error)
             vecbuf_1_rxp(1,1:nao) = vecbuf_1_rxp(1,1:nao) - vecbuf_1_d(1,1:nao)
          ENDIF
          CALL cp_fm_set_submatrix(psi_rxp_istate,vecbuf_1_rxp,1,1,nao,1,&
               &                   transpose=.TRUE.,error=error)
          !
          ! psi1_p_iiB_istate and psi1_p_iiiB_istate
          CALL cp_fm_get_submatrix(psi1_p(ispin,iiB)%matrix,&
               &                   vecbuf_1_p(iiB)%array,1,istate,nao,1,&
               &                   transpose=.TRUE.,error=error)
          CALL cp_fm_set_submatrix(psi_p1_istate,vecbuf_1_p(iiB)%array,1,1,&
               &                   nao,1,transpose=.TRUE.,error=error)
          CALL cp_fm_get_submatrix(psi1_p(ispin,iiiB)%matrix,&
               &                   vecbuf_1_p(iiiB)%array,1,istate,nao,1,&
               &                   transpose=.TRUE.,error=error)
          CALL cp_fm_set_submatrix(psi_p2_istate,vecbuf_1_p(iiiB)%array,1,1,&
               &                   nao,1,transpose=.TRUE.,error=error)
          !
          ! Multuply left and right by the appropriate coefficients and sum into the
          ! correct component of the chi tensor using the appropriate multiplicative factor
          ! (don't forget the occupation number)
          ! Loop over the cartesian components of the tensor
          ! The loop over the components of the external field is external, thereby
          ! only one column of the chi tensor is computed here
          DO idir = 1,3
             ! get ii and iii
             CALL set_vecp(idir,ii,iii)
             !
             ! term: 2[C0| (r-dk)_ii |d_iii(C1(rxp-D))]-2[C0| (r-dk)_iii |d_ii(C1(rxp-D))]
             ! the factor 2 should be already included in the matrix elements
             contrib = 0.0_dp
             CALL cp_sm_fm_multiply(op_mom_der_ao(ii,iii)%matrix,psi_rxp_istate,&
                  &                 momxpsi_istate,ncol=1,alpha=1.e0_dp,error=error)
             CALL cp_fm_trace(psi0_istate,momxpsi_istate,contrib,error=error)
             chi(idir) = chi(idir) + 2.0_dp * contrib
             
             contrib = 0.0_dp
             CALL cp_sm_fm_multiply(op_mom_der_ao(iii,ii)%matrix,psi_rxp_istate,&
                  &                 momxpsi_istate,ncol=1,alpha=1.e0_dp,error=error)
             CALL cp_fm_trace(psi0_istate,momxpsi_istate,contrib,error=error)
             chi(idir) = chi(idir) - 2.0_dp *contrib
             !
             ! correction: dk_ii*2[C0| d_iii(C1(rxp-D))] - dk_iii*2[C0| d_ii(C1(rxp-D))]
             ! factor 2 not included in the matrix elements
             contrib = 0.0_dp
             CALL cp_sm_fm_multiply(op_p_ao(iii)%matrix,psi_rxp_istate,&
                  &                 momxpsi_istate,ncol=1,alpha=1.e0_dp,error=error)
             CALL cp_fm_trace(psi0_istate,momxpsi_istate,contrib,error=error)
             chi(idir) = chi(idir) + 2.0_dp*dk(ii)*contrib
             !
             contrib2 = 0.0_dp
             CALL cp_sm_fm_multiply(op_p_ao(ii)%matrix,psi_rxp_istate,&
                  &                 momxpsi_istate,ncol=1,alpha=1.e0_dp,error=error)
             CALL cp_fm_trace(psi0_istate,momxpsi_istate,contrib2,error=error)
             chi(idir) = chi(idir) - 2.0_dp*dk(iii)*contrib2
             !
             ! term: -2[C0| (r-dk)_ii  (r-dk)_iiB | d_iii(C1(piiiB))] \
             !       +2[C0| (r-dk)_iii (r-dk)_iiB | d_ii(C1(piiiB))]
             ! the factor 2 should be already included in the matrix elements
             contrib = 0.0_dp
             idir2 = ind_m2(ii,iiB)
             CALL cp_sm_fm_multiply(op_mom_der_ao(idir2,iii)%matrix,psi_p2_istate,&
                  &                 momxpsi_istate,ncol=1,alpha=1.e0_dp,error=error)
             CALL cp_fm_trace(psi0_istate,momxpsi_istate,contrib,error=error)
             chi(idir) = chi(idir) - 2.0_dp * contrib
             contrib2 = 0.0_dp
             IF(iiB==iii) THEN
                CALL cp_sm_fm_multiply(op_mom_ao(ii)%matrix,psi_p2_istate,&
                     &                 momxpsi_istate,ncol=1,alpha=1.e0_dp,error=error)
                CALL cp_fm_trace(psi0_istate,momxpsi_istate,contrib2,error=error)
                chi(idir) = chi(idir) - contrib2
             ENDIF

             contrib = 0.0_dp
             idir2 = ind_m2(iii,iiB)
             CALL cp_sm_fm_multiply(op_mom_der_ao(idir2,ii)%matrix,psi_p2_istate,&
                  &                 momxpsi_istate,ncol=1,alpha=1.e0_dp,error=error)
             CALL cp_fm_trace(psi0_istate,momxpsi_istate,contrib,error=error)
             chi(idir) = chi(idir) + 2.0_dp * contrib
             contrib2 = 0.0_dp
             IF(iiB==ii) THEN
                CALL cp_sm_fm_multiply(op_mom_ao(iii)%matrix,psi_p2_istate,&
                     &                 momxpsi_istate,ncol=1,alpha=1.e0_dp,error=error)
                CALL cp_fm_trace(psi0_istate,momxpsi_istate,contrib2,error=error)
                chi(idir) = chi(idir) + contrib2
             ENDIF
             !
             ! correction: -dk_ii * 2[C0|(r-dk)_iiB | d_iii(C1(piiiB))] \
             !             +dk_iii * 2[C0|(r-dk)_iiB | d_ii(C1(piiiB))]
             ! the factor 2 should be already included in the matrix elements
             ! no additional correction terms because of the orthogonality between C0 and C1
             contrib = 0.0_dp
             CALL cp_sm_fm_multiply(op_mom_der_ao(iiB,iii)%matrix,psi_p2_istate,&
                  &                 momxpsi_istate,ncol=1,alpha=1.e0_dp,error=error)
             CALL cp_fm_trace(psi0_istate,momxpsi_istate,contrib,error=error)
             chi(idir) = chi(idir) - 2.0_dp*dk(ii)*contrib
             !
             contrib2 = 0.0_dp
             CALL cp_sm_fm_multiply(op_mom_der_ao(iiB,ii)%matrix,psi_p2_istate,&
                  &                 momxpsi_istate,ncol=1,alpha=1.e0_dp,error=error)
             CALL cp_fm_trace(psi0_istate,momxpsi_istate,contrib2,error=error)
             chi(idir) = chi(idir) + 2.0_dp*dk(iii)*contrib2
             !
             ! term: +2[C0| (r-dk)_ii  (r-dk)_iiiB | d_iii(C1(piiB))] \
             !       -2[C0| (r-dk)_iii (r-dk)_iiiB | d_ii(C1(piiB))]
             ! the factor 2 should be already included in the matrix elements
             contrib = 0.0_dp
             idir2 = ind_m2(ii,iiiB)
             CALL cp_sm_fm_multiply(op_mom_der_ao(idir2,iii)%matrix,psi_p1_istate,&
                  &                 momxpsi_istate,ncol=1,alpha=1.e0_dp,error=error)
             CALL cp_fm_trace(psi0_istate,momxpsi_istate,contrib,error=error)
             chi(idir) = chi(idir) + 2.0_dp*contrib
             contrib2 = 0.0_dp
             IF(iiiB==iii) THEN
                CALL cp_sm_fm_multiply(op_mom_ao(ii)%matrix,psi_p1_istate,&
                     &                 momxpsi_istate,ncol=1,alpha=1.e0_dp,error=error)
                CALL cp_fm_trace(psi0_istate,momxpsi_istate,contrib2,error=error)
                chi(idir) = chi(idir) + contrib2
             ENDIF
             
             contrib = 0.0_dp
             idir2 = ind_m2(iii,iiiB)
             CALL cp_sm_fm_multiply(op_mom_der_ao(idir2,ii)%matrix,psi_p1_istate,&
                  &                 momxpsi_istate,ncol=1,alpha=1.e0_dp,error=error)
             CALL cp_fm_trace(psi0_istate,momxpsi_istate,contrib,error=error)
             chi(idir) = chi(idir) - 2.0_dp * contrib
             contrib2 = 0.0_dp
             IF(iiiB==ii) THEN
                CALL cp_sm_fm_multiply(op_mom_ao(iii)%matrix,psi_p1_istate,&
                     &                 momxpsi_istate,ncol=1,alpha=1.e0_dp,error=error)
                CALL cp_fm_trace(psi0_istate,momxpsi_istate,contrib2,error=error)
                chi(idir) = chi(idir) - contrib2
             ENDIF
             !
             ! correction: +dk_ii * 2[C0|(r-dk)_iiiB | d_iii(C1(piiB))] +\
             !             -dk_iii * 2[C0|(r-dk)_iiiB | d_ii(C1(piiB))]
             ! the factor 2 should be already included in the matrix elements
             contrib = 0.0_dp
             CALL cp_sm_fm_multiply(op_mom_der_ao(iiiB,iii)%matrix,psi_p1_istate,&
                  &                 momxpsi_istate,ncol=1,alpha=1.e0_dp,error=error)
             CALL cp_fm_trace(psi0_istate,momxpsi_istate,contrib,error=error)
             chi(idir) = chi(idir) + 2.0_dp*dk(ii)*contrib
             !
             contrib2 = 0.0_dp
             CALL cp_sm_fm_multiply(op_mom_der_ao(iiiB,ii)%matrix,psi_p1_istate,&
                  &                 momxpsi_istate,ncol=1,alpha=1.e0_dp,error=error)
             CALL cp_fm_trace(psi0_istate,momxpsi_istate,contrib2,error=error)
             chi(idir) = chi(idir) - 2.0_dp*dk(iii)*contrib2
             
          ENDDO  ! idir

       ENDDO ! istate

    ENDDO  ! icenter
    !
    !IF (output_unit>0) WRITE(output_unit,*) "CHI TENSOR : iB =", iB
    DO idir = 1,3
       nmr_env%chi_tensor(idir,iB) = nmr_env%chi_tensor(idir,iB)+&
            &                        chi(idir)*nmr_env%chi_factor
       !IF (output_unit>0) WRITE(output_unit,*) idir, nmr_env%chi_tensor(idir,iB)
    ENDDO
    !
    ! deallocate the sparse matrices
    CALL deallocate_matrix_set(op_mom_ao,error=error)
    CALL deallocate_matrix_set(op_mom_der_ao,error=error)
    CALL cp_fm_release (psi0_istate, error=error )
    CALL cp_fm_release (psi_rxp_istate, error=error )
    CALL cp_fm_release (psi_p1_istate, error=error )
    CALL cp_fm_release (psi_p2_istate, error=error )
    CALL cp_fm_release (momxpsi_istate, error=error )
    DEALLOCATE(vecbuf_0,vecbuf_1_rxp,vecbuf_1_d,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DO idir = 1,3
       DEALLOCATE(vecbuf_1_p(idir)%array,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ENDDO

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE chi_analytic_new

! *****************************************************************************
!> \brief interpolate the shift calculated on the PW grid in order to ger
!>       the value on arbitrary points in real space
!> \param nme_env to get the shift tensor and the list of additional points
!> \param particle_set for the atomic position
!> \param cell to take into account the pbs, and to have the volume
!> \param shift_pw_rspace specific component of the shift tensor on the pw grid
!> \param i_B component of the magnetic field for which the shift is calculated (row)
!> \param idir component of the vector \int_{r}[ ((r-r') x j(r))/|r-r'|^3 ] = Bind(r')
!> \author MI
! *****************************************************************************
  SUBROUTINE interpolate_shift_pwgrid(nmr_env,pw_env,particle_set,cell,shift_pw_rspace,&
       i_B,idir,nmr_section,error)

    TYPE(nmr_env_type)                       :: nmr_env
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(pw_p_type)                          :: shift_pw_rspace
    INTEGER, INTENT(IN)                      :: i_B, idir
    TYPE(section_vals_type), POINTER         :: nmr_section
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'interpolate_shift_pwgrid', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: aint_precond, iat, iatom, &
                                                max_iter, n_nics, natom, &
                                                precond_kind
    INTEGER, DIMENSION(:), POINTER           :: cs_atom_list
    LOGICAL                                  :: failure, success
    REAL(dp)                                 :: eps_r, eps_x, R_iatom(3), &
                                                ra(3), shift_val
    REAL(dp), DIMENSION(:, :), POINTER       :: r_nics
    REAL(dp), DIMENSION(:, :, :), POINTER    :: chemical_shift, &
                                                chemical_shift_nics
    TYPE(pw_p_type)                          :: shiftspl
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(pw_spline_precond_type), POINTER    :: precond
    TYPE(section_vals_type), POINTER         :: interp_section

    failure = .FALSE.
#ifdef __parallel
    CALL stop_program(routineP,"Doesnt support parallel yet.")
#endif
    NULLIFY (interp_section)
    NULLIFY (auxbas_pw_pool,precond)
    NULLIFY (cs_atom_list,chemical_shift,chemical_shift_nics,r_nics)

    CPPrecondition(ASSOCIATED(shift_pw_rspace%pw),cp_failure_level,routineP,error,failure)
    IF(.NOT. failure) THEN

      interp_section =>  section_vals_get_subs_vals(nmr_section,&
                        "INTERPOLATOR",error=error)
      CALL section_vals_val_get(interp_section,"aint_precond", &
                                i_val=aint_precond, error=error)
      CALL section_vals_val_get(interp_section,"precond",i_val=precond_kind, error=error)
      CALL section_vals_val_get(interp_section,"max_iter",i_val=max_iter, error=error)
      CALL section_vals_val_get(interp_section,"eps_r",r_val=eps_r, error=error)
      CALL section_vals_val_get(interp_section,"eps_x",r_val=eps_x, error=error)

      ! calculate spline coefficients
      CALL pw_env_get(pw_env,  auxbas_pw_pool=auxbas_pw_pool,error=error)
      CALL pw_pool_create_pw(auxbas_pw_pool,shiftspl%pw, &
                             use_data=REALDATA3D,in_space=REALSPACE,&
                             error=error)

      CALL pw_spline_precond_create(precond,precond_kind=aint_precond,&
          pool=auxbas_pw_pool,pbc=.TRUE.,transpose=.FALSE.,error=error)
      CALL pw_spline_do_precond(precond,shift_pw_rspace%pw,shiftspl%pw,error=error)
      CALL pw_spline_precond_set_kind(precond,precond_kind,error=error)
      success=find_coeffs(values=shift_pw_rspace%pw,coeffs=shiftspl%pw,&
             linOp=spl3_pbc,preconditioner=precond,pool=auxbas_pw_pool, &
             eps_r=eps_r,eps_x=eps_x,max_iter=max_iter, &
             error=error)
      CPPostconditionNoFail(success,cp_warning_level,routineP,error)
      CALL pw_spline_precond_release(precond,error=error)

      CALL get_nmr_env(nmr_env=nmr_env, cs_atom_list=cs_atom_list,&
           chemical_shift=chemical_shift,&
           chemical_shift_nics=chemical_shift_nics,&
           n_nics=n_nics,r_nics=r_nics,error=error)

      IF(ASSOCIATED(cs_atom_list)) THEN
        natom = SIZE(cs_atom_list,1)
      ELSE
        natom = -1
      ENDIF

      DO iat = 1,natom
         iatom = cs_atom_list(iat)
         R_iatom = pbc(particle_set(iatom)%r,cell)
         shift_val = Eval_Interp_Spl3_pbc(R_iatom,shiftspl%pw,error=error)
         chemical_shift(idir,i_B,iatom)= chemical_shift(idir,i_B,iatom)+&
                                         nmr_env%shift_factor * twopi**2 * shift_val
      END DO

      IF(nmr_env%do_nics) THEN
        DO iatom = 1,n_nics
          ra(1:3) = r_nics(1:3,iatom)
          R_iatom = pbc(ra,cell)
          shift_val = Eval_Interp_Spl3_pbc(R_iatom,shiftspl%pw,error)
          chemical_shift_nics(idir,i_B,iatom)= chemical_shift_nics(idir,i_B,iatom)+&
                                               nmr_env%shift_factor * twopi**2 * shift_val
        END DO
      END IF

      CALL pw_pool_give_back_pw(auxbas_pw_pool,shiftspl%pw,error=error)

    END IF

  END SUBROUTINE interpolate_shift_pwgrid

! *****************************************************************************
  SUBROUTINE gsum_shift_pwgrid(nmr_env,pw_env,particle_set,cell,shift_pw_gspace,&
       &                       i_B,idir,error)
    TYPE(nmr_env_type)                       :: nmr_env
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(pw_p_type)                          :: shift_pw_gspace
    INTEGER, INTENT(IN)                      :: i_B, idir
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'gsum_shift_pwgrid', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iat, iatom, n_nics, &
                                                natom
    INTEGER, DIMENSION(:), POINTER           :: cs_atom_list
    LOGICAL                                  :: failure
    REAL(dp)                                 :: R_iatom(3), ra(3)
    REAL(dp), DIMENSION(:, :), POINTER       :: r_nics
    REAL(dp), DIMENSION(:, :, :), POINTER    :: chemical_shift, &
                                                chemical_shift_nics

!
!

    failure = .FALSE.
    !
    CALL timeset(routineN,"I"," ",handle)
    !
    NULLIFY(cs_atom_list,chemical_shift,chemical_shift_nics,r_nics)
    CPPrecondition(ASSOCIATED(shift_pw_gspace%pw),cp_failure_level,routineP,error,failure)
    IF(.NOT. failure) THEN
       !
       CALL get_nmr_env(nmr_env=nmr_env, cs_atom_list=cs_atom_list,&
            &           chemical_shift=chemical_shift,&
            &           chemical_shift_nics=chemical_shift_nics,&
            &           n_nics=n_nics,r_nics=r_nics,error=error)
       !
       IF(ASSOCIATED(cs_atom_list)) THEN
          natom = SIZE(cs_atom_list,1)
       ELSE
          natom = -1
       ENDIF
       !
       ! compute the chemical shift
       DO iat = 1,natom
          iatom = cs_atom_list(iat)
          R_iatom = pbc(particle_set(iatom)%r,cell)
          chemical_shift(idir,i_B,iatom)= chemical_shift(idir,i_B,iatom)+&
               & nmr_env%shift_factor * twopi**2 * REAL(gsumr(R_iatom,shift_pw_gspace%pw),dp)
       ENDDO
       !
       ! compute nics
       IF(nmr_env%do_nics) THEN
          DO iat = 1,n_nics
             ra = pbc(r_nics(:,iat),cell)
             chemical_shift_nics(idir,i_B,iat)= chemical_shift_nics(idir,i_B,iat)+&
                  & nmr_env%shift_factor * twopi**2 * REAL(gsumr(ra,shift_pw_gspace%pw),dp)
          ENDDO
       ENDIF
    ENDIF

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE gsum_shift_pwgrid

! *****************************************************************************
  FUNCTION gsumr(r,pw)
    REAL(dp), INTENT(IN)                     :: r(3)
    TYPE(pw_type), POINTER                   :: pw
    COMPLEX(dp)                              :: gsumr

    COMPLEX(dp)                              :: rg
    INTEGER                                  :: ig
    TYPE(pw_grid_type), POINTER              :: grid

    grid => pw%pw_grid
    gsumr = CMPLX(0.0_dp,0.0_dp,KIND=dp)
    DO ig = grid%first_gne0,grid%ngpts_cut_local
       rg = ( grid%g(1,ig)*r(1) + grid%g(2,ig)*r(2) + grid%g(3,ig)*r(3) ) * gaussi
       gsumr = gsumr + pw%cc(ig)*EXP(rg)
    ENDDO
    IF(grid%have_g0) gsumr = gsumr + pw%cc(1)
    CALL mp_sum(gsumr,grid%para%group)
  END FUNCTION gsumr

! *****************************************************************************
!> \brief Shielding tensor and Chi are printed into a file
!>       if required from input
!>       It is possible to print only for a subset of atoms or
!>       or points in non-ionic positions
!> \author MI
! *****************************************************************************
  SUBROUTINE nmr_shift_print(nmr_env,qs_env,error)
    TYPE(nmr_env_type)                       :: nmr_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'nmr_shift_print', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=2)                         :: element_symbol
    CHARACTER(LEN=default_string_length)     :: name, title
    INTEGER                                  :: i, iatom, ir, istat, j, &
                                                n_nics, nat_print, natom, &
                                                output_unit, unit_atoms, &
                                                unit_nics
    INTEGER, DIMENSION(:), POINTER           :: cs_atom_list
    LOGICAL                                  :: failure, gapw, ionode
    REAL(dp) :: chi_aniso, chi_iso, chi_sym(3,3), chi_sym_tot(3,3), &
      chi_tensor(3,3), chi_tensor_loc(3,3), chi_tmp(3,3), eig(3), eig_tmp(3), &
      rpos(3), shift_aniso, shift_iso, shift_nosym(3,3), shift_sym(3,3), &
      shift_tmp(3,3)
    REAL(dp), DIMENSION(:, :), POINTER       :: r_nics
    REAL(dp), DIMENSION(:, :, :), POINTER :: chemical_shift, &
      chemical_shift_loc, chemical_shift_nics, chemical_shift_nics_loc, &
      chemical_shift_nics_tot, chemical_shift_tot
    REAL(dp), EXTERNAL                       :: DDOT
    TYPE(atomic_kind_type), POINTER          :: atom_kind
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(section_vals_type), POINTER         :: lr_section, nmr_section

    failure = .FALSE.

    NULLIFY(chemical_shift, chemical_shift_nics, r_nics, chemical_shift_loc, &
         &  chemical_shift_nics_loc, logger, particle_set,atom_kind,dft_control)

    logger => cp_error_get_logger(error)
    ionode = logger%para_env%mepos==logger%para_env%source
    output_unit= cp_logger_get_default_unit_nr(logger)

    lr_section => section_vals_get_subs_vals(qs_env%input,"PROPERTIES%LINRES",error=error)
    nmr_section => section_vals_get_subs_vals(qs_env%input, &
                                                 "PROPERTIES%LINRES%NMR",error=error)

    CALL get_nmr_env(nmr_env=nmr_env, &
         &           chemical_shift=chemical_shift, &
         &           chemical_shift_nics=chemical_shift_nics, &
         &           chemical_shift_loc=chemical_shift_loc, &
         &           chemical_shift_nics_loc=chemical_shift_nics_loc, &
         &           chi_tensor=chi_tensor, &
         &           chi_tensor_loc=chi_tensor_loc, &
         &           cs_atom_list=cs_atom_list, &
         &           n_nics=n_nics, &
         &           r_nics=r_nics, &
         &           error=error)

    CALL get_qs_env(qs_env=qs_env, &
         &          dft_control=dft_control, &
         &          particle_set=particle_set, &
         &          error=error)

    natom = SIZE(particle_set,1)
    gapw = dft_control%qs_control%gapw
    nat_print=SIZE(cs_atom_list,1)

    ALLOCATE(chemical_shift_tot(3,3,nat_print),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    IF(nmr_env%do_nics) THEN
       ALLOCATE(chemical_shift_nics_tot(3,3,n_nics),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ENDIF
    ! Finalize Chi calculation
    IF(gapw) THEN
       ! Symmetrize
       DO i = 1,3
          DO j = 1,3
             chi_sym(i,j) = (chi_tensor(i,j)+chi_tensor(j,i))/2.0_dp
             chi_sym_tot(i,j) = (chi_tensor(i,j)+chi_tensor(j,i))/2.0_dp
             chi_sym_tot(i,j) = (chi_tensor(i,j)+chi_tensor(j,i)+&
                                 chi_tensor_loc(i,j)+chi_tensor_loc(j,i))/2.0_dp
          END DO
       END DO
    ELSE
       ! Symmetrize
       DO i = 1,3
          DO j = 1,3
             chi_sym(i,j) = (chi_tensor(i,j)+chi_tensor(j,i))/2.0_dp
          END DO
       END DO
       CALL  dcopy(9,chi_sym,1,chi_sym_tot,1)
    END IF
    CALL  dcopy(9,chi_sym_tot,1,chi_tmp,1)
    CALL  diamat_all(chi_tmp,eig)
    chi_iso   = (eig(1)+ eig(2)+eig(3))/3.e0_dp
    chi_aniso =  eig(3)-(eig(2)+eig(1))/2.e0_dp
    !
    IF(output_unit > 0) THEN
       WRITE(output_unit,'(T2,A,E24.16)') 'CheckSum Chi =',&
            & SQRT(DDOT(9,chi_tensor(1,1),1,chi_tensor(1,1),1))
    ENDIF
    !
    IF (BTEST(cp_print_key_should_output(logger%iter_info,nmr_section,&
              "PRINT%CHI_TENSOR",error=error),cp_p_file)) THEN

       unit_atoms=cp_print_key_unit_nr(logger,nmr_section,"PRINT%CHI_TENSOR",&
                  extension=".data",middle_name="CHI",log_filename=.FALSE.,&
                  error=error)

       WRITE(title,'(A)') "Magnetic Susceptibility Tensor "
       IF(unit_atoms > 0) THEN
          WRITE(unit_atoms,'(T2,A)') title
          IF(gapw) THEN
             WRITE(unit_atoms,'(T30,A)') " CHI from SOFT J in 10^-30 J/T^2 units"
             WRITE(unit_atoms,'(3(A,f15.6))') ' XX=',chi_tensor(1,1),&
                  &                           ' XY=',chi_tensor(1,2),&
                  &                           ' XZ=',chi_tensor(1,3)
             WRITE(unit_atoms,'(3(A,f15.6))') ' YX=',chi_tensor(2,1),&
                  &                           ' YY=',chi_tensor(2,2),&
                  &                           ' YZ=',chi_tensor(2,3)
             WRITE(unit_atoms,'(3(A,f15.6))') ' ZX=',chi_tensor(3,1),&
                  &                           ' ZY=',chi_tensor(3,2),&
                  &                           ' ZZ=',chi_tensor(3,3)
             WRITE(unit_atoms,'(T30,A)') " CHI from LOCAL J in 10^-30 J/T^2 units"
             WRITE(unit_atoms,'(3(A,f15.6))') ' XX=',chi_tensor_loc(1,1),&
                  ' XY=',chi_tensor_loc(1,2),' XZ=',chi_tensor_loc(1,3)
             WRITE(unit_atoms,'(3(A,f15.6))') ' YX=',chi_tensor_loc(2,1),&
                  ' YY=',chi_tensor_loc(2,2),' YZ=',chi_tensor_loc(2,3)
             WRITE(unit_atoms,'(3(A,f15.6))') ' ZX=',chi_tensor_loc(3,1),&
                 ' ZY=',chi_tensor_loc(3,2),' ZZ=',chi_tensor_loc(3,3)
          ELSE
             WRITE(unit_atoms,'(T30,A)') " CHI from J in atomic units"
             WRITE(unit_atoms,'(3(A,f15.6))')' XX=',chi_tensor(1,1),&
                  &                          ' XY=',chi_tensor(1,2),&
                  &                          ' XZ=',chi_tensor(1,3)
             WRITE(unit_atoms,'(3(A,f15.6))')' YX=',chi_tensor(2,1),&
                  &                          ' YY=',chi_tensor(2,2),&
                  &                          ' YZ=',chi_tensor(2,3)
             WRITE(unit_atoms,'(3(A,f15.6))')' ZX=',chi_tensor(3,1),&
                  &                          ' ZY=',chi_tensor(3,2),&
                  &                          ' ZZ=',chi_tensor(3,3)
          ENDIF
          WRITE(unit_atoms,'(/T2,3(A,F15.5))')' PV1=',eig(1),&
               &                              ' PV2=',eig(2),&
               &                              ' PV3=',eig(3)
          WRITE (unit_atoms,'(T2,A,F15.5,10X,A,F15.5)')&
               '  iso=',chi_iso,  'aniso=',chi_aniso
          WRITE(unit_atoms,'(T30,A)') " Total CHI in atomic units"
          WRITE(unit_atoms,'(3(A,f15.6))')' XX=',chi_sym_tot(1,1),&
               &                          ' XY=',chi_sym_tot(1,2),&
               &                          ' XZ=',chi_sym_tot(1,3)
          WRITE(unit_atoms,'(3(A,f15.6))')' YX=',chi_sym_tot(2,1),&
               &                          ' YY=',chi_sym_tot(2,2),&
               &                          ' YZ=',chi_sym_tot(2,3)
          WRITE(unit_atoms,'(3(A,f15.6))')' ZX=',chi_sym_tot(3,1),&
               &                          ' ZY=',chi_sym_tot(3,2),&
               &                          ' ZZ=',chi_sym_tot(3,3)
          CALL dscal(3*3,nmr_env%chi_SI2ppmcgs,chi_sym_tot,1)
          WRITE(unit_atoms,'(T30,A)') " Total CHI in ppm-cgs units"
          WRITE(unit_atoms,'(3(A,f15.6))')' XX=',chi_sym_tot(1,1),&
               &                          ' XY=',chi_sym_tot(1,2),&
               &                          ' XZ=',chi_sym_tot(1,3)
          WRITE(unit_atoms,'(3(A,f15.6))')' YX=',chi_sym_tot(2,1),&
               &                          ' YY=',chi_sym_tot(2,2),&
               &                          ' YZ=',chi_sym_tot(2,3)
          WRITE(unit_atoms,'(3(A,f15.6))')' ZX=',chi_sym_tot(3,1),&
               &                          ' ZY=',chi_sym_tot(3,2),&
               &                          ' ZZ=',chi_sym_tot(3,3)
          WRITE(unit_atoms, '(/T2,3(A,F15.5))') &
               ' PV1=',nmr_env%chi_SI2ppmcgs*eig(1),&
               ' PV2=',nmr_env%chi_SI2ppmcgs*eig(2),&
               ' PV3=',nmr_env%chi_SI2ppmcgs*eig(3)
          WRITE(unit_atoms,'(T2,A,F15.5,10X,A,F15.5)')&
               '  iso=',nmr_env%chi_SI2ppmcgs*chi_iso, &
               'aniso=',nmr_env%chi_SI2ppmcgs*chi_aniso
       ENDIF
       
       CALL cp_print_key_finished_output(unit_atoms, logger,nmr_section,&
            &                            "PRINT%CHI_TENSOR", error=error)
    ENDIF ! print chi
    !
    ! Add the chi part to the shifts
    chemical_shift_tot = 0.0_dp
    DO ir = 1,nat_print
       iatom = cs_atom_list(ir)
       rpos(1:3) = particle_set(iatom)%r(1:3)
       atom_kind => particle_set(iatom)%atomic_kind
       CALL get_atomic_kind(atom_kind,name=name,element_symbol=element_symbol)
       DO i = 1,3
          DO j = 1,3
             chemical_shift_tot(i,j,ir) = chi_tensor(i,j) * nmr_env%chi_SI2shiftppm &
                  & + chemical_shift(i,j,iatom)
             IF(gapw) THEN
                chemical_shift_tot(i,j,ir) = chemical_shift_tot(i,j,ir)+chemical_shift_loc(i,j,iatom)
             ENDIF
          ENDDO
       ENDDO
    END DO  ! ir
    IF(output_unit > 0) THEN
       WRITE(output_unit,'(T2,A,E24.16)') 'CheckSum Shifts =',&
            & SQRT(DDOT(9*SIZE(chemical_shift_tot,3),chemical_shift_tot(1,1,1),1,&
            &           chemical_shift_tot(1,1,1),1))
    ENDIF
    !
    ! print shifts
    IF (BTEST(cp_print_key_should_output(logger%iter_info,nmr_section,&
              "PRINT%SHIELDING_TENSOR",error=error),cp_p_file)) THEN

       unit_atoms=cp_print_key_unit_nr(logger,nmr_section,"PRINT%SHIELDING_TENSOR",&
            &                          extension=".data",middle_name="SHIFT",&
            &                          log_filename=.FALSE.,error=error)

       nat_print = SIZE(cs_atom_list,1)
       IF(unit_atoms > 0) THEN
          WRITE(title,'(A,1X,I5)') "Shielding atom at atomic positions. # tensors printed ", nat_print
          WRITE(unit_atoms,'(T2,A)') title
          DO ir = 1,nat_print
             iatom = cs_atom_list(ir)
             rpos(1:3) = particle_set(iatom)%r(1:3)
             atom_kind => particle_set(iatom)%atomic_kind
             CALL get_atomic_kind(atom_kind,name=name,element_symbol=element_symbol)
             DO i = 1,3
                DO j = 1,3
                   shift_nosym(i,j) = chi_tensor(i,j) * nmr_env%chi_SI2shiftppm &
                        &           + chemical_shift(i,j,iatom)
                   IF(gapw) shift_nosym(i,j) = shift_nosym(i,j)+chemical_shift_loc(i,j,iatom)
                ENDDO
             ENDDO
             shift_sym = 0.5_dp*(shift_nosym+TRANSPOSE(shift_nosym))
             CALL  dcopy(9,shift_sym(1,1),1,shift_tmp(1,1),1)
             CALL  diamat_all(shift_tmp,eig_tmp)
             shift_iso   = (shift_nosym(1,1)+shift_nosym(2,2)+shift_nosym(3,3))/3.0_dp
             shift_aniso = eig_tmp(3)-0.5_dp*(eig_tmp(1)+eig_tmp(2))
             WRITE(unit_atoms,'(T2,I5,A,2X,A2,2X,3f15.6)') iatom,TRIM(name),element_symbol,rpos(1:3)
             !
             IF(gapw) THEN
                WRITE(unit_atoms,'(T1,A)') " SIGMA from SOFT J"
                WRITE(unit_atoms,'(A,3f15.6)') 'S_x = ',chemical_shift(1,:,iatom)
                WRITE(unit_atoms,'(A,3f15.6)') 'S_y = ',chemical_shift(2,:,iatom)
                WRITE(unit_atoms,'(A,3f15.6)') 'S_z = ',chemical_shift(3,:,iatom)
                WRITE(unit_atoms,'(T1,A)') " SIGMA from LOCAL J"
                WRITE(unit_atoms,'(A,3f15.6)') 'S_x = ',chemical_shift_loc(1,:,iatom)
                WRITE(unit_atoms,'(A,3f15.6)') 'S_y = ',chemical_shift_loc(2,:,iatom)
                WRITE(unit_atoms,'(A,3f15.6)') 'S_z = ',chemical_shift_loc(3,:,iatom)
                WRITE(unit_atoms,'(T1,2(A,F13.4))') " SIGMA TOTAL ISO = ",shift_iso,&
                     & "   ANISO = ",shift_aniso
                WRITE(unit_atoms,'(A,3f15.6)') 'S_x = ',shift_nosym(1, 1:3)
                WRITE(unit_atoms,'(A,3f15.6)') 'S_y = ',shift_nosym(2, 1:3)
                WRITE(unit_atoms,'(A,3f15.6)') 'S_z = ',shift_nosym(3, 1:3)
             ELSE
                WRITE(unit_atoms,'(T1,2(A,F13.4))') " SIGMA TOTAL ISO = ",shift_iso,&
                     & "   ANISO = ",shift_aniso
                WRITE(unit_atoms,'(A,3f15.6)') 'S_x = ',shift_nosym(1, 1:3)
                WRITE(unit_atoms,'(A,3f15.6)') 'S_y = ',shift_nosym(2, 1:3)
                WRITE(unit_atoms,'(A,3f15.6)') 'S_z = ',shift_nosym(3, 1:3)
             ENDIF
          ENDDO  ! ir
       ENDIF
       CALL cp_print_key_finished_output(unit_atoms,logger,nmr_section,&
            &                            "PRINT%SHIELDING_TENSOR",error=error)

       IF(nmr_env%do_nics) THEN
          !
          ! Add the chi part to the nics
          chemical_shift_nics_tot = 0.0_dp
          DO ir = 1,n_nics
             DO i = 1,3
                DO j = 1,3
                   ! here we need to multiply chi and the shift by some factors
                   chemical_shift_nics_tot(i,j,ir) = chi_tensor(i,j) * nmr_env%chi_SI2shiftppm &
                        &                          + chemical_shift_nics(i,j,ir)
                   IF(gapw) THEN
                      CALL stop_program(routineP,"GAPW and NICS doesnt work yet.")
                   ENDIF
                ENDDO
             ENDDO
          END DO  ! ir
          IF(output_unit > 0) THEN
             WRITE(output_unit,'(T2,A,E24.16)') 'CheckSum NICS =',&
                  & SQRT(DDOT(9*SIZE(chemical_shift_nics_tot,3),chemical_shift_nics_tot(1,1,1),&
                  &           1,chemical_shift_nics_tot(1,1,1),1))
          ENDIF
          !
          unit_nics=cp_print_key_unit_nr(logger,nmr_section,"PRINT%SHIELDING_TENSOR",&
               &                         extension=".data",middle_name="NICS",&
               &                         log_filename=.FALSE.,error=error)
          IF(unit_nics > 0) THEN
             WRITE(title,'(A,1X,I5)') "Shielding atom at atomic positions. # tensors printed ", n_nics
             WRITE(unit_nics,'(T2,A)') title
             DO ir = 1,n_nics
                WRITE(unit_nics,'(T2,I5,2X,3f15.6)') ir,r_nics(1:3,ir)
                IF(gapw) THEN
                   CALL stop_program(routineP,"GAPW and NICS doesnt work yet.")
                   !WRITE(unit_nics,'(3f15.6,5x,3f15.6)') chemical_shift_nics(1, 1:3, iatom),&
                   !     chemical_shift_nics_loc(1,1:3, iatom)
                   !WRITE(unit_nics,'(3f15.6,5x,3f15.6)') chemical_shift_nics(2, 1:3, iatom),&
                   !     chemical_shift_nics_loc(2,1:3, iatom)
                   !WRITE(unit_nics,'(3f15.6,5x,3f15.6)') chemical_shift_nics(3, 1:3, iatom),&
                   !     chemical_shift_nics_loc(3,1:3, iatom)
                ELSE
                   WRITE(unit_nics,'(3f15.6)') chemical_shift_nics_tot(1, 1:3, ir)
                   WRITE(unit_nics,'(3f15.6)') chemical_shift_nics_tot(2, 1:3, ir)
                   WRITE(unit_nics,'(3f15.6)') chemical_shift_nics_tot(3, 1:3, ir)
                ENDIF
             ENDDO
          ENDIF
          CALL cp_print_key_finished_output(unit_nics,logger,nmr_section,&
               &                            "PRINT%SHIELDING_TENSOR",error=error)
       ENDIF
    ENDIF  ! print shift
    !
    ! clean up
    DEALLOCATE(chemical_shift_tot, STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    IF(nmr_env%do_nics) THEN
       DEALLOCATE(chemical_shift_nics_tot,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ENDIF
    !
  END SUBROUTINE nmr_shift_print

END MODULE qs_linres_nmr_shift

