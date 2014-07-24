!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
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
  USE cell_methods,                    ONLY: cell_type,&
                                             pbc
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE f77_blas
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
                                             pw_transfer,&
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
                                             set_vecp,&
                                             set_vecp_rev
  USE qs_linres_types,                 ONLY: current_env_type,&
                                             get_current_env,&
                                             get_nmr_env,&
                                             jrho_atom_type,&
                                             nmr_env_type
  USE realspace_grid_types,            ONLY: realspace_grid_desc_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: get_limit
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  ! *** Public subroutines ***
  PUBLIC :: nmr_shift_print,&
            nmr_shift

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_linres_nmr_shift'

! **************

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param nmr_env ...
!> \param current_env ...
!> \param qs_env ...
!> \param iB ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE nmr_shift(nmr_env,current_env,qs_env,iB,error)

    TYPE(nmr_env_type)                       :: nmr_env
    TYPE(current_env_type)                   :: current_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    INTEGER, INTENT(IN)                      :: iB
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'nmr_shift', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, idir, idir2, idir3, &
                                                iiB, iiiB, ispin, istat, &
                                                natom, nspins
    LOGICAL                                  :: failure, gapw, &
                                                interpolate_shift
    REAL(dp)                                 :: scale_fac
    REAL(dp), DIMENSION(:, :, :), POINTER    :: chemical_shift, &
                                                chemical_shift_loc, &
                                                chemical_shift_nics, &
                                                chemical_shift_nics_loc
    TYPE(cell_type), POINTER                 :: cell
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type)                          :: pw_gspace_work, &
                                                shift_pw_rspace
    TYPE(pw_p_type), DIMENSION(:, :), &
      POINTER                                :: shift_pw_gspace
    TYPE(pw_p_type), POINTER                 :: jrho_gspace
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(realspace_grid_desc_type), POINTER  :: auxbas_rs_desc
    TYPE(section_vals_type), POINTER         :: nmr_section

    failure = .FALSE.
    CALL timeset(routineN,handle)

    NULLIFY(chemical_shift,chemical_shift_loc,chemical_shift_nics,chemical_shift_nics_loc,&
            cell,dft_control,pw_env,auxbas_rs_desc,auxbas_pw_pool,jrho_gspace,&
            pw_pools,particle_set)

    CALL get_qs_env(qs_env=qs_env,cell=cell,dft_control=dft_control,&
                    particle_set=particle_set,error=error)

    gapw = dft_control%qs_control%gapw
    natom = SIZE(particle_set,1)
    nspins = dft_control%nspins

    CALL get_nmr_env(nmr_env=nmr_env,chemical_shift=chemical_shift,&
                     chemical_shift_loc=chemical_shift_loc,&
                     chemical_shift_nics=chemical_shift_nics,&
                     chemical_shift_nics_loc=chemical_shift_nics_loc,&
                     interpolate_shift=interpolate_shift,&
                     error=error)

    CALL get_qs_env(qs_env=qs_env,pw_env=pw_env,error=error)
    CALL pw_env_get(pw_env,auxbas_rs_desc=auxbas_rs_desc,&
                    auxbas_pw_pool=auxbas_pw_pool,pw_pools=pw_pools,&
                    error=error)
    !
    !
    nmr_section => section_vals_get_subs_vals(qs_env%input, &
         & "PROPERTIES%LINRES%NMR",error=error)
    !
    ! Initialize
    ! Allocate grids for the calculation of jrho and the shift
    ALLOCATE(shift_pw_gspace(3,nspins),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DO ispin = 1,nspins
       DO idir = 1,3
          CALL pw_pool_create_pw(auxbas_pw_pool,shift_pw_gspace(idir,ispin)%pw,&
                                 use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,&
                                 error=error)
          CALL pw_zero(shift_pw_gspace(idir,ispin)%pw,error=error)
       ENDDO
    ENDDO
    !
    !
    CALL set_vecp(iB,iiB,iiiB)
    !
    CALL pw_pool_create_pw(auxbas_pw_pool,pw_gspace_work%pw,&
                           use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,&
                           error=error)
    CALL pw_zero(pw_gspace_work%pw,error=error)
    DO ispin = 1,nspins
       !
       DO idir = 1,3
          jrho_gspace => current_env%jrho1_set(idir)%rho%rho_g(ispin)
          ! Field gradient
          ! loop over the Gvec  components: x,y,z
          DO idir2 = 1,3
             IF(idir /= idir2) THEN
                ! in reciprocal space multiply (G_idir2(i)/G(i)^2)J_(idir)(G(i))
                CALL mult_G_ov_G2_grid(cell,auxbas_pw_pool,jrho_gspace,&
                                       pw_gspace_work,idir2,0.0_dp,error=error)
                !
                ! scale and add to the correct component of the shift column
                CALL set_vecp_rev(idir,idir2,idir3)
                scale_fac=fac_vecp(idir3,idir2,idir)
                CALL pw_axpy(pw_gspace_work%pw,shift_pw_gspace(idir3,ispin)%pw,scale_fac,error=error)
             ENDIF
          ENDDO
          !
       ENDDO ! idir
    ENDDO ! ispin
    !
    CALL pw_pool_give_back_pw(auxbas_pw_pool,pw_gspace_work%pw,error=error)
    !
    ! compute shildings
    IF(interpolate_shift) THEN
       CALL pw_pool_create_pw(auxbas_pw_pool,shift_pw_rspace%pw,&
                              use_data=REALDATA3D,in_space=REALSPACE,error=error)
       DO ispin = 1,nspins
          DO idir = 1,3
             ! Here first G->R and then interpolation to get the shifts.
             ! The interpolation doesnt work in parallel yet.
             ! The choice between both methods should be left to the user.
             CALL pw_transfer(shift_pw_gspace(idir,ispin)%pw,shift_pw_rspace%pw,error=error)
             CALL interpolate_shift_pwgrid(nmr_env,pw_env,particle_set,cell,shift_pw_rspace,&
                                           iB,idir,nmr_section,error=error)
          ENDDO
       ENDDO
       CALL pw_pool_give_back_pw(auxbas_pw_pool,shift_pw_rspace%pw,error=error)
    ELSE
       DO ispin = 1,nspins
          DO idir = 1,3
             ! Here the shifts are computed from summation of the coeff on the G-grip .
             CALL gsum_shift_pwgrid(nmr_env,particle_set,cell,&
                                    shift_pw_gspace(idir,ispin),iB,idir,&
                                    error=error)
          ENDDO
       ENDDO
    ENDIF
    !
    IF(gapw) THEN
       DO idir = 1,3
          ! Finally the radial functions are multiplied by the YLM and properly summed
          ! The resulting array is J on the local grid. One array per atom.
          ! Local contributions by numerical integration over the spherical grids
          CALL nmr_shift_gapw(nmr_env,current_env,qs_env,iB,idir,error=error)
       ENDDO ! idir
    ENDIF
    !
    ! Dellocate grids for the calculation of jrho and the shift
    DO ispin = 1,nspins
       DO idir = 1,3
          CALL pw_pool_give_back_pw(auxbas_pw_pool,shift_pw_gspace(idir,ispin)%pw,error=error)
       ENDDO
    ENDDO
    DEALLOCATE(shift_pw_gspace,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    !
    ! Finalize
    CALL timestop(handle)
    !
  END SUBROUTINE nmr_shift


! *****************************************************************************
!> \brief ...
!> \param nmr_env ...
!> \param current_env ...
!> \param qs_env ...
!> \param iB ...
!> \param idir ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE nmr_shift_gapw(nmr_env,current_env,qs_env,iB,idir,error)

    TYPE(nmr_env_type)                       :: nmr_env
    TYPE(current_env_type)                   :: current_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    INTEGER, INTENT(IN)                      :: IB, idir
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'nmr_shift_gapw', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, ia, iat, iatom, idir2_1, idir3_1, ikind, ir, ira, &
      ispin, istat, j, jatom, mepos, n_nics, na, natom, natom_local, &
      natom_tot, nkind, nnics_local, nr, nra, nspins, num_pe, output_unit
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: list_j, list_nics_j
    INTEGER, DIMENSION(2)                    :: bo
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: do_nics, failure, paw_atom
    REAL(dp) :: ddiff, dist, dum, itegrated_jrho, r_jatom(3), rdiff(3), &
      rij(3), s_1, s_2, scale_fac_1, shift_gapw_radius
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: j_grid
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: cs_loc_tmp, cs_nics_loc_tmp, &
                                                dist_ij, dist_nics_ij, r_grid
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

    CALL timeset(routineN,handle)
    !
    failure =.FALSE.
    NULLIFY(atomic_kind_set,cell,dft_control,para_env,particle_set,&
            chemical_shift_loc,chemical_shift_nics_loc,jrho1_atom_set,&
            jrho1_atom,r_nics,jrho_h_grid,jrho_s_grid,atom_kind,&
            atom_list,grid_atom,harmonics,logger)
    !
    logger => cp_error_get_logger(error)
    output_unit= cp_logger_get_default_io_unit(logger)
    !
    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    cell=cell,&
                    dft_control=dft_control,&
                    para_env=para_env,&
                    particle_set=particle_set,&
                    error=error)

    CALL get_nmr_env(nmr_env=nmr_env,&
                     chemical_shift_loc=chemical_shift_loc,&
                     chemical_shift_nics_loc=chemical_shift_nics_loc,&
                     shift_gapw_radius=shift_gapw_radius,&
                     n_nics=n_nics,&
                     r_nics=r_nics,&
                     do_nics=do_nics,&
                     error=error)

    CALL get_current_env(current_env=current_env,&
                         jrho1_atom_set=jrho1_atom_set,&
                         error=error)
    !
    nkind = SIZE(atomic_kind_set,1)
    natom_tot = SIZE(particle_set,1)
    nspins = dft_control%nspins
    itegrated_jrho = 0.0_dp
    !
    idir2_1 = MODULO(idir,3)+1
    idir3_1 = MODULO(idir+1,3)+1
    scale_fac_1 = fac_vecp(idir3_1,idir2_1,idir)
    !
    ALLOCATE(cs_loc_tmp(3,natom_tot),list_j(natom_tot),&
             dist_ij(3,natom_tot),STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    cs_loc_tmp = 0.0_dp
    IF(do_nics) THEN
       ALLOCATE(cs_nics_loc_tmp(3,n_nics),list_nics_j(n_nics),&
                dist_nics_ij(3,n_nics),STAT=istat)
       CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
       cs_nics_loc_tmp = 0.0_dp
    ENDIF
    !
    ! Loop over atoms to collocate the current on each atomic grid, JA
    ! Per each JA, loop over the points where the shift needs to be computed
    DO ikind = 1,nkind

       NULLIFY(atom_kind,atom_list,grid_atom,harmonics)
       atom_kind => atomic_kind_set(ikind)
       CALL get_atomic_kind(atomic_kind=atom_kind,&
                            atom_list=atom_list,&
                            natom=natom,&
                            paw_atom=paw_atom,&
                            harmonics=harmonics,&
                            grid_atom=grid_atom)
       !
       na = grid_atom%ng_sphere
       nr = grid_atom%nr
       nra = nr*na
       ALLOCATE(r_grid(3,nra),j_grid(nra))
       ira = 1
       DO ia = 1,na
          DO ir = 1,nr
             r_grid(:,ira) = grid_atom%rad(ir)*harmonics%a(:,ia)
             ira = ira  + 1
          ENDDO
       ENDDO
       !
       ! Quick cycle if needed
       IF(paw_atom) THEN
          !
          ! Distribute the atoms of this kind
          num_pe = para_env%num_pe
          mepos  = para_env%mepos
          bo = get_limit( natom, num_pe, mepos )
          !
          DO iat = bo(1),bo(2)
             iatom = atom_list(iat)
             !
             ! find all the atoms within the radius
             natom_local = 0
             DO jatom = 1,natom_tot
                rij(:) = pbc(particle_set(iatom)%r,particle_set(jatom)%r,cell)
                dist = SQRT(rij(1)*rij(1)+rij(2)*rij(2)+rij(3)*rij(3))
                IF(dist.LE.shift_gapw_radius) THEN
                   natom_local = natom_local + 1
                   list_j(natom_local) = jatom
                   dist_ij(:,natom_local) = rij(:)
                ENDIF
             ENDDO
             !
             ! ... also for nics
             IF(do_nics) THEN
                nnics_local = 0
                DO jatom = 1,n_nics
                   r_jatom(:) = r_nics(:,jatom)
                   rij(:) = pbc(particle_set(iatom)%r,r_jatom,cell)
                   dist = SQRT(rij(1)*rij(1)+rij(2)*rij(2)+rij(3)*rij(3))
                   IF(dist.LE.shift_gapw_radius) THEN
                      nnics_local = nnics_local + 1
                      list_nics_j(nnics_local) = jatom
                      dist_nics_ij(:,nnics_local) = rij(:)
                   ENDIF
                ENDDO
             ENDIF
             !
             NULLIFY(jrho1_atom,jrho_h_grid,jrho_s_grid)
             jrho1_atom => jrho1_atom_set(iatom)
             !
             DO ispin = 1,nspins
                jrho_h_grid => jrho1_atom%jrho_vec_rad_h(idir,ispin)%r_coef
                jrho_s_grid => jrho1_atom%jrho_vec_rad_s(idir,ispin)%r_coef
                !
                ! loop over the atoms neighbors of iatom in terms of the current density
                ! for each compute the contribution to the shift coming from the
                ! local current density at iatom
                ira = 1
                DO ia = 1,na
                   DO ir = 1,nr
                      j_grid(ira) = (jrho_h_grid(ir,ia) - jrho_s_grid(ir,ia)) * grid_atom%weight(ia,ir)
                      itegrated_jrho = itegrated_jrho + j_grid(ira)
                      ira = ira + 1
                   ENDDO
                ENDDO
                !
                DO j = 1,natom_local
                   jatom = list_j(j)
                   rij(:) = dist_ij(:,j)
                   !
                   s_1 = 0.0_dp
                   s_2 = 0.0_dp
                   DO ira = 1,nra
                      !
                      rdiff(:) = rij(:) - r_grid(:,ira)
                      ddiff = SQRT( rdiff(1)*rdiff(1) + rdiff(2)*rdiff(2) + rdiff(3)*rdiff(3) )
                      IF(ddiff.GT.1.0E-12_dp) THEN
                         dum = scale_fac_1 * j_grid(ira) / (ddiff*ddiff*ddiff)
                         s_1 = s_1 + rdiff(idir2_1) * dum
                         s_2 = s_2 + rdiff(idir3_1) * dum
                      ENDIF ! ddiff
                   ENDDO ! ira
                   cs_loc_tmp(idir3_1,jatom) = cs_loc_tmp(idir3_1,jatom) + s_1
                   cs_loc_tmp(idir2_1,jatom) = cs_loc_tmp(idir2_1,jatom) - s_2
                ENDDO ! j
                !
                IF(do_nics) THEN
                   DO j = 1,nnics_local
                      jatom = list_nics_j(j)
                      rij(:) = dist_nics_ij(:,j)
                      !
                      s_1 = 0.0_dp
                      s_2 = 0.0_dp
                      DO ira = 1,nra
                         !
                         rdiff(:) = rij(:) - r_grid(:,ira)
                         ddiff = SQRT( rdiff(1)*rdiff(1) + rdiff(2)*rdiff(2) + rdiff(3)*rdiff(3) )
                         IF(ddiff.GT.1.0E-12_dp) THEN
                            dum = scale_fac_1 * j_grid(ira) / (ddiff*ddiff*ddiff)
                            s_1 = s_1 + rdiff(idir2_1) * dum
                            s_2 = s_2 + rdiff(idir3_1) * dum
                         ENDIF ! ddiff
                      ENDDO ! ira
                      cs_nics_loc_tmp(idir3_1,jatom) = cs_nics_loc_tmp(idir3_1,jatom) + s_1
                      cs_nics_loc_tmp(idir2_1,jatom) = cs_nics_loc_tmp(idir2_1,jatom) - s_2
                   ENDDO ! j
                ENDIF ! do_nics
             ENDDO ! ispin
          ENDDO ! iat
       ENDIF
       DEALLOCATE(r_grid,j_grid)
    ENDDO ! ikind
    !
    !
    CALL mp_sum(itegrated_jrho,para_env%group)
    IF(output_unit>0) THEN
       WRITE(output_unit,'(T2,A,E24.16)') 'Integrated local j_'&
            &//ACHAR(idir+119)//ACHAR(iB+119)//'(r)=',itegrated_jrho
    ENDIF
    !
    CALL mp_sum(cs_loc_tmp,para_env%group)
    chemical_shift_loc(:,iB,:) = chemical_shift_loc(:,iB,:) &
         & - nmr_env%shift_factor_gapw * cs_loc_tmp(:,:) / 2.0_dp
    !
    DEALLOCATE(cs_loc_tmp,list_j,dist_ij,STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    !
    IF(do_nics) THEN
       CALL mp_sum(cs_nics_loc_tmp,para_env%group)
       chemical_shift_nics_loc(:,iB,:) = chemical_shift_nics_loc(:,iB,:) &
            & - nmr_env%shift_factor_gapw * cs_nics_loc_tmp(:,:) / 2.0_dp
       !
       DEALLOCATE(cs_nics_loc_tmp,list_nics_j,dist_nics_ij,STAT=istat)
       CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    ENDIF
    !
    CALL timestop(handle)
    !
  END SUBROUTINE nmr_shift_gapw

! *****************************************************************************
!> \brief interpolate the shift calculated on the PW grid in order to ger
!>       the value on arbitrary points in real space
!> \param nmr_env ...
!> \param pw_env ...
!> \param particle_set for the atomic position
!> \param cell to take into account the pbs, and to have the volume
!> \param shift_pw_rspace specific component of the shift tensor on the pw grid
!> \param i_B component of the magnetic field for which the shift is calculated (row)
!> \param idir component of the vector \int_{r}[ ((r-r') x j(r))/|r-r'|^3 ] = Bind(r')
!> \param nmr_section ...
!> \param error ...
!> \param nme_env to get the shift tensor and the list of additional points UNMATCHED_PROCEDURE_ARGUMENT: please check 
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

    INTEGER                                  :: aint_precond, handle, iat, &
                                                iatom, max_iter, n_nics, &
                                                natom, precond_kind
    INTEGER, DIMENSION(:), POINTER           :: cs_atom_list
    LOGICAL                                  :: do_nics, failure, success
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
    !
    CALL timeset(routineN,handle)
    !
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
           n_nics=n_nics,r_nics=r_nics,&
           do_nics=do_nics,error=error)

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

      IF(do_nics) THEN
         DO iatom = 1,n_nics
            ra(:) = r_nics(:,iatom)
            R_iatom = pbc(ra,cell)
            shift_val = Eval_Interp_Spl3_pbc(R_iatom,shiftspl%pw,error)
            chemical_shift_nics(idir,i_B,iatom)= chemical_shift_nics(idir,i_B,iatom)+&
                                                 nmr_env%shift_factor * twopi**2 * shift_val
         END DO
      END IF

      CALL pw_pool_give_back_pw(auxbas_pw_pool,shiftspl%pw,error=error)

    END IF
    !
    CALL timestop(handle)
    !
  END SUBROUTINE interpolate_shift_pwgrid

! *****************************************************************************
!> \brief ...
!> \param nmr_env ...
!> \param particle_set ...
!> \param cell ...
!> \param shift_pw_gspace ...
!> \param i_B ...
!> \param idir ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE gsum_shift_pwgrid(nmr_env,particle_set,cell,shift_pw_gspace,&
                               i_B,idir,error)
    TYPE(nmr_env_type)                       :: nmr_env
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(pw_p_type)                          :: shift_pw_gspace
    INTEGER, INTENT(IN)                      :: i_B, idir
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'gsum_shift_pwgrid', &
      routineP = moduleN//':'//routineN

    COMPLEX(dp)                              :: cplx
    INTEGER                                  :: handle, iat, iatom, n_nics, &
                                                natom
    INTEGER, DIMENSION(:), POINTER           :: cs_atom_list
    LOGICAL                                  :: do_nics, failure
    REAL(dp)                                 :: R_iatom(3), ra(3)
    REAL(dp), DIMENSION(:, :), POINTER       :: r_nics
    REAL(dp), DIMENSION(:, :, :), POINTER    :: chemical_shift, &
                                                chemical_shift_nics

!

    failure = .FALSE.
    !
    CALL timeset(routineN,handle)
    !
    NULLIFY(cs_atom_list,chemical_shift,chemical_shift_nics,r_nics)
    CPPrecondition(ASSOCIATED(shift_pw_gspace%pw),cp_failure_level,routineP,error,failure)
    IF(.NOT. failure) THEN
       !
       CALL get_nmr_env(nmr_env=nmr_env, cs_atom_list=cs_atom_list,&
                        chemical_shift=chemical_shift,&
                        chemical_shift_nics=chemical_shift_nics,&
                        n_nics=n_nics,r_nics=r_nics,do_nics=do_nics,&
                        error=error)
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
          CALL gsumr(R_iatom,shift_pw_gspace%pw,cplx)
          chemical_shift(idir,i_B,iatom)= chemical_shift(idir,i_B,iatom)+&
                 nmr_env%shift_factor * twopi**2 * REAL(cplx,dp)
       ENDDO
       !
       ! compute nics
       IF(do_nics) THEN
          DO iat = 1,n_nics
             ra = pbc(r_nics(:,iat),cell)
             CALL gsumr(ra,shift_pw_gspace%pw,cplx)
             chemical_shift_nics(idir,i_B,iat)= chemical_shift_nics(idir,i_B,iat)+&
                    nmr_env%shift_factor * twopi**2 * REAL(cplx,dp)
          ENDDO
       ENDIF
    ENDIF

    CALL timestop(handle)

  END SUBROUTINE gsum_shift_pwgrid

! *****************************************************************************
!> \brief ...
!> \param r ...
!> \param pw ...
!> \param cplx ...
! *****************************************************************************
  SUBROUTINE gsumr(r,pw,cplx)
    REAL(dp), INTENT(IN)                     :: r(3)
    TYPE(pw_type), POINTER                   :: pw
    COMPLEX(dp)                              :: cplx

    COMPLEX(dp)                              :: rg
    INTEGER                                  :: ig
    TYPE(pw_grid_type), POINTER              :: grid

    grid => pw%pw_grid
    cplx = CMPLX(0.0_dp,0.0_dp,KIND=dp)
    DO ig = grid%first_gne0,grid%ngpts_cut_local
       rg = ( grid%g(1,ig)*r(1) + grid%g(2,ig)*r(2) + grid%g(3,ig)*r(3) ) * gaussi
       cplx = cplx + pw%cc(ig)*EXP(rg)
    ENDDO
    IF(grid%have_g0) cplx = cplx + pw%cc(1)
    CALL mp_sum(cplx,grid%para%group)
  END SUBROUTINE gsumr

! *****************************************************************************
!> \brief Shielding tensor and Chi are printed into a file
!>       if required from input
!>       It is possible to print only for a subset of atoms or
!>       or points in non-ionic positions
!> \param nmr_env ...
!> \param current_env ...
!> \param qs_env ...
!> \param error ...
!> \author MI
! *****************************************************************************
  SUBROUTINE nmr_shift_print(nmr_env,current_env,qs_env,error)
    TYPE(nmr_env_type)                       :: nmr_env
    TYPE(current_env_type)                   :: current_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'nmr_shift_print', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=2)                         :: element_symbol
    CHARACTER(LEN=default_string_length)     :: name, title
    INTEGER                                  :: iatom, ir, istat, n_nics, &
                                                nat_print, natom, &
                                                output_unit, unit_atoms, &
                                                unit_nics
    INTEGER, DIMENSION(:), POINTER           :: cs_atom_list
    LOGICAL                                  :: do_nics, failure, gapw
    REAL(dp) :: chi_aniso, chi_iso, chi_sym_tot(3,3), chi_tensor(3,3,2), &
      chi_tensor_loc(3,3,2), chi_tensor_loc_tmp(3,3), chi_tensor_tmp(3,3), &
      chi_tmp(3,3), eig(3), rpos(3), shift_aniso, shift_iso, &
      shift_sym_tot(3,3)
    REAL(dp), DIMENSION(:, :), POINTER       :: r_nics
    REAL(dp), DIMENSION(:, :, :), POINTER    :: cs, cs_loc, cs_nics, &
                                                cs_nics_loc, cs_nics_tot, &
                                                cs_tot
    REAL(dp), EXTERNAL                       :: DDOT
    TYPE(atomic_kind_type), POINTER          :: atom_kind
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(section_vals_type), POINTER         :: nmr_section

    failure = .FALSE.

    NULLIFY(cs,cs_nics,r_nics,cs_loc,cs_nics_loc,logger,particle_set,atom_kind,dft_control)

    logger => cp_error_get_logger(error)
    output_unit= cp_logger_get_default_io_unit(logger)

    nmr_section => section_vals_get_subs_vals(qs_env%input, &
                                                 "PROPERTIES%LINRES%NMR",error=error)

    CALL get_nmr_env(nmr_env=nmr_env, &
                     chemical_shift=cs, &
                     chemical_shift_nics=cs_nics, &
                     chemical_shift_loc=cs_loc, &
                     chemical_shift_nics_loc=cs_nics_loc, &
                     cs_atom_list=cs_atom_list, &
                     n_nics=n_nics, &
                     r_nics=r_nics, &
                     do_nics=do_nics, &
                     error=error)
    !
    CALL get_current_env(current_env=current_env,&
                         chi_tensor=chi_tensor,&
                         chi_tensor_loc=chi_tensor_loc,&
                         error=error)
    !
    ! multiply by the appropriate factor
    chi_tensor_tmp(:,:)     = 0.0_dp
    chi_tensor_loc_tmp(:,:) = 0.0_dp
    chi_tensor_tmp(:,:) = ( chi_tensor(:,:,1) + chi_tensor(:,:,2) ) * nmr_env%chi_factor
    !chi_tensor_loc_tmp(:,:) = (chi_tensor_loc(:,:,1) + chi_tensor_loc(:,:,2)) * here there is another factor
    !
    CALL get_qs_env(qs_env=qs_env, &
                    dft_control=dft_control, &
                    particle_set=particle_set, &
                    error=error)

    natom = SIZE(particle_set,1)
    gapw = dft_control%qs_control%gapw
    nat_print=SIZE(cs_atom_list,1)

    ALLOCATE(cs_tot(3,3,nat_print),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    IF(do_nics) THEN
       ALLOCATE(cs_nics_tot(3,3,n_nics),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ENDIF
    ! Finalize Chi calculation
    ! Symmetrize
    chi_sym_tot(:,:) = ( chi_tensor_tmp(:,:) + TRANSPOSE(chi_tensor_tmp(:,:)) ) / 2.0_dp
    IF(gapw) THEN
       chi_sym_tot(:,:) = chi_sym_tot(:,:) &
            & + ( chi_tensor_loc_tmp(:,:) + TRANSPOSE(chi_tensor_loc_tmp(:,:)) ) / 2.0_dp
    ENDIF
    chi_tmp(:,:) = chi_sym_tot(:,:)
    CALL diamat_all(chi_tmp,eig,error=error)
    chi_iso   = (eig(1)+ eig(2)+eig(3))/3.0_dp
    chi_aniso =  eig(3)-(eig(2)+eig(1))/2.0_dp
    !
    IF(output_unit > 0) THEN
       WRITE(output_unit,'(T2,A,E14.6)') 'CheckSum Chi =',&
              SQRT(DDOT(9,chi_tensor_tmp(1,1),1,chi_tensor_tmp(1,1),1))
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
          WRITE(unit_atoms,'(T1,A)') " CHI from SOFT J in 10^-30 J/T^2 units"
          WRITE(unit_atoms,'(3(A,f10.4))') '  XX = ',chi_tensor_tmp(1,1),&
               &                           '  XY = ',chi_tensor_tmp(1,2),&
               &                           '  XZ = ',chi_tensor_tmp(1,3)
          WRITE(unit_atoms,'(3(A,f10.4))') '  YX = ',chi_tensor_tmp(2,1),&
               &                           '  YY = ',chi_tensor_tmp(2,2),&
               &                           '  YZ = ',chi_tensor_tmp(2,3)
          WRITE(unit_atoms,'(3(A,f10.4))') '  ZX = ',chi_tensor_tmp(3,1),&
               &                           '  ZY = ',chi_tensor_tmp(3,2),&
               &                           '  ZZ = ',chi_tensor_tmp(3,3)
          IF(gapw) THEN
             WRITE(unit_atoms,'(T1,A)') " CHI from LOCAL J in 10^-30 J/T^2 units"
             WRITE(unit_atoms,'(3(A,f10.4))') '  XX = ',chi_tensor_loc_tmp(1,1),&
                  &                           '  XY = ',chi_tensor_loc_tmp(1,2),&
                  &                           '  XZ = ',chi_tensor_loc_tmp(1,3)
             WRITE(unit_atoms,'(3(A,f10.4))') '  YX = ',chi_tensor_loc_tmp(2,1),&
                  &                           '  YY = ',chi_tensor_loc_tmp(2,2),&
                  &                           '  YZ = ',chi_tensor_loc_tmp(2,3)
             WRITE(unit_atoms,'(3(A,f10.4))') '  ZX = ',chi_tensor_loc_tmp(3,1),&
                  &                           '  ZY = ',chi_tensor_loc_tmp(3,2),&
                  &                           '  ZZ = ',chi_tensor_loc_tmp(3,3)
          ENDIF
          WRITE(unit_atoms,'(T1,A)') " Total CHI in 10^-30 J/T^2 units"
          WRITE(unit_atoms,'(3(A,f10.4))')'  XX = ',chi_sym_tot(1,1),&
               &                          '  XY = ',chi_sym_tot(1,2),&
               &                          '  XZ = ',chi_sym_tot(1,3)
          WRITE(unit_atoms,'(3(A,f10.4))')'  YX = ',chi_sym_tot(2,1),&
               &                          '  YY = ',chi_sym_tot(2,2),&
               &                          '  YZ = ',chi_sym_tot(2,3)
          WRITE(unit_atoms,'(3(A,f10.4))')'  ZX = ',chi_sym_tot(3,1),&
               &                          '  ZY = ',chi_sym_tot(3,2),&
               &                          '  ZZ = ',chi_sym_tot(3,3)
          chi_sym_tot(:,:) = chi_sym_tot(:,:) * nmr_env%chi_SI2ppmcgs
          WRITE(unit_atoms,'(T1,A)') " Total CHI in ppm-cgs units"
          WRITE(unit_atoms,'(3(A,f10.4))')'  XX = ',chi_sym_tot(1,1),&
               &                          '  XY = ',chi_sym_tot(1,2),&
               &                          '  XZ = ',chi_sym_tot(1,3)
          WRITE(unit_atoms,'(3(A,f10.4))')'  YX = ',chi_sym_tot(2,1),&
               &                          '  YY = ',chi_sym_tot(2,2),&
               &                          '  YZ = ',chi_sym_tot(2,3)
          WRITE(unit_atoms,'(3(A,f10.4))')'  ZX = ',chi_sym_tot(3,1),&
               &                          '  ZY = ',chi_sym_tot(3,2),&
               &                          '  ZZ = ',chi_sym_tot(3,3)
          WRITE(unit_atoms, '(/T1,3(A,f10.4))') &
               '  PV1=',nmr_env%chi_SI2ppmcgs*eig(1),&
               '  PV2=',nmr_env%chi_SI2ppmcgs*eig(2),&
               '  PV3=',nmr_env%chi_SI2ppmcgs*eig(3)
          WRITE(unit_atoms,'(T1,A,F10.4,10X,A,F10.4)')&
               '  ISO=',nmr_env%chi_SI2ppmcgs*chi_iso, &
               'ANISO=',nmr_env%chi_SI2ppmcgs*chi_aniso
       ENDIF

       CALL cp_print_key_finished_output(unit_atoms, logger,nmr_section,&
            &                            "PRINT%CHI_TENSOR", error=error)
    ENDIF ! print chi
    !
    ! Add the chi part to the shifts
    cs_tot(:,:,:) = 0.0_dp
    DO ir = 1,nat_print
       iatom = cs_atom_list(ir)
       rpos(1:3) = particle_set(iatom)%r(1:3)
       atom_kind => particle_set(iatom)%atomic_kind
       CALL get_atomic_kind(atom_kind,name=name,element_symbol=element_symbol)
       cs_tot(:,:,ir) = chi_tensor_tmp(:,:) * nmr_env%chi_SI2shiftppm + cs(:,:,iatom)
       IF(gapw) cs_tot(:,:,ir) = cs_tot(:,:,ir) + cs_loc(:,:,iatom)
    ENDDO  ! ir
    IF(output_unit > 0) THEN
       WRITE(output_unit,'(T2,A,E14.6)') 'CheckSum Shifts =',&
              SQRT(DDOT(9*SIZE(cs_tot,3),cs_tot(1,1,1),1,cs_tot(1,1,1),1))
    ENDIF
    !
    ! print shifts
    IF (BTEST(cp_print_key_should_output(logger%iter_info,nmr_section,&
              "PRINT%SHIELDING_TENSOR",error=error),cp_p_file)) THEN

       unit_atoms=cp_print_key_unit_nr(logger,nmr_section,"PRINT%SHIELDING_TENSOR",&
                                       extension=".data",middle_name="SHIFT",&
                                       log_filename=.FALSE.,error=error)

       nat_print = SIZE(cs_atom_list,1)
       IF(unit_atoms > 0) THEN
          WRITE(title,'(A,1X,I5)') "Shielding atom at atomic positions. # tensors printed ", nat_print
          WRITE(unit_atoms,'(T2,A)') title
          DO ir = 1,nat_print
             iatom = cs_atom_list(ir)
             rpos(1:3) = particle_set(iatom)%r(1:3)
             atom_kind => particle_set(iatom)%atomic_kind
             CALL get_atomic_kind(atom_kind,name=name,element_symbol=element_symbol)
             shift_sym_tot(:,:) = 0.5_dp * (cs_tot(:,:,ir) + TRANSPOSE(cs_tot(:,:,ir)))
             CALL  diamat_all(shift_sym_tot,eig,error=error)
             shift_iso   = (eig(1) + eig(2) + eig(3)) / 3.0_dp
             shift_aniso =  eig(3) - (eig(2) + eig(1)) / 2.0_dp
             !
             WRITE(unit_atoms,'(T2,I5,A,2X,A2,2X,3f15.6)') iatom,TRIM(name),element_symbol,rpos(1:3)
             !
             IF(gapw) THEN
                WRITE(unit_atoms,'(T1,A)') " SIGMA from SOFT J"
                WRITE(unit_atoms,'(3(A,f10.4))') '  XX = ',cs(1,1,iatom),&
                     &                           '  XY = ',cs(1,2,iatom),&
                     &                           '  XZ = ',cs(1,3,iatom)
                WRITE(unit_atoms,'(3(A,f10.4))') '  YX = ',cs(2,1,iatom),&
                     &                           '  YY = ',cs(2,2,iatom),&
                     &                           '  YZ = ',cs(2,3,iatom)
                WRITE(unit_atoms,'(3(A,f10.4))') '  ZX = ',cs(3,1,iatom),&
                     &                           '  ZY = ',cs(3,2,iatom),&
                     &                           '  ZZ = ',cs(3,3,iatom)
                WRITE(unit_atoms,'(T1,A)') " SIGMA from LOCAL J"
                WRITE(unit_atoms,'(3(A,f10.4))') '  XX = ',cs_loc(1,1,iatom),&
                     &                           '  XY = ',cs_loc(1,2,iatom),&
                     &                           '  XZ = ',cs_loc(1,3,iatom)
                WRITE(unit_atoms,'(3(A,f10.4))') '  YX = ',cs_loc(2,1,iatom),&
                     &                           '  YY = ',cs_loc(2,2,iatom),&
                     &                           '  YZ = ',cs_loc(2,3,iatom)
                WRITE(unit_atoms,'(3(A,f10.4))') '  ZX = ',cs_loc(3,1,iatom),&
                     &                           '  ZY = ',cs_loc(3,2,iatom),&
                     &                           '  ZZ = ',cs_loc(3,3,iatom)
             ENDIF
             WRITE(unit_atoms,'(T1,A)') " SIGMA TOTAL"
             WRITE(unit_atoms,'(3(A,f10.4))') '  XX = ',cs_tot(1,1,ir),&
                  &                           '  XY = ',cs_tot(1,2,ir),&
                  &                           '  XZ = ',cs_tot(1,3,ir)
             WRITE(unit_atoms,'(3(A,f10.4))') '  YX = ',cs_tot(2,1,ir),&
                  &                           '  YY = ',cs_tot(2,2,ir),&
                  &                           '  YZ = ',cs_tot(2,3,ir)
             WRITE(unit_atoms,'(3(A,f10.4))') '  ZX = ',cs_tot(3,1,ir),&
                  &                           '  ZY = ',cs_tot(3,2,ir),&
                  &                           '  ZZ = ',cs_tot(3,3,ir)
             WRITE(unit_atoms,'(T1,2(A,f12.4))') '  ISOTROPY = ',shift_iso,&
                  &                            '  ANISOTROPY = ',shift_aniso
          ENDDO ! ir
       ENDIF
       CALL cp_print_key_finished_output(unit_atoms,logger,nmr_section,&
            &                            "PRINT%SHIELDING_TENSOR",error=error)

       IF(do_nics) THEN
          !
          ! Add the chi part to the nics
          cs_nics_tot(:,:,:) = 0.0_dp
          DO ir = 1,n_nics
             cs_nics_tot(:,:,ir) = chi_tensor_tmp(:,:) * nmr_env%chi_SI2shiftppm + cs_nics(:,:,ir)
             IF(gapw) cs_nics_tot(:,:,ir) = cs_nics_tot(:,:,ir) + cs_nics_loc(:,:,ir)
          ENDDO  ! ir
          IF(output_unit > 0) THEN
             WRITE(output_unit,'(T2,A,E14.6)') 'CheckSum NICS =',&
                    SQRT(DDOT(9*SIZE(cs_nics_tot,3),cs_nics_tot(1,1,1),1,cs_nics_tot(1,1,1),1))
          ENDIF
          !
          unit_nics=cp_print_key_unit_nr(logger,nmr_section,"PRINT%SHIELDING_TENSOR",&
                                         extension=".data",middle_name="NICS",&
                                         log_filename=.FALSE.,error=error)
          IF(unit_nics > 0) THEN
             WRITE(title,'(A,1X,I5)') "Shielding at nics positions. # tensors printed ", n_nics
             WRITE(unit_nics,'(T2,A)') title
             DO ir = 1,n_nics
                shift_sym_tot(:,:) = 0.5_dp * (cs_nics_tot(:,:,ir) + TRANSPOSE(cs_nics_tot(:,:,ir)))
                CALL  diamat_all(shift_sym_tot,eig,error=error)
                shift_iso   = (eig(1) + eig(2) + eig(3)) / 3.0_dp
                shift_aniso =  eig(3) - (eig(2) + eig(1)) / 2.0_dp
                !
                WRITE(unit_nics,'(T2,I5,2X,3f15.6)') ir,r_nics(1:3,ir)
                !
                IF(gapw) THEN
                   WRITE(unit_nics,'(T1,A)') " SIGMA from SOFT J"
                   WRITE(unit_nics,'(3(A,f10.4))') '  XX = ',cs_nics(1,1,ir),&
                        &                          '  XY = ',cs_nics(1,2,ir),&
                        &                          '  XZ = ',cs_nics(1,3,ir)
                   WRITE(unit_nics,'(3(A,f10.4))') '  YX = ',cs_nics(2,1,ir),&
                        &                          '  YY = ',cs_nics(2,2,ir),&
                        &                          '  YZ = ',cs_nics(2,3,ir)
                   WRITE(unit_nics,'(3(A,f10.4))') '  ZX = ',cs_nics(3,1,ir),&
                        &                          '  ZY = ',cs_nics(3,2,ir),&
                        &                          '  ZZ = ',cs_nics(3,3,ir)
                   WRITE(unit_nics,'(T1,A)') " SIGMA from LOCAL J"
                   WRITE(unit_nics,'(3(A,f10.4))') '  XX = ',cs_nics_loc(1,1,ir),&
                        &                          '  XY = ',cs_nics_loc(1,2,ir),&
                        &                          '  XZ = ',cs_nics_loc(1,3,ir)
                   WRITE(unit_nics,'(3(A,f10.4))') '  YX = ',cs_nics_loc(2,1,ir),&
                        &                          '  YY = ',cs_nics_loc(2,2,ir),&
                        &                          '  YZ = ',cs_nics_loc(2,3,ir)
                   WRITE(unit_nics,'(3(A,f10.4))') '  ZX = ',cs_nics_loc(3,1,ir),&
                        &                          '  ZY = ',cs_nics_loc(3,2,ir),&
                        &                          '  ZZ = ',cs_nics_loc(3,3,ir)
                ENDIF
                WRITE(unit_nics,'(T1,A)') " SIGMA TOTAL"
                WRITE(unit_nics,'(3(A,f10.4))') '  XX = ',cs_nics_tot(1,1,ir),&
                     &                          '  XY = ',cs_nics_tot(1,2,ir),&
                     &                          '  XZ = ',cs_nics_tot(1,3,ir)
                WRITE(unit_nics,'(3(A,f10.4))') '  YX = ',cs_nics_tot(2,1,ir),&
                     &                          '  YY = ',cs_nics_tot(2,2,ir),&
                     &                          '  YZ = ',cs_nics_tot(2,3,ir)
                WRITE(unit_nics,'(3(A,f10.4))') '  ZX = ',cs_nics_tot(3,1,ir),&
                     &                          '  ZY = ',cs_nics_tot(3,2,ir),&
                     &                          '  ZZ = ',cs_nics_tot(3,3,ir)
                WRITE(unit_nics,'(T1,2(A,f12.4))') '  ISOTROPY = ',shift_iso,&
                     &                           '  ANISOTROPY = ',shift_aniso
             ENDDO
          ENDIF
          CALL cp_print_key_finished_output(unit_nics,logger,nmr_section,&
               &                            "PRINT%SHIELDING_TENSOR",error=error)
       ENDIF
    ENDIF ! print shift
    !
    ! clean up
    DEALLOCATE(cs_tot,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    IF(do_nics) THEN
       DEALLOCATE(cs_nics_tot,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ENDIF
    !
!100 FORMAT(A,1X,I5)
!101 FORMAT(T2,A)
!102 FORMAT(T2,I5,A,2X,A2,2X,3f15.6)
!103 FORMAT(T2,I5,2X,3f15.6)
!104 FORMAT(T1,A)
!105 FORMAT(3(A,f10.4))
!106 FORMAT(T1,2(A,f12.4))
  END SUBROUTINE nmr_shift_print

END MODULE qs_linres_nmr_shift

