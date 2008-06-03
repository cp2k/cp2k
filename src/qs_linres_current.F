!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief given the response wavefunctions obtained by the application
!>      of the (rxp), p, and ((dk-dl)xp) operators,
!>      here the current density vector (jx, jy, jz)
!>      is computed for the 3 directions of the magnetic field (Bx, By, Bz)
!> \par History
!>      created 02-2006 [MI]
!> \author MI
! *****************************************************************************
MODULE qs_linres_current

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE cp_array_i_utils,                ONLY: cp_2d_i_p_type
  USE cp_array_r_utils,                ONLY: cp_2d_r_p_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale,&
                                             cp_fm_scale_and_add,&
                                             cp_fm_trace
  USE cp_fm_pool_types,                ONLY: cp_fm_pool_p_type,&
                                             fm_pools_create_fm_vect,&
                                             fm_pools_give_back_fm_vect
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_submatrix,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_set_all,&
                                             cp_fm_set_submatrix,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_rs_pool_types,                ONLY: cp_rs_pool_p_type,&
                                             cp_rs_pool_type,&
                                             rs_pool_create_rs,&
                                             rs_pool_give_back_rs,&
                                             rs_pools_create_rs_vect
  USE cp_sm_fm_interactions,           ONLY: cp_sm_fm_multiply,&
                                             cp_sm_plus_fm_fm_t
  USE cube_utils,                      ONLY: cube_info_type
  USE gaussian_gridlevels,             ONLY: gaussian_gridlevel,&
                                             gridlevel_info_type
  USE input_section_types,             ONLY: section_get_ivals,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type
  USE kinds,                           ONLY: default_path_length,&
                                             dp,&
                                             int_8
  USE mathconstants,                   ONLY: twopi
  USE memory_utilities,                ONLY: reallocate
  USE orbital_pointers,                ONLY: ncoset
  USE particle_types,                  ONLY: particle_type
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_methods,                      ONLY: pw_copy,&
                                             pw_integrate_function,&
                                             pw_scale,&
                                             pw_zero
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_p_type,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_p_type
  USE qs_collocate_density,            ONLY: collocate_pgf_product_rspace,&
                                             density_rs2pw,&
                                             lgrid_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_linres_atom_current,          ONLY: calculate_jrho_atom,&
                                             calculate_jrho_atom_coeff,&
                                             calculate_jrho_atom_rad
  USE qs_linres_op,                    ONLY: ind_m2,&
                                             set_vecp
  USE qs_linres_types,                 ONLY: current_env_type,&
                                             get_current_env
  USE qs_matrix_pools,                 ONLY: mpools_get,&
                                             qs_matrix_pools_type
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
  USE qs_modify_pab_block,             ONLY: FUNC_ADBmDAB,&
                                             FUNC_ARDBmDARB
  USE qs_neighbor_list_types,          ONLY: &
       first_list, first_node, get_neighbor_list, get_neighbor_list_set, &
       get_neighbor_node, neighbor_list_set_p_type, neighbor_list_type, &
       neighbor_node_type, next
  USE qs_operators_ao,                 ONLY: rRc_xyz_der_ao,&
                                             set_up_op_sm
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE realspace_grid_cube,             ONLY: pw_to_cube
  USE realspace_grid_types,            ONLY: realspace_grid_p_type,&
                                             realspace_grid_type,&
                                             rs_find_node,&
                                             rs_grid_zero
  USE realspace_task_selection,        ONLY: distribute_matrix,&
                                             int2pair,&
                                             pair2int,&
                                             rs_get_my_tasks
  USE sparse_matrix_types,             ONLY: &
       add_block_node, allocate_matrix, allocate_matrix_set, &
       deallocate_matrix, deallocate_matrix_set, get_block_node, &
       real_matrix_p_type, real_matrix_type, replicate_matrix_structure, &
       set_matrix
  USE task_list_methods,               ONLY: estimate_task_cost,&
                                             get_task_bounds
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  ! *** Public subroutines ***
  PUBLIC :: current_build_current,calculate_jrho_resp,current_build_chi

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_linres_current'

CONTAINS

! *****************************************************************************
!> \brief First calculate the density matrixes, for each component of the current
!>       they are 3 because of the r dependent terms
!>       Next it collocates on the grid to have J(r)
!>       In the GAPW case one need to collocate on the PW grid only the soft part
!>       while the rest goes on Lebedev grids
!>       The contributions to the shift and to the susceptibility will be
!>       calulated separately and added only at the end
!>       The calculation of the shift tensor is performed on the position of the atoms
!>       and on other selected points in real space summing up the contributions
!>       from the PW grid current density and the local densities
!>       Spline interpolation is used
!> \param psi 1, p_psi1 : scratch MOS coefficients
!> \note
!>       The susceptibility is needed to compute the G=0 term of the shift
!>       in reciprocal space. \chi_{ij} = \int (r x Jj)_i
!>       (where Jj id the current density generated by the field in direction j)
!>       To calculate the susceptibility on the PW grids it is necessary to apply
!>       the position operator yet another time.
!>       This cannot be done on directly on the full J(r) because it is not localized
!>       Therefore it is done state by state (see linres_nmr_shift)
!> \author MI
! *****************************************************************************
  SUBROUTINE current_build_current(current_env,qs_env,iB,error)
    !
    TYPE(current_env_type)                   :: current_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    INTEGER, INTENT(IN)                      :: iB
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'current_build_current', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: ext, filename
    INTEGER :: handle, homo, iao, idir, idir3, iiB, iiiB, ispin, istat, &
      istate, j, jstate, nao, natom, nmo, nspins, output_unit, unit_nr
    LOGICAL                                  :: failure, gapw, ionode, uni_occ
    REAL(dp)                                 :: alpha, dk(3), &
                                                jrho_tot_G(3,3), &
                                                jrho_tot_R(3,3), maxocc, &
                                                rmu(3), rmu_dk(3), scale_fac
    REAL(dp), DIMENSION(:), POINTER          :: occupation
    REAL(dp), DIMENSION(:, :), POINTER       :: vecbuf_psi1
    REAL(dp), EXTERNAL                       :: DDOT
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_2d_i_p_type), DIMENSION(:), &
      POINTER                                :: center_list
    TYPE(cp_2d_r_p_type), DIMENSION(3)       :: vecbuf
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: p_psi1, psi1
    TYPE(cp_fm_p_type), DIMENSION(:, :), &
      POINTER                                :: psi1_D, psi1_p, psi1_rxp
    TYPE(cp_fm_pool_p_type), DIMENSION(:), &
      POINTER                                :: ao_mo_fm_pools
    TYPE(cp_fm_struct_type), POINTER         :: tmp_fm_struct
    TYPE(cp_fm_type), POINTER                :: fm_work1, mo_coeff, psi0_fk, &
                                                psi_a_iB
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_rs_pool_type), POINTER           :: auxbas_rs_pool
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type)                          :: pw_gspace_work, &
                                                shift_pw_rspace, wf_r
    TYPE(pw_p_type), DIMENSION(:, :), &
      POINTER                                :: shift_pw_gspace
    TYPE(pw_p_type), POINTER                 :: rho_gspace, rho_rspace
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(qs_matrix_pools_type), POINTER      :: mpools
    TYPE(qs_rho_type), POINTER               :: rho_struct
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: density_matrix, &
                                                density_matrix0, &
                                                density_matrix_ii, &
                                                density_matrix_iii
    TYPE(realspace_grid_type), POINTER       :: rs
    TYPE(section_vals_type), POINTER         :: current_section

    failure = .FALSE.
    CALL timeset(routineN,"I"," ",handle)
    !
    NULLIFY(rho_rspace,rho_gspace,logger,current_section,density_matrix0,density_matrix,&
         &  density_matrix_ii,density_matrix_iii,cell,dft_control,mos,rho_struct,particle_set,&
         &  pw_env,auxbas_rs_pool,auxbas_pw_pool,pw_pools,shift_pw_gspace,para_env,center_list,&
         &  mo_coeff,psi0_fk,fm_work1,psi_a_iB,psi1,p_psi1,psi1_p,psi1_D,psi1_rxp)

    logger => cp_error_get_logger(error)
    ionode = logger%para_env%mepos==logger%para_env%source
    output_unit= cp_logger_get_default_unit_nr(logger)
    !
    !
    CALL get_current_env(current_env=current_env,&
         &               center_list=center_list,&
         &               psi1_rxp=psi1_rxp,&
         &               psi1_D=psi1_D,&
         &               psi1_p=psi1_p,&
         &               error=error)
    !
    !
    CALL get_qs_env(qs_env=qs_env,&
         &          rho=rho_struct,&
         &          cell=cell,&
         &          dft_control=dft_control,&
         &          mos=mos,&
         &          mpools=mpools,&
         &          para_env=para_env,&
         &          particle_set=particle_set,&
         &          error=error)

    gapw  = dft_control%qs_control%gapw
    nspins = dft_control%nspins
    natom = SIZE(particle_set,1)
    CPPrecondition(ASSOCIATED(rho_struct),cp_failure_level,routineP,error,failure)
    CPPrecondition(rho_struct%ref_count>0,cp_failure_level,routineP,error,failure)
    !
    !
    CALL mpools_get(mpools, ao_mo_fm_pools=ao_mo_fm_pools,error=error)
    CALL fm_pools_create_fm_vect(ao_mo_fm_pools, psi1, name=routineP//":psi1",error=error)
    CALL fm_pools_create_fm_vect(ao_mo_fm_pools, p_psi1, name=routineP//":p_psi1",error=error)
    !
    !
    ALLOCATE(density_matrix0(nspins),density_matrix(nspins),&
         &   density_matrix_ii(nspins),density_matrix_iii(nspins),&
         &   STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DO ispin = 1,nspins
       !
       !density_matrix0
       NULLIFY(density_matrix0(ispin)%matrix)
       CALL replicate_matrix_structure(current_env%op_rmd_ao(1)%matrix, &
            &                          density_matrix0(ispin)%matrix,"density_matrix0"//&
            &                          "-"//TRIM(ADJUSTL(cp_to_string(ispin))),&
            &                          target_symmetry="none",error=error)
       !
       !density_matrix
       NULLIFY(density_matrix(ispin)%matrix)
       CALL replicate_matrix_structure(current_env%op_rmd_ao(1)%matrix, &
            &                          density_matrix(ispin)%matrix,"density_matrix"//&
            &                          "-"//TRIM(ADJUSTL(cp_to_string(ispin))),&
            &                          target_symmetry="none",error=error)
       !
       !density_matrix_ii
       NULLIFY(density_matrix_ii(ispin)%matrix)
       CALL replicate_matrix_structure(current_env%op_rmd_ao(1)%matrix, &
            &                          density_matrix_ii(ispin)%matrix,"density_matrix_ii"//&
            &                          "-"//TRIM(ADJUSTL(cp_to_string(ispin))),&
            &                          target_symmetry="none",error=error)
       !
       !density_matrix_iii
       NULLIFY(density_matrix_iii(ispin)%matrix)
       CALL replicate_matrix_structure(current_env%op_rmd_ao(1)%matrix, &
            &                          density_matrix_iii(ispin)%matrix,"density_matrix_iii"//&
            &                          "-"//TRIM(ADJUSTL(cp_to_string(ispin))),&
            &                          target_symmetry="none",error=error)
    ENDDO
    !

    current_section => section_vals_get_subs_vals(qs_env%input,"PROPERTIES%LINRES%CURRENT",&
         &                                        error=error)
    !
    IF (.NOT. failure) THEN
       CALL get_qs_env(qs_env=qs_env,pw_env=pw_env,error=error)
       CALL pw_env_get(pw_env, auxbas_rs_pool=auxbas_rs_pool,&
            &          auxbas_pw_pool=auxbas_pw_pool,pw_pools=pw_pools,error=error)
       IF (BTEST(cp_print_key_should_output(logger%iter_info,current_section,&
            &    "PRINT%CURRENT_CUBES",error=error),cp_p_file)) THEN
          CALL rs_pool_create_rs(auxbas_rs_pool,rs,error=error)
          CALL pw_pool_create_pw(auxbas_pw_pool,wf_r%pw,use_data=REALDATA3D,&
               &                  in_space=REALSPACE,error=error)
       ENDIF
       !
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
            &                  use_data=REALDATA3D,in_space=REALSPACE,error=error)
       CALL pw_zero(shift_pw_rspace%pw)
       CALL pw_pool_create_pw(auxbas_pw_pool,pw_gspace_work%pw,&
            &                  use_data=COMPLEXDATA1D,&
            &                  in_space=RECIPROCALSPACE,error=error)
       CALL pw_zero(pw_gspace_work%pw)
       !
       jrho_tot_G = 0.0_dp
       jrho_tot_R = 0.0_dp
       !
       ! built density matrix
       DO ispin = 1,nspins
          CALL get_mo_set(mo_set=mos(ispin)%mo_set,mo_coeff=mo_coeff,&
               &          occupation_numbers=occupation,homo=homo,nao=nao,nmo=nmo,&
               &          uniform_occupation=uni_occ,maxocc=maxocc)
          !
          ! Multiply by the occupation number for the density matrix
          IF ( .NOT. uni_occ ) THEN
             alpha = 1.0_dp
             CALL cp_fm_to_fm(mo_coeff,psi0_fk,error=error)
             CALL cp_fm_column_scale(psi0_fk,occupation(1:homo))
          ELSE
             alpha = maxocc
             psi0_fk => mo_coeff
          ENDIF
          !
          ! Build the first density matrix
          CALL set_matrix(density_matrix0(ispin)%matrix,0.0_dp)
          CALL cp_sm_plus_fm_fm_t(sparse_matrix=density_matrix0(ispin)%matrix ,&
               &                  matrix_v=psi0_fk,matrix_g=psi0_fk,&
               &                  ncol=homo,alpha=alpha,error=error)
       ENDDO
       !
       ! Lets go!
       DO ispin = 1,nspins
          CALL cp_fm_set_all(psi1(ispin)%matrix,0.0_dp,error=error)
          CALL cp_fm_set_all(p_psi1(ispin)%matrix,0.0_dp,error=error)
       ENDDO

       CALL set_vecp(iB,iiB,iiiB)
       DO ispin = 1,nspins
          CALL get_mo_set(mo_set=mos(ispin)%mo_set,mo_coeff=mo_coeff,&
               &          occupation_numbers=occupation,homo=homo,nao=nao,nmo=nmo,&
               &          uniform_occupation=uni_occ,maxocc=maxocc)
          NULLIFY(fm_work1)
          !
          ! create a new matrix
          NULLIFY(tmp_fm_struct)
          CALL cp_fm_struct_create(tmp_fm_struct,nrow_global=nao,&
               &                   ncol_global=nmo,para_env=para_env,&
               &                   context=mo_coeff%matrix_struct%context,&
               &                   error=error)
          IF(ASSOCIATED(fm_work1))THEN
             CALL cp_fm_release(fm_work1,error=error)
          ENDIF
          CALL cp_fm_create(fm_work1,tmp_fm_struct,error=error)
          CALL cp_fm_set_all(fm_work1,0.0_dp,error=error)
          IF ( .NOT. uni_occ ) THEN
             CALL cp_fm_create(psi0_fk,tmp_fm_struct,error=error)
             CALL cp_fm_set_all(psi0_fk,0.0_dp,error=error)
          ELSE
             NULLIFY(psi0_fk)
          END IF
          CALL cp_fm_struct_release(tmp_fm_struct,error=error)
          !
          ! Allocate buffer vectors
          DO idir = 1,3
             ALLOCATE(vecbuf(idir)%array(1,nao),STAT=istat)
             CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          END DO
          ALLOCATE(vecbuf_psi1(1,nao),STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          !
          ! Construct the 3 density matrices for the field in direction iB
          !
          ! First the full matrix psi_a_iB
          !   psi1_rxp - psi1_D - (R_\nu-dk)_ii psi1_piiiB + (R_\nu-dk)_iii psi1_piiB
          psi_a_iB => psi1(ispin)%matrix
          CALL cp_fm_set_all(psi_a_iB,0.0_dp,error=error)
          !
          ! contributions from the response psi1_p_ii and psi1_p_iii
          DO istate = 1,current_env%nbr_center(ispin)
             dk(1:3) = current_env%centers_set(ispin)%array(1:3,istate)
             !
             ! Copy the vector in the full matrix psi1
             DO j = center_list(ispin)%array(1,istate),center_list(ispin)%array(1,istate+1)-1
                jstate = center_list(ispin)%array(2,j)
                CALL cp_fm_get_submatrix(psi1_p(ispin,iiB)%matrix,&
                     &                   vecbuf(iiB)%array,1,jstate,nao,1,&
                     &                   transpose=.TRUE.,error=error)
                CALL cp_fm_get_submatrix(psi1_p(ispin,iiiB)%matrix,&
                     &                   vecbuf(iiiB)%array,1,jstate,nao,1,&
                     &                   transpose=.TRUE.,error=error)
                DO iao = 1,nao
                   Rmu(1:3) = current_env%basisfun_center(1:3,iao)
                   Rmu_dk = pbc(dk,Rmu,cell)
                   vecbuf_psi1(1,iao) = Rmu_dk(iiB)*vecbuf(iiiB)%array(1,iao)-&
                        &               Rmu_dk(iiiB)*vecbuf(iiB)%array(1,iao)
                ENDDO  ! iao
                CALL cp_fm_set_submatrix(psi1(ispin)%matrix,vecbuf_psi1,&
                     &                   1,jstate,nao,1,transpose=.TRUE.,&
                     &                   error=error)
             ENDDO
          ENDDO ! istate
          !
          ! contribution from the response psi1_rxp
          CALL cp_fm_scale_and_add(-1.0_dp,psi_a_iB,1.0_dp,&
               &                   psi1_rxp(ispin,iB)%matrix,&
               &                   error=error)
          IF(current_env%full) THEN
             ! contribution from the response psi1_D
             CALL  cp_fm_scale_and_add(1.0_dp,psi_a_iB,-1.0_dp,&
                  &                    psi1_D(ispin,iB)%matrix,&
                  &                    error=error)
          ENDIF
          !
          ! Multiply by the occupation number for the density matrix
          IF ( .NOT. uni_occ ) THEN
             alpha = 1.0_dp
             CALL cp_fm_to_fm(mo_coeff,psi0_fk,error=error)
             CALL cp_fm_column_scale(psi0_fk,occupation(1:homo))
          ELSE
             alpha = maxocc
             psi0_fk => mo_coeff
          ENDIF
          !
          ! Build the first density matrix
          CALL set_matrix(density_matrix(ispin)%matrix ,0.0_dp)
          CALL cp_sm_plus_fm_fm_t(sparse_matrix=density_matrix(ispin)%matrix ,&
               &                  matrix_v=psi0_fk,matrix_g=psi_a_iB,&
               &                  ncol=homo,alpha=alpha,error=error)
          !
          ! Build the second density matrix
          CALL set_matrix(density_matrix_iii(ispin)%matrix,0.0_dp)
          CALL cp_sm_plus_fm_fm_t(sparse_matrix=density_matrix_iii(ispin)%matrix ,&
               &                  matrix_v=psi0_fk,matrix_g=psi1_p(ispin,iiiB)%matrix,&
               &                  ncol=homo,alpha=alpha,error=error)
          !
          ! Build the third density matrix
          CALL set_matrix(density_matrix_ii(ispin)%matrix ,0.0_dp)
          CALL cp_sm_plus_fm_fm_t(sparse_matrix=density_matrix_ii(ispin)%matrix ,&
               &                  matrix_v=psi0_fk,matrix_g=psi1_p(ispin,iiB)%matrix,&
               &                  ncol=homo,alpha=alpha,error=error)
          DO idir3=1,3
             ! set to zero for the calculation of the shift
             CALL pw_zero(shift_pw_gspace(idir3,ispin)%pw)
          ENDDO
          DO idir = 1,3
             !
             ! Calculate the current density on the pw grid (only soft if GAPW)
             ! idir is the cartesian component of the response current density
             ! generated by the magnetic field pointing in cartesian direction iB
             ! Use the qs_rho_type already  used for rho during the scf
             rho_rspace => current_env%jrho1_set(idir)%rho%rho_r(ispin)
             rho_gspace => current_env%jrho1_set(idir)%rho%rho_g(ispin)
             CALL pw_zero(rho_rspace%pw)
             CALL pw_zero(rho_gspace%pw)
             CALL calculate_jrho_resp(density_matrix(ispin)%matrix, &
                  &                   density_matrix_ii(ispin)%matrix,&
                  &                   density_matrix_iii(ispin)%matrix, &
                  &                   iB,idir,rho_rspace,rho_gspace,qs_env,&
                  &                   gapw,error=error)

             scale_fac = cell%deth / twopi
             CALL pw_scale(rho_rspace%pw,scale_fac)
             CALL pw_scale(rho_gspace%pw,scale_fac)
             
             jrho_tot_G(idir,iB) = pw_integrate_function(rho_gspace%pw,isign=-1)
             jrho_tot_R(idir,iB) = pw_integrate_function(rho_rspace%pw,isign=-1)

             IF(output_unit>0) THEN
                WRITE(output_unit,'(T2,2(A,E24.16))') 'Integrated j_'&
                     &//ACHAR(idir+119)//ACHAR(iB+119)//'(r): G-space=',&
                     &jrho_tot_G(idir,iB),' R-space=',jrho_tot_R(idir,iB)
             ENDIF
             !
          ENDDO ! idir
          !
          ! Deallocate buffer vectors
          DO idir = 1,3
             DEALLOCATE(vecbuf(idir)%array,STAT=istat)
             CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          ENDDO
          DEALLOCATE(vecbuf_psi1,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

          CALL cp_fm_release(fm_work1,error=error)
          IF(.NOT.uni_occ) THEN
             CALL cp_fm_release(psi0_fk,error=error)
          ELSE
             NULLIFY(psi0_fk)
          ENDIF
       ENDDO !  ispin
       
       IF(gapw) THEN
          DO idir = 1,3
             !
             ! compute the atomic response current densities on the spherical grids
             ! First the sparse matrices are multiplied by the expansion coefficients
             ! this is the equivalent of the CPC for the charge density
             CALL calculate_jrho_atom_coeff(qs_env,current_env,density_matrix0,density_matrix,&
                  &                           density_matrix_ii,density_matrix_iii,&
                  &                           iB,idir,error=error)
             !
             ! Then the radial parts are computed on the local radial grid, atom by atom
             ! 8 functions are computed for each atom, per grid point
             ! and per LM angular momentum. The multiplication by the Clebsh-Gordon
             ! coefficients or they correspondent for the derivatives, is also done here
             CALL calculate_jrho_atom_rad(qs_env,current_env,iB,idir,error=error)
             !
             ! The current on the atomic grids
             CALL calculate_jrho_atom(current_env,qs_env,iB,idir,error=error)
          ENDDO ! idir
       ENDIF
       !
       ! Cube files 
       IF (BTEST(cp_print_key_should_output(logger%iter_info,current_section,&
            &    "PRINT%CURRENT_CUBES",error=error),cp_p_file)) THEN
          IF(output_unit>0) THEN
             WRITE(output_unit,*) 'WARNING WARNING WARNING WARNING WARNING WARNING WARNING'
             WRITE(output_unit,*) 'Are you sure you want to punch current density cube files?'
             WRITE(output_unit,*) 'WARNING WARNING WARNING WARNING WARNING WARNING WARNING'
          ENDIF
          DO idir = 1,3
             CALL pw_zero(wf_r%pw)
             !!WARNING
             ! if nspins=2 here ist should be add instead of copy, or separated cubefiles
             ! for the different spin should be written
             !!WARNING
             DO ispin =1 ,nspins
                CALL pw_copy(current_env%jrho1_set(idir)%rho%rho_r(ispin)%pw,wf_r%pw)
             ENDDO
             IF(gapw) THEN
                ! Add the local hard and soft contributions
                ! This can be done atom by atom by a spline extrapolation of the  values
                ! on the spherical grid to the grid points.
                CALL stop_program(routineP,"GAPW needs to be finalized") 
             ENDIF
             filename="jresp"
             WRITE(ext,'(a2,I1,a2,I1,a5)')  "iB",iB,"_d",idir,".cube"
             WRITE(ext,'(a2,a1,a2,a1,a5)')  "iB",ACHAR(iB+119),"_d",ACHAR(idir+119),".cube"
             unit_nr=cp_print_key_unit_nr(logger,current_section,"PRINT%CURRENT_CUBES",&
                  &                       extension=TRIM(ext),middle_name=TRIM(filename),&
                  &                       log_filename=.FALSE.,file_position="REWIND",&
                  &                       error=error)
             
             CALL pw_to_cube(wf_r%pw,unit_nr,ionode,"RESPONSE CURRENT DENSITY ",&
                  &             stride=section_get_ivals(current_section,"PRINT%CURRENT_CUBES%STRIDE",&
                  &             error=error),error=error)
             CALL cp_print_key_finished_output(unit_nr,logger,current_section,&
                  &                            "PRINT%CURRENT_CUBES",error=error)
          ENDDO
          
       ENDIF  ! current cube
       !
       ! Integrated current response checksum
       IF(output_unit>0) THEN
          WRITE(output_unit,'(T2,A,E24.16)') 'CheckSum R-integrated j=',&
               &            SQRT(DDOT(9,jrho_tot_R(1,1),1,jrho_tot_R(1,1),1))
       ENDIF
       !
       IF (BTEST(cp_print_key_should_output(logger%iter_info,current_section,&
            &    "PRINT%CURRENT_CUBES",error=error),cp_p_file)) THEN
          CALL pw_pool_give_back_pw(auxbas_pw_pool,wf_r%pw,error=error)
          CALL rs_pool_give_back_rs(auxbas_rs_pool,rs, error=error)
       ENDIF
       !
       ! Dellocate grids for the calculation of jrho and the shift
       CALL fm_pools_give_back_fm_vect(ao_mo_fm_pools, psi1,error=error)
       CALL fm_pools_give_back_fm_vect(ao_mo_fm_pools, p_psi1,error=error)
       !
       !
       CALL pw_pool_give_back_pw(auxbas_pw_pool,pw_gspace_work%pw,error=error)
       DO ispin = 1,nspins
          DO idir = 1,3
             CALL pw_pool_give_back_pw(auxbas_pw_pool,shift_pw_gspace(idir,ispin)%pw,error=error)
          ENDDO
       ENDDO
       DEALLOCATE(shift_pw_gspace,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       CALL pw_pool_give_back_pw(auxbas_pw_pool,shift_pw_rspace%pw,error=error)
       !
       !
       DO ispin = 1,nspins
          CALL deallocate_matrix(density_matrix0(ispin)%matrix,error=error)
          CALL deallocate_matrix(density_matrix(ispin)%matrix,error=error)
          CALL deallocate_matrix(density_matrix_ii(ispin)%matrix,error=error)
          CALL deallocate_matrix(density_matrix_iii(ispin)%matrix,error=error)
       ENDDO
       DEALLOCATE(density_matrix0,density_matrix,density_matrix_ii,density_matrix_iii,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ENDIF !  failure
    !
    ! Finalize
    CALL timestop(0.0_dp,handle)
    !
  END SUBROUTINE current_build_current

! *****************************************************************************
!> \brief Calculation of the idir component of the response current density
!>       in the presence of a constant magnetic field in direction iB
!>       the current density is collocated on the pw grid in real space
!> \note
!>       The collocate is done in three parts, one for each density matrix
!>       In all cases the density matrices and therefore the collocation
!>       are not symmetric, that means that all the pairs (ab and ba) have
!>       to be considered separately
!> 
!>       mat_jp_{\mu\nu} is multiplied by
!>           f_{\mu\nu} = \phi_{\mu} (d\phi_{\nu}/dr)_{idir} -
!>                        (d\phi_{\mu}/dr)_{idir} \phi_{\nu}
!> 
!>       mat_jp_rii_{\mu\nu} is multiplied by
!>           f_{\mu\nu} = \phi_{\mu} (r - R_{\nu})_{iiiB} (d\phi_{\nu}/dr)_{idir} -
!>                        (d\phi_{\mu}/dr)_{idir} (r - R_{\nu})_{iiiB} \phi_{\nu} +
!>                         \phi_{\mu} \phi_{\nu}  (last term only if iiiB=idir)
!> 
!>       mat_jp_riii_{\mu\nu} is multiplied by
!>                             (be careful: change in sign with respect to previous)
!>           f_{\mu\nu} = -\phi_{\mu} (r - R_{\nu})_{iiB} (d\phi_{\nu}/dr)_{idir} +
!>                        (d\phi_{\mu}/dr)_{idir} (r - R_{\nu})_{iiB} \phi_{\nu} -
!>                         \phi_{\mu} \phi_{\nu}  (last term only if iiB=idir)
!> 
!>       All the terms sum up to the same grid
! *****************************************************************************
  SUBROUTINE calculate_jrho_resp(mat_jp,mat_jp_rii,mat_jp_riii,iB,idir,&
       rho_rs, rho_gs, qs_env, soft_valid, error)

    TYPE(real_matrix_type), POINTER          :: mat_jp, mat_jp_rii, &
                                                mat_jp_riii
    INTEGER, INTENT(IN)                      :: iB, idir
    TYPE(pw_p_type), INTENT(INOUT)           :: rho_rs, rho_gs
    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(IN), OPTIONAL            :: soft_valid
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_jrho_resp', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: add_tasks = 1000, &
                                                max_tasks = 2000
    REAL(kind=dp), PARAMETER                 :: mult_tasks = 2.0_dp

    INTEGER :: ab, added_tasks, bcol, brow, cost, curr_tasks, dest, dir, &
      handle, i, iatom, iatom_old, igrid_level, iiB, iiiB, ikind, ikind_old, &
      ilist, inode, ipgf, iset, iset_old, istat, itask, ithread, j, jatom, &
      jatom_old, jkind, jkind_old, jpgf, jset, jset_old, lb_cube(3), maxco, &
      maxpgf, maxset, maxsgf, maxsgf_set, na1, na2, natom, nb1, nb2, ncoa, &
      ncob, nkind, nlist, nnode, nseta, nsetb, ntasks, nthread, sgfa, sgfb, &
      stat, tp(3), ub_cube(3)
    INTEGER(kind=int_8), DIMENSION(:), &
      POINTER                                :: atom_pair_recv, atom_pair_send
    INTEGER(kind=int_8), DIMENSION(:, :), &
      POINTER                                :: tasks
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, lb, lb_max, &
                                                lb_min, npgfa, npgfb, nsgfa, &
                                                nsgfb
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb
    LOGICAL :: atom_pair_changed, distributed_rs_grids, failure, &
      map_consistent, my_minimum_image, my_soft
    REAL(KIND=dp)                            :: dab, eps_rho_rspace, &
                                                kind_radius_a, kind_radius_b, &
                                                Lxo2, Lyo2, Lzo2, rab2, &
                                                scale, zetp
    REAL(KIND=dp), DIMENSION(3)              :: ra, rab, rb, rp
    REAL(KIND=dp), DIMENSION(:), POINTER     :: set_radius_a, set_radius_b
    REAL(KIND=dp), DIMENSION(:, :), POINTER :: dist_ab, jp_block, &
      jp_block_rii, jp_block_riii, jpab, jpab_ii, jpab_iii, jpblock, &
      jpblock_rii, jpblock_riii, rpgfa, rpgfb, sphi_a, sphi_b, work, zeta, &
      zetb
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: jpabt, jpabt_ii, jpabt_iii, &
                                                workt
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_rs_pool_p_type), DIMENSION(:), &
      POINTER                                :: rs_pools
    TYPE(cube_info_type), DIMENSION(:), &
      POINTER                                :: cube_info
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(gridlevel_info_type), POINTER       :: gridlevel_info
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(lgrid_type)                         :: lgrid
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb
    TYPE(neighbor_list_type), POINTER        :: sab_orb_neighbor_list
    TYPE(neighbor_node_type), POINTER        :: sab_orb_neighbor_node
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(real_matrix_type), POINTER          :: deltajp, deltajp_rii, &
                                                deltajp_riii
    TYPE(realspace_grid_p_type), &
      DIMENSION(:), POINTER                  :: rs_rho
    TYPE(section_vals_type), POINTER         :: input, interp_section

!   ---------------------------------------------------------------------------

    failure=.FALSE.
    NULLIFY(atomic_kind,cell,dft_control,orb_basis_set,sab_orb_neighbor_list,&
         sab_orb_neighbor_node,atomic_kind_set,sab_orb,particle_set,&
         rs_rho,pw_env,rs_pools,para_env,dist_ab,&
         set_radius_a,set_radius_b,la_max,la_min,&
         lb_max,lb_min,npgfa,npgfb,nsgfa,nsgfb,&
         rpgfa,rpgfb,sphi_a,sphi_b,zeta,zetb,first_sgfa,first_sgfb,&
         dist_ab,tasks,workt)
    NULLIFY(deltajp,deltajp_rii,deltajp_riii)
    NULLIFY(jp_block,jp_block_rii,jp_block_riii)
    NULLIFY(jpblock,jpblock_rii,jpblock_riii)
    NULLIFY(jpabt,jpabt_ii,jpabt_iii)
    NULLIFY(lgrid%r)

    !    debug_count=debug_count+1

    CALL timeset(routineN,"I"," ",handle)

    CALL get_qs_env(qs_env=qs_env,&
         atomic_kind_set=atomic_kind_set,&
         cell=cell,&
         dft_control=dft_control,&
         particle_set=particle_set,&
         sab_all=sab_orb,&
         para_env=para_env,&
         input=input,&
         pw_env=pw_env,error=error)

    ! Component of appearing in the vector product rxp, iiB and iiiB
    CALL set_vecp(iB,iiB,iiiB)

    ! *** assign from pw_env
    gridlevel_info=>pw_env%gridlevel_info
    cube_info=>pw_env%cube_info

    interp_section => section_vals_get_subs_vals(input,"DFT%MGRID%INTERPOLATOR",&
         error=error)
    !    CALL section_vals_val_get(interp_section,"KIND",i_val=interp_kind,error=error)

    !   Check that the neighbor list with all the pairs is associated
    CPPrecondition(ASSOCIATED(sab_orb),cp_failure_level,routineP,error,failure)
    ! *** set up the pw multi-grids
    CPPrecondition(ASSOCIATED(pw_env),cp_failure_level,routineP,error,failure)
    CALL pw_env_get(pw_env, rs_pools=rs_pools, error=error)

    ! *** set up the rs multi-grids
    distributed_rs_grids=.FALSE.
    CALL rs_pools_create_rs_vect(rs_pools, rs_rho, error=error)
    DO igrid_level=1,gridlevel_info%ngrid_levels
       CALL rs_grid_zero(rs_rho(igrid_level)%rs_grid)
       IF ( .NOT. ALL (rs_rho(igrid_level)%rs_grid%perd == 1 ) ) THEN
          distributed_rs_grids=.TRUE.
       ENDIF
    END DO

    eps_rho_rspace = dft_control%qs_control%eps_rho_rspace
    map_consistent = dft_control%qs_control%map_consistent
    nthread = 1

    !   *** Allocate work storage ***

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
         maxco=maxco,&
         maxsgf=maxsgf,&
         maxsgf_set=maxsgf_set)
    my_minimum_image = .TRUE.
    !    IF(PRESENT(minimum_image)) THEN
    !       my_minimum_image=minimum_image
    Lxo2 = SQRT ( SUM ( cell % hmat ( :, 1 ) ** 2 ) )/2.0_dp
    Lyo2 = SQRT ( SUM ( cell % hmat ( :, 2 ) ** 2 ) )/2.0_dp
    Lzo2 = SQRT ( SUM ( cell % hmat ( :, 3 ) ** 2 ) )/2.0_dp
    !    END IF

    my_soft=.FALSE.
    IF (PRESENT(soft_valid)) my_soft = soft_valid

    ntasks = 0
    IF ( nthread > 1 ) THEN
       DO igrid_level = 1,gridlevel_info%ngrid_levels
          ntasks = MAX(ntasks,rs_rho(igrid_level)%rs_grid%ngpts_local)
       END DO
       ntasks = ntasks*nthread
       CALL reallocate(lgrid%r,1,ntasks)
    END IF

    nkind = SIZE(atomic_kind_set)

    CALL reallocate(jpabt,1,maxco,1,maxco,0,nthread-1)
    CALL reallocate(jpabt_ii,1,maxco,1,maxco,0,nthread-1)
    CALL reallocate(jpabt_iii,1,maxco,1,maxco,0,nthread-1)
    CALL reallocate(workt,1,maxco,1,maxsgf_set,0,nthread-1)
    CALL reallocate(tasks,1,5,1,max_tasks)
    CALL reallocate(dist_ab,1,3,1,max_tasks)

    tasks = 0
    ntasks = 0
    curr_tasks = SIZE(tasks,2)

    !   get maximum numbers
    natom = SIZE( particle_set )
    maxset=0
    maxpgf=0

    DO ikind=1,nkind
       atomic_kind => atomic_kind_set(ikind)

       CALL get_atomic_kind(atomic_kind=atomic_kind,&
            orb_basis_set=orb_basis_set)

       IF (.NOT.ASSOCIATED(orb_basis_set)) CYCLE

       CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
            npgf=npgfa, nset=nseta )

       maxset=MAX(nseta,maxset)
       maxpgf=MAX(MAXVAL(npgfa),maxpgf)
    END DO

    !   *** Initialize working density matrix ***

    ! distributed rs grids require a matrix that will be changed (rs_get_my_tasks)
    ! whereas this is not the case for replicated grids
    IF (distributed_rs_grids) THEN
       CALL allocate_matrix(matrix=deltajp,&
            nrow=mat_jp%nrow,&
            ncol=mat_jp%ncol,&
            nblock_row=mat_jp%nblock_row,&
            nblock_col=mat_jp%nblock_col,&
            first_row=mat_jp%first_row(:),&
            last_row=mat_jp%last_row(:),&
            first_col=mat_jp%first_col(:),&
            last_col=mat_jp%last_col(:),&
            matrix_name="DeltaP",&
            sparsity_id=-1, &   ! basically unknown sparsity in parallel
            matrix_symmetry=mat_jp%symmetry,error=error)
       CALL allocate_matrix(matrix=deltajp_rii,&
            nrow=mat_jp%nrow,&
            ncol=mat_jp%ncol,&
            nblock_row=mat_jp%nblock_row,&
            nblock_col=mat_jp%nblock_col,&
            first_row=mat_jp%first_row(:),&
            last_row=mat_jp%last_row(:),&
            first_col=mat_jp%first_col(:),&
            last_col=mat_jp%last_col(:),&
            matrix_name="DeltaP",&
            sparsity_id=-1, &   ! basically unknown sparsity in parallel
            matrix_symmetry=mat_jp%symmetry,error=error)
       CALL allocate_matrix(matrix=deltajp_riii,&
            nrow=mat_jp%nrow,&
            ncol=mat_jp%ncol,&
            nblock_row=mat_jp%nblock_row,&
            nblock_col=mat_jp%nblock_col,&
            first_row=mat_jp%first_row(:),&
            last_row=mat_jp%last_row(:),&
            first_col=mat_jp%first_col(:),&
            last_col=mat_jp%last_col(:),&
            matrix_name="DeltaP",&
            sparsity_id=-1, &   ! basically unknown sparsity in parallel
            matrix_symmetry=mat_jp%symmetry,error=error)
    ELSE
       deltajp=>mat_jp
       deltajp_rii=>mat_jp_rii
       deltajp_riii=>mat_jp_riii
    ENDIF

    loop_ikind: DO ikind=1,nkind

       atomic_kind => atomic_kind_set(ikind)

       CALL get_atomic_kind(atomic_kind=atomic_kind,&
            softb = my_soft, &
            orb_basis_set=orb_basis_set)

       IF (.NOT.ASSOCIATED(orb_basis_set)) CYCLE

       CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
            first_sgf=first_sgfa,&
            kind_radius=kind_radius_a,&
            lmax=la_max,&
            lmin=la_min,&
            npgf=npgfa,&
            nset=nseta,&
            nsgf_set=nsgfa,&
            pgf_radius=rpgfa,&
            set_radius=set_radius_a,&
            sphi=sphi_a,&
            zet=zeta)

       loop_jkind: DO jkind=1,nkind

          atomic_kind => atomic_kind_set(jkind)

          CALL get_atomic_kind(atomic_kind=atomic_kind,&
               softb = my_soft, &
               orb_basis_set=orb_basis_set)

          IF (.NOT.ASSOCIATED(orb_basis_set)) CYCLE

          CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
               first_sgf=first_sgfb,&
               kind_radius=kind_radius_b,&
               lmax=lb_max,&
               lmin=lb_min,&
               npgf=npgfb,&
               nset=nsetb,&
               nsgf_set=nsgfb,&
               pgf_radius=rpgfb,&
               set_radius=set_radius_b,&
               sphi=sphi_b,&
               zet=zetb)

          ab = ikind + nkind*(jkind - 1)

          IF (ASSOCIATED(sab_orb(ab)%neighbor_list_set)) THEN

             CALL get_neighbor_list_set(neighbor_list_set=&
                  sab_orb(ab)%neighbor_list_set,&
                  nlist=nlist)
             sab_orb_neighbor_list => first_list(sab_orb(ab)%neighbor_list_set)
          ELSE
             nlist=0
          END IF

          loop_ilist:DO ilist = 1, nlist

             CALL get_neighbor_list(neighbor_list=sab_orb_neighbor_list,&
                  atom=iatom,nnode=nnode)

             ra(:) = pbc(particle_set(iatom)%r,cell)

             sab_orb_neighbor_node => first_node(sab_orb_neighbor_list)

             loop_inode: DO inode = 1, nnode

                CALL get_neighbor_node(neighbor_node=sab_orb_neighbor_node,&
                     neighbor=jatom,&
                     r=rab(:))

                IF(my_minimum_image) THEN
                   IF(ABS(rab(1)) > Lxo2 .OR. ABS(rab(2)) > Lyo2 .OR. ABS(rab(3)) > Lzo2) THEN
                      sab_orb_neighbor_node => next(sab_orb_neighbor_node)
                      CYCLE
                   END IF
                END IF

                brow = iatom
                bcol = jatom

                CALL get_block_node(matrix=mat_jp,&
                     block_row=brow,&
                     block_col=bcol,&
                     BLOCK=jp_block)
                CALL get_block_node(matrix=mat_jp_rii,&
                     block_row=brow,&
                     block_col=bcol,&
                     BLOCK=jp_block_rii)
                CALL get_block_node(matrix=mat_jp_riii,&
                     block_row=brow,&
                     block_col=bcol,&
                     BLOCK=jp_block_riii)

                IF (.NOT.ASSOCIATED(jp_block)) THEN
                   sab_orb_neighbor_node => next(sab_orb_neighbor_node)
                   CYCLE
                END IF

                IF (distributed_rs_grids) THEN
!                   NULLIFY (jpblock,jpblock_rii,jp_block_riii )
                   NULLIFY (jpblock,jpblock_rii,jpblock_riii )
                   CALL add_block_node ( deltajp, brow, bcol, jpblock ,error=error)
                   jpblock = jp_block
                   CALL add_block_node ( deltajp_rii, brow, bcol, jpblock_rii ,error=error)
                   jpblock_rii = jp_block_rii
                   CALL add_block_node ( deltajp_riii, brow, bcol, jpblock_riii ,error=error)
                   jpblock_riii = jp_block_riii
                ELSE
                   jpblock => jp_block
                   jpblock_rii => jp_block_rii
                   jpblock_riii => jp_block_riii
                ENDIF

                IF (.NOT. map_consistent) THEN
                   IF ( ALL ( 100.0_dp*ABS(jpblock)      < eps_rho_rspace) .AND. & 
                        ALL ( 100.0_dp*ABS(jpblock_rii)  < eps_rho_rspace) .AND. & 
                        ALL ( 100.0_dp*ABS(jpblock_riii) < eps_rho_rspace) ) THEN 
                      sab_orb_neighbor_node => next(sab_orb_neighbor_node)
                      CYCLE
                   END IF
                END IF

                rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
                dab = SQRT(rab2)

                loop_iset : DO iset=1,nseta
                   IF (set_radius_a(iset) + kind_radius_b < dab) CYCLE

                   loop_jset: DO jset = 1,nsetb
                      IF (set_radius_a(iset) + set_radius_b(jset) < dab) CYCLE

                      loop_ipgf: DO ipgf=1,npgfa(iset)
                         IF (rpgfa(ipgf,iset) + set_radius_b(jset) < dab) CYCLE

                         loop_jpgf: DO jpgf=1,npgfb(jset)
                            IF (rpgfa(ipgf,iset) + rpgfb(jpgf,jset) < dab) CYCLE

                            zetp = zeta(ipgf,iset) + zetb(jpgf,jset)

                            IF (dab.lt.0.1E0_dp .AND. dft_control%qs_control%map_paa) THEN
                               igrid_level = 1
                            ELSE
                               igrid_level = gaussian_gridlevel(gridlevel_info,zetp)
                            ENDIF

                            ntasks = ntasks + 1
                            IF ( ntasks > curr_tasks ) THEN
                               curr_tasks = curr_tasks*mult_tasks
                               CALL reallocate(tasks,1,5,1,curr_tasks)
                               CALL reallocate(dist_ab,1,3,1,curr_tasks)
                            END IF

                            rp(:) = ra(:) + zetb(jpgf,jset)/zetp*rab(:)
                            rp(:) = pbc(rp,cell)
                            DO dir = 1,3
                               tp(dir) = FLOOR(DOT_PRODUCT(cell%h_inv(dir,:),rp)*rs_rho(igrid_level)%rs_grid%npts(dir))
                               tp(dir) = MODULO ( tp(dir), rs_rho(igrid_level)%rs_grid%npts(dir) )
                               tp(dir) = tp(dir) + rs_rho(igrid_level)%rs_grid%lb(dir)
                            END DO

                            IF(rs_rho(igrid_level)%rs_grid%distributed) THEN
                               ! find the bounds required for proper collocation
                               CALL get_task_bounds(rs_rho(igrid_level)%rs_grid,cube_info(igrid_level),&
                                    la_max(iset),zeta(ipgf,iset),la_min(iset),&
                                    lb_max(jset),zetb(jpgf,jset),lb_min(jset),&
                                    ra,rab,rab2,dft_control%qs_control%eps_rho_rspace,&
                                    lb_cube,ub_cube)
                               ! finds the node(s) that need to process this task
                               ! on exit tasks(4,:) is 1 for distributed tasks and 2 for generalised tasks
                               CALL rs_find_node(rs_rho(igrid_level)%rs_grid,tp,dest,&
                                    ntasks=ntasks,tasks=tasks,lb_cube=lb_cube,ub_cube=ub_cube,added_tasks=added_tasks)
                               IF ( ntasks > curr_tasks ) THEN
                                  curr_tasks = (curr_tasks+add_tasks)*mult_tasks
                                  CALL reallocate(dist_ab,1,3,1,curr_tasks)
                               END IF
                            ELSE
                               dest =  rs_rho(igrid_level)%rs_grid%my_pos
                               tasks (1,ntasks) = dest
                               tasks(4,ntasks) = 0
                               added_tasks = 1
                            ENDIF

                            ! cost_factor for this task
                            CALL estimate_task_cost(cost,rs_rho(igrid_level)%rs_grid,cube_info(igrid_level),&
                                 la_max(iset),zeta(ipgf,iset),la_min(iset),&
                                 lb_max(jset),zetb(jpgf,jset),lb_min(jset),&
                                 ra,rab,rab2,dft_control%qs_control%eps_rho_rspace)

                            DO j= 1, added_tasks
                               tasks (2,ntasks-added_tasks+j) = rs_rho(igrid_level)% rs_grid % my_pos
                               tasks (5,ntasks-added_tasks+j) = cost
                               !encode the atom pairs and basis info as a single long integer
                               CALL pair2int(tasks(3,ntasks-added_tasks+j),igrid_level,&
                                    iatom,jatom,iset,jset,ipgf,jpgf,natom,maxset,maxpgf)
                               dist_ab (1,ntasks-added_tasks+j) = rab(1)
                               dist_ab (2,ntasks-added_tasks+j) = rab(2)
                               dist_ab (3,ntasks-added_tasks+j) = rab(3)
                            ENDDO
                           
                         END DO loop_jpgf

                      END DO loop_ipgf

                   END DO loop_jset

                END DO loop_iset

                sab_orb_neighbor_node => next(sab_orb_neighbor_node)

             END DO loop_inode

             sab_orb_neighbor_list => next(sab_orb_neighbor_list)

          END DO loop_ilist

       END DO loop_jkind

    END DO loop_ikind

    ! sorts / redistributes the task list
    CALL rs_get_my_tasks ( rs_rho(1)%rs_grid, distributed_rs_grids, ntasks, natom,&
            maxset, maxpgf, tasks, dist_ab, atom_pair_send, atom_pair_recv,&
            .FALSE. , error)

    IF (distributed_rs_grids) THEN
        CALL distribute_matrix (rs_rho(1)%rs_grid, deltajp, atom_pair_send, atom_pair_recv, natom, scatter=.TRUE., error=error)
        CALL distribute_matrix (rs_rho(1)%rs_grid, deltajp_rii, atom_pair_send, atom_pair_recv, natom, scatter=.TRUE., error=error)
        CALL distribute_matrix (rs_rho(1)%rs_grid, deltajp_riii, atom_pair_send, atom_pair_recv, natom, scatter=.TRUE., error=error)
    ENDIF

    ithread = 0
    jpab => jpabt(:,:,ithread)
    jpab_ii => jpabt_ii(:,:,ithread)
    jpab_iii => jpabt_iii(:,:,ithread)
    work => workt(:,:,ithread)

    iatom_old = -1 ; jatom_old = -1 ; iset_old = -1 ; jset_old = -1 
    ikind_old = -1 ; jkind_old = -1

    loop_tasks: DO itask = 1,ntasks

       CALL int2pair(tasks(3,itask),igrid_level,iatom,jatom,iset,jset,ipgf,jpgf,natom,maxset,maxpgf)

       IF (iatom .NE. iatom_old .OR. jatom .NE. jatom_old) THEN 

          ikind = particle_set(iatom)%atomic_kind%kind_number
          jkind = particle_set(jatom)%atomic_kind%kind_number

          IF(iatom .NE. iatom_old ) ra(:) = pbc(particle_set(iatom)%r,cell)

          brow = iatom
          bcol = jatom

          IF (ikind .NE. ikind_old ) THEN
             CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
                  softb = my_soft, &
                  orb_basis_set=orb_basis_set)

             CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                  first_sgf=first_sgfa,&
                  lmax=la_max,&
                  lmin=la_min,&
                  npgf=npgfa,&
                  nset=nseta,&
                  nsgf_set=nsgfa,&
                  pgf_radius=rpgfa,&
                  set_radius=set_radius_a,&
                  sphi=sphi_a,&
                  zet=zeta)
          ENDIF

          IF (jkind .NE. jkind_old ) THEN

             CALL get_atomic_kind(atomic_kind=particle_set(jatom)%atomic_kind,&
                  softb = my_soft, &
                  orb_basis_set=orb_basis_set)

             CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                  first_sgf=first_sgfb,&
                  kind_radius=kind_radius_b,&
                  lmax=lb_max,&
                  lmin=lb_min,&
                  npgf=npgfb,&
                  nset=nsetb,&
                  nsgf_set=nsgfb,&
                  pgf_radius=rpgfb,&
                  set_radius=set_radius_b,&
                  sphi=sphi_b,&
                  zet=zetb)          

          ENDIF

          CALL get_block_node(matrix=deltajp,&
               block_row=brow,&
               block_col=bcol,&
               BLOCK=jp_block)
          CALL get_block_node(matrix=deltajp_rii,&
               block_row=brow,&
               block_col=bcol,&
               BLOCK=jp_block_rii)
          CALL get_block_node(matrix=deltajp_riii,&
               block_row=brow,&
               block_col=bcol,&
               BLOCK=jp_block_riii)
          IF (.NOT.ASSOCIATED(jp_block)) &
               CALL stop_program(routineP,"p_block not associated in deltap")

          iatom_old = iatom
          jatom_old = jatom
          ikind_old = ikind
          jkind_old = jkind

          atom_pair_changed = .TRUE.

       ELSE

          atom_pair_changed = .FALSE.

       ENDIF

       IF (atom_pair_changed .OR. iset_old .NE. iset .OR. jset_old .NE. jset) THEN

          ncoa = npgfa(iset)*ncoset(la_max(iset))
          sgfa = first_sgfa(1,iset)
          ncob = npgfb(jset)*ncoset(lb_max(jset))
          sgfb = first_sgfb(1,jset)
          ! Decontraction step for the selected blocks of the 3 density matrices
          CALL dgemm("N","N",ncoa,nsgfb(jset),nsgfa(iset),&
               1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
               jp_block(sgfa,sgfb),SIZE(jp_block,1),&
               0.0_dp,work(1,1),maxco)
          CALL dgemm("N","T",ncoa,ncob,nsgfb(jset),&
               1.0_dp,work(1,1),maxco,&
               sphi_b(1,sgfb),SIZE(sphi_b,1),&
               0.0_dp,jpab(1,1),maxco)
          CALL dgemm("N","N",ncoa,nsgfb(jset),nsgfa(iset),&
               1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
               jp_block_rii(sgfa,sgfb),SIZE(jp_block_rii,1),&
               0.0_dp,work(1,1),maxco)
          CALL dgemm("N","T",ncoa,ncob,nsgfb(jset),&
               1.0_dp,work(1,1),maxco,&
               sphi_b(1,sgfb),SIZE(sphi_b,1),&
               0.0_dp,jpab_ii(1,1),maxco)
          CALL dgemm("N","N",ncoa,nsgfb(jset),nsgfa(iset),&
               1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
               jp_block_riii(sgfa,sgfb),SIZE(jp_block_riii,1),&
               0.0_dp,work(1,1),maxco)
          CALL dgemm("N","T",ncoa,ncob,nsgfb(jset),&
               1.0_dp,work(1,1),maxco,&
               sphi_b(1,sgfb),SIZE(sphi_b,1),&
               0.0_dp,jpab_iii(1,1),maxco)

          iset_old = iset
          jset_old = jset

       ENDIF

       rab(:) = dist_ab (:,itask)
       rab2  = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
       rb(:) = ra(:) + rab(:)
       zetp = zeta(ipgf,iset) + zetb(jpgf,jset)

       na1 = (ipgf - 1)*ncoset(la_max(iset)) + 1
       na2 = ipgf*ncoset(la_max(iset))
       nb1 = (jpgf - 1)*ncoset(lb_max(jset)) + 1
       nb2 = jpgf*ncoset(lb_max(jset))

       scale = 1.0_dp

       ! Four calls to the general collocate density, to multply the correct function
       ! to each density matrix

       ! here the decontracted mat_jp_{ab} is multiplied by
       !     f_{ab} = g_{a} (dg_{b}/dr)_{idir} - (dg_{a}/dr)_{idir} g_{b}
       CALL collocate_pgf_product_rspace(la_max(iset),zeta(ipgf,iset),&
            la_min(iset),lb_max(jset),zetb(jpgf,jset),lb_min(jset),&
            ra,rab,rab2,scale,jpab,na1-1,nb1-1,&
            rs_rho(igrid_level)%rs_grid,cell,cube_info(igrid_level),&
            eps_rho_rspace,&
            ga_gb_function=FUNC_ADBmDAB,&
            idir=idir,&
            map_consistent=map_consistent,error=error)
       ! here the decontracted mat_jp_rii{ab} is multiplied by
       !     f_{ab} = g_{a} (r - R_{b})_{iiB} (dg_{b}/dr)_{idir} -
       !             (dg_{a}/dr)_{idir} (r - R_{b})_{iiB} g_{b}
       CALL collocate_pgf_product_rspace(la_max(iset),zeta(ipgf,iset),&
            la_min(iset),lb_max(jset),zetb(jpgf,jset),lb_min(jset),&
            ra,rab,rab2,scale,jpab_ii,na1-1,nb1-1,&
            rs_rho(igrid_level)%rs_grid,cell,cube_info(igrid_level),&
            eps_rho_rspace,&
            ga_gb_function=FUNC_ARDBmDARB,&
            idir=idir,ir=iiiB,&
            map_consistent=map_consistent,error=error)
       ! here the decontracted mat_jp_riii{ab} is multiplied by
       !     f_{ab} = -g_{a} (r - R_{b})_{iiB} (dg_{b}/dr)_{idir} +
       !             (dg_{a}/dr)_{idir} (r - R_{b})_{iiB} g_{b}
       scale = -1.0_dp
       CALL collocate_pgf_product_rspace(la_max(iset),zeta(ipgf,iset),&
            la_min(iset),lb_max(jset),zetb(jpgf,jset),lb_min(jset),&
            ra,rab,rab2,scale,jpab_iii,na1-1,nb1-1,&
            rs_rho(igrid_level)%rs_grid,cell,cube_info(igrid_level),&
            eps_rho_rspace,&
            ga_gb_function=FUNC_ARDBmDARB,&
            idir=idir,ir=iiB,&
            map_consistent=map_consistent,error=error)

    END DO loop_tasks

    !   *** Release work storage ***

    IF (distributed_rs_grids) THEN
       CALL deallocate_matrix ( deltajp ,error=error)
       CALL deallocate_matrix ( deltajp_rii ,error=error)
       CALL deallocate_matrix ( deltajp_riii ,error=error)
    END IF

    IF ( nthread > 1 ) THEN
       DEALLOCATE (lgrid%r,STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineP,"lgrid%r")
    END IF

    DEALLOCATE (jpabt,jpabt_ii,jpabt_iii,workt,tasks,dist_ab,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"jpabt,workt,tasks,dist_ab")

    IF (distributed_rs_grids) THEN
       DEALLOCATE(atom_pair_send,STAT=stat)
       IF (stat/=0) CALL stop_memory("collocate_density","atom_pair_send")        
       DEALLOCATE(atom_pair_recv,STAT=stat)
       IF (stat/=0) CALL stop_memory("collocate_density","atom_pair_recv")
    ENDIF

    CALL density_rs2pw(pw_env,rs_rho,rho_rs,rho_gs,&
         interp_section=interp_section,error=error)

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE calculate_jrho_resp

  SUBROUTINE current_build_chi(current_env,qs_env,iB,error)
    !
    TYPE(current_env_type)                   :: current_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    INTEGER, INTENT(IN)                      :: iB
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'current_build_chi', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, i, icenter, idir, idir2, ii, iiB, iii, iiiB, ispin, &
      istat, istate, j, nao, nbr_center(2), nspins, output_unit
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

    nspins = dft_control%nspins
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
    ! Allocate full matrices for only one vector
    CALL get_mo_set(mo_set=mos(1)%mo_set,mo_coeff=mo_coeff)
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
    ! get iiB and iiiB
    CALL set_vecp(iB,iiB,iiiB)
    DO ispin = 1,nspins
       !
       ! get ground state MOS
       CALL get_mo_set(mo_set=mos(ispin)%mo_set,mo_coeff=mo_coeff,&
            &          occupation_numbers=occupation)
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
       DO idir = 1,3
          current_env%chi_tensor(idir,iB,ispin) = current_env%chi_tensor(idir,iB,ispin)+&
               &                                  chi(idir)
       ENDDO
    ENDDO ! ispin
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

  END SUBROUTINE current_build_chi

END MODULE qs_linres_current

