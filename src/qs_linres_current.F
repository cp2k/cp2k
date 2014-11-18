!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
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
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_p_type,&
                                             gto_basis_set_type
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE cp_array_i_utils,                ONLY: cp_2d_i_p_type
  USE cp_array_r_utils,                ONLY: cp_2d_r_p_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
  USE cp_dbcsr_interface,              ONLY: &
       array_i1d_obj, array_new, array_nullify, array_release, &
       convert_offsets_to_sizes, cp_dbcsr_add_block_node, &
       cp_dbcsr_allocate_matrix_set, cp_dbcsr_col_block_sizes, cp_dbcsr_copy, &
       cp_dbcsr_create, cp_dbcsr_deallocate_matrix, &
       cp_dbcsr_deallocate_matrix_set, cp_dbcsr_distribution, &
       cp_dbcsr_finalize, cp_dbcsr_get_block_p, cp_dbcsr_get_data_size, &
       cp_dbcsr_get_matrix_type, cp_dbcsr_get_num_blocks, cp_dbcsr_init, &
       cp_dbcsr_p_type, cp_dbcsr_row_block_sizes, cp_dbcsr_set, &
       cp_dbcsr_type, dbcsr_distribution_obj, dbcsr_type_antisymmetric, &
       dbcsr_type_no_symmetry
  USE cp_dbcsr_operations,             ONLY: cp_dbcsr_plus_fm_fm_t,&
                                             cp_dbcsr_sm_fm_multiply
  USE cp_fm_basic_linalg,              ONLY: cp_fm_scale_and_add,&
                                             cp_fm_trace
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_set_all,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_realspace_grid_cube,          ONLY: cp_pw_to_cube
  USE cube_utils,                      ONLY: cube_info_type
  USE gaussian_gridlevels,             ONLY: gridlevel_info_type
  USE input_constants,                 ONLY: current_gauge_atom
  USE input_section_types,             ONLY: section_get_ivals,&
                                             section_get_lval,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type
  USE kinds,                           ONLY: default_path_length,&
                                             dp,&
                                             int_8
  USE lgrid_types,                     ONLY: lgrid_type
  USE mathconstants,                   ONLY: twopi
  USE memory_utilities,                ONLY: reallocate
  USE orbital_pointers,                ONLY: ncoset
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_methods,                ONLY: get_particle_set
  USE particle_types,                  ONLY: particle_type
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_methods,                      ONLY: pw_axpy,&
                                             pw_integrate_function,&
                                             pw_scale,&
                                             pw_zero
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_type
  USE pw_types,                        ONLY: REALDATA3D,&
                                             REALSPACE,&
                                             pw_p_type
  USE qs_collocate_density,            ONLY: collocate_pgf_product_rspace
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_kind_types,                   ONLY: get_qs_kind,&
                                             get_qs_kind_set,&
                                             qs_kind_type
  USE qs_linres_atom_current,          ONLY: calculate_jrho_atom,&
                                             calculate_jrho_atom_coeff,&
                                             calculate_jrho_atom_rad
  USE qs_linres_op,                    ONLY: fac_vecp,&
                                             fm_scale_by_pbc_AC,&
                                             ind_m2,&
                                             set_vecp,&
                                             set_vecp_rev
  USE qs_linres_types,                 ONLY: current_env_type,&
                                             get_current_env,&
                                             realspaces_grid_p_type
  USE qs_matrix_pools,                 ONLY: qs_matrix_pools_type
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
  USE qs_modify_pab_block,             ONLY: FUNC_AB,&
                                             FUNC_ADBmDAB,&
                                             FUNC_ARDBmDARB
  USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                             neighbor_list_iterate,&
                                             neighbor_list_iterator_create,&
                                             neighbor_list_iterator_p_type,&
                                             neighbor_list_iterator_release,&
                                             neighbor_list_set_p_type
  USE qs_operators_ao,                 ONLY: build_lin_mom_matrix,&
                                             rRc_xyz_der_ao
  USE qs_rho_types,                    ONLY: qs_rho_get
  USE qs_subsys_types,                 ONLY: qs_subsys_get,&
                                             qs_subsys_type
  USE realspace_grid_types,            ONLY: &
       realspace_grid_desc_p_type, realspace_grid_desc_type, &
       realspace_grid_p_type, realspace_grid_type, rs_grid_create, &
       rs_grid_mult_and_add, rs_grid_release, rs_grid_retain, rs_grid_zero
  USE rs_pw_interface,                 ONLY: density_rs2pw
  USE task_list_methods,               ONLY: distribute_tasks,&
                                             int2pair,&
                                             rs_distribute_matrix,&
                                             task_list_inner_loop
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  ! *** Public subroutines ***
  PUBLIC :: current_build_current, current_build_chi

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_linres_current'

  TYPE box_type
     INTEGER :: n
     REAL(dp), POINTER, DIMENSION(:,:) :: r
  END TYPE box_type
 REAL(dp), DIMENSION(3, 3, 3),PARAMETER  :: Levi_Civita = RESHAPE((/&
       0.0_dp, 0.0_dp,0.0_dp,0.0_dp,0.0_dp,-1.0_dp, 0.0_dp,1.0_dp,0.0_dp,&
       0.0_dp, 0.0_dp,1.0_dp,0.0_dp,0.0_dp, 0.0_dp,-1.0_dp,0.0_dp,0.0_dp,&
       0.0_dp,-1.0_dp,0.0_dp,1.0_dp,0.0_dp, 0.0_dp, 0.0_dp,0.0_dp,0.0_dp/),(/3,3,3/))


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
!> \param current_env ...
!> \param qs_env ...
!> \param iB ...
!> \param error ...
!> \author MI
!> \note
!>       The susceptibility is needed to compute the G=0 term of the shift
!>       in reciprocal space. \chi_{ij} = \int (r x Jj)_i
!>       (where Jj id the current density generated by the field in direction j)
!>       To calculate the susceptibility on the PW grids it is necessary to apply
!>       the position operator yet another time.
!>       This cannot be done on directly on the full J(r) because it is not localized
!>       Therefore it is done state by state (see linres_nmr_shift)
! *****************************************************************************
  SUBROUTINE current_build_current(current_env,qs_env,iB,error)
    !
    TYPE(current_env_type)                   :: current_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    INTEGER, INTENT(IN)                      :: iB
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'current_build_current', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: ext, filename, my_pos
    INTEGER :: handle, idir, iiB, iiiB, ispin, istat, istate, j, jstate, nao, &
      natom, nmo, nspins, nstates(2), output_unit, unit_nr
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: first_sgf, last_sgf
    INTEGER, DIMENSION(:), POINTER           :: rbs
    LOGICAL                                  :: append_cube, failure, gapw
    REAL(dp)                                 :: dk(3), jrho_tot_G(3,3), &
                                                jrho_tot_R(3,3), maxocc, &
                                                scale_fac
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: ddk
    REAL(dp), EXTERNAL                       :: DDOT
    TYPE(array_i1d_obj)                      :: row_blk_sizes
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_2d_i_p_type), DIMENSION(:), &
      POINTER                                :: center_list
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: density_matrix0, &
                                                density_matrix_a, &
                                                density_matrix_ii, &
                                                density_matrix_iii
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: p_psi1, psi0_order, psi1
    TYPE(cp_fm_p_type), DIMENSION(:, :), &
      POINTER                                :: psi1_D, psi1_p, psi1_rxp
    TYPE(cp_fm_type), POINTER                :: mo_coeff, psi_a_iB, psi_buf
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dbcsr_distribution_obj), POINTER    :: dbcsr_dist
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_all
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type)                          :: wf_r
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: jrho1_g, jrho1_r
    TYPE(pw_p_type), POINTER                 :: jrho_gspace, jrho_rspace
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(qs_matrix_pools_type), POINTER      :: mpools
    TYPE(qs_subsys_type), POINTER            :: subsys
    TYPE(realspace_grid_desc_type), POINTER  :: auxbas_rs_desc
    TYPE(section_vals_type), POINTER         :: current_section

    failure = .FALSE.
    CALL timeset(routineN,handle)
    !
    NULLIFY(jrho_rspace,jrho_gspace,logger,current_section,density_matrix0,density_matrix_a,&
            density_matrix_ii,density_matrix_iii,cell,dft_control,mos,&
            particle_set,pw_env,auxbas_rs_desc,auxbas_pw_pool,&
            para_env,center_list,mo_coeff,psi_a_iB,jrho1_r, jrho1_g,&
            psi1,p_psi1,psi1_p,psi1_D,psi1_rxp,psi0_order,sab_all, qs_kind_set)

    logger => cp_error_get_logger(error)
    output_unit= cp_logger_get_default_io_unit(logger)
    !
    !
    CALL get_current_env(current_env=current_env,&
                         center_list=center_list,&
                         psi1_rxp=psi1_rxp,&
                         psi1_D=psi1_D,&
                         psi1_p=psi1_p,&
                         psi0_order=psi0_order,&
                         nstates=nstates,&
                         nao=nao,&
                         error=error)
    !
    !
    CALL get_qs_env(qs_env=qs_env,&
                    cell=cell,&
                    dft_control=dft_control,&
                    mos=mos,&
                    mpools=mpools,&
                    pw_env=pw_env,&
                    para_env=para_env,&
                    subsys=subsys,&
                    sab_all=sab_all,&
                    particle_set=particle_set,&
                    qs_kind_set=qs_kind_set,&
                    dbcsr_dist=dbcsr_dist,&
                    error=error)

    CALL qs_subsys_get(subsys,particles=particles,error=error)

    gapw  = dft_control%qs_control%gapw
    nspins = dft_control%nspins
    natom = SIZE(particle_set,1)
    !
    ! allocate temporary arrays
    ALLOCATE(psi1(nspins),p_psi1(nspins),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DO ispin = 1,nspins
       CALL cp_fm_create(  psi1(ispin)%matrix,psi0_order(ispin)%matrix%matrix_struct,error=error)
       CALL cp_fm_create(p_psi1(ispin)%matrix,psi0_order(ispin)%matrix%matrix_struct,error=error)
       CALL cp_fm_set_all(  psi1(ispin)%matrix,0.0_dp,error=error)
       CALL cp_fm_set_all(p_psi1(ispin)%matrix,0.0_dp,error=error)
    ENDDO
    !
    !
    CALL cp_dbcsr_allocate_matrix_set(density_matrix0,nspins,error=error)
    CALL cp_dbcsr_allocate_matrix_set(density_matrix_a,nspins,error=error)
    CALL cp_dbcsr_allocate_matrix_set(density_matrix_ii,nspins,error=error)
    CALL cp_dbcsr_allocate_matrix_set(density_matrix_iii,nspins,error=error)
    !
    ! prepare for allocation
    ALLOCATE (first_sgf(natom),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (last_sgf(natom),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    CALL get_particle_set(particle_set, qs_kind_set,&
                          first_sgf=first_sgf,&
                          last_sgf=last_sgf,error=error)
    ALLOCATE (rbs(natom), STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    CALL convert_offsets_to_sizes (first_sgf, rbs, last_sgf)
    CALL array_nullify (row_blk_sizes)
    CALL array_new (row_blk_sizes, rbs, gift=.TRUE.)
    DEALLOCATE (first_sgf,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (last_sgf,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    !
    !
    DO ispin = 1,nspins
       !
       !density_matrix0
       ALLOCATE(density_matrix0(ispin)%matrix)
       CALL cp_dbcsr_init(density_matrix0(ispin)%matrix,error=error)
       CALL cp_dbcsr_create(matrix=density_matrix0(ispin)%matrix, &
            name="density_matrix0", &
            dist=dbcsr_dist, matrix_type=dbcsr_type_no_symmetry,&
            row_blk_size=row_blk_sizes, col_blk_size=row_blk_sizes, &
            nblks=0, nze=0, mutable_work=.TRUE., &
            error=error)
       CALL cp_dbcsr_alloc_block_from_nbl(density_matrix0(ispin)%matrix,sab_all,error=error)
       !
       !density_matrix_a
       ALLOCATE(density_matrix_a(ispin)%matrix)
       CALL cp_dbcsr_init(density_matrix_a(ispin)%matrix,error=error)
       CALL cp_dbcsr_copy(density_matrix_a(ispin)%matrix,density_matrix0(ispin)%matrix,&
            name="density_matrix_a",error=error)
       !
       !density_matrix_ii
       ALLOCATE(density_matrix_ii(ispin)%matrix)
       CALL cp_dbcsr_init(density_matrix_ii(ispin)%matrix,error=error)
       CALL cp_dbcsr_copy(density_matrix_ii(ispin)%matrix,density_matrix0(ispin)%matrix,&
            name="density_matrix_ii",error=error)
       !
       !density_matrix_iii
       ALLOCATE(density_matrix_iii(ispin)%matrix)
       CALL cp_dbcsr_init(density_matrix_iii(ispin)%matrix,error=error)
       CALL cp_dbcsr_copy(density_matrix_iii(ispin)%matrix,density_matrix0(ispin)%matrix,&
            name="density_matrix_iii",error=error)
    ENDDO
    !
    CALL array_release (row_blk_sizes)
    !
    !
    current_section => section_vals_get_subs_vals(qs_env%input,"PROPERTIES%LINRES%CURRENT",error=error)
    !
    IF(.NOT. failure) THEN
       !
       !
       jrho_tot_G = 0.0_dp
       jrho_tot_R = 0.0_dp
       !
       ! Lets go!
       CALL set_vecp(iB,iiB,iiiB)
       DO ispin = 1,nspins
          nmo = nstates(ispin)
          mo_coeff => psi0_order(ispin)%matrix
          !maxocc = max_occ(ispin)
          !
          CALL get_mo_set(mo_set=mos(ispin)%mo_set,maxocc=maxocc)
          !
          !
          ! Build the first density matrix
          CALL cp_dbcsr_set(density_matrix0(ispin)%matrix,0.0_dp,error=error)
          CALL cp_dbcsr_plus_fm_fm_t(sparse_matrix=density_matrix0(ispin)%matrix ,&
                                  matrix_v=mo_coeff,matrix_g=mo_coeff,&
                                  ncol=nmo,alpha=maxocc,error=error)
          !
          ! Allocate buffer vectors
          ALLOCATE(ddk(3,nmo),STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          !
          ! Construct the 3 density matrices for the field in direction iB
          !
          ! First the full matrix psi_a_iB
          psi_a_iB => psi1(ispin)%matrix
          psi_buf => p_psi1(ispin)%matrix
          CALL cp_fm_set_all(psi_a_iB,0.0_dp,error=error)
          CALL cp_fm_set_all(psi_buf,0.0_dp,error=error)
          ! psi_a_iB = - (R_\nu-dk)_ii psi1_piiiB + (R_\nu-dk)_iii psi1_piiB
          !
          ! contributions from the response psi1_p_ii and psi1_p_iii
          DO istate = 1,current_env%nbr_center(ispin)
             dk(1:3) = current_env%centers_set(ispin)%array(1:3,istate)
             !
             ! Copy the vector in the full matrix psi1
             !nstate_loc = center_list(ispin)%array(1,icenter+1)-center_list(ispin)%array(1,icenter)
             DO j = center_list(ispin)%array(1,istate),center_list(ispin)%array(1,istate+1)-1
                jstate = center_list(ispin)%array(2,j)
                CALL cp_fm_to_fm(psi1_p(ispin,iiB)%matrix,psi_a_iB,1,jstate,jstate)
                CALL cp_fm_to_fm(psi1_p(ispin,iiiB)%matrix,psi_buf,1,jstate,jstate)
                ddk(:,jstate) = dk(1:3)
             ENDDO
          ENDDO ! istate
          CALL fm_scale_by_pbc_AC(psi_a_iB,current_env%basisfun_center,ddk,cell,iiiB,error)
          CALL fm_scale_by_pbc_AC(psi_buf,current_env%basisfun_center,ddk,cell,iiB,error)
          CALL cp_fm_scale_and_add(-1.0_dp,psi_a_iB,1.0_dp,psi_buf,error=error)
          !
          !psi_a_iB = psi_a_iB + psi1_rxp
          !
          ! contribution from the response psi1_rxp
          CALL cp_fm_scale_and_add(-1.0_dp,psi_a_iB,1.0_dp,psi1_rxp(ispin,iB)%matrix,error=error)
          !
          !psi_a_iB = psi_a_iB - psi1_D
          IF(current_env%full) THEN
             !
             ! contribution from the response psi1_D
             CALL  cp_fm_scale_and_add(1.0_dp,psi_a_iB,-1.0_dp,psi1_D(ispin,iB)%matrix,error=error)
          ENDIF
          !
          ! Multiply by the occupation number for the density matrix
          !
          ! Build the first density matrix
          CALL cp_dbcsr_set(density_matrix_a(ispin)%matrix ,0.0_dp,error=error)
          CALL cp_dbcsr_plus_fm_fm_t(sparse_matrix=density_matrix_a(ispin)%matrix ,&
                                  matrix_v=mo_coeff,matrix_g=psi_a_iB,&
                                  ncol=nmo,alpha=maxocc,error=error)
          !
          ! Build the second density matrix
          CALL cp_dbcsr_set(density_matrix_iii(ispin)%matrix,0.0_dp,error=error)
          CALL cp_dbcsr_plus_fm_fm_t(sparse_matrix=density_matrix_iii(ispin)%matrix ,&
                                  matrix_v=mo_coeff,matrix_g=psi1_p(ispin,iiiB)%matrix,&
                                  ncol=nmo,alpha=maxocc,error=error)
          !
          ! Build the third density matrix
          CALL cp_dbcsr_set(density_matrix_ii(ispin)%matrix,0.0_dp,error=error)
          CALL cp_dbcsr_plus_fm_fm_t(sparse_matrix=density_matrix_ii(ispin)%matrix ,&
                                  matrix_v=mo_coeff,matrix_g=psi1_p(ispin,iiB)%matrix,&
                                  ncol=nmo,alpha=maxocc,error=error)
          DO idir = 1,3
             !
             ! Calculate the current density on the pw grid (only soft if GAPW)
             ! idir is the cartesian component of the response current density
             ! generated by the magnetic field pointing in cartesian direction iB
             ! Use the qs_rho_type already  used for rho during the scf
             CALL qs_rho_get(current_env%jrho1_set(idir)%rho, rho_r=jrho1_r, error=error)
             CALL qs_rho_get(current_env%jrho1_set(idir)%rho, rho_g=jrho1_g, error=error)
             jrho_rspace => jrho1_r(ispin)
             jrho_gspace => jrho1_g(ispin)
             CALL pw_zero(jrho_rspace%pw,error=error)
             CALL pw_zero(jrho_gspace%pw,error=error)
             CALL calculate_jrho_resp(density_matrix0(ispin)%matrix, &
                                      density_matrix_a(ispin)%matrix, &
                                      density_matrix_ii(ispin)%matrix, &
                                      density_matrix_iii(ispin)%matrix, &
                                      iB,idir,jrho_rspace,jrho_gspace,qs_env, &
                                      current_env,gapw,error=error)

             scale_fac = cell%deth / twopi
             CALL pw_scale(jrho_rspace%pw,scale_fac,error=error)
             CALL pw_scale(jrho_gspace%pw,scale_fac,error=error)

             jrho_tot_G(idir,iB) = pw_integrate_function(jrho_gspace%pw,isign=-1,error=error)
             jrho_tot_R(idir,iB) = pw_integrate_function(jrho_rspace%pw,isign=-1,error=error)

             IF(output_unit>0) THEN
                WRITE(output_unit,'(T2,2(A,E24.16))') 'Integrated j_'&
                     &//ACHAR(idir+119)//ACHAR(iB+119)//'(r): G-space=',&
                      jrho_tot_G(idir,iB),' R-space=',jrho_tot_R(idir,iB)
             ENDIF
             !
          ENDDO ! idir
          !
          ! Deallocate buffer vectors
          DEALLOCATE(ddk,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          !
       ENDDO ! ispin

       IF(gapw) THEN
          DO idir = 1,3
             !
             ! compute the atomic response current densities on the spherical grids
             ! First the sparse matrices are multiplied by the expansion coefficients
             ! this is the equivalent of the CPC for the charge density
             CALL calculate_jrho_atom_coeff(qs_env,current_env, &
                                            density_matrix0, &
                                            density_matrix_a, &
                                            density_matrix_ii, &
                                            density_matrix_iii, &
                                            iB,idir,error=error)
             !
             ! Then the radial parts are computed on the local radial grid, atom by atom
             ! 8 functions are computed for each atom, per grid point
             ! and per LM angular momentum. The multiplication by the Clebsh-Gordon
             ! coefficients or they correspondent for the derivatives, is also done here
             CALL calculate_jrho_atom_rad(qs_env,current_env,idir,error=error)
             !
             ! The current on the atomic grids
             CALL calculate_jrho_atom(current_env,qs_env,iB,idir,error=error)
          ENDDO ! idir
       ENDIF
       !
       ! Cube files
       IF(BTEST(cp_print_key_should_output(logger%iter_info,current_section,&
            &   "PRINT%CURRENT_CUBES",error=error),cp_p_file)) THEN
          append_cube = section_get_lval(current_section,"PRINT%CURRENT_CUBES%APPEND",error=error)
          my_pos = "REWIND"
          IF(append_cube) THEN
            my_pos = "APPEND"
          END IF
          !
          CALL pw_env_get(pw_env, auxbas_rs_desc=auxbas_rs_desc,&
                          auxbas_pw_pool=auxbas_pw_pool,&
                          error=error)
          !
          CALL pw_pool_create_pw(auxbas_pw_pool,wf_r%pw,use_data=REALDATA3D,&
                                  in_space=REALSPACE,error=error)
          !
          DO idir = 1,3
             CALL pw_zero(wf_r%pw,error=error)
             CALL qs_rho_get(current_env%jrho1_set(idir)%rho, rho_r=jrho1_r, error=error)
             DO ispin =1 ,nspins
                CALL pw_axpy(jrho1_r(ispin)%pw,wf_r%pw,1.0_dp, error=error)
             ENDDO
             !
             IF(gapw) THEN
                ! Add the local hard and soft contributions
                ! This can be done atom by atom by a spline extrapolation of the  values
                ! on the spherical grid to the grid points.
                CALL stop_program(routineN,moduleN,__LINE__,"GAPW needs to be finalized")
             ENDIF
             filename="jresp"
             WRITE(ext,'(a2,I1,a2,I1,a5)')  "iB",iB,"_d",idir,".cube"
             WRITE(ext,'(a2,a1,a2,a1,a5)')  "iB",ACHAR(iB+119),"_d",ACHAR(idir+119),".cube"
             unit_nr=cp_print_key_unit_nr(logger,current_section,"PRINT%CURRENT_CUBES",&
                                          extension=TRIM(ext),middle_name=TRIM(filename),&
                                          log_filename=.FALSE.,file_position=my_pos,&
                                          error=error)

             CALL cp_pw_to_cube(wf_r%pw,unit_nr,"RESPONSE CURRENT DENSITY ",&
                             particles=particles,&
                             stride=section_get_ivals(current_section,"PRINT%CURRENT_CUBES%STRIDE",&
                             error=error),error=error)
             CALL cp_print_key_finished_output(unit_nr,logger,current_section,&
                  &                            "PRINT%CURRENT_CUBES",error=error)
          ENDDO
          !
          CALL pw_pool_give_back_pw(auxbas_pw_pool,wf_r%pw,error=error)
       ENDIF ! current cube
       !
       ! Integrated current response checksum
       IF(output_unit>0) THEN
          WRITE(output_unit,'(T2,A,E24.16)') 'CheckSum R-integrated j=',&
                            SQRT(DDOT(9,jrho_tot_R(1,1),1,jrho_tot_R(1,1),1))
       ENDIF
       !
       !
       ! Dellocate grids for the calculation of jrho and the shift
       DO ispin = 1,nspins
          CALL cp_fm_release(  psi1(ispin)%matrix,error)
          CALL cp_fm_release(p_psi1(ispin)%matrix,error)
       ENDDO
       DEALLOCATE(psi1,p_psi1,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

       CALL cp_dbcsr_deallocate_matrix_set(density_matrix0,error=error)
       CALL cp_dbcsr_deallocate_matrix_set(density_matrix_a,error=error)
       CALL cp_dbcsr_deallocate_matrix_set(density_matrix_ii,error=error)
       CALL cp_dbcsr_deallocate_matrix_set(density_matrix_iii,error=error)

       !
    ENDIF ! failure
    !
    ! Finalize
    CALL timestop(handle)
    !
  END SUBROUTINE current_build_current

! *****************************************************************************
!> \brief Calculation of the idir component of the response current density
!>       in the presence of a constant magnetic field in direction iB
!>       the current density is collocated on the pw grid in real space
!> \param mat_d0 ...
!> \param mat_jp ...
!> \param mat_jp_rii ...
!> \param mat_jp_riii ...
!> \param iB ...
!> \param idir ...
!> \param current_rs ...
!> \param current_gs ...
!> \param qs_env ...
!> \param current_env ...
!> \param soft_valid ...
!> \param error ...
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
  SUBROUTINE calculate_jrho_resp(mat_d0,mat_jp,mat_jp_rii,mat_jp_riii,iB,idir,&
       current_rs, current_gs, qs_env, current_env, soft_valid, error)

    TYPE(cp_dbcsr_type), POINTER             :: mat_d0, mat_jp, mat_jp_rii, &
                                                mat_jp_riii
    INTEGER, INTENT(IN)                      :: iB, idir
    TYPE(pw_p_type), INTENT(INOUT)           :: current_rs, current_gs
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(current_env_type)                   :: current_env
    LOGICAL, INTENT(IN), OPTIONAL            :: soft_valid
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_jrho_resp', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: add_tasks = 1000, &
                                                max_tasks = 2000
    REAL(kind=dp), PARAMETER                 :: mult_tasks = 2.0_dp

    INTEGER :: bcol, brow, curr_tasks, handle, i, iatom, iatom_old, idir2, &
      igrid_level, iiB, iiiB, ikind, ikind_old, ipgf, iset, iset_old, istat, &
      itask, ithread, jatom, jatom_old, jkind, jkind_old, jpgf, jset, &
      jset_old, maxco, maxpgf, maxset, maxsgf, maxsgf_set, na1, na2, natom, &
      nb1, nb2, ncoa, ncob, ngpts, nkind, nseta, nsetb, ntasks, nthread, &
      sgfa, sgfb
    INTEGER(kind=int_8), DIMENSION(:), &
      POINTER                                :: atom_pair_recv, atom_pair_send
    INTEGER(kind=int_8), DIMENSION(:, :), &
      POINTER                                :: tasks
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, lb_max, &
                                                lb_min, npgfa, npgfb, nsgfa, &
                                                nsgfb
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb
    LOGICAL :: atom_pair_changed, den_found, den_found_a, &
      distributed_rs_grids, do_igaim, failure, map_consistent, &
      my_minimum_image, my_soft
    REAL(dp), DIMENSION(:, :, :), POINTER    :: my_current, my_gauge, my_rho
    REAL(KIND=dp)                            :: eps_rho_rspace, &
                                                kind_radius_a, kind_radius_b, &
                                                Lxo2, Lyo2, Lzo2, rab2, &
                                                scale, scale2, zetp
    REAL(KIND=dp), DIMENSION(3)              :: ra, rab, rb
    REAL(KIND=dp), DIMENSION(:), POINTER     :: set_radius_a, set_radius_b
    REAL(KIND=dp), DIMENSION(:, :), POINTER :: dist_ab, jp_block_a, &
      jp_block_b, jp_block_c, jp_block_d, jpab_a, jpab_b, jpab_c, jpab_d, &
      jpblock_a, jpblock_b, jpblock_c, jpblock_d, rpgfa, rpgfb, sphi_a, &
      sphi_b, work, zeta, zetb
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: jpabt_a, jpabt_b, jpabt_c, &
                                                jpabt_d, workt
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_dbcsr_type), POINTER             :: deltajp_a, deltajp_b, &
                                                deltajp_c, deltajp_d, mat_a, &
                                                mat_b, mat_c, mat_d
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cube_info_type), DIMENSION(:), &
      POINTER                                :: cube_info
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(gridlevel_info_type), POINTER       :: gridlevel_info
    TYPE(gto_basis_set_p_type), &
      DIMENSION(:), POINTER                  :: basis_set_list
    TYPE(gto_basis_set_type), POINTER        :: basis_set_a, basis_set_b, &
                                                orb_basis_set
    TYPE(lgrid_type)                         :: lgrid
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: nl_iterator
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(qs_kind_type), POINTER              :: qs_kind
    TYPE(realspace_grid_desc_p_type), &
      DIMENSION(:), POINTER                  :: rs_descs
    TYPE(realspace_grid_p_type), &
      DIMENSION(:), POINTER                  :: rs_current, rs_rho
    TYPE(realspaces_grid_p_type), &
      DIMENSION(:), POINTER                  :: rs_gauge

    failure=.FALSE.
    NULLIFY(qs_kind, cell, dft_control, orb_basis_set, rs_rho, &
         qs_kind_set, sab_orb, particle_set, rs_current, pw_env, &
         rs_descs, para_env, dist_ab, set_radius_a, set_radius_b, la_max, &
         la_min, lb_max, lb_min, npgfa, npgfb, nsgfa, nsgfb, rpgfa, rpgfb, &
         sphi_a, sphi_b, zeta, zetb, first_sgfa, first_sgfb, dist_ab, &
         tasks, workt, mat_a, mat_b, mat_c, mat_d, rs_gauge)
    NULLIFY(deltajp_a,deltajp_b,deltajp_c,deltajp_d)
    NULLIFY(jp_block_a,jp_block_b,jp_block_c,jp_block_d)
    NULLIFY(jpblock_a,jpblock_b,jpblock_c,jpblock_d)
    NULLIFY(jpabt_a,jpabt_b,jpabt_c,jpabt_d)
    NULLIFY(lgrid%r)

    CALL timeset(routineN,handle)

    !
    ! Set pointers for the different gauge
    do_igaim = current_env%gauge.EQ.current_gauge_atom

    IF(do_igaim) THEN
       mat_a => mat_jp
       mat_b => mat_jp_rii
       mat_c => mat_jp_riii
       mat_d => mat_d0
    ELSE
       mat_a => mat_jp
       mat_b => mat_jp_rii
       mat_c => mat_jp_riii
    ENDIF

    CALL get_qs_env(qs_env=qs_env,&
         qs_kind_set=qs_kind_set,&
         cell=cell,&
         dft_control=dft_control,&
         particle_set=particle_set,&
         sab_all=sab_orb,&
         para_env=para_env,&
         pw_env=pw_env,error=error)

    CALL get_current_env(current_env=current_env,&
                         rs_gauge=rs_gauge,&
                         error=error)

    ! Component of appearing in the vector product rxp, iiB and iiiB
    CALL set_vecp(iB,iiB,iiiB)
    !
    !
    scale2 = 0.0_dp
    idir2 = 1
    IF(idir.NE.iB) THEN
       CALL set_vecp_rev(idir,iB,idir2)
       scale2 = fac_vecp(idir,iB,idir2)
    ENDIF
    !
    ! *** assign from pw_env
    gridlevel_info=>pw_env%gridlevel_info
    cube_info=>pw_env%cube_info

    !   Check that the neighbor list with all the pairs is associated
    CPPrecondition(ASSOCIATED(sab_orb),cp_failure_level,routineP,error,failure)
    ! *** set up the pw multi-grids
    CPPrecondition(ASSOCIATED(pw_env),cp_failure_level,routineP,error,failure)
    CALL pw_env_get(pw_env, rs_descs=rs_descs, rs_grids=rs_rho, error=error)

    distributed_rs_grids=.FALSE.
    DO igrid_level=1,gridlevel_info%ngrid_levels
       IF(.NOT.ALL(rs_descs(igrid_level)%rs_desc%perd == 1)) THEN
          distributed_rs_grids = .TRUE.
       ENDIF
    ENDDO
    eps_rho_rspace = dft_control%qs_control%eps_rho_rspace
    map_consistent = dft_control%qs_control%map_consistent
    nthread = 1

    !   *** Allocate work storage ***
    CALL get_qs_kind_set(qs_kind_set=qs_kind_set,&
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

    nkind = SIZE(qs_kind_set)

    CALL reallocate(jpabt_a,1,maxco,1,maxco,0,nthread-1)
    CALL reallocate(jpabt_b,1,maxco,1,maxco,0,nthread-1)
    CALL reallocate(jpabt_c,1,maxco,1,maxco,0,nthread-1)
    CALL reallocate(jpabt_d,1,maxco,1,maxco,0,nthread-1)
    CALL reallocate(workt,1,maxco,1,maxsgf_set,0,nthread-1)
    CALL reallocate(tasks,1,6,1,max_tasks)
    CALL reallocate(dist_ab,1,3,1,max_tasks)

    tasks = 0
    ntasks = 0
    curr_tasks = SIZE(tasks,2)

    !   get maximum numbers
    natom = SIZE( particle_set )
    maxset=0
    maxpgf=0

    DO ikind=1,nkind
       qs_kind => qs_kind_set(ikind)

       CALL get_qs_kind(qs_kind=qs_kind,&
            orb_basis_set=orb_basis_set)

       IF (.NOT.ASSOCIATED(orb_basis_set)) CYCLE

       CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
            npgf=npgfa, nset=nseta )

       maxset=MAX(nseta,maxset)
       maxpgf=MAX(MAXVAL(npgfa),maxpgf)
    END DO

    !   *** Initialize working density matrix ***

    ! distributed rs grids require a matrix that will be changed (distribute_tasks)
    ! whereas this is not the case for replicated grids
    IF (distributed_rs_grids) THEN
       ALLOCATE(deltajp_a,deltajp_b,deltajp_c)
       CALL cp_dbcsr_init(deltajp_a, error=error)
       CALL cp_dbcsr_init(deltajp_b, error=error)
       CALL cp_dbcsr_init(deltajp_c, error=error)
       IF(do_igaim) THEN
          ALLOCATE(deltajp_d)
          CALL cp_dbcsr_init(deltajp_d, error=error)
       ENDIF

       CALL cp_dbcsr_create(deltajp_a, ' deltajp_a ', &
            cp_dbcsr_distribution (mat_a), cp_dbcsr_get_matrix_type (mat_a),&
            cp_dbcsr_row_block_sizes(mat_a),&
            cp_dbcsr_col_block_sizes(mat_a), cp_dbcsr_get_num_blocks(mat_a), cp_dbcsr_get_data_size(mat_a),&
            error=error)
       CALL cp_dbcsr_create(deltajp_b, ' deltajp_b ', &
            cp_dbcsr_distribution (mat_a), cp_dbcsr_get_matrix_type (mat_a),&
            cp_dbcsr_row_block_sizes(mat_a),&
            cp_dbcsr_col_block_sizes(mat_a), cp_dbcsr_get_num_blocks(mat_a), cp_dbcsr_get_data_size(mat_a),&
            error=error)
       CALL cp_dbcsr_create(deltajp_c, ' deltajp_c ', &
            cp_dbcsr_distribution (mat_a), cp_dbcsr_get_matrix_type (mat_a),&
            cp_dbcsr_row_block_sizes(mat_a),&
            cp_dbcsr_col_block_sizes(mat_a), cp_dbcsr_get_num_blocks(mat_a), cp_dbcsr_get_data_size(mat_a),&
            error=error)
       IF(do_igaim) CALL cp_dbcsr_create(deltajp_d, ' deltajp_d ', &
            cp_dbcsr_distribution (mat_a), cp_dbcsr_get_matrix_type (mat_a),&
            cp_dbcsr_row_block_sizes(mat_a),&
            cp_dbcsr_col_block_sizes(mat_a), cp_dbcsr_get_num_blocks(mat_a), cp_dbcsr_get_data_size(mat_a),&
            error=error)
    ELSE
       deltajp_a => mat_a !mat_jp
       deltajp_b => mat_b !mat_jp_rii
       deltajp_c => mat_c !mat_jp_riii
       IF(do_igaim) deltajp_d => mat_d !mat_d0
    ENDIF

    ALLOCATE (basis_set_list(nkind),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DO ikind=1,nkind
      qs_kind => qs_kind_set(ikind)
      CALL get_qs_kind(qs_kind=qs_kind,softb=my_soft,orb_basis_set=basis_set_a)
      IF (ASSOCIATED(basis_set_a)) THEN
        basis_set_list(ikind)%gto_basis_set => basis_set_a
      ELSE
        NULLIFY(basis_set_list(ikind)%gto_basis_set)
      END IF
    END DO
    CALL neighbor_list_iterator_create(nl_iterator,sab_orb)
    DO WHILE (neighbor_list_iterate(nl_iterator)==0)
       CALL get_iterator_info(nl_iterator,ikind=ikind,jkind=jkind,iatom=iatom,jatom=jatom,r=rab)
       basis_set_a => basis_set_list(ikind)%gto_basis_set
       IF (.NOT.ASSOCIATED(basis_set_a)) CYCLE
       basis_set_b => basis_set_list(jkind)%gto_basis_set
       IF (.NOT.ASSOCIATED(basis_set_b)) CYCLE
       ra(:) = pbc(particle_set(iatom)%r,cell)
       ! basis ikind
       first_sgfa   =>  basis_set_a%first_sgf
       la_max       =>  basis_set_a%lmax
       la_min       =>  basis_set_a%lmin
       npgfa        =>  basis_set_a%npgf
       nseta        =   basis_set_a%nset
       nsgfa        =>  basis_set_a%nsgf_set
       rpgfa        =>  basis_set_a%pgf_radius
       set_radius_a =>  basis_set_a%set_radius
       kind_radius_a=   basis_set_a%kind_radius
       sphi_a       =>  basis_set_a%sphi
       zeta         =>  basis_set_a%zet
       ! basis jkind
       first_sgfb   =>  basis_set_b%first_sgf
       lb_max       =>  basis_set_b%lmax
       lb_min       =>  basis_set_b%lmin
       npgfb        =>  basis_set_b%npgf
       nsetb        =   basis_set_b%nset
       nsgfb        =>  basis_set_b%nsgf_set
       rpgfb        =>  basis_set_b%pgf_radius
       set_radius_b =>  basis_set_b%set_radius
       kind_radius_b=   basis_set_b%kind_radius
       sphi_b       =>  basis_set_b%sphi
       zetb         =>  basis_set_b%zet

       IF(my_minimum_image) THEN
          IF(ABS(rab(1)) > Lxo2 .OR. ABS(rab(2)) > Lyo2 .OR. ABS(rab(3)) > Lzo2) THEN
             CYCLE
          END IF
       END IF

       brow = iatom
       bcol = jatom

       CALL cp_dbcsr_get_block_p(matrix=mat_a,row=brow,col=bcol,&
            block=jp_block_a,found=den_found_a)
       CALL cp_dbcsr_get_block_p(matrix=mat_b,row=brow,col=bcol,&
            block=jp_block_b,found=den_found)
       CALL cp_dbcsr_get_block_p(matrix=mat_c,row=brow,col=bcol,&
            block=jp_block_c,found=den_found)
       IF(do_igaim) CALL cp_dbcsr_get_block_p(matrix=mat_d,row=brow,col=bcol,&
            block=jp_block_d,found=den_found)

       IF (.NOT.ASSOCIATED(jp_block_a)) CYCLE

       IF (distributed_rs_grids) THEN
          NULLIFY (jpblock_a,jpblock_b,jpblock_c,jpblock_d)
          CALL cp_dbcsr_add_block_node ( deltajp_a, brow, bcol, jpblock_a ,error=error)
          jpblock_a = jp_block_a
          CALL cp_dbcsr_add_block_node ( deltajp_b, brow, bcol, jpblock_b ,error=error)
          jpblock_b = jp_block_b
          CALL cp_dbcsr_add_block_node ( deltajp_c, brow, bcol, jpblock_c ,error=error)
          jpblock_c = jp_block_c
          IF(do_igaim) CALL cp_dbcsr_add_block_node ( deltajp_d, brow, bcol, jpblock_d ,error=error)
          IF(do_igaim) jpblock_d = jp_block_d
       ELSE
          jpblock_a => jp_block_a
          jpblock_b => jp_block_b
          jpblock_c => jp_block_c
          IF(do_igaim) jpblock_d => jp_block_d
       ENDIF

       IF (.NOT. map_consistent) THEN
          IF ( ALL ( 100.0_dp*ABS(jpblock_a) < eps_rho_rspace ) .AND. &
               ALL ( 100.0_dp*ABS(jpblock_b) < eps_rho_rspace ) .AND. &
               ALL ( 100.0_dp*ABS(jpblock_c) < eps_rho_rspace ) ) THEN
             CYCLE
          END IF
       END IF

       CALL task_list_inner_loop(tasks, dist_ab, ntasks, curr_tasks, rs_descs,dft_control,cube_info,gridlevel_info,cell,&
                                 iatom,jatom,rpgfa,rpgfb,zeta,zetb,kind_radius_b,set_radius_a,set_radius_b,ra,rab,&
                                 la_max,la_min,lb_max,lb_min,npgfa,npgfb,maxpgf,maxset,natom,nseta,nsetb,error)

    END DO
    CALL neighbor_list_iterator_release(nl_iterator)

    DEALLOCATE (basis_set_list,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    IF (distributed_rs_grids) THEN
       CALL cp_dbcsr_finalize(deltajp_a, error=error)
       CALL cp_dbcsr_finalize(deltajp_b, error=error)
       CALL cp_dbcsr_finalize(deltajp_c, error=error)
       IF(do_igaim) CALL cp_dbcsr_finalize(deltajp_d, error=error)
    ENDIF

    ! sorts / redistributes the task list
    CALL distribute_tasks ( rs_descs, ntasks, natom,&
            maxset, maxpgf, tasks, dist_ab, atom_pair_send, atom_pair_recv,&
            symmetric=.FALSE., reorder_rs_grid_ranks=.TRUE., &
            skip_load_balance_distributed=.FALSE., error=error)

    ALLOCATE(rs_current(gridlevel_info%ngrid_levels))
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    DO igrid_level=1,gridlevel_info%ngrid_levels
       ! Here we need to reallocate the distributed rs_grids, which may have been reordered
       ! by distribute_tasks
       IF (rs_descs(igrid_level)%rs_desc%distributed) THEN
         CALL rs_grid_release(rs_rho(igrid_level)%rs_grid, error=error)
         NULLIFY (rs_rho(igrid_level)%rs_grid)
         CALL rs_grid_create(rs_rho(igrid_level)%rs_grid, rs_descs(igrid_level)%rs_desc, error=error)
       ELSE
         CALL rs_grid_retain(rs_rho(igrid_level)%rs_grid, error=error)
       ENDIF
       CALL rs_grid_zero(rs_rho(igrid_level)%rs_grid)
       CALL rs_grid_create(rs_current(igrid_level)%rs_grid, rs_descs(igrid_level)%rs_desc,&
                           error=error)
       CALL rs_grid_zero(rs_current(igrid_level)%rs_grid)
    ENDDO

    ngpts = 0
    IF ( nthread > 1 ) THEN
       DO igrid_level = 1,gridlevel_info%ngrid_levels
          ngpts = MAX(ngpts,rs_rho(igrid_level)%rs_grid%ngpts_local)
       END DO
       ngpts = ngpts*nthread
       CALL reallocate(lgrid%r,1,ngpts)
    END IF

    !
    ! we need to build the gauge here
    IF(.NOT.current_env%gauge_init) THEN
       CALL current_set_gauge(current_env, qs_env, error=error)
       current_env%gauge_init = .TRUE.
    ENDIF
    !
    ! for any case double check the bounds !
    DO igrid_level=1,gridlevel_info%ngrid_levels
       my_rho=>rs_rho(igrid_level)%rs_grid%r
       my_current=>rs_current(igrid_level)%rs_grid%r
       IF(LBOUND( my_rho, 3 ).NE.LBOUND( my_current, 3 ).OR.&
          LBOUND( my_rho, 2 ).NE.LBOUND( my_current, 2 ).OR.&
          LBOUND( my_rho, 1 ).NE.LBOUND( my_current, 1 ).OR.&
          UBOUND( my_rho, 3 ).NE.UBOUND( my_current, 3 ).OR.&
          UBOUND( my_rho, 2 ).NE.UBOUND( my_current, 2 ).OR.&
          UBOUND( my_rho, 1 ).NE.UBOUND( my_current, 1 ) )THEN
          WRITE(*,*) 'LBOUND(my_rho,3),LBOUND(my_current,3)',LBOUND(my_rho,3),LBOUND(my_current,3)
          WRITE(*,*) 'LBOUND(my_rho,2),LBOUND(my_current,2)',LBOUND(my_rho,2),LBOUND(my_current,2)
          WRITE(*,*) 'LBOUND(my_rho,1),LBOUND(my_current,1)',LBOUND(my_rho,1),LBOUND(my_current,1)
          WRITE(*,*) 'UBOUND(my_rho,3),UBOUND(my_current,3)',UBOUND(my_rho,3),UBOUND(my_current,3)
          WRITE(*,*) 'UBOUND(my_rho,2),UBOUND(my_current,2)',UBOUND(my_rho,2),UBOUND(my_current,2)
          WRITE(*,*) 'UBOUND(my_rho,1),UBOUND(my_current,1)',UBOUND(my_rho,1),UBOUND(my_current,1)
          CALL stop_program(routineN,moduleN,__LINE__,"Bug")
       ENDIF

       IF(do_igaim) THEN
          my_gauge=>rs_gauge(1)%rs(igrid_level)%rs_grid%r
          IF(LBOUND( my_rho, 3 ).NE.LBOUND( my_gauge, 3 ).OR.&
             LBOUND( my_rho, 2 ).NE.LBOUND( my_gauge, 2 ).OR.&
             LBOUND( my_rho, 1 ).NE.LBOUND( my_gauge, 1 ).OR.&
             UBOUND( my_rho, 3 ).NE.UBOUND( my_gauge, 3 ).OR.&
             UBOUND( my_rho, 2 ).NE.UBOUND( my_gauge, 2 ).OR.&
             UBOUND( my_rho, 1 ).NE.UBOUND( my_gauge, 1 )) THEN
             WRITE(*,*) 'LBOUND(my_rho,3),LBOUND(my_gauge,3)',LBOUND(my_rho,3),LBOUND(my_gauge,3)
             WRITE(*,*) 'LBOUND(my_rho,2),LBOUND(my_gauge,2)',LBOUND(my_rho,2),LBOUND(my_gauge,2)
             WRITE(*,*) 'LBOUND(my_rho,1),LBOUND(my_gauge,1)',LBOUND(my_rho,1),LBOUND(my_gauge,1)
             WRITE(*,*) 'UBOUND(my_rho,3),UbOUND(my_gauge,3)',UBOUND(my_rho,3),UBOUND(my_gauge,3)
             WRITE(*,*) 'UBOUND(my_rho,2),UBOUND(my_gauge,2)',UBOUND(my_rho,2),UBOUND(my_gauge,2)
             WRITE(*,*) 'UBOUND(my_rho,1),UBOUND(my_gauge,1)',UBOUND(my_rho,1),UBOUND(my_gauge,1)
             CALL stop_program(routineN,moduleN,__LINE__,"Bug")
          ENDIF
       ENDIF
    ENDDO
    !
    !-------------------------------------------------------------

    IF (distributed_rs_grids) THEN
        CALL rs_distribute_matrix (rs_descs, deltajp_a, atom_pair_send, atom_pair_recv, &
             natom, scatter=.TRUE., error=error)
        CALL rs_distribute_matrix (rs_descs, deltajp_b, atom_pair_send, atom_pair_recv, &
             natom, scatter=.TRUE., error=error)
        CALL rs_distribute_matrix (rs_descs, deltajp_c, atom_pair_send, atom_pair_recv, &
             natom, scatter=.TRUE., error=error)
        IF(do_igaim) CALL rs_distribute_matrix (rs_descs, deltajp_d, atom_pair_send, atom_pair_recv, &
             natom, scatter=.TRUE., error=error)
    ENDIF

    ithread = 0
    jpab_a => jpabt_a(:,:,ithread)
    jpab_b => jpabt_b(:,:,ithread)
    jpab_c => jpabt_c(:,:,ithread)
    IF(do_igaim) jpab_d => jpabt_d(:,:,ithread)
    work => workt(:,:,ithread)

    iatom_old = -1 ; jatom_old = -1 ; iset_old = -1 ; jset_old = -1
    ikind_old = -1 ; jkind_old = -1

    loop_tasks: DO itask = 1,ntasks

       CALL int2pair(tasks(3,itask),igrid_level,iatom,jatom,iset,jset,ipgf,jpgf,natom,maxset,maxpgf)

       ! apparently generalised collocation not implemented correctly yet
       CPPostcondition(tasks(4,itask).NE.2,cp_failure_level,routineP,error,failure)

       IF (iatom .NE. iatom_old .OR. jatom .NE. jatom_old) THEN

          ikind = particle_set(iatom)%atomic_kind%kind_number
          jkind = particle_set(jatom)%atomic_kind%kind_number

          IF(iatom .NE. iatom_old ) ra(:) = pbc(particle_set(iatom)%r,cell)

          brow = iatom
          bcol = jatom

          IF (ikind .NE. ikind_old ) THEN
             CALL get_qs_kind(qs_kind_set(ikind),&
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

             CALL get_qs_kind(qs_kind_set(jkind),&
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

          CALL cp_dbcsr_get_block_p(matrix=deltajp_a,row=brow,col=bcol,&
               block=jp_block_a,found=den_found)
          CALL cp_dbcsr_get_block_p(matrix=deltajp_b,row=brow,col=bcol,&
               block=jp_block_b,found=den_found)
          CALL cp_dbcsr_get_block_p(matrix=deltajp_c,row=brow,col=bcol,&
               block=jp_block_c,found=den_found)
          IF(do_igaim) CALL cp_dbcsr_get_block_p(matrix=deltajp_d,row=brow,col=bcol,&
               block=jp_block_d,found=den_found)

          IF (.NOT.ASSOCIATED(jp_block_a)) &
               CALL stop_program(routineN,moduleN,__LINE__,&
                                 "p_block not associated in deltap")

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
               jp_block_a(sgfa,sgfb),SIZE(jp_block_a,1),&
               0.0_dp,work(1,1),maxco)
          CALL dgemm("N","T",ncoa,ncob,nsgfb(jset),&
               1.0_dp,work(1,1),maxco,&
               sphi_b(1,sgfb),SIZE(sphi_b,1),&
               0.0_dp,jpab_a(1,1),maxco)

          CALL dgemm("N","N",ncoa,nsgfb(jset),nsgfa(iset),&
               1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
               jp_block_b(sgfa,sgfb),SIZE(jp_block_b,1),&
               0.0_dp,work(1,1),maxco)
          CALL dgemm("N","T",ncoa,ncob,nsgfb(jset),&
               1.0_dp,work(1,1),maxco,&
               sphi_b(1,sgfb),SIZE(sphi_b,1),&
               0.0_dp,jpab_b(1,1),maxco)

          CALL dgemm("N","N",ncoa,nsgfb(jset),nsgfa(iset),&
               1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
               jp_block_c(sgfa,sgfb),SIZE(jp_block_c,1),&
               0.0_dp,work(1,1),maxco)
          CALL dgemm("N","T",ncoa,ncob,nsgfb(jset),&
               1.0_dp,work(1,1),maxco,&
               sphi_b(1,sgfb),SIZE(sphi_b,1),&
               0.0_dp,jpab_c(1,1),maxco)

          IF(do_igaim) THEN
             CALL dgemm("N","N",ncoa,nsgfb(jset),nsgfa(iset),&
                  1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                  jp_block_d(sgfa,sgfb),SIZE(jp_block_d,1),&
                  0.0_dp,work(1,1),maxco)
             CALL dgemm("N","T",ncoa,ncob,nsgfb(jset),&
                  1.0_dp,work(1,1),maxco,&
                  sphi_b(1,sgfb),SIZE(sphi_b,1),&
                  0.0_dp,jpab_d(1,1),maxco)
          ENDIF

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

       ! Four calls to the general collocate density, to multply the correct function
       ! to each density matrix
       IF(do_igaim) THEN
          !
          ! here the decontracted mat_jp_{ab} is multiplied by
          !     f_{ab} = g_{a} (dg_{b}/dr)_{idir} - (dg_{a}/dr)_{idir} g_{b}
          scale = 1.0_dp
          CALL collocate_pgf_product_rspace(la_max(iset),zeta(ipgf,iset),&
               la_min(iset),lb_max(jset),zetb(jpgf,jset),lb_min(jset),&
               ra,rab,rab2,scale,jpab_a,na1-1,nb1-1,&
               rs_current(igrid_level)%rs_grid,cell,cube_info(igrid_level),&
               eps_rho_rspace,&
               ga_gb_function=FUNC_ADBmDAB,&
               idir=idir,&
               map_consistent=map_consistent,error=error)
          ! here the decontracted mat_jb_{ab} is multiplied by
          !     f_{ab} = g_{a} * g_{b} ! THIS GOES OUTSIDE THE LOOP !
          IF(scale2.NE.0.0_dp) THEN
             CALL collocate_pgf_product_rspace(la_max(iset),zeta(ipgf,iset),&
                  la_min(iset),lb_max(jset),zetb(jpgf,jset),lb_min(jset),&
                  ra,rab,rab2,scale2,jpab_d,na1-1,nb1-1,&
                  rs_rho(igrid_level)%rs_grid,cell,cube_info(igrid_level),&
                  eps_rho_rspace,&
                  ga_gb_function=FUNC_AB,&
                  map_consistent=map_consistent,error=error)
          ENDIF!rm
          ! here the decontracted mat_jp_rii{ab} is multiplied by
          !     f_{ab} = g_{a} (d(r) - R_{b})_{iiB} (dg_{b}/dr)_{idir} -
          !             (dg_{a}/dr)_{idir} (d(r) - R_{b})_{iiB} g_{b}
          scale = 1.0_dp
          CALL collocate_pgf_product_rspace(la_max(iset),zeta(ipgf,iset),&
               la_min(iset),lb_max(jset),zetb(jpgf,jset),lb_min(jset),&
               ra,rab,rab2,scale,jpab_b,na1-1,nb1-1,&
               rs_current(igrid_level)%rs_grid,cell,cube_info(igrid_level),&
               eps_rho_rspace,&
               ga_gb_function=FUNC_ADBmDAB,&
               idir=idir,ir=iiiB,&
               rsgauge=rs_gauge(iiiB)%rs(igrid_level)%rs_grid,&
               rsbuf=current_env%rs_buf(igrid_level)%rs_grid,&
               map_consistent=map_consistent,error=error)
          ! here the decontracted mat_jp_riii{ab} is multiplied by
          !     f_{ab} = -g_{a} (d(r) - R_{b})_{iiB} (dg_{b}/dr)_{idir} +
          !             (dg_{a}/dr)_{idir} (d(r) - R_{b})_{iiB} g_{b}
          scale = -1.0_dp
          CALL collocate_pgf_product_rspace(la_max(iset),zeta(ipgf,iset),&
               la_min(iset),lb_max(jset),zetb(jpgf,jset),lb_min(jset),&
               ra,rab,rab2,scale,jpab_c,na1-1,nb1-1,&
               rs_current(igrid_level)%rs_grid,cell,cube_info(igrid_level),&
               eps_rho_rspace,&
               ga_gb_function=FUNC_ADBmDAB,&
               idir=idir,ir=iiB,&
               rsgauge=rs_gauge(iiB)%rs(igrid_level)%rs_grid,&
               rsbuf=current_env%rs_buf(igrid_level)%rs_grid,&
               map_consistent=map_consistent,error=error)
       ELSE
          ! here the decontracted mat_jp_{ab} is multiplied by
          !     f_{ab} = g_{a} (dg_{b}/dr)_{idir} - (dg_{a}/dr)_{idir} g_{b}
          scale = 1.0_dp
          CALL collocate_pgf_product_rspace(la_max(iset),zeta(ipgf,iset),&
               la_min(iset),lb_max(jset),zetb(jpgf,jset),lb_min(jset),&
               ra,rab,rab2,scale,jpab_a,na1-1,nb1-1,&
               rs_current(igrid_level)%rs_grid,cell,cube_info(igrid_level),&
               eps_rho_rspace,&
               ga_gb_function=FUNC_ADBmDAB,&
               idir=idir,&
               map_consistent=map_consistent,error=error)
          ! here the decontracted mat_jp_rii{ab} is multiplied by
          !     f_{ab} = g_{a} (r - R_{b})_{iiB} (dg_{b}/dr)_{idir} -
          !             (dg_{a}/dr)_{idir} (r - R_{b})_{iiB} g_{b}
          scale = 1.0_dp
          CALL collocate_pgf_product_rspace(la_max(iset),zeta(ipgf,iset),&
               la_min(iset),lb_max(jset),zetb(jpgf,jset),lb_min(jset),&
               ra,rab,rab2,scale,jpab_b,na1-1,nb1-1,&
               rs_current(igrid_level)%rs_grid,cell,cube_info(igrid_level),&
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
               ra,rab,rab2,scale,jpab_c,na1-1,nb1-1,&
               rs_current(igrid_level)%rs_grid,cell,cube_info(igrid_level),&
               eps_rho_rspace,&
               ga_gb_function=FUNC_ARDBmDARB,&
               idir=idir,ir=iiB,&
               map_consistent=map_consistent,error=error)
       ENDIF

    END DO loop_tasks
    !
    ! Scale the density with the gauge rho * ( r - d(r) ) if needed
    IF(do_igaim) THEN
       DO igrid_level=1,gridlevel_info%ngrid_levels
          CALL rs_grid_mult_and_add(rs_current(igrid_level)%rs_grid,rs_rho(igrid_level)%rs_grid,&
                                    rs_gauge(idir2)%rs(igrid_level)%rs_grid,1.0_dp)
       ENDDO
    ENDIF
    !   *** Release work storage ***

    IF (distributed_rs_grids) THEN
       CALL cp_dbcsr_deallocate_matrix ( deltajp_a ,error=error)
       CALL cp_dbcsr_deallocate_matrix ( deltajp_b ,error=error)
       CALL cp_dbcsr_deallocate_matrix ( deltajp_c ,error=error)
       IF(do_igaim) CALL cp_dbcsr_deallocate_matrix ( deltajp_d ,error=error)
    END IF

    IF ( nthread > 1 ) THEN
       DEALLOCATE (lgrid%r,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    END IF

    DEALLOCATE (jpabt_a,jpabt_b,jpabt_c,jpabt_d,workt,tasks,dist_ab,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    IF (distributed_rs_grids) THEN
       DEALLOCATE(atom_pair_send,atom_pair_recv,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ENDIF

    CALL density_rs2pw(pw_env,rs_current,current_rs,current_gs,error=error)

    IF (ASSOCIATED(rs_rho)) THEN
      DO i=1, SIZE(rs_rho)
        CALL rs_grid_release(rs_rho(i)%rs_grid, error=error)
      END DO
    END IF

    ! Free the array of grids (grids themselves are released in density_rs2pw)
    DEALLOCATE (rs_current,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE calculate_jrho_resp


! *****************************************************************************
!> \brief ...
!> \param current_env ...
!> \param qs_env ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE current_set_gauge(current_env,qs_env,error)
    !
    TYPE(current_env_type)                   :: current_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'current_set_gauge', &
      routineP = moduleN//':'//routineN

    REAL(dp)                                 :: dbox(3)
    REAL(dp), ALLOCATABLE, DIMENSION(:,:)    :: box_data
    INTEGER                                  :: handle, igrid_level, nbox(3), gauge, istat
    INTEGER, ALLOCATABLE, DIMENSION(:, :, :) :: box_ptr
    TYPE(realspace_grid_desc_p_type), DIMENSION(:), &
      POINTER                                :: rs_descs
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(realspaces_grid_p_type), DIMENSION(:), POINTER :: rs_gauge

    TYPE(box_type), DIMENSION(:,:,:), POINTER :: box
    LOGICAL :: use_old_gauge_atom, failure

    NULLIFY(rs_gauge,box)

    CALL timeset(routineN,handle)

    CALL get_current_env(current_env=current_env,&
                         use_old_gauge_atom=use_old_gauge_atom,&
                         rs_gauge=rs_gauge,&
                         gauge=gauge,&
                         error=error)

    IF(gauge.EQ.current_gauge_atom) THEN
       CALL get_qs_env(qs_env=qs_env,&
            pw_env=pw_env,&
            error=error)
       CALL pw_env_get(pw_env=pw_env,&
            rs_descs=rs_descs,&
            error=error)
       !
       ! box the atoms
       IF(use_old_gauge_atom) THEN
          CALL box_atoms(qs_env,error)
       ELSE
          CALL  box_atoms_new(current_env,qs_env,box,error)
       ENDIF
       !
       ! allocate and build the gauge
       ALLOCATE (rs_gauge(1)%rs(pw_env%gridlevel_info%ngrid_levels),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE (rs_gauge(2)%rs(pw_env%gridlevel_info%ngrid_levels),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE (rs_gauge(3)%rs(pw_env%gridlevel_info%ngrid_levels),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       DO igrid_level = pw_env%gridlevel_info%ngrid_levels,1,-1

          CALL rs_grid_create(rs_gauge(1)%rs(igrid_level)%rs_grid, rs_descs(igrid_level)%rs_desc, error=error)
          CALL rs_grid_create(rs_gauge(2)%rs(igrid_level)%rs_grid, rs_descs(igrid_level)%rs_desc, error=error)
          CALL rs_grid_create(rs_gauge(3)%rs(igrid_level)%rs_grid, rs_descs(igrid_level)%rs_desc, error=error)

          IF(use_old_gauge_atom) THEN
             CALL collocate_gauge(current_env,qs_env,&
                                  rs_gauge(1)%rs(igrid_level)%rs_grid,&
                                  rs_gauge(2)%rs(igrid_level)%rs_grid,&
                                  rs_gauge(3)%rs(igrid_level)%rs_grid,&
                                  error)
          ELSE
             CALL collocate_gauge_new(current_env,qs_env,&
                                      rs_gauge(1)%rs(igrid_level)%rs_grid,&
                                      rs_gauge(2)%rs(igrid_level)%rs_grid,&
                                      rs_gauge(3)%rs(igrid_level)%rs_grid,&
                                      box,error)
          ENDIF
       ENDDO
       !
       ! allocate the buf
       ALLOCATE (current_env%rs_buf(pw_env%gridlevel_info%ngrid_levels),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       DO igrid_level = 1,pw_env%gridlevel_info%ngrid_levels
         CALL rs_grid_create(current_env%rs_buf(igrid_level)%rs_grid, rs_descs(igrid_level)%rs_desc, error=error)
       END DO
       !
       DEALLOCATE(box_ptr,box_data)
       CALL deallocate_box(box)
    ENDIF

    CALL timestop(handle)

  CONTAINS

! *****************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param error ...
! *****************************************************************************
    SUBROUTINE box_atoms(qs_env,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    REAL(kind=dp), PARAMETER                 :: box_size_guess = 5.0_dp

    INTEGER                                  :: i, iatom, ibox, ii, jbox, &
                                                kbox, natms
    REAL(dp)                                 :: offset(3)
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: ratom
    TYPE(cell_type), POINTER                 :: cell
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set

      CALL get_qs_env(qs_env=qs_env,&
           qs_kind_set=qs_kind_set,&
           cell=cell,&
           particle_set=particle_set,&
           error=error)

      natms = SIZE(particle_set,1)
      ALLOCATE(ratom(3,natms))
      !
      ! box the atoms
      nbox(1) = CEILING(cell%hmat(1,1) / box_size_guess)
      nbox(2) = CEILING(cell%hmat(2,2) / box_size_guess)
      nbox(3) = CEILING(cell%hmat(3,3) / box_size_guess)
      !write(*,*) 'nbox',nbox
      dbox(1) = cell%hmat(1,1) / REAL(nbox(1),dp)
      dbox(2) = cell%hmat(2,2) / REAL(nbox(2),dp)
      dbox(3) = cell%hmat(3,3) / REAL(nbox(3),dp)
      !write(*,*) 'dbox',dbox
      ALLOCATE(box_ptr(0:nbox(1),0:nbox(2)-1,0:nbox(3)-1),box_data(3,natms))
      box_data(:,:) = HUGE(0.0_dp)
      box_ptr(:,:,:) = HUGE(0)
      !
      offset(1) = cell%hmat(1,1)*0.5_dp
      offset(2) = cell%hmat(2,2)*0.5_dp
      offset(3) = cell%hmat(3,3)*0.5_dp
      DO iatom = 1,natms
         ratom(:,iatom) = pbc(particle_set(iatom)%r(:),cell)+offset(:)
      ENDDO
      !
      i = 1
      DO kbox = 0,nbox(3)-1
      DO jbox = 0,nbox(2)-1
         box_ptr(0,jbox,kbox) = i
      DO ibox = 0,nbox(1)-1
         ii = 0
         DO iatom = 1,natms
            IF(INT(ratom(1,iatom)/dbox(1)).EQ.ibox.AND.&
               INT(ratom(2,iatom)/dbox(2)).EQ.jbox.AND.&
               INT(ratom(3,iatom)/dbox(3)).EQ.kbox) THEN
               box_data(:,i) = ratom(:,iatom)-offset(:)
               i = i + 1
               ii = ii + 1
            ENDIF
         ENDDO
         box_ptr(ibox+1,jbox,kbox) = box_ptr(ibox,jbox,kbox) + ii
      ENDDO
      ENDDO
      ENDDO
      !
      IF(.FALSE.) THEN
         DO kbox = 0,nbox(3)-1
         DO jbox = 0,nbox(2)-1
         DO ibox = 0,nbox(1)-1
            WRITE(*,*) 'box=',ibox,jbox,kbox
            WRITE(*,*) 'nbr atom=',box_ptr(ibox+1,jbox,kbox)-box_ptr(ibox,jbox,kbox)
            DO iatom = box_ptr(ibox,jbox,kbox),box_ptr(ibox+1,jbox,kbox)-1
               WRITE(*,*) 'iatom=',iatom
               WRITE(*,'(A,3E14.6)') 'coor=',box_data(:,iatom)
            ENDDO
         ENDDO
         ENDDO
         ENDDO
      ENDIF
      DEALLOCATE(ratom)
    END SUBROUTINE box_atoms


! *****************************************************************************
!> \brief ...
!> \param current_env ...
!> \param qs_env ...
!> \param rs_grid_x ...
!> \param rs_grid_y ...
!> \param rs_grid_z ...
!> \param error ...
! *****************************************************************************
    SUBROUTINE collocate_gauge(current_env,qs_env,rs_grid_x,rs_grid_y,rs_grid_z,error)
      !
    TYPE(current_env_type)                   :: current_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(realspace_grid_type), POINTER       :: rs_grid_x, rs_grid_y, &
                                                rs_grid_z
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    INTEGER :: i, iatom, ibeg, ibox, iend, imax, imin, j, jatom, jbox, jmax, &
      jmin, k, kbox, kmax, kmin, lb(3), lb_local(3), natms, natms_local, ng(3)
    REAL(KIND=dp) :: ab, buf_tmp, dist, dr(3), gauge_atom_radius, offset(3), &
      pa, pb, point(3), pra(3), r(3), res(3), summe, tmp, x, y, z
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:) :: buf, nrm_atms_pnt
    REAL(kind=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: atms_pnt, ratom
    REAL(kind=dp), DIMENSION(:, :, :), &
      POINTER                                :: grid_x, grid_y, grid_z
    TYPE(cell_type), POINTER                 :: cell
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set

!

      CALL get_current_env(current_env=current_env,&
           gauge_atom_radius=gauge_atom_radius,&
           error=error)
      !
      CALL get_qs_env(qs_env=qs_env,&
           qs_kind_set=qs_kind_set,&
           cell=cell,&
           particle_set=particle_set,&
           error=error)
      !
      natms = SIZE(particle_set,1)
      dr(1) = rs_grid_x%desc%dh(1,1)
      dr(2) = rs_grid_x%desc%dh(2,2)
      dr(3) = rs_grid_x%desc%dh(3,3)
      lb(:) = rs_grid_x%desc%lb(:)
      lb_local(:) = rs_grid_x%lb_local(:)
      grid_x => rs_grid_x%r(:,:,:)
      grid_y => rs_grid_y%r(:,:,:)
      grid_z => rs_grid_z%r(:,:,:)
      ng(:) = UBOUND(grid_x)
      offset(1) = cell%hmat(1,1)*0.5_dp
      offset(2) = cell%hmat(2,2)*0.5_dp
      offset(3) = cell%hmat(3,3)*0.5_dp
      ALLOCATE(buf(natms),ratom(3,natms),atms_pnt(3,natms),nrm_atms_pnt(natms))
      !
      ! go over the grid
      DO k = 1,ng(3)
         DO j = 1,ng(2)
            DO i = 1,ng(1)
               !
               point(3) = REAL(k-1+lb_local(3)-lb(3),dp)*dr(3)
               point(2) = REAL(j-1+lb_local(2)-lb(2),dp)*dr(2)
               point(1) = REAL(i-1+lb_local(1)-lb(1),dp)*dr(1)
               point = pbc(point,cell)
               !
               ! run over the overlaping boxes
               natms_local = 0
               kmin = INT((point(3)+offset(3)-gauge_atom_radius)/dbox(3))
               kmax = INT((point(3)+offset(3)+gauge_atom_radius)/dbox(3))
               IF(kmax-kmin+1.GT.nbox(3)) THEN
                  kmin = 0
                  kmax = nbox(3)-1
               ENDIF
               DO kbox = kmin,kmax
                  jmin = INT((point(2)+offset(2)-gauge_atom_radius)/dbox(2))
                  jmax = INT((point(2)+offset(2)+gauge_atom_radius)/dbox(2))
                  IF(jmax-jmin+1.GT.nbox(2)) THEN
                     jmin = 0
                     jmax = nbox(2)-1
                  ENDIF
               DO jbox = jmin,jmax
                  imin = INT((point(1)+offset(1)-gauge_atom_radius)/dbox(1))
                  imax = INT((point(1)+offset(1)+gauge_atom_radius)/dbox(1))
                  IF(imax-imin+1.GT.nbox(1)) THEN
                     imin = 0
                     imax = nbox(1)-1
                  ENDIF
               DO ibox = imin,imax
                  ibeg = box_ptr(MODULO(ibox,nbox(1))  ,MODULO(jbox,nbox(2)),MODULO(kbox,nbox(3)))
                  iend = box_ptr(MODULO(ibox,nbox(1))+1,MODULO(jbox,nbox(2)),MODULO(kbox,nbox(3)))-1
                  DO iatom = ibeg,iend
                     r(:) = pbc(box_data(:,iatom)-point(:),cell) + point(:)
                     dist = (r(1)-point(1))**2 + (r(2)-point(2))**2 + (r(3)-point(3))**2
                     IF(dist.LT.gauge_atom_radius**2) THEN
                        natms_local = natms_local + 1
                        ratom(:,natms_local) = r(:)
                        !
                        ! compute the distance atoms-point
                        x = point(1) - r(1)
                        y = point(2) - r(2)
                        z = point(3) - r(3)
                        atms_pnt(1,natms_local) = x
                        atms_pnt(2,natms_local) = y
                        atms_pnt(3,natms_local) = z
                        nrm_atms_pnt(natms_local) = SQRT( x*x + y*y + z*z )
                     ENDIF
                  ENDDO
               ENDDO
               ENDDO
               ENDDO
               !
               IF(natms_local.GT.0) THEN
                  !
                  !
                  DO iatom = 1,natms_local
                     buf_tmp = 1.0_dp
                     pra(1) = atms_pnt(1,iatom)
                     pra(2) = atms_pnt(2,iatom)
                     pra(3) = atms_pnt(3,iatom)
                     pa     = nrm_atms_pnt(iatom)
                     DO jatom = 1,natms_local
                        IF(iatom.EQ.jatom) CYCLE
                        pb = nrm_atms_pnt(jatom)
                        x  = pra(1) - atms_pnt(1,jatom)
                        y  = pra(2) - atms_pnt(2,jatom)
                        z  = pra(3) - atms_pnt(3,jatom)
                        ab  = SQRT( x*x + y*y + z*z )
                        !
                        tmp = ( pa - pb ) / ab
                        tmp = 0.5_dp * ( 3.0_dp - tmp * tmp ) * tmp
                        tmp = 0.5_dp * ( 3.0_dp - tmp * tmp ) * tmp
                        tmp = 0.5_dp * ( 3.0_dp - tmp * tmp ) * tmp
                        buf_tmp = buf_tmp * 0.5_dp * ( 1.0_dp - tmp )
                     ENDDO
                     buf(iatom) = buf_tmp
                  ENDDO
                  res(1) = 0.0_dp
                  res(2) = 0.0_dp
                  res(3) = 0.0_dp
                  summe  = 0.0_dp
                  DO iatom = 1,natms_local
                     res(1) = res(1) + ratom(1,iatom) * buf(iatom)
                     res(2) = res(2) + ratom(2,iatom) * buf(iatom)
                     res(3) = res(3) + ratom(3,iatom) * buf(iatom)
                     summe  = summe + buf(iatom)
                  ENDDO
                  res(1) = res(1) / summe
                  res(2) = res(2) / summe
                  res(3) = res(3) / summe
                  grid_x(i,j,k) = point(1) - res(1)
                  grid_y(i,j,k) = point(2) - res(2)
                  grid_z(i,j,k) = point(3) - res(3)
               ELSE
                  grid_x(i,j,k) = 0.0_dp
                  grid_y(i,j,k) = 0.0_dp
                  grid_z(i,j,k) = 0.0_dp
               ENDIF
            ENDDO
         ENDDO
      ENDDO

      DEALLOCATE(buf,ratom,atms_pnt,nrm_atms_pnt)

    END SUBROUTINE collocate_gauge

! *****************************************************************************
!> \brief ...
!> \param current_env ...
!> \param qs_env ...
!> \param box ...
!> \param error ...
! *****************************************************************************
    SUBROUTINE box_atoms_new(current_env,qs_env,box,error)
    TYPE(current_env_type)                   :: current_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(box_type), DIMENSION(:, :, :), &
      POINTER                                :: box
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'box_atoms_new', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, i, iatom, ibeg, ibox, iend, ifind, ii, imax, imin, j, &
      jatom, jbox, jmax, jmin, k, kbox, kmax, kmin, natms, natms_local
    REAL(dp)                                 :: gauge_atom_radius, offset(3), &
                                                scale
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: ratom
    REAL(dp), DIMENSION(:, :), POINTER       :: r_ptr
    REAL(kind=dp)                            :: box_center(3), &
                                                box_center_wrap(3), &
                                                box_size_guess, r(3)
    TYPE(cell_type), POINTER                 :: cell
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set

      CALL timeset(routineN,handle)

      CALL get_qs_env(qs_env=qs_env,&
           qs_kind_set=qs_kind_set,&
           cell=cell,&
           particle_set=particle_set,&
           error=error)

      CALL get_current_env(current_env=current_env,&
           gauge_atom_radius=gauge_atom_radius,&
           error=error)

      scale = 2.0_dp

      box_size_guess = gauge_atom_radius/scale

      natms = SIZE(particle_set,1)
      ALLOCATE(ratom(3,natms))

      !
      ! box the atoms
      nbox(1) = CEILING(cell%hmat(1,1) / box_size_guess)
      nbox(2) = CEILING(cell%hmat(2,2) / box_size_guess)
      nbox(3) = CEILING(cell%hmat(3,3) / box_size_guess)
      dbox(1) = cell%hmat(1,1) / REAL(nbox(1),dp)
      dbox(2) = cell%hmat(2,2) / REAL(nbox(2),dp)
      dbox(3) = cell%hmat(3,3) / REAL(nbox(3),dp)
      ALLOCATE(box_ptr(0:nbox(1),0:nbox(2)-1,0:nbox(3)-1),box_data(3,natms))
      box_data(:,:)  = HUGE(0.0_dp)
      box_ptr(:,:,:) = HUGE(0)
      !
      offset(1) = cell%hmat(1,1)*0.5_dp
      offset(2) = cell%hmat(2,2)*0.5_dp
      offset(3) = cell%hmat(3,3)*0.5_dp
      DO iatom = 1,natms
         ratom(:,iatom) = pbc(particle_set(iatom)%r(:),cell)
      ENDDO
      !
      i = 1
      DO kbox = 0,nbox(3)-1
      DO jbox = 0,nbox(2)-1
         box_ptr(0,jbox,kbox) = i
      DO ibox = 0,nbox(1)-1
         ii = 0
         DO iatom = 1,natms
            IF(MODULO(FLOOR(ratom(1,iatom)/dbox(1)),nbox(1)).EQ.ibox.AND.&
               MODULO(FLOOR(ratom(2,iatom)/dbox(2)),nbox(2)).EQ.jbox.AND.&
               MODULO(FLOOR(ratom(3,iatom)/dbox(3)),nbox(3)).EQ.kbox) THEN
               box_data(:,i) = ratom(:,iatom)
               i = i + 1
               ii = ii + 1
            ENDIF
         ENDDO
         box_ptr(ibox+1,jbox,kbox) = box_ptr(ibox,jbox,kbox) + ii
      ENDDO
      ENDDO
      ENDDO
      !
      IF(.FALSE.) THEN
         DO kbox = 0,nbox(3)-1
         DO jbox = 0,nbox(2)-1
         DO ibox = 0,nbox(1)-1
            IF(box_ptr(ibox+1,jbox,kbox)-box_ptr(ibox,jbox,kbox).GT.0) THEN
               WRITE(*,*) 'box=',ibox,jbox,kbox
               WRITE(*,*) 'nbr atom=',box_ptr(ibox+1,jbox,kbox)-box_ptr(ibox,jbox,kbox)
               DO iatom = box_ptr(ibox,jbox,kbox),box_ptr(ibox+1,jbox,kbox)-1
                  WRITE(*,'(A,I3,3E14.6)') 'coor=',iatom,box_data(:,iatom)
               ENDDO
            ENDIF
         ENDDO
         ENDDO
         ENDDO
      ENDIF
      !
      NULLIFY(box)
      ALLOCATE(box(0:nbox(1)-1,0:nbox(2)-1,0:nbox(3)-1))
      !
      ! build the list
      DO k = 0,nbox(3)-1
      DO j = 0,nbox(2)-1
      DO i = 0,nbox(1)-1
         !
         box_center(1) = ( REAL(i,dp) + 0.5_dp ) * dbox(1)
         box_center(2) = ( REAL(j,dp) + 0.5_dp ) * dbox(2)
         box_center(3) = ( REAL(k,dp) + 0.5_dp ) * dbox(3)
         box_center_wrap = pbc(box_center,cell)
         !
         ! find the atoms that are in the overlaping boxes
         natms_local = 0
         kmin = FLOOR( ( box_center(3) - gauge_atom_radius ) / dbox(3) )
         kmax = FLOOR( ( box_center(3) + gauge_atom_radius ) / dbox(3) )
         IF(kmax-kmin+1.GT.nbox(3)) THEN
            kmin = 0
            kmax = nbox(3)-1
         ENDIF
         DO kbox = kmin,kmax
            jmin = FLOOR( ( box_center(2) - gauge_atom_radius ) / dbox(2) )
            jmax = FLOOR( ( box_center(2) + gauge_atom_radius ) / dbox(2) )
            IF(jmax-jmin+1.GT.nbox(2)) THEN
               jmin = 0
               jmax = nbox(2)-1
            ENDIF
         DO jbox = jmin,jmax
            imin = FLOOR( ( box_center(1) - gauge_atom_radius ) / dbox(1) )
            imax = FLOOR( ( box_center(1) + gauge_atom_radius ) / dbox(1) )
            IF(imax-imin+1.GT.nbox(1)) THEN
               imin = 0
               imax = nbox(1)-1
            ENDIF
         DO ibox = imin,imax
            ibeg = box_ptr(MODULO(ibox,nbox(1))  ,MODULO(jbox,nbox(2)),MODULO(kbox,nbox(3)))
            iend = box_ptr(MODULO(ibox,nbox(1))+1,MODULO(jbox,nbox(2)),MODULO(kbox,nbox(3)))-1
            DO iatom = ibeg,iend
               r = pbc(box_center_wrap(:)-box_data(:,iatom),cell)
               IF(ABS(r(1)).LE.(scale+0.5_dp)*dbox(1).AND.&
                  ABS(r(2)).LE.(scale+0.5_dp)*dbox(2).AND.&
                  ABS(r(3)).LE.(scale+0.5_dp)*dbox(3)) THEN
                  natms_local = natms_local + 1
                  ratom(:,natms_local) = box_data(:,iatom)
               ENDIF
            ENDDO
         ENDDO ! box
         ENDDO
         ENDDO
         !
         ! set the list
         box(i,j,k)%n = natms_local
         NULLIFY(box(i,j,k)%r)
         IF(natms_local.GT.0) THEN
            ALLOCATE(box(i,j,k)%r(3,natms_local))
            r_ptr => box(i,j,k)%r
            CALL dcopy(3*natms_local,ratom(1,1),1,r_ptr(1,1),1)
         ENDIF
      ENDDO ! list
      ENDDO
      ENDDO

      IF(.FALSE.) THEN
         DO k = 0,nbox(3)-1
         DO j = 0,nbox(2)-1
         DO i = 0,nbox(1)-1
            IF(box(i,j,k)%n.GT.0) THEN
               WRITE(*,*)
               WRITE(*,*) 'box=',i,j,k
               box_center(1) = ( REAL(i,dp)+0.5_dp ) * dbox(1)
               box_center(2) = ( REAL(j,dp)+0.5_dp ) * dbox(2)
               box_center(3) = ( REAL(k,dp)+0.5_dp ) * dbox(3)
               box_center = pbc(box_center,cell)
               WRITE(*,'(A,3E14.6)') 'box_center=',box_center
               WRITE(*,*) 'nbr atom=',box(i,j,k)%n
               r_ptr => box(i,j,k)%r
               DO iatom = 1,box(i,j,k)%n
                  WRITE(*,'(A,I3,3E14.6)') 'coor=',iatom,r_ptr(:,iatom)
                  r(:) = pbc(box_center(:)-r_ptr(:,iatom),cell)
                  IF(ABS(r(1)).GT.(scale+0.5_dp)*dbox(1).OR.&
                     ABS(r(2)).GT.(scale+0.5_dp)*dbox(2).OR.&
                     ABS(r(3)).GT.(scale+0.5_dp)*dbox(3)) THEN
                     WRITE(*,*) 'error too many atoms'
                     WRITE(*,*) 'dist=',ABS(r(:))
                     WRITE(*,*) 'large_dist=',(scale+0.5_dp)*dbox
                     STOP
                  ENDIF
               ENDDO
            ENDIF
         ENDDO ! list
         ENDDO
         ENDDO
      ENDIF

      IF(.TRUE.) THEN
         DO k = 0,nbox(3)-1
         DO j = 0,nbox(2)-1
         DO i = 0,nbox(1)-1
            box_center(1) = ( REAL(i,dp)+0.5_dp ) * dbox(1)
            box_center(2) = ( REAL(j,dp)+0.5_dp ) * dbox(2)
            box_center(3) = ( REAL(k,dp)+0.5_dp ) * dbox(3)
            box_center = pbc(box_center,cell)
            r_ptr => box(i,j,k)%r
            DO iatom = 1,natms
               r(:) = pbc(box_center(:)-ratom(:,iatom),cell)
               ifind = 0
               DO jatom = 1,box(i,j,k)%n
                  IF(SUM(ABS(ratom(:,iatom)-r_ptr(:,jatom))).LT.1E-10_dp) ifind = 1
               ENDDO

               IF(ifind.EQ.0) THEN
                  IF(SQRT(DOT_PRODUCT(r,r)).LT.gauge_atom_radius) THEN
                     WRITE(*,*) 'error atom too close'
                     WRITE(*,*) 'iatom',iatom
                     WRITE(*,*) 'box_center',box_center
                     WRITE(*,*) 'ratom',ratom(:,iatom)
                     WRITE(*,*) 'gauge_atom_radius',gauge_atom_radius
                     STOP
                  ENDIF
               ENDIF
            ENDDO
         ENDDO ! list
         ENDDO
         ENDDO
      ENDIF

      DEALLOCATE(ratom)

      CALL timestop(handle)

    END SUBROUTINE box_atoms_new

! *****************************************************************************
!> \brief ...
!> \param current_env ...
!> \param qs_env ...
!> \param rs_grid_x ...
!> \param rs_grid_y ...
!> \param rs_grid_z ...
!> \param box ...
!> \param error ...
! *****************************************************************************
    SUBROUTINE collocate_gauge_new(current_env,qs_env,rs_grid_x,rs_grid_y,rs_grid_z,box,error)
      !
    TYPE(current_env_type)                   :: current_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(realspace_grid_type), POINTER       :: rs_grid_x, rs_grid_y, &
                                                rs_grid_z
    TYPE(box_type), DIMENSION(:, :, :), &
      POINTER                                :: box
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'collocate_gauge_new', &
      routineP = moduleN//':'//routineN

    INTEGER :: delta_lb(3), handle, i, iatom, ib, ibe, ibox, ibs, ie, is, j, &
      jatom, jb, jbe, jbox, jbs, je, js, k, kb, kbe, kbox, kbs, ke, ks, &
      lb(3), lb_local(3), natms, natms_local0, natms_local1, ng(3)
    REAL(dp), DIMENSION(:, :), POINTER       :: r_ptr
    REAL(KIND=dp) :: ab, box_center(3), buf_tmp, dist, dr(3), &
      gauge_atom_radius, offset(3), pa, pb, point(3), pra(3), r(3), res(3), &
      summe, tmp, x, y, z
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:) :: buf, nrm_atms_pnt
    REAL(kind=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: atms_pnt, ratom
    REAL(kind=dp), DIMENSION(:, :, :), &
      POINTER                                :: grid_x, grid_y, grid_z
    TYPE(cell_type), POINTER                 :: cell
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set

      CALL timeset(routineN,handle)

!
      CALL get_current_env(current_env=current_env,&
           gauge_atom_radius=gauge_atom_radius,&
           error=error)
      !
      CALL get_qs_env(qs_env=qs_env,&
           qs_kind_set=qs_kind_set,&
           cell=cell,&
           particle_set=particle_set,&
           error=error)
      !
      natms = SIZE(particle_set,1)
      dr(1) = rs_grid_x%desc%dh(1,1)
      dr(2) = rs_grid_x%desc%dh(2,2)
      dr(3) = rs_grid_x%desc%dh(3,3)
      lb(:) = rs_grid_x%desc%lb(:)
      lb_local(:) = rs_grid_x%lb_local(:)
      grid_x => rs_grid_x%r(:,:,:)
      grid_y => rs_grid_y%r(:,:,:)
      grid_z => rs_grid_z%r(:,:,:)
      ng(:) = UBOUND(grid_x)
      delta_lb(:) = lb_local(:)-lb(:)
      offset(1) = cell%hmat(1,1) * 0.5_dp
      offset(2) = cell%hmat(2,2) * 0.5_dp
      offset(3) = cell%hmat(3,3) * 0.5_dp
      ALLOCATE(buf(natms),ratom(3,natms),atms_pnt(3,natms),nrm_atms_pnt(natms))
      !
      ! find the boxes that match the grid
      ibs = FLOOR( REAL(         delta_lb(1), dp ) * dr(1) / dbox(1) )
      ibe = FLOOR( REAL( ng(1)-1+delta_lb(1), dp ) * dr(1) / dbox(1) )
      jbs = FLOOR( REAL(         delta_lb(2), dp ) * dr(2) / dbox(2) )
      jbe = FLOOR( REAL( ng(2)-1+delta_lb(2), dp ) * dr(2) / dbox(2) )
      kbs = FLOOR( REAL(         delta_lb(3), dp ) * dr(3) / dbox(3) )
      kbe = FLOOR( REAL( ng(3)-1+delta_lb(3), dp ) * dr(3) / dbox(3) )
      !
      ! go over the box-list
      DO kb = kbs,kbe
      DO jb = jbs,jbe
      DO ib = ibs,ibe
         ibox = MODULO(ib,nbox(1))
         jbox = MODULO(jb,nbox(2))
         kbox = MODULO(kb,nbox(3))
         !
         is = MAX(CEILING(REAL(ib  ,dp)*dbox(1)/dr(1)),        delta_lb(1))-delta_lb(1)+1
         ie = MIN(  FLOOR(REAL(ib+1,dp)*dbox(1)/dr(1)),ng(1)-1+delta_lb(1))-delta_lb(1)+1
         js = MAX(CEILING(REAL(jb  ,dp)*dbox(2)/dr(2)),        delta_lb(2))-delta_lb(2)+1
         je = MIN(  FLOOR(REAL(jb+1,dp)*dbox(2)/dr(2)),ng(2)-1+delta_lb(2))-delta_lb(2)+1
         ks = MAX(CEILING(REAL(kb  ,dp)*dbox(3)/dr(3)),        delta_lb(3))-delta_lb(3)+1
         ke = MIN(  FLOOR(REAL(kb+1,dp)*dbox(3)/dr(3)),ng(3)-1+delta_lb(3))-delta_lb(3)+1
         !
         ! sanity checks
         IF(.TRUE.) THEN
            IF(REAL(ks-1+delta_lb(3),dp)*dr(3).LT.REAL(kb,dp)*dbox(3).OR.&
               REAL(ke-1+delta_lb(3),dp)*dr(3).GT.REAL(kb+1,dp)*dbox(3)) THEN
               WRITE(*,*) 'box_k',REAL(kb,dp)*dbox(3),REAL(kb+1,dp)*dbox(3)
               WRITE(*,*) 'point_k',REAL(ks-1+delta_lb(3),dp)*dr(3),REAL(ke-1+delta_lb(3),dp)*dr(3)
               WRITE(*,*) 'ibox',ibox,'jbox',jbox,'kbox',kbox
               WRITE(*,*) 'is,ie',is,ie,' js,je',js,je,' ks,ke',ks,ke
               WRITE(*,*) 'ibs,ibe',ibs,ibe,' jbs,jbe',jbs,jbe,' kbs,kbe',kbs,kbe
               STOP 'we stop_k'
            ENDIF
            IF(REAL(js-1+delta_lb(2),dp)*dr(2).LT.REAL(jb,dp)*dbox(2).OR.&
               REAL(je-1+delta_lb(2),dp)*dr(2).GT.REAL(jb+1,dp)*dbox(2)) THEN
               WRITE(*,*) 'box_j',REAL(jb,dp)*dbox(2),REAL(jb+1,dp)*dbox(2)
               WRITE(*,*) 'point_j',REAL(js-1+delta_lb(2),dp)*dr(2),REAL(je-1+delta_lb(2),dp)*dr(2)
               WRITE(*,*) 'is,ie',is,ie,' js,je',js,je,' ks,ke',ks,ke
               WRITE(*,*) 'ibs,ibe',ibs,ibe,' jbs,jbe',jbs,jbe,' kbs,kbe',kbs,kbe
               STOP 'we stop_j'
            ENDIF
            IF(REAL(is-1+delta_lb(1),dp)*dr(1).LT.REAL(ib,dp)*dbox(1).OR.&
               REAL(ie-1+delta_lb(1),dp)*dr(1).GT.REAL(ib+1,dp)*dbox(1)) THEN
               WRITE(*,*) 'box_i',REAL(ib,dp)*dbox(1),REAL(ib+1,dp)*dbox(1)
               WRITE(*,*) 'point_i',REAL(is-1+delta_lb(1),dp)*dr(1),REAL(ie-1+delta_lb(1),dp)*dr(1)
               WRITE(*,*) 'is,ie',is,ie,' js,je',js,je,' ks,ke',ks,ke
               WRITE(*,*) 'ibs,ibe',ibs,ibe,' jbs,jbe',jbs,jbe,' kbs,kbe',kbs,kbe
               STOP 'we stop_i'
            ENDIF
         ENDIF
         !
         ! the center of the box
         box_center(1) = ( REAL(ibox,dp) + 0.5_dp ) * dbox(1)
         box_center(2) = ( REAL(jbox,dp) + 0.5_dp ) * dbox(2)
         box_center(3) = ( REAL(kbox,dp) + 0.5_dp ) * dbox(3)
         !
         ! find the atoms that are in the overlaping boxes
         natms_local0 = box(ibox,jbox,kbox)%n
         r_ptr => box(ibox,jbox,kbox)%r
         !
         ! go over the grid inside the box
         IF(natms_local0.GT.0) THEN
            !
            ! here there are some atoms...
            DO k = ks,ke
            DO j = js,je
            DO i = is,ie
               point(1) = REAL(i-1+delta_lb(1),dp)*dr(1)
               point(2) = REAL(j-1+delta_lb(2),dp)*dr(2)
               point(3) = REAL(k-1+delta_lb(3),dp)*dr(3)
               point = pbc(point,cell)
               !
               ! compute atom-point distances
               natms_local1 = 0
               DO iatom = 1,natms_local0
                  r(:) = pbc(r_ptr(:,iatom)-point(:),cell) + point(:)!needed?
                  dist = (r(1)-point(1))**2 + (r(2)-point(2))**2 + (r(3)-point(3))**2
                  IF(dist.LT.gauge_atom_radius**2) THEN
                     natms_local1 = natms_local1 + 1
                     ratom(:,natms_local1) = r(:)
                     !
                     ! compute the distance atoms-point
                     x = point(1) - r(1)
                     y = point(2) - r(2)
                     z = point(3) - r(3)
                     atms_pnt(1,natms_local1) = x
                     atms_pnt(2,natms_local1) = y
                     atms_pnt(3,natms_local1) = z
                     nrm_atms_pnt(natms_local1) = SQRT( x*x + y*y + z*z )
                  ENDIF
               ENDDO
               !
               !
               IF(natms_local1.GT.0) THEN
                  !
                  ! build the step
                  DO iatom = 1,natms_local1
                     buf_tmp = 1.0_dp
                     pra(1) = atms_pnt(1,iatom)
                     pra(2) = atms_pnt(2,iatom)
                     pra(3) = atms_pnt(3,iatom)
                     pa     = nrm_atms_pnt(iatom)
                     DO jatom = 1,natms_local1
                        IF(iatom.EQ.jatom) CYCLE
                        pb = nrm_atms_pnt(jatom)
                        x  = pra(1) - atms_pnt(1,jatom)
                        y  = pra(2) - atms_pnt(2,jatom)
                        z  = pra(3) - atms_pnt(3,jatom)
                        ab  = SQRT( x*x + y*y + z*z )
                        !
                        tmp = ( pa - pb ) / ab
                        tmp = 0.5_dp * ( 3.0_dp - tmp * tmp ) * tmp
                        tmp = 0.5_dp * ( 3.0_dp - tmp * tmp ) * tmp
                        tmp = 0.5_dp * ( 3.0_dp - tmp * tmp ) * tmp
                        buf_tmp = buf_tmp * 0.5_dp * ( 1.0_dp - tmp )
                     ENDDO
                     buf(iatom) = buf_tmp
                  ENDDO
                  res(1) = 0.0_dp
                  res(2) = 0.0_dp
                  res(3) = 0.0_dp
                  summe  = 0.0_dp
                  DO iatom = 1,natms_local1
                     res(1) = res(1) + ratom(1,iatom) * buf(iatom)
                     res(2) = res(2) + ratom(2,iatom) * buf(iatom)
                     res(3) = res(3) + ratom(3,iatom) * buf(iatom)
                     summe  = summe + buf(iatom)
                  ENDDO
                  res(1) = res(1) / summe
                  res(2) = res(2) / summe
                  res(3) = res(3) / summe
                  grid_x(i,j,k) = point(1) - res(1)
                  grid_y(i,j,k) = point(2) - res(2)
                  grid_z(i,j,k) = point(3) - res(3)
               ELSE
                  grid_x(i,j,k) = 0.0_dp
                  grid_y(i,j,k) = 0.0_dp
                  grid_z(i,j,k) = 0.0_dp
               ENDIF
            ENDDO ! grid
            ENDDO
            ENDDO
            !
         ELSE
            !
            ! here there is no atom
            DO k = ks,ke
            DO j = js,je
            DO i = is,ie
               grid_x(i,j,k) = 0.0_dp
               grid_y(i,j,k) = 0.0_dp
               grid_z(i,j,k) = 0.0_dp
            ENDDO ! grid
            ENDDO
            ENDDO
            !
         ENDIF
         !
      ENDDO ! list
      ENDDO
      ENDDO

      DEALLOCATE(buf,ratom,atms_pnt,nrm_atms_pnt)

      CALL timestop(handle)

    END SUBROUTINE collocate_gauge_new

! *****************************************************************************
!> \brief ...
!> \param box ...
! *****************************************************************************
    SUBROUTINE deallocate_box(box)
    TYPE(box_type), DIMENSION(:, :, :), &
      POINTER                                :: box

    INTEGER                                  :: i, j, k

      IF(ASSOCIATED(box)) THEN
         DO k = LBOUND(box,3),UBOUND(box,3)
         DO j = LBOUND(box,2),UBOUND(box,2)
         DO i = LBOUND(box,1),UBOUND(box,1)
            IF(ASSOCIATED(box(i,j,k)%r)) THEN
               DEALLOCATE(box(i,j,k)%r)
            ENDIF
         ENDDO
         ENDDO
         ENDDO
         DEALLOCATE(box)
      ENDIF
    END SUBROUTINE deallocate_box
  END SUBROUTINE current_set_gauge


! *****************************************************************************
!> \brief ...
!> \param current_env ...
!> \param qs_env ...
!> \param iB ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE current_build_chi(current_env,qs_env,iB,error)
    !
    TYPE(current_env_type)                   :: current_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    INTEGER, INTENT(IN)                      :: iB
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    IF(current_env%full) THEN
       CALL current_build_chi_many_centers(current_env,qs_env,iB,error)
    ELSE
       CALL current_build_chi_one_center(current_env,qs_env,iB,error)
    ENDIF

  END SUBROUTINE current_build_chi


! *****************************************************************************
!> \brief ...
!> \param current_env ...
!> \param qs_env ...
!> \param iB ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE current_build_chi_many_centers(current_env,qs_env,iB,error)
    !
    TYPE(current_env_type)                   :: current_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    INTEGER, INTENT(IN)                      :: iB
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'current_build_chi_many_centers', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, icenter, idir, idir2, ii, iiB, iii, iiiB, ispin, &
      istat, istate, j, jstate, max_states, nao, natom, nbr_center(2), nmo, &
      nspins, nstate_loc, nstates(2), output_unit
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: first_sgf, last_sgf
    INTEGER, DIMENSION(:), POINTER           :: rbs
    LOGICAL                                  :: chi_pbc, failure, gapw
    REAL(dp)                                 :: chi(3), chi_tmp, contrib, &
                                                contrib2, dk(3), &
                                                int_current(3), &
                                                int_current_tmp, maxocc
    TYPE(array_i1d_obj)                      :: row_blk_sizes
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_2d_i_p_type), DIMENSION(:), &
      POINTER                                :: center_list
    TYPE(cp_2d_r_p_type), DIMENSION(:), &
      POINTER                                :: centers_set
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: op_mom_ao, op_p_ao
    TYPE(cp_dbcsr_p_type), DIMENSION(:, :), &
      POINTER                                :: op_mom_der_ao
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: p_rxp, psi0_order, r_p1, r_p2
    TYPE(cp_fm_p_type), DIMENSION(:, :), &
      POINTER                                :: psi1_D, psi1_p, psi1_rxp, &
                                                rr_p1, rr_p2, rr_rxp
    TYPE(cp_fm_struct_type), POINTER         :: tmp_fm_struct
    TYPE(cp_fm_type), POINTER                :: mo_coeff, psi0, psi_D, &
                                                psi_p1, psi_p2, psi_rxp
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dbcsr_distribution_obj), POINTER    :: dbcsr_dist
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_all, sab_orb
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set

    failure = .FALSE.
    !
    CALL timeset(routineN,handle)
    !
    NULLIFY(dft_control,mos,para_env,mo_coeff,op_mom_ao,&
            op_mom_der_ao,center_list,centers_set,psi0,psi_rxp,psi_D,psi_p1,psi_p2,&
            op_p_ao,psi1_p,psi1_rxp,psi1_D,p_rxp,r_p1,r_p2,rr_rxp,rr_p1,rr_p2,&
            cell,psi0_order,particle_set,qs_kind_set)

    logger => cp_error_get_logger(error)
    output_unit= cp_logger_get_default_io_unit(logger)

    CALL get_qs_env(qs_env=qs_env,&
                    dft_control=dft_control,&
                    mos=mos,&
                    para_env=para_env,&
                    cell=cell,&
                    dbcsr_dist=dbcsr_dist,&
                    particle_set=particle_set,&
                    qs_kind_set=qs_kind_set,&
                    sab_all=sab_all,&
                    sab_orb=sab_orb,&
                    error=error)

    nspins = dft_control%nspins
    gapw = dft_control%qs_control%gapw

    CALL get_current_env(current_env=current_env,&
                         chi_pbc=chi_pbc,&
                         nao=nao,&
                         nbr_center=nbr_center,&
                         center_list=center_list,&
                         centers_set=centers_set,&
                         psi1_p=psi1_p,&
                         psi1_rxp=psi1_rxp,&
                         psi1_D=psi1_D,&
                         nstates=nstates,&
                         psi0_order=psi0_order,&
                         error=error)
    !
    ! get max nbr of states per center
    max_states = 0
    DO ispin = 1,nspins
       DO icenter = 1,nbr_center(ispin)
          max_states = MAX(max_states,center_list(ispin)%array(1,icenter+1)&
               &                     -center_list(ispin)%array(1,icenter))
       ENDDO
    ENDDO
    !
    ! Allocate sparse matrices for dipole, quadrupole and their derivatives => 9x3
    ! Remember the derivatives are antisymmetric
    CALL cp_dbcsr_allocate_matrix_set(op_mom_ao,9,error=error)
    CALL cp_dbcsr_allocate_matrix_set(op_mom_der_ao,9,3,error=error)
    !
    ! prepare for allocation
    natom = SIZE(particle_set,1)
    ALLOCATE (first_sgf(natom),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (last_sgf(natom),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    CALL get_particle_set(particle_set, qs_kind_set,&
                          first_sgf=first_sgf,&
                          last_sgf=last_sgf,error=error)
    ALLOCATE (rbs(natom), STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    CALL convert_offsets_to_sizes (first_sgf, rbs, last_sgf)
    CALL array_nullify (row_blk_sizes)
    CALL array_new (row_blk_sizes, rbs, gift=.TRUE.)
    DEALLOCATE (first_sgf,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (last_sgf,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    !
    !
    ALLOCATE(op_mom_ao(1)%matrix)
    CALL cp_dbcsr_init(op_mom_ao(1)%matrix,error=error)
    CALL cp_dbcsr_create(matrix=op_mom_ao(1)%matrix, &
         name="op_mom", &
         dist=dbcsr_dist, matrix_type=dbcsr_type_no_symmetry,&
         row_blk_size=row_blk_sizes, col_blk_size=row_blk_sizes, &
         nblks=0, nze=0, mutable_work=.TRUE., &
         error=error)
    CALL cp_dbcsr_alloc_block_from_nbl(op_mom_ao(1)%matrix,sab_all,error=error)

    DO idir2=1,3
       ALLOCATE(op_mom_der_ao(1,idir2)%matrix)
       CALL cp_dbcsr_init(op_mom_der_ao(1,idir2)%matrix, error=error)
       CALL cp_dbcsr_copy(op_mom_der_ao(1,idir2)%matrix,op_mom_ao(1)%matrix,&
            "op_mom_der_ao"//"-"//TRIM(ADJUSTL(cp_to_string(idir2))),error)
    ENDDO

    DO idir = 2,SIZE(op_mom_ao,1)
       ALLOCATE(op_mom_ao(idir)%matrix)
       CALL cp_dbcsr_init(op_mom_ao(idir)%matrix, error=error)
       CALL cp_dbcsr_copy(op_mom_ao(idir)%matrix,op_mom_ao(1)%matrix,&
            "op_mom_ao"//"-"//TRIM(ADJUSTL(cp_to_string(idir))),error)
       DO idir2=1,3
          ALLOCATE(op_mom_der_ao(idir,idir2)%matrix)
          CALL cp_dbcsr_init(op_mom_der_ao(idir,idir2)%matrix, error=error)
          CALL cp_dbcsr_copy(op_mom_der_ao(idir,idir2)%matrix,op_mom_ao(1)%matrix,&
               "op_mom_der_ao"//"-"//TRIM(ADJUSTL(cp_to_string(idir*idir2))),error)
       ENDDO
    ENDDO
    !
    CALL cp_dbcsr_allocate_matrix_set(op_p_ao,3,error=error)
    ALLOCATE(op_p_ao(1)%matrix)
    CALL cp_dbcsr_init(op_p_ao(1)%matrix,error=error)
    CALL cp_dbcsr_create(matrix=op_p_ao(1)%matrix, &
         name="op_p_ao", &
         dist=dbcsr_dist, matrix_type=dbcsr_type_antisymmetric,&
         row_blk_size=row_blk_sizes, col_blk_size=row_blk_sizes, &
         nblks=0, nze=0, mutable_work=.TRUE., &
         error=error)
    CALL cp_dbcsr_alloc_block_from_nbl(op_p_ao(1)%matrix, sab_orb, error=error)

    DO idir=2,3
       ALLOCATE(op_p_ao(idir)%matrix)
       CALL cp_dbcsr_init(op_p_ao(idir)%matrix, error=error)
       CALL cp_dbcsr_copy(op_p_ao(idir)%matrix,op_p_ao(1)%matrix,&
            "op_p_ao"//"-"//TRIM(ADJUSTL(cp_to_string(idir))),error)
    ENDDO
    !
    !
    CALL array_release (row_blk_sizes)
    !
    !
    ! Allocate full matrices for only one vector
    mo_coeff => psi0_order(1)%matrix
    NULLIFY(tmp_fm_struct)
    CALL cp_fm_struct_create(tmp_fm_struct,nrow_global=nao,&
                             ncol_global=max_states,para_env=para_env,&
                             context=mo_coeff%matrix_struct%context,&
                             error=error)
    CALL cp_fm_create(psi0,tmp_fm_struct,error=error)
    CALL cp_fm_create(psi_D,tmp_fm_struct,error=error)
    CALL cp_fm_create(psi_rxp,tmp_fm_struct,error=error)
    CALL cp_fm_create(psi_p1,tmp_fm_struct,error=error)
    CALL cp_fm_create(psi_p2,tmp_fm_struct,error=error)
    CALL cp_fm_struct_release(tmp_fm_struct,error=error)
    !
    ALLOCATE(p_rxp(3),r_p1(3),r_p2(3),rr_rxp(9,3),rr_p1(9,3),rr_p2(9,3),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    CALL cp_fm_struct_create(tmp_fm_struct,nrow_global=nao,&
                             ncol_global=max_states,para_env=para_env,&
                             context=mo_coeff%matrix_struct%context,&
                             error=error)
    DO idir = 1,3
       CALL cp_fm_create(p_rxp(idir)%matrix,tmp_fm_struct,error=error)
       CALL cp_fm_create(r_p1(idir)%matrix,tmp_fm_struct,error=error)
       CALL cp_fm_create(r_p2(idir)%matrix,tmp_fm_struct,error=error)
       DO idir2 = 1,9
          CALL cp_fm_create(rr_rxp(idir2,idir)%matrix,tmp_fm_struct,error=error)
          CALL cp_fm_create(rr_p1(idir2,idir)%matrix,tmp_fm_struct,error=error)
          CALL cp_fm_create(rr_p2(idir2,idir)%matrix,tmp_fm_struct,error=error)
       ENDDO
    ENDDO
    CALL cp_fm_struct_release(tmp_fm_struct,error=error)
    !
    !
    !
    ! recompute the linear momentum matrices
    CALL build_lin_mom_matrix(qs_env,op_p_ao,error)
    !CALL p_xyz_ao(op_p_ao,qs_env,minimum_image=.FALSE.,error=error)
    !
    !
    ! get iiB and iiiB
    CALL set_vecp(iB,iiB,iiiB)
    DO ispin = 1,nspins
       !
       ! get ground state MOS
       nmo = nstates(ispin)
       mo_coeff => psi0_order(ispin)%matrix
       CALL get_mo_set(mo_set=mos(ispin)%mo_set,maxocc=maxocc)
       !
       ! Initialize the temporary vector chi
       chi = 0.0_dp
       int_current = 0.0_dp
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
             CALL cp_dbcsr_set(op_mom_ao(idir)%matrix,0.0_dp,error=error)
             DO idir2 = 1,3
                CALL cp_dbcsr_set(op_mom_der_ao(idir,idir2)%matrix,0.0_dp,error=error)
             ENDDO
          ENDDO
          CALL rRc_xyz_der_ao(op_mom_ao,op_mom_der_ao,qs_env,dk,order=2,&
                              minimum_image=.FALSE.,soft=gapw,error=error)
          !
          ! collecte the states that belong to a given center
          CALL cp_fm_set_all(psi0,0.0_dp,error=error)
          CALL cp_fm_set_all(psi_rxp,0.0_dp,error=error)
          CALL cp_fm_set_all(psi_D,0.0_dp,error=error)
          CALL cp_fm_set_all(psi_p1,0.0_dp,error=error)
          CALL cp_fm_set_all(psi_p2,0.0_dp,error=error)
          nstate_loc = center_list(ispin)%array(1,icenter+1)-center_list(ispin)%array(1,icenter)
          jstate = 1
          DO j = center_list(ispin)%array(1,icenter),center_list(ispin)%array(1,icenter+1)-1
             istate = center_list(ispin)%array(2,j)
             !
             ! block the states that belong to this center
             CALL cp_fm_to_fm(mo_coeff,psi0,1,istate,jstate)
             !
             CALL cp_fm_to_fm(psi1_rxp(ispin,iB)%matrix,psi_rxp,1,istate,jstate)
             IF(current_env%full) CALL cp_fm_to_fm(psi1_D(ispin,iB)%matrix,psi_D,1,istate,jstate)
             !
             ! psi1_p_iiB_istate and psi1_p_iiiB_istate
             CALL cp_fm_to_fm(psi1_p(ispin, iiB)%matrix,psi_p1,1,istate,jstate)
             CALL cp_fm_to_fm(psi1_p(ispin,iiiB)%matrix,psi_p2,1,istate,jstate)
             !
             jstate = jstate+1
          ENDDO ! istate
          !
          ! scale the ordered mos
          IF(current_env%full) CALL cp_fm_scale_and_add(1.0_dp,psi_rxp,-1.0_dp,psi_D,error=error)
          !
          DO idir = 1,3
             CALL set_vecp(idir,ii,iii)
             CALL cp_dbcsr_sm_fm_multiply(op_p_ao(idir)%matrix,psi_rxp,&
                                    p_rxp(idir)%matrix,ncol=nstate_loc,alpha=1.e0_dp,&
                                    error=error)
             IF(iiiB.EQ.iii.OR.iiiB.EQ.ii) THEN
                CALL cp_dbcsr_sm_fm_multiply(op_mom_ao(idir)%matrix,psi_p1,&
                                       r_p1(idir)%matrix,ncol=nstate_loc,alpha=1.e0_dp,&
                                       error=error)
             ENDIF
             IF(iiB.EQ.iii.OR.iiB.EQ.ii) THEN
                CALL cp_dbcsr_sm_fm_multiply(op_mom_ao(idir)%matrix,psi_p2,&
                                       r_p2(idir)%matrix,ncol=nstate_loc,alpha=1.e0_dp,&
                                       error=error)
             ENDIF
             DO idir2 = 1,9
                IF(idir2.EQ.ii.OR.idir2.EQ.iii) THEN
                   CALL cp_dbcsr_sm_fm_multiply(op_mom_der_ao(idir2,idir)%matrix,psi_rxp,&
                                          rr_rxp(idir2,idir)%matrix,ncol=nstate_loc,alpha=1.e0_dp,&
                                          error=error)
                ENDIF
                !
                IF(idir2.EQ.ind_m2(ii,iiiB).OR.idir2.EQ.ind_m2(iii,iiiB).OR.idir2.EQ.iiiB) THEN
                   CALL cp_dbcsr_sm_fm_multiply(op_mom_der_ao(idir2,idir)%matrix,psi_p1,&
                                          rr_p1(idir2,idir)%matrix,ncol=nstate_loc,alpha=1.e0_dp,&
                                          error=error)
                ENDIF
                !
                IF(idir2.EQ.ind_m2(ii,iiB).OR.idir2.EQ.ind_m2(iii,iiB).OR.idir2.EQ.iiB) THEN
                   CALL cp_dbcsr_sm_fm_multiply(op_mom_der_ao(idir2,idir)%matrix,psi_p2,&
                                          rr_p2(idir2,idir)%matrix,ncol=nstate_loc,alpha=1.e0_dp,&
                                          error=error)
                ENDIF
             ENDDO
          ENDDO
          !
          ! Multuply left and right by the appropriate coefficients and sum into the
          ! correct component of the chi tensor using the appropriate multiplicative factor
          ! (don't forget the occupation number)
          ! Loop over the cartesian components of the tensor
          ! The loop over the components of the external field is external, thereby
          ! only one column of the chi tensor is computed here
          DO idir = 1,3
             chi_tmp = 0.0_dp
             int_current_tmp = 0.0_dp
             !
             ! get ii and iii
             CALL set_vecp(idir,ii,iii)
             !
             ! term: 2[C0| (r-dk)_ii |d_iii(C1(rxp-D))]-2[C0| (r-dk)_iii |d_ii(C1(rxp-D))]
             ! the factor 2 should be already included in the matrix elements
             contrib = 0.0_dp
             CALL cp_fm_trace(psi0,rr_rxp(ii,iii)%matrix,contrib,error=error)
             chi_tmp = chi_tmp + 2.0_dp * contrib
             !
             contrib = 0.0_dp
             CALL cp_fm_trace(psi0,rr_rxp(iii,ii)%matrix,contrib,error=error)
             chi_tmp = chi_tmp - 2.0_dp * contrib
             !
             ! correction: dk_ii*2[C0| d_iii(C1(rxp-D))] - dk_iii*2[C0| d_ii(C1(rxp-D))]
             ! factor 2 not included in the matrix elements
             contrib = 0.0_dp
             CALL cp_fm_trace(psi0,p_rxp(iii)%matrix,contrib,error=error)
             IF(.NOT.chi_pbc) chi_tmp = chi_tmp + 2.0_dp * dk(ii) * contrib
             int_current_tmp = int_current_tmp + 2.0_dp * contrib
             !
             contrib2 = 0.0_dp
             CALL cp_fm_trace(psi0,p_rxp(ii)%matrix,contrib2,error=error)
             IF(.NOT.chi_pbc) chi_tmp = chi_tmp - 2.0_dp * dk(iii) * contrib2
             !
             ! term: -2[C0| (r-dk)_ii  (r-dk)_iiB | d_iii(C1(piiiB))] \
             !       +2[C0| (r-dk)_iii (r-dk)_iiB | d_ii(C1(piiiB))]
             ! the factor 2 should be already included in the matrix elements
             contrib = 0.0_dp
             idir2 = ind_m2(ii,iiB)
             CALL cp_fm_trace(psi0,rr_p2(idir2,iii)%matrix,contrib,error=error)
             chi_tmp = chi_tmp - 2.0_dp * contrib
             contrib2 = 0.0_dp
             IF(iiB==iii) THEN
                CALL cp_fm_trace(psi0,r_p2(ii)%matrix,contrib2,error=error)
                chi_tmp = chi_tmp - contrib2
             ENDIF
             !
             contrib = 0.0_dp
             idir2 = ind_m2(iii,iiB)
             CALL cp_fm_trace(psi0,rr_p2(idir2,ii)%matrix,contrib,error=error)
             chi_tmp = chi_tmp + 2.0_dp * contrib
             contrib2 = 0.0_dp
             IF(iiB==ii) THEN
                CALL cp_fm_trace(psi0,r_p2(iii)%matrix,contrib2,error=error)
                chi_tmp = chi_tmp + contrib2
             ENDIF
             !
             ! correction: -dk_ii * 2[C0|(r-dk)_iiB | d_iii(C1(piiiB))] \
             !             +dk_iii * 2[C0|(r-dk)_iiB | d_ii(C1(piiiB))]
             ! the factor 2 should be already included in the matrix elements
             ! no additional correction terms because of the orthogonality between C0 and C1
             contrib = 0.0_dp
             CALL cp_fm_trace(psi0,rr_p2(iiB,iii)%matrix,contrib,error=error)
             IF(.NOT.chi_pbc) chi_tmp = chi_tmp - 2.0_dp * dk(ii) * contrib
             int_current_tmp = int_current_tmp - 2.0_dp * contrib
             !
             contrib2 = 0.0_dp
             CALL cp_fm_trace(psi0,rr_p2(iiB,ii)%matrix,contrib2,error=error)
             IF(.NOT.chi_pbc) chi_tmp = chi_tmp + 2.0_dp * dk(iii) * contrib2
             !
             ! term: +2[C0| (r-dk)_ii  (r-dk)_iiiB | d_iii(C1(piiB))] \
             !       -2[C0| (r-dk)_iii (r-dk)_iiiB | d_ii(C1(piiB))]
             ! the factor 2 should be already included in the matrix elements
             contrib = 0.0_dp
             idir2 = ind_m2(ii,iiiB)
             CALL cp_fm_trace(psi0,rr_p1(idir2,iii)%matrix,contrib,error=error)
             chi_tmp = chi_tmp + 2.0_dp * contrib
             contrib2 = 0.0_dp
             IF(iiiB==iii) THEN
                CALL cp_fm_trace(psi0,r_p1(ii)%matrix,contrib2,error=error)
                chi_tmp = chi_tmp + contrib2
             ENDIF
             !
             contrib = 0.0_dp
             idir2 = ind_m2(iii,iiiB)
             CALL cp_fm_trace(psi0,rr_p1(idir2,ii)%matrix,contrib,error=error)
             chi_tmp = chi_tmp - 2.0_dp * contrib
             contrib2 = 0.0_dp
             IF(iiiB==ii) THEN
                CALL cp_fm_trace(psi0,r_p1(iii)%matrix,contrib2,error=error)
                chi_tmp = chi_tmp - contrib2
             ENDIF
             !
             ! correction: +dk_ii * 2[C0|(r-dk)_iiiB | d_iii(C1(piiB))] +\
             !             -dk_iii * 2[C0|(r-dk)_iiiB | d_ii(C1(piiB))]
             ! the factor 2 should be already included in the matrix elements
             contrib = 0.0_dp
             CALL cp_fm_trace(psi0,rr_p1(iiiB,iii)%matrix,contrib,error=error)
             IF(.NOT.chi_pbc) chi_tmp = chi_tmp + 2.0_dp * dk(ii) * contrib
             int_current_tmp = int_current_tmp + 2.0_dp * contrib
             !
             contrib2 = 0.0_dp
             CALL cp_fm_trace(psi0,rr_p1(iiiB,ii)%matrix,contrib2,error=error)
             IF(.NOT.chi_pbc) chi_tmp = chi_tmp - 2.0_dp * dk(iii) * contrib2
             !
             ! accumulate
             chi(idir) = chi(idir) + maxocc * chi_tmp
             int_current(iii) = int_current(iii) + int_current_tmp
          ENDDO  ! idir

       ENDDO  ! icenter
       !
       DO idir = 1,3
          current_env%chi_tensor(idir,iB,ispin) = current_env%chi_tensor(idir,iB,ispin)+&
                                                  chi(idir)
          IF (output_unit>0) THEN
             !WRITE(output_unit,'(A,E12.6)') ' chi_'//ACHAR(119+idir)//ACHAR(119+iB)//&
             !     &                         ' = ',chi(idir)
             !WRITE(output_unit,'(A,E12.6)') ' analytic \int j_'//ACHAR(119+idir)//ACHAR(119+iB)//&
             !     &                         '(r) d^3r = ',int_current(idir)
          ENDIF
       ENDDO
       !
    ENDDO ! ispin
    !
    ! deallocate the sparse matrices
    CALL cp_dbcsr_deallocate_matrix_set(op_mom_ao,error=error)
    CALL cp_dbcsr_deallocate_matrix_set(op_mom_der_ao,error=error)
    CALL cp_dbcsr_deallocate_matrix_set(op_p_ao,error=error)

    CALL cp_fm_release(psi0,error=error)
    CALL cp_fm_release(psi_rxp,error=error)
    CALL cp_fm_release(psi_D,error=error)
    CALL cp_fm_release(psi_p1,error=error)
    CALL cp_fm_release(psi_p2,error=error)
    DO idir = 1,3
       CALL cp_fm_release(p_rxp(idir)%matrix,error=error)
       CALL cp_fm_release(r_p1(idir)%matrix,error=error)
       CALL cp_fm_release(r_p2(idir)%matrix,error=error)
       DO idir2 = 1,9
          CALL cp_fm_release(rr_rxp(idir2,idir)%matrix,error=error)
          CALL cp_fm_release(rr_p1(idir2,idir)%matrix,error=error)
          CALL cp_fm_release(rr_p2(idir2,idir)%matrix,error=error)
       ENDDO
    ENDDO
    DEALLOCATE(p_rxp,r_p1,r_p2,rr_rxp,rr_p1,rr_p2,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE current_build_chi_many_centers


! *****************************************************************************
!> \brief ...
!> \param current_env ...
!> \param qs_env ...
!> \param iB ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE current_build_chi_one_center(current_env,qs_env,iB,error)
    !
    TYPE(current_env_type)                   :: current_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    INTEGER, INTENT(IN)                      :: iB
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'current_build_chi_one_center', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, idir, idir2, iiB, iiiB, ispin, istat, jdir, jjdir, &
      kdir, max_states, nao, natom, nbr_center(2), nmo, nspins, nstates(2), &
      output_unit
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: first_sgf, last_sgf
    INTEGER, DIMENSION(:), POINTER           :: rbs
    LOGICAL                                  :: chi_pbc, failure, gapw
    REAL(dp)                                 :: chi(3), contrib, dk(3), &
                                                int_current(3), maxocc
    TYPE(array_i1d_obj)                      :: row_blk_sizes
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_2d_i_p_type), DIMENSION(:), &
      POINTER                                :: center_list
    TYPE(cp_2d_r_p_type), DIMENSION(:), &
      POINTER                                :: centers_set
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: op_mom_ao, op_p_ao
    TYPE(cp_dbcsr_p_type), DIMENSION(:, :), &
      POINTER                                :: op_mom_der_ao
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: psi0_order
    TYPE(cp_fm_p_type), DIMENSION(:, :), &
      POINTER                                :: psi1_D, psi1_p, psi1_rxp
    TYPE(cp_fm_struct_type), POINTER         :: tmp_fm_struct
    TYPE(cp_fm_type), POINTER                :: buf, mo_coeff
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dbcsr_distribution_obj), POINTER    :: dbcsr_dist
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_all, sab_orb
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set

    failure = .FALSE.
    !
    CALL timeset(routineN,handle)
    !
    NULLIFY(dft_control,mos,para_env,mo_coeff,op_mom_ao,&
            op_mom_der_ao,center_list,centers_set,&
            op_p_ao,psi1_p,psi1_rxp,psi1_D,buf,cell,psi0_order)

    logger => cp_error_get_logger(error)
    output_unit= cp_logger_get_default_io_unit(logger)

    CALL get_qs_env(qs_env=qs_env,&
                    dft_control=dft_control,&
                    mos=mos,&
                    para_env=para_env,&
                    cell=cell,&
                    particle_set=particle_set,&
                    qs_kind_set=qs_kind_set,&
                    sab_all=sab_all,&
                    sab_orb=sab_orb,&
                    dbcsr_dist=dbcsr_dist,&
                    error=error)

    nspins = dft_control%nspins
    gapw = dft_control%qs_control%gapw

    CALL get_current_env(current_env=current_env,&
                         chi_pbc=chi_pbc,&
                         nao=nao,&
                         nbr_center=nbr_center,&
                         center_list=center_list,&
                         centers_set=centers_set,&
                         psi1_p=psi1_p,&
                         psi1_rxp=psi1_rxp,&
                         psi1_D=psi1_D,&
                         nstates=nstates,&
                         psi0_order=psi0_order,&
                         error=error)
    !
    max_states = MAXVAL(nstates(1:nspins))
    !
    ! Allocate sparse matrices for dipole, quadrupole and their derivatives => 9x3
    ! Remember the derivatives are antisymmetric
    CALL cp_dbcsr_allocate_matrix_set(op_mom_ao,9,error=error)
    CALL cp_dbcsr_allocate_matrix_set(op_mom_der_ao,9,3,error=error)
    !
    ! prepare for allocation
    natom = SIZE(particle_set,1)
    ALLOCATE (first_sgf(natom),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (last_sgf(natom),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    CALL get_particle_set(particle_set, qs_kind_set,&
                          first_sgf=first_sgf,&
                          last_sgf=last_sgf,error=error)
    ALLOCATE (rbs(natom), STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    CALL convert_offsets_to_sizes (first_sgf, rbs, last_sgf)
    CALL array_nullify (row_blk_sizes)
    CALL array_new (row_blk_sizes, rbs, gift=.TRUE.)
    DEALLOCATE (first_sgf,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (last_sgf,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    !
    !
    ALLOCATE(op_mom_ao(1)%matrix)
    CALL cp_dbcsr_init(op_mom_ao(1)%matrix,error=error)
    CALL cp_dbcsr_create(matrix=op_mom_ao(1)%matrix, &
         name="op_mom", &
         dist=dbcsr_dist, matrix_type=dbcsr_type_no_symmetry,&
         row_blk_size=row_blk_sizes, col_blk_size=row_blk_sizes, &
         nblks=0, nze=0, mutable_work=.TRUE., &
         error=error)
    CALL cp_dbcsr_alloc_block_from_nbl(op_mom_ao(1)%matrix,sab_all,error=error)

    DO idir2=1,3
       ALLOCATE(op_mom_der_ao(1,idir2)%matrix)
       CALL cp_dbcsr_init(op_mom_der_ao(1,idir2)%matrix, error=error)
       CALL cp_dbcsr_copy(op_mom_der_ao(1,idir2)%matrix,op_mom_ao(1)%matrix,&
            "op_mom_der_ao"//"-"//TRIM(ADJUSTL(cp_to_string(idir2))),error)
    ENDDO

    DO idir = 2,SIZE(op_mom_ao,1)
       ALLOCATE(op_mom_ao(idir)%matrix)
       CALL cp_dbcsr_init(op_mom_ao(idir)%matrix, error=error)
       CALL cp_dbcsr_copy(op_mom_ao(idir)%matrix,op_mom_ao(1)%matrix,&
            "op_mom_ao"//"-"//TRIM(ADJUSTL(cp_to_string(idir))),error)
       DO idir2=1,3
          ALLOCATE(op_mom_der_ao(idir,idir2)%matrix)
          CALL cp_dbcsr_init(op_mom_der_ao(idir,idir2)%matrix, error=error)
          CALL cp_dbcsr_copy(op_mom_der_ao(idir,idir2)%matrix,op_mom_ao(1)%matrix,&
               "op_mom_der_ao"//"-"//TRIM(ADJUSTL(cp_to_string(idir*idir2))),error)
       ENDDO
    ENDDO
    !
    CALL cp_dbcsr_allocate_matrix_set(op_p_ao,3,error=error)
    ALLOCATE(op_p_ao(1)%matrix)
    CALL cp_dbcsr_init(op_p_ao(1)%matrix,error=error)
    CALL cp_dbcsr_create(matrix=op_p_ao(1)%matrix, &
         name="op_p_ao", &
         dist=dbcsr_dist, matrix_type=dbcsr_type_antisymmetric,&
         row_blk_size=row_blk_sizes, col_blk_size=row_blk_sizes, &
         nblks=0, nze=0, mutable_work=.TRUE., &
         error=error)
    CALL cp_dbcsr_alloc_block_from_nbl(op_p_ao(1)%matrix, sab_orb, error=error)

    DO idir=2,3
       ALLOCATE(op_p_ao(idir)%matrix)
       CALL cp_dbcsr_init(op_p_ao(idir)%matrix, error=error)
       CALL cp_dbcsr_copy(op_p_ao(idir)%matrix,op_p_ao(1)%matrix,&
            "op_p_ao"//"-"//TRIM(ADJUSTL(cp_to_string(idir))),error)
    ENDDO
    !
    !
    CALL array_release (row_blk_sizes)
    !
    !
    ! Allocate full matrices for only one vector
    mo_coeff => psi0_order(1)%matrix
    NULLIFY(tmp_fm_struct)
    CALL cp_fm_struct_create(tmp_fm_struct,nrow_global=nao,&
                             ncol_global=max_states,para_env=para_env,&
                             context=mo_coeff%matrix_struct%context,&
                             error=error)
    CALL cp_fm_create(buf,tmp_fm_struct,error=error)
    CALL cp_fm_struct_release(tmp_fm_struct,error=error)
    !
    !
    !
    ! recompute the linear momentum matrices
    CALL build_lin_mom_matrix(qs_env,op_p_ao,error)
    !CALL p_xyz_ao(op_p_ao,qs_env,minimum_image=.FALSE.,error=error)
    !
    !
    ! get iiB and iiiB
    CALL set_vecp(iB,iiB,iiiB)
    DO ispin = 1,nspins
       !
       CPPostcondition(nbr_center(ispin)==1,cp_failure_level,routineP,error,failure)
       !
       ! get ground state MOS
       nmo = nstates(ispin)
       mo_coeff => psi0_order(ispin)%matrix
       CALL get_mo_set(mo_set=mos(ispin)%mo_set,maxocc=maxocc)
       !
       ! Initialize the temporary vector chi
       chi = 0.0_dp
       int_current = 0.0_dp
       !
       !
       ! Get the Wannier center of the istate-th ground state orbital
       dk(1:3) = centers_set(ispin)%array(1:3,1)
       !
       ! Compute the multipole integrals for the state istate,
       ! using as reference center the corresponding Wannier center
       DO idir = 1,9
          CALL cp_dbcsr_set(op_mom_ao(idir)%matrix,0.0_dp,error=error)
          DO idir2 = 1,3
             CALL cp_dbcsr_set(op_mom_der_ao(idir,idir2)%matrix,0.0_dp,error=error)
          ENDDO
       ENDDO
       CALL rRc_xyz_der_ao(op_mom_ao,op_mom_der_ao,qs_env,dk,order=2,&
                           minimum_image=.FALSE.,soft=gapw,error=error)
       !
       !
       ! Multuply left and right by the appropriate coefficients and sum into the
       ! correct component of the chi tensor using the appropriate multiplicative factor
       ! (don't forget the occupation number)
       ! Loop over the cartesian components of the tensor
       ! The loop over the components of the external field is external, thereby
       ! only one column of the chi tensor is computed here
       DO idir = 1,3
          !
          !
          !
          ! term: dk_ii*2[C0| d_iii(C1(rxp-D))] - dk_iii*2[C0| d_ii(C1(rxp-D))]
          IF(.NOT.chi_pbc) THEN
             CALL cp_dbcsr_sm_fm_multiply(op_p_ao(idir)%matrix,mo_coeff,&
                                    buf,ncol=nmo,alpha=1.e0_dp,error=error)
             DO jdir = 1,3
                DO kdir = 1,3
                   IF(Levi_Civita(kdir,jdir,idir).EQ.0.0_dp) CYCLE
                   CALL cp_fm_trace(buf,psi1_rxp(ispin,iB)%matrix,contrib,error=error)
                   chi(kdir) = chi(kdir) - Levi_Civita(kdir,jdir,idir) * 2.0_dp * dk(jdir) * contrib
                ENDDO
             ENDDO
          ENDIF
          !
          !
          !
          ! term: 2[C0| (r-dk)_ii |d_iii(C1(rxp-D))]-2[C0| (r-dk)_iii |d_ii(C1(rxp-D))]
          ! and
          ! term: -dk_ii * 2[C0|(r-dk)_iiB | d_iii(C1(piiiB))] +
          !       +dk_iii * 2[C0|(r-dk)_iiB | d_ii(C1(piiiB))]
          ! and
          ! term: +dk_ii * 2[C0|(r-dk)_iiiB | d_iii(C1(piiB))] +
          !       -dk_iii * 2[C0|(r-dk)_iiiB | d_ii(C1(piiB))]
          DO jdir = 1,3
             CALL cp_dbcsr_sm_fm_multiply(op_mom_der_ao(jdir,idir)%matrix,mo_coeff,&
                                    buf,ncol=nmo,alpha=1.e0_dp,error=error)
             DO kdir = 1,3
                IF(Levi_Civita(kdir,jdir,idir).EQ.0.0_dp) CYCLE
                CALL cp_fm_trace(buf,psi1_rxp(ispin,iB)%matrix,contrib,error=error)
                chi(kdir) = chi(kdir) - Levi_Civita(kdir,jdir,idir) * 2.0_dp * contrib
             ENDDO
             !
             IF(.NOT.chi_pbc) THEN
                IF(jdir.EQ.iiB) THEN
                   DO jjdir = 1,3
                      DO kdir = 1,3
                         IF(Levi_Civita(kdir,jjdir,idir).EQ.0.0_dp) CYCLE
                         CALL cp_fm_trace(buf,psi1_p(ispin,iiiB)%matrix,contrib,error=error)
                         chi(kdir) = chi(kdir) + Levi_Civita(kdir,jjdir,idir) * 2.0_dp * dk(jjdir) * contrib
                      ENDDO
                   ENDDO
                ENDIF
                !
                IF(jdir.EQ.iiiB) THEN
                   DO jjdir = 1,3
                      DO kdir = 1,3
                         IF(Levi_Civita(kdir,jjdir,idir).EQ.0.0_dp) CYCLE
                         CALL cp_fm_trace(buf,psi1_p(ispin,iiB)%matrix,contrib,error=error)
                         chi(kdir) = chi(kdir) - Levi_Civita(kdir,jjdir,idir) * 2.0_dp * dk(jjdir) * contrib
                      ENDDO
                   ENDDO
                ENDIF
             ENDIF
          ENDDO
          !
          !
          !
          ! term1: -2[C0| (r-dk)_ii  (r-dk)_iiB | d_iii(C1(piiiB))] +
          !        +2[C0| (r-dk)_iii (r-dk)_iiB | d_ii(C1(piiiB))]
          ! and
          ! term1: +2[C0| (r-dk)_ii  (r-dk)_iiiB | d_iii(C1(piiB))] +
          !        -2[C0| (r-dk)_iii (r-dk)_iiiB | d_ii(C1(piiB))]
          ! HERE THERE IS ONE EXTRA MULTIPLY
          DO jdir = 1,3
             CALL cp_dbcsr_sm_fm_multiply(op_mom_der_ao(ind_m2(jdir,iiB),idir)%matrix,mo_coeff,&
                                    buf,ncol=nmo,alpha=1.e0_dp,error=error)
             DO kdir = 1,3
                IF(Levi_Civita(kdir,jdir,idir).EQ.0.0_dp) CYCLE
                CALL cp_fm_trace(buf,psi1_p(ispin,iiiB)%matrix,contrib,error=error)
                chi(kdir) = chi(kdir) + Levi_Civita(kdir,jdir,idir) * 2.0_dp * contrib
             ENDDO
             !
             CALL cp_dbcsr_sm_fm_multiply(op_mom_der_ao(ind_m2(jdir,iiiB),idir)%matrix,mo_coeff,&
                                    buf,ncol=nmo,alpha=1.e0_dp,error=error)
             DO kdir = 1,3
                IF(Levi_Civita(kdir,jdir,idir).EQ.0.0_dp) CYCLE
                CALL cp_fm_trace(buf,psi1_p(ispin,iiB)%matrix,contrib,error=error)
                chi(kdir) = chi(kdir) - Levi_Civita(kdir,jdir,idir) * 2.0_dp * contrib
             ENDDO
          ENDDO
          !
          !
          !
          ! term2: -2[C0| (r-dk)_ii  (r-dk)_iiB | d_iii(C1(piiiB))] +
          !        +2[C0| (r-dk)_iii (r-dk)_iiB | d_ii(C1(piiiB))]
          ! and
          ! term2: +2[C0| (r-dk)_ii  (r-dk)_iiiB | d_iii(C1(piiB))] +
          !        -2[C0| (r-dk)_iii (r-dk)_iiiB | d_ii(C1(piiB))]
          CALL cp_dbcsr_sm_fm_multiply(op_mom_ao(idir)%matrix,mo_coeff,&
                                 buf,ncol=nmo,alpha=1.e0_dp,error=error)
          DO jdir = 1,3
             DO kdir = 1,3
                IF(Levi_Civita(kdir,idir,jdir).EQ.0.0_dp) CYCLE
                IF(iiB==jdir) THEN
                   CALL cp_fm_trace(buf,psi1_p(ispin,iiiB)%matrix,contrib,error=error)
                   chi(kdir) = chi(kdir) + Levi_Civita(kdir,idir,jdir) * contrib
                ENDIF
             ENDDO
          ENDDO
          !
          DO jdir = 1,3
             DO kdir = 1,3
                IF(Levi_Civita(kdir,idir,jdir).EQ.0.0_dp) CYCLE
                IF(iiiB==jdir) THEN
                   CALL cp_fm_trace(buf,psi1_p(ispin,iiB)%matrix,contrib,error=error)
                   chi(kdir) = chi(kdir) - Levi_Civita(kdir,idir,jdir) * contrib
                ENDIF
                !
             ENDDO
          ENDDO
          !
          !
          !
          !
       ENDDO ! idir
       !
       DO idir = 1,3
          current_env%chi_tensor(idir,iB,ispin) = current_env%chi_tensor(idir,iB,ispin)+&
                                                  maxocc * chi(idir)
          IF (output_unit>0) THEN
             !WRITE(output_unit,'(A,E12.6)') ' chi_'//ACHAR(119+idir)//ACHAR(119+iB)//&
             !     &                         ' = ',maxocc * chi(idir)
          ENDIF
       ENDDO
       !
    ENDDO ! ispin
    !
    ! deallocate the sparse matrices
    CALL cp_dbcsr_deallocate_matrix_set(op_mom_ao,error=error)
    CALL cp_dbcsr_deallocate_matrix_set(op_mom_der_ao,error=error)
    CALL cp_dbcsr_deallocate_matrix_set(op_p_ao,error=error)
    CALL cp_fm_release(buf,error=error)

    CALL timestop(handle)

  END SUBROUTINE current_build_chi_one_center

END MODULE qs_linres_current
