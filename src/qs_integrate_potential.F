!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2003 CP2K developers group                           !
!-----------------------------------------------------------------------------!
!!****** cp2k/qs_integrate_potential [1.0] *
!!
!!   NAME
!!     qs_integrate_potential
!!
!!   FUNCTION
!!     Build up the plane wave density by collocating the primitive Gaussian
!!     functions (pgf).
!!
!!   AUTHOR
!!     Matthias Krack (03.04.2001)
!!
!!   MODIFICATION HISTORY
!!     Joost VandeVondele (02.2002)
!!           1) rewrote collocate_pgf for increased accuracy and speed
!!           2) collocate_core hack for PGI compiler
!!           3) added multiple grid feature
!!           4) new way to go over the grid
!!     Joost VandeVondele (05.2002)
!!           1) prelim. introduction of the real space grid type
!!     JGH [30.08.02] multigrid arrays independent from potential
!!     JGH [17.07.03] distributed real space code
!!     JGH [23.11.03] refactoring and new loop ordering
!!     JGH [04.12.03] OpneMP parallelization of main loops
!!     Joost VandeVondele (12.2003)
!!          1) modified to compute tau
!!     Joost removed incremental build feature
!!
!!   SOURCE
!***********************************************************************

MODULE qs_integrate_potential

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE coefficient_types,               ONLY: coeff_copy,&
                                             coeff_transform_space,&
                                             coeff_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_rs_pool_types,                ONLY: cp_rs_pool_p_type,&
                                             cp_rs_pool_type,&
                                             rs_pool_create_rs,&
                                             rs_pool_give_back_rs,&
                                             rs_pools_create_rs_vect,&
                                             rs_pools_give_back_rs_vect
  USE cube_utils,                      ONLY: cube_info_type,&
                                             return_cube
  USE gaussian_gridlevels,             ONLY: gaussian_gridlevel,&
                                             gridlevel_info_type
  USE kinds,                           ONLY: dp,&
                                             int_size
  USE l_utils,                         ONLY: l_info_type,&
                                             return_l_info
  USE memory_utilities,                ONLY: reallocate
  USE message_passing,                 ONLY: mp_max,&
                                             mp_shift
  USE orbital_pointers,                ONLY: coset,&
                                             ncoset
  USE particle_types,                  ONLY: particle_type
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_pool_types,                   ONLY: pw_pool_p_type,&
                                             pw_pools_give_back_coeffs,&
                                             pw_pools_init_coeffs
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_interactions,                 ONLY: exp_radius_very_extended
  USE qs_neighbor_list_types,          ONLY: first_list,&
                                             first_node,&
                                             get_neighbor_list,&
                                             get_neighbor_list_set,&
                                             get_neighbor_node,&
                                             neighbor_list_set_p_type,&
                                             neighbor_list_type,&
                                             neighbor_node_type,&
                                             next
  USE realspace_grid_types,            ONLY: realspace_grid_p_type,&
                                             realspace_grid_type,&
                                             rs_get_loop_vars,&
                                             rs_get_my_tasks,&
                                             rs_pw_transfer
  USE simulation_cell,                 ONLY: cell_type,&
                                             pbc
  USE sparse_matrix_types,             ONLY: add_block_node,&
                                             allocate_matrix,&
                                             cp_sm_get_block_list,&
                                             cp_sm_scale_and_add,&
                                             deallocate_matrix,&
                                             get_block_node,&
                                             real_matrix_p_type,&
                                             real_matrix_type,&
                                             set_matrix
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: get_limit

  IMPLICIT NONE

  PRIVATE

! *** Public subroutines ***

  INTERFACE integrate_v_rspace
      MODULE PROCEDURE integrate_v_rspace_light, integrate_v_rspace_full
  END INTERFACE

  PUBLIC :: integrate_v_rspace,&
            integrate_v_core_rspace,&
            integrate_pgf_product_rspace

!!***
! **********************************************************************

CONTAINS

! **********************************************************************

  SUBROUTINE integrate_v_core_rspace(v_rspace,cube_info,l_info,qs_env)

    TYPE(coeff_type), INTENT(INOUT)          :: v_rspace
    TYPE(cube_info_type), INTENT(IN)         :: cube_info
    TYPE(l_info_type), INTENT(IN)            :: l_info
    TYPE(qs_environment_type), INTENT(INOUT) :: qs_env

    CHARACTER(LEN=*), PARAMETER :: routine = &
      "SUBROUTINE integrate_v_core_rspace (MODULE qs_integrate_potential)"

    INTEGER                                  :: atom_a, dir, handle, &
                                                iatom, ierr, ikind, j, &
                                                natom_of_kind, npme
    INTEGER, DIMENSION(2)                    :: bo
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    INTEGER, DIMENSION(:, :), POINTER        :: tasks
    REAL(KIND=dp)                            :: alpha_core_charge, &
                                                ccore_charge, eps_rho_rspace
    REAL(KIND=dp), DIMENSION(3)              :: force_a, force_b, ra
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: hab, pab
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_rs_pool_type), POINTER           :: auxbas_rs_pool
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(realspace_grid_type), POINTER       :: rs_v

!   --------------------------------------------------------------------

    CALL timeset("integrate_v_core_rspace","I","",handle)

    NULLIFY(pw_env,auxbas_rs_pool,hab,pab,tasks)
    ALLOCATE(hab(1,1),pab(1,1),STAT=ierr)
    IF(ierr/=0) CALL stop_memory("integrate_v_core_rspace","hab,pab",1)

    CALL get_qs_env(qs_env=qs_env,pw_env=pw_env)
    CALL pw_env_get(pw_env=pw_env,auxbas_rs_pool=auxbas_rs_pool)
    CALL rs_pool_create_rs(auxbas_rs_pool,rs_v)

    CALL rs_pw_transfer(rs_v,v_rspace%pw,"BACKWARD")

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    cell=cell,&
                    dft_control=dft_control,&
                    particle_set=particle_set,&
                    para_env=para_env,pw_env=pw_env,&
                    force=force)

    eps_rho_rspace = dft_control%qs_control%eps_rho_rspace

    DO ikind=1,SIZE(atomic_kind_set)

      atomic_kind => atomic_kind_set(ikind)

      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           natom=natom_of_kind,&
                           atom_list=atom_list,&
                           alpha_core_charge=alpha_core_charge,&
                           ccore_charge=ccore_charge)

      pab(1,1) = -ccore_charge

      bo = get_limit ( natom_of_kind, para_env%num_pe, para_env%mepos )

      npme = bo(2) - bo(1) + 1
      CALL reallocate(tasks,1,2,1,npme)

      dir = rs_v%direction

      DO iatom=bo(1),bo(2)

        j = iatom - bo(1) + 1
        tasks ( 1, j ) = iatom
        IF ( dir /= 0) THEN
          atom_a = atom_list(iatom)
          ra(:) = pbc(particle_set(atom_a)%r,cell)
          tasks ( 2, j ) = FLOOR(ra(dir)/rs_v%dr(dir))
          tasks ( 2, j ) = MODULO ( tasks ( 2, j ), rs_v%npts(dir) )
          tasks ( 2, j ) = tasks ( 2, j ) + rs_v%lb(dir)
        END IF

      END DO

      CALL rs_get_my_tasks ( rs_v, tasks, npme )

      DO j=1,npme

        iatom = tasks(1,j)
        atom_a = atom_list(iatom)
        ra(:) = pbc(particle_set(atom_a)%r,cell)
        hab(1,1) = 0.0_dp
        force_a(:) = 0.0_dp
        force_b(:) = 0.0_dp

        CALL integrate_pgf_product_rspace(0,alpha_core_charge,0,&
             0,0.0_dp,0,ra,(/0.0_dp,0.0_dp,0.0_dp/),0.0_dp,&
             rs_v,cube_info,l_info,hab,pab=pab,o1=0,o2=0,&
             eps_gvg_rspace=eps_rho_rspace,&
             calculate_forces=.TRUE.,force_a=force_a,&
             force_b=force_b)

        force(ikind)%rho_core(:,iatom) =&
          force(ikind)%rho_core(:,iatom) + force_a(:)

      END DO

    END DO

    CALL rs_pool_give_back_rs(auxbas_rs_pool,rs_v)

    DEALLOCATE(hab,pab,tasks,STAT=ierr)
    IF(ierr /= 0) CALL stop_memory("integrate_v_core_rspace",&
                  "hab,pab,tasks")

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE integrate_v_core_rspace

! **********************************************************************
! The simplest possible interface to integrate_v_rspace 
! real space pw object in, sparse matrix out
! gets all the needed info from qs_env
! I guess one could use this style for the more general function as well
! **********************************************************************
  SUBROUTINE integrate_v_rspace_light(v_rspace,matrix_v,qs_env)

    TYPE(coeff_type)                         :: v_rspace
    TYPE(real_matrix_p_type), INTENT(INOUT)  :: matrix_v
    TYPE(qs_environment_type), INTENT(IN)    :: qs_env

    TYPE(coeff_type), DIMENSION(:), POINTER  :: mgrid_gspace, mgrid_rspace
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools

    NULLIFY(mgrid_gspace,mgrid_rspace,pw_pools,pw_env)

    CALL get_qs_env(qs_env=qs_env,pw_env=pw_env)
    ! set up all the necessary extra grids
    CALL pw_env_get(pw_env,pw_pools=pw_pools)
    ALLOCATE(mgrid_rspace(SIZE(pw_pools)), mgrid_gspace(SIZE(pw_pools)))
    CALL pw_pools_init_coeffs(pw_pools,mgrid_gspace,&
            use_data = COMPLEXDATA1D,&
            in_space = RECIPROCALSPACE)
    CALL pw_pools_init_coeffs(pw_pools,mgrid_rspace,&
            use_data = REALDATA3D,&
            in_space = REALSPACE)
    CALL set_matrix(matrix_v%matrix,0.0_dp)
    CALL integrate_v_rspace_full(v_rspace=v_rspace,&
                  mgrid_rspace=mgrid_rspace,&
                  mgrid_gspace=mgrid_gspace,&
                  auxbas_grid=pw_env%auxbas_grid,&
                  gridlevel_info=pw_env%gridlevel_info,&
                  cube_info=pw_env%cube_info,&
                  l_info=pw_env%l_info,&
                  h=matrix_v,&
                  qs_env=qs_env,&
                  calculate_forces=.FALSE.)

    CALL pw_pools_give_back_coeffs(pw_pools,mgrid_rspace)
    CALL pw_pools_give_back_coeffs(pw_pools,mgrid_gspace)
    DEALLOCATE(mgrid_rspace,mgrid_gspace)

  END SUBROUTINE integrate_v_rspace_light

! **********************************************************************
! integrates a given potential (or other object on a real 
! space grid) = v_rspace using a multi grid technique (mgrid_*)
! over the basis set producing a number for every element of h 
! (should have the same sparsity structure of S)
! additional screening is available using the magnitude of the 
! elements in p (? I'm not sure this is a very good idea)
! this argument is optional
! derivatives of these matrix elements with respect to the ionic 
! coordinates can be computed as well
! ***********************************************************************
 
  SUBROUTINE integrate_v_rspace_full(v_rspace,mgrid_rspace,mgrid_gspace,&
                                auxbas_grid,gridlevel_info,&
                                cube_info,l_info,&
                                p,h,qs_env,calculate_forces,compute_tau,error)

    TYPE(coeff_type)                         :: v_rspace
    TYPE(coeff_type), DIMENSION(:), POINTER  :: mgrid_rspace, mgrid_gspace
    INTEGER, INTENT(IN)                      :: auxbas_grid
    TYPE(gridlevel_info_type), INTENT(IN)    :: gridlevel_info
    TYPE(cube_info_type), DIMENSION(:), &
      POINTER                                :: cube_info
    TYPE(l_info_type), INTENT(IN)            :: l_info
    TYPE(real_matrix_p_type), INTENT(IN), &
      OPTIONAL                               :: p
    TYPE(real_matrix_p_type), INTENT(INOUT)  :: h
    TYPE(qs_environment_type), INTENT(IN)    :: qs_env
    LOGICAL, INTENT(IN)                      :: calculate_forces
    LOGICAL, INTENT(IN), OPTIONAL            :: compute_tau
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routine = &
      "SUBROUTINE integrate_v_rspace_full (MODULE qs_integrate_potential)"
    INTEGER, PARAMETER                       :: add_tasks = 1000, &
                                                max_tasks = 3000

    INTEGER :: ab, atom_a, atom_b, bcol, brow, curr_tasks, dir, handle, i, &
      iatom, igrid_level, ijatoms, ijsets, ikind, ilist, inode, ipgf, iset, &
      istat, itask, ithread, jatom, jkind, jpgf, jset, maxco, maxsgf_set, n, &
      na1, na2, nalimit(2), natom, natom_pairs, nb1, nb2, ncoa, ncob, nkind, &
      nlist, nnode, npme, nseta, nsetb, nthread, sgfa, sgfb, tp
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, lb_max, &
                                                lb_min, npgfa, npgfb, nsgfa, &
                                                nsgfb, ntasks
    INTEGER, DIMENSION(:, :), POINTER        :: asets, atasks, first_sgfa, &
                                                first_sgfb, ival, latom, &
                                                tasks_local
    INTEGER, DIMENSION(:, :, :), POINTER     :: tasks
    LOGICAL                                  :: h_duplicated, my_compute_tau, &
                                                p_duplicated
    REAL(KIND=dp)                            :: dab, eps_gvg_rspace, &
                                                kind_radius_b, maxvdiff, &
                                                maxvnew, rab2, scale,  &
                                                zetp
    REAL(KIND=dp), DIMENSION(3)              :: force_a, force_b, ra, rab, &
                                                rb, rp
    REAL(KIND=dp), DIMENSION(:), POINTER     :: set_radius_a, set_radius_b
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: dab_local, dp_block, h_block, &
                                                hab, p_block, pab, rpgfa, &
                                                rpgfb, sphi_a, sphi_b, work, &
                                                zeta, zetb
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: dist_ab, habt, pabt, workt
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_rs_pool_p_type), DIMENSION(:), &
      POINTER                                :: rs_pools
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb
    TYPE(neighbor_list_type), POINTER        :: sab_orb_neighbor_list
    TYPE(neighbor_node_type), POINTER        :: sab_orb_neighbor_node
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(real_matrix_type), POINTER          :: deltap, dh
    TYPE(realspace_grid_p_type), &
      DIMENSION(:), POINTER                  :: rs_v
    INTEGER :: omp_get_thread_num, omp_get_max_threads

!   ---------------------------------------------------------------------------
! default, if compute_tau is not present we assume it is false

    IF (PRESENT(compute_tau)) THEN
        my_compute_tau = compute_tau
    ELSE
        my_compute_tau = .FALSE.
    ENDIF

    IF (my_compute_tau) THEN
      CALL timeset("integrate_v_tau","I","",handle)
    ELSE
      CALL timeset("integrate_v_rspace","I","",handle)
    END IF

    NULLIFY(pw_env,rs_pools,ntasks,tasks,dist_ab,tasks_local,&
            ival,latom,dab_local,asets,atasks)

!   *** Get difference potential ***
    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    cell=cell,&
                    dft_control=dft_control,&
                    particle_set=particle_set,&
                    sab_orb=sab_orb,&
                    para_env=para_env,&
                    force=force,pw_env=pw_env)

!   *** Get the potential on the subgrids in real space, via fft ***

    CALL coeff_transform_space(v_rspace,mgrid_gspace(auxbas_grid))

    DO igrid_level=1,gridlevel_info%ngrid_levels

      IF ( igrid_level /= auxbas_grid )  &
           CALL coeff_copy(mgrid_gspace(auxbas_grid),mgrid_gspace(igrid_level))
      CALL coeff_transform_space(mgrid_gspace(igrid_level),&
                                 mgrid_rspace(igrid_level))

!     *** Multiply by the grid volume element ratio ***
      IF ( igrid_level /= auxbas_grid ) THEN
         scale = mgrid_rspace(igrid_level)%pw%pw_grid%dvol/&
                 mgrid_rspace(auxbas_grid)%pw%pw_grid%dvol
         mgrid_rspace(igrid_level)%pw%cr3d = &
                                   scale*mgrid_rspace(igrid_level)%pw%cr3d
      END IF

    END DO

    CALL pw_env_get(pw_env, rs_pools=rs_pools)
    CALL rs_pools_create_rs_vect(rs_pools, rs_v)
    DO igrid_level=1,gridlevel_info%ngrid_levels
       CALL rs_pw_transfer(rs_v(igrid_level)%rs_grid,&
                           mgrid_rspace(igrid_level)%pw,"BACKWARD")
    ENDDO

    nkind = SIZE(atomic_kind_set)
    natom = SIZE(particle_set)

    IF (calculate_forces) THEN
      ALLOCATE (atom_of_kind(natom),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"atom_of_kind",natom*int_size)
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                               atom_of_kind=atom_of_kind)
    END IF

    eps_gvg_rspace = dft_control%qs_control%eps_gvg_rspace


    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxco=maxco,&
                             maxsgf_set=maxsgf_set)

    h_duplicated = .FALSE.
    dh => h%matrix
    DO igrid_level = 1, gridlevel_info%ngrid_levels
      IF ( rs_v(igrid_level)%rs_grid%direction /= 0 ) THEN
        NULLIFY ( dh )
        CALL allocate_matrix(matrix=dh,&
                         nrow=h%matrix%nrow,&
                         ncol=h%matrix%ncol,&
                         nblock_row=h%matrix%nblock_row,&
                         nblock_col=h%matrix%nblock_col,&
                         first_row=h%matrix%first_row(:),&
                         last_row=h%matrix%last_row(:),&
                         first_col=h%matrix%first_col(:),&
                         last_col=h%matrix%last_col(:),&
                         matrix_name="LocalH",&
                         matrix_symmetry=h%matrix%symmetry)
        h_duplicated = .TRUE.
        EXIT
      END IF
    END DO
    p_duplicated = .FALSE.
    IF ( PRESENT ( p ) ) THEN
       deltap => p%matrix
       DO igrid_level = 1, gridlevel_info%ngrid_levels
         IF ( rs_v(igrid_level)%rs_grid%direction /= 0 ) THEN
           NULLIFY ( deltap )
           CALL allocate_matrix(matrix=deltap,&
                         nrow=p%matrix%nrow,&
                         ncol=p%matrix%ncol,&
                         nblock_row=p%matrix%nblock_row,&
                         nblock_col=p%matrix%nblock_col,&
                         first_row=p%matrix%first_row(:),&
                         last_row=p%matrix%last_row(:),&
                         first_col=p%matrix%first_col(:),&
                         last_col=p%matrix%last_col(:),&
                         matrix_name="LocalP",&
                         matrix_symmetry=p%matrix%symmetry)
           p_duplicated = .TRUE.
           EXIT
         END IF
       END DO
    END IF

    nthread = 1
!$  nthread = omp_get_max_threads()

!   *** Allocate work storage ***

    NULLIFY ( pabt, habt, workt )
    CALL reallocate(habt,1,maxco,1,maxco,0,nthread)
    CALL reallocate(workt,1,maxco,1,maxsgf_set,0,nthread)
    IF (PRESENT(p)) THEN
      CALL reallocate(pabt,1,maxco,1,maxco,0,nthread)
    ELSE
      IF (calculate_forces) CALL stop_program("integrate_v_rspace",&
                                        "need p for forces")
    ENDIF

    CALL reallocate(ntasks,1,gridlevel_info%ngrid_levels)
    CALL reallocate(tasks,1,8,1,max_tasks,1,gridlevel_info%ngrid_levels)
    CALL reallocate(dist_ab,1,3,1,max_tasks,1,gridlevel_info%ngrid_levels)
    CALL reallocate(tasks_local,1,2,1,max_tasks)
    CALL reallocate(ival,1,6,1,max_tasks)
    CALL reallocate(latom,1,2,1,max_tasks)
    CALL reallocate(dab_local,1,3,1,max_tasks)
    CALL reallocate(atasks,1,2,1,max_tasks)
    CALL reallocate(asets,1,2,1,max_tasks)
    curr_tasks = max_tasks

    DO ikind=1,nkind

      atomic_kind => atomic_kind_set(ikind)

      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           orb_basis_set=orb_basis_set)

      IF (.NOT.ASSOCIATED(orb_basis_set)) CYCLE

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

      DO jkind=1,nkind

        atomic_kind => atomic_kind_set(jkind)

        CALL get_atomic_kind(atomic_kind=atomic_kind,&
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
           nlist = 0
        END IF

        tasks = 0
        ntasks = 0

        DO ilist=1,nlist

          CALL get_neighbor_list(neighbor_list=sab_orb_neighbor_list,&
                                 atom=iatom,&
                                 nnode=nnode)

          ra(:) = pbc(particle_set(iatom)%r,cell)

          sab_orb_neighbor_node => first_node(sab_orb_neighbor_list)

          DO inode=1,nnode

            CALL get_neighbor_node(neighbor_node=sab_orb_neighbor_node,&
                                   neighbor=jatom,&
                                   r=rab(:))

            IF (iatom <= jatom) THEN
              CALL get_block_node(matrix=h%matrix,&
                                  block_row=iatom,&
                                  block_col=jatom,&
                                  block=h_block)
            ELSE
              CALL get_block_node(matrix=h%matrix,&
                                  block_row=jatom,&
                                  block_col=iatom,&
                                  block=h_block)
            END IF

!           *** Check, if the atomic block has to be ***
!           *** calculated by the current processor  ***

            IF (.NOT.ASSOCIATED(h_block)) THEN
              sab_orb_neighbor_node => next(sab_orb_neighbor_node)
              CYCLE
            END IF

            IF ( p_duplicated ) THEN
              IF (iatom <= jatom) THEN
                 NULLIFY ( dp_block )
                 CALL add_block_node ( deltap, iatom, jatom, dp_block )
                 CALL get_block_node(matrix=p%matrix,&
                                     block_row=iatom,&
                                     block_col=jatom,&
                                     block=p_block)
               ELSE
                 NULLIFY ( dp_block )
                 CALL add_block_node ( deltap, jatom, iatom, dp_block )
                 CALL get_block_node(matrix=p%matrix,&
                                     block_row=jatom,&
                                     block_col=iatom,&
                                     block=p_block)
               END IF
               dp_block = p_block
            ENDIF
 
             rb(:) = ra(:) + rab(:)
             rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
             dab = SQRT(rab2)
 
             DO iset=1,nseta
 
               IF (set_radius_a(iset) + kind_radius_b < dab) CYCLE
 
               ncoa = npgfa(iset)*ncoset(la_max(iset))
               sgfa = first_sgfa(1,iset)
 
               DO jset=1,nsetb
 
                 IF (set_radius_a(iset) + set_radius_b(jset) < dab) CYCLE
 
                 ncob = npgfb(jset)*ncoset(lb_max(jset))
                 sgfb = first_sgfb(1,jset)
    
                 DO ipgf=1,npgfa(iset)
 
                   IF (rpgfa(ipgf,iset) + set_radius_b(jset) < dab) CYCLE
 
                   na1 = (ipgf - 1)*ncoset(la_max(iset)) + 1
                   na2 = ipgf*ncoset(la_max(iset))
 
                   DO jpgf=1,npgfb(jset)
 
                     IF (rpgfa(ipgf,iset) + rpgfb(jpgf,jset) < dab) CYCLE
 
                     nb1 = (jpgf - 1)*ncoset(lb_max(jset)) + 1
                     nb2 = jpgf*ncoset(lb_max(jset))
 
                     zetp = zeta(ipgf,iset) + zetb(jpgf,jset)
                     igrid_level = gaussian_gridlevel(gridlevel_info,zetp)

                     ntasks(igrid_level) = ntasks(igrid_level) + 1
                     n = ntasks(igrid_level)
                     IF ( n > curr_tasks ) THEN
                       curr_tasks = curr_tasks + add_tasks
                       CALL reallocate ( tasks, 1, 8, 1, curr_tasks, &
                                         1, gridlevel_info%ngrid_levels )
                       CALL reallocate ( dist_ab, 1, 3, 1, curr_tasks, &
                                         1, gridlevel_info%ngrid_levels )
                     END IF
                     dir = rs_v(igrid_level)%rs_grid%direction
                     tasks (1,n,igrid_level) = n
                     IF ( dir /= 0) THEN
                       rp(:) = ra(:) + zetb(jpgf,jset)/zetp*rab(:)
                       rp(:) = pbc(rp,cell)
                       tp = FLOOR(rp(dir)/rs_v(igrid_level)%rs_grid%dr(dir))
                       tp = MODULO ( tp, rs_v(igrid_level)%rs_grid%npts(dir) )
                       tasks (2,n,igrid_level) = tp + &
                          rs_v(igrid_level)%rs_grid%lb(dir)
                     END IF
                     tasks (3,n,igrid_level) = iatom
                     tasks (4,n,igrid_level) = jatom
                     tasks (5,n,igrid_level) = iset
                     tasks (6,n,igrid_level) = jset
                     tasks (7,n,igrid_level) = ipgf
                     tasks (8,n,igrid_level) = jpgf
                     dist_ab (:,n,igrid_level) = rab(:)

                   END DO
                 END DO
               END DO
             END DO

             sab_orb_neighbor_node => next(sab_orb_neighbor_node)

          END DO

          sab_orb_neighbor_list => next(sab_orb_neighbor_list)

        END DO

        DO igrid_level = 1, gridlevel_info%ngrid_levels
          n = ntasks ( igrid_level )
          IF ( n > SIZE ( tasks_local, 2 ) ) &
            CALL reallocate(tasks_local,1,2,1,n)
          IF ( n > SIZE ( ival, 2 ) ) &
            CALL reallocate(ival,1,6,1,n)
          IF ( n > SIZE ( dab_local, 2 ) ) &
            CALL reallocate(dab_local,1,3,1,n)
!$OMP parallel do private(i)
          DO i=1,n
            tasks_local(1,i) = tasks(1,i,igrid_level)
            tasks_local(2,i) = tasks(2,i,igrid_level)
            ival(1,i) = tasks(3,i,igrid_level)
            ival(2,i) = tasks(4,i,igrid_level)
            ival(3,i) = tasks(5,i,igrid_level)
            ival(4,i) = tasks(6,i,igrid_level)
            ival(5,i) = tasks(7,i,igrid_level)
            ival(6,i) = tasks(8,i,igrid_level)
            dab_local(1,i) = dist_ab(1,i,igrid_level)
            dab_local(2,i) = dist_ab(2,i,igrid_level)
            dab_local(3,i) = dist_ab(3,i,igrid_level)
          END DO
!$OMP parallel do private(i)
          DO i=n+1,SIZE(tasks_local,2)
            tasks_local(1,i)=0
            tasks_local(2,i)=0
          END DO

          npme = 0
          IF (PRESENT(P)) THEN
            IF ( n > SIZE ( latom, 2 ) ) CALL reallocate(latom,1,2,1,n)
!$OMP parallel do private(i)
            DO i=1,n
              latom(1,i) = tasks(3,i,igrid_level)
              latom(2,i) = tasks(4,i,igrid_level)
            END DO
            CALL rs_get_my_tasks ( rs_v(igrid_level)%rs_grid, tasks_local, &
                                   npme, ival=ival, rval=dab_local, &
                                   pmat=deltap, pcor=latom )
          ELSE
            CALL rs_get_my_tasks ( rs_v(igrid_level)%rs_grid, tasks_local, &
                                   npme, ival=ival, rval=dab_local )
          END IF
          CALL rs_get_loop_vars ( npme, ival, natom_pairs, asets, atasks )

!$OMP parallel &
!$OMP default(none) &
!$OMP private(ijatoms,ithread,nalimit,itask,iatom,jatom,ra,brow,bcol) &
!$OMP private(h_block,p_block,atom_a,atom_b,ijsets,iset,jset) &
!$OMP private(ncoa,ncob,sgfa,sgfb,work,pab,hab,rab,rb,rab2) &
!$OMP private(ipgf,jpgf,na1,na2,nb1,nb2,force_a,force_b,istat) &
!$OMP shared(natom_pairs,nthread,asets,atasks,particle_set,cell,dh,p) &
!$OMP shared(deltap,calculate_forces,atom_of_kind,ival,npgfa,npgfb,la_max,lb_max) &
!$OMP shared(first_sgfa,first_sgfb,sphi_a,sphi_b,maxco,maxsgf_set,dab_local) &
!$OMP shared(ncoset,zeta,zetb,la_min,lb_min,rs_v,cube_info,l_info,eps_gvg_rspace) &
!$OMP shared(force,ikind,jkind,igrid_level,nsgfa,nsgfb,pabt,habt,workt,my_compute_tau)
          ithread = 0
!$        ithread = omp_get_thread_num()
          nalimit = get_limit(natom_pairs,nthread,ithread)
          work => workt(:,:,ithread)
          hab => habt(:,:,ithread)
          IF (PRESENT(p)) THEN
            pab => pabt(:,:,ithread)
          ENDIF
!         DO ijatoms = nalimit(1),nalimit(2)
!$OMP do
          DO ijatoms = 1,natom_pairs
            itask = atasks(1,asets(1,ijatoms))
            iatom  = ival (1,itask)
            jatom  = ival (2,itask)
            ra(:) = pbc(particle_set(iatom)%r,cell)
            IF (iatom <= jatom) THEN
              brow = iatom
              bcol = jatom
            ELSE
              brow = jatom
              bcol = iatom
            END IF
            CALL get_block_node(matrix=dh,&
                                block_row=brow,&
                                block_col=bcol,&
                                block=h_block)
            IF (.NOT.ASSOCIATED(h_block)) &
               CALL add_block_node ( dh, brow, bcol, h_block )
            IF (PRESENT(P)) THEN
              CALL get_block_node(matrix=deltap,&
                                  block_row=brow,&
                                  block_col=bcol,&
                                  BLOCK=p_block)
              IF (.NOT.ASSOCIATED(p_block)) &
                 CALL stop_program(routine,"p_block not associated in deltap")
            END IF
            IF (calculate_forces) THEN
               atom_a = atom_of_kind(iatom)
               atom_b = atom_of_kind(jatom)
               force_a(:) = 0.0_dp
               force_b(:) = 0.0_dp
            END IF
            DO ijsets = asets(1,ijatoms), asets(2,ijatoms)
              itask = atasks(1,ijsets)
              iset   = ival (3,itask)
              jset   = ival (4,itask)
              ncoa = npgfa(iset)*ncoset(la_max(iset))
              sgfa = first_sgfa(1,iset)
              ncob = npgfb(jset)*ncoset(lb_max(jset))
              sgfb = first_sgfb(1,jset)
              IF (PRESENT(P)) THEN
                IF (iatom <= jatom) THEN
                  CALL dgemm("N","N",ncoa,nsgfb(jset),nsgfa(iset),&
                              1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                              p_block(sgfa,sgfb),SIZE(p_block,1),&
                              0.0_dp,work(1,1),SIZE(work,1))
                  CALL dgemm("N","T",ncoa,ncob,nsgfb(jset),&
                              1.0_dp,work(1,1),SIZE(work,1),&
                              sphi_b(1,sgfb),SIZE(sphi_b,1),&
                              0.0_dp,pab(1,1),SIZE(pab,1))
                ELSE
                  CALL dgemm("N","N",ncob,nsgfa(iset),nsgfb(jset),&
                             1.0_dp,sphi_b(1,sgfb),SIZE(sphi_b,1),&
                             p_block(sgfb,sgfa),SIZE(p_block,1),&
                             0.0_dp,work(1,1),SIZE(work,1))
                  CALL dgemm("N","T",ncob,ncoa,nsgfa(iset),&
                             1.0_dp,work(1,1),SIZE(work,1),&
                             sphi_a(1,sgfa),SIZE(sphi_a,1),&
                             0.0_dp,pab(1,1),SIZE(pab,1))
                END IF
              END IF
              hab = 0._dp

              DO itask = atasks(1,ijsets),atasks(2,ijsets)

                rab(1) = dab_local (1,itask)
                rab(2) = dab_local (2,itask)
                rab(3) = dab_local (3,itask)
                rab2  = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
                rb(1) = ra(1) + rab(1)
                rb(2) = ra(2) + rab(2)
                rb(3) = ra(3) + rab(3)
                ipgf   = ival (5,itask)
                jpgf   = ival (6,itask)
                na1 = (ipgf - 1)*ncoset(la_max(iset)) + 1
                na2 = ipgf*ncoset(la_max(iset))
                nb1 = (jpgf - 1)*ncoset(lb_max(jset)) + 1
                nb2 = jpgf*ncoset(lb_max(jset))
                IF (PRESENT(P)) THEN
                   IF (iatom <= jatom) THEN
                     CALL integrate_pgf_product_rspace(&
                          la_max(iset),zeta(ipgf,iset),la_min(iset),&
                          lb_max(jset),zetb(jpgf,jset),lb_min(jset),&
                          ra,rab,rab2,rs_v(igrid_level)%rs_grid,&
                          cube_info(igrid_level),l_info,&
                          hab,pab=pab,o1=na1-1,o2=nb1-1, &
                          eps_gvg_rspace=eps_gvg_rspace,&
                          calculate_forces=calculate_forces,&
                          force_a=force_a,force_b=force_b,ithread=ithread,compute_tau=my_compute_tau)
                   ELSE
                     CALL integrate_pgf_product_rspace(&
                          lb_max(jset),zetb(jpgf,jset),lb_min(jset),&
                          la_max(iset),zeta(ipgf,iset),la_min(iset),&
                          rb,-rab,rab2,rs_v(igrid_level)%rs_grid,&
                          cube_info(igrid_level),l_info,&
                          hab,pab=pab,o1=nb1-1,o2=na1-1, &
                          eps_gvg_rspace=eps_gvg_rspace,&
                          calculate_forces=calculate_forces,&
                          force_a=force_b,force_b=force_a,ithread=ithread,compute_tau=my_compute_tau)
                   END IF
                ELSE
                   IF (iatom <= jatom) THEN
                     CALL integrate_pgf_product_rspace(&
                          la_max(iset),zeta(ipgf,iset),la_min(iset),&
                          lb_max(jset),zetb(jpgf,jset),lb_min(jset),&
                          ra,rab,rab2,rs_v(igrid_level)%rs_grid,&
                          cube_info(igrid_level),l_info,&
                          hab,o1=na1-1,o2=nb1-1,&
                          eps_gvg_rspace=eps_gvg_rspace,&
                          calculate_forces=calculate_forces,&
                          force_a=force_a,force_b=force_b,ithread=ithread,compute_tau=my_compute_tau)
                   ELSE
                     CALL integrate_pgf_product_rspace(&
                          lb_max(jset),zetb(jpgf,jset),lb_min(jset),&
                          la_max(iset),zeta(ipgf,iset),la_min(iset),&
                          rb,-rab,rab2,rs_v(igrid_level)%rs_grid,&
                          cube_info(igrid_level),l_info,&
                          hab,o1=nb1-1,o2=na1-1,&
                          eps_gvg_rspace=eps_gvg_rspace,&
                          calculate_forces=calculate_forces,&
                          force_a=force_b,force_b=force_a,ithread=ithread,compute_tau=my_compute_tau)
                   END IF
                END IF

              END DO

              IF (iatom <= jatom) THEN
                CALL dgemm("N","N",ncoa,nsgfb(jset),ncob,&
                           1.0_dp,hab(1,1),SIZE(hab,1),&
                           sphi_b(1,sgfb),SIZE(sphi_b,1),&
                           0.0_dp,work(1,1),SIZE(work,1))
                CALL dgemm("T","N",nsgfa(iset),nsgfb(jset),ncoa,&
                           1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                           work(1,1),SIZE(work,1),&
                           1.0_dp,h_block(sgfa,sgfb),SIZE(h_block,1))
              ELSE
                CALL dgemm("N","N",ncob,nsgfa(iset),ncoa,&
                           1.0_dp,hab(1,1),SIZE(hab,1),&
                           sphi_a(1,sgfa),SIZE(sphi_a,1),&
                           0.0_dp,work(1,1),SIZE(work,1))
                CALL dgemm("T","N",nsgfb(jset),nsgfa(iset),ncob,&
                           1.0_dp,sphi_b(1,sgfb),SIZE(sphi_b,1),&
                           work(1,1),SIZE(work,1),&
                           1.0_dp,h_block(sgfb,sgfa),SIZE(h_block,1))
              END IF

            END DO

            IF (calculate_forces) THEN
!$OMP critical (qs_integrate_force)
               force(ikind)%rho_elec(:,atom_a) =&
                  force(ikind)%rho_elec(:,atom_a) + 2.0_dp*force_a(:)
               IF (iatom /= jatom) THEN
                  force(jkind)%rho_elec(:,atom_b) =&
                    force(jkind)%rho_elec(:,atom_b) + 2.0_dp*force_b(:)
               END IF
!$OMP end critical (qs_integrate_force)
            END IF

          END DO
!$OMP end parallel

        END DO
    
      END DO

    END DO

!   *** Release work storage ***

    IF ( h_duplicated ) THEN
      ! Reconstruct H matrix 
      CALL redistribute_matrix ( dh, h%matrix, para_env )
      CALL deallocate_matrix ( dh )
    ELSE
      NULLIFY ( dh )
    END IF

    IF ( p_duplicated ) THEN
      CALL deallocate_matrix ( deltap )
    ELSE
      NULLIFY ( deltap )
    END IF

    DEALLOCATE (habt,workt,ntasks,tasks,tasks_local,ival,latom,&
                dist_ab,dab_local,asets,atasks,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"habt,workt,ntasks,tasks,"//&
       "tasks_local,ival,latom,dist_ab,dab_local,asets,atasks")
    IF ( PRESENT(p) ) THEN
      DEALLOCATE (pabt,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"pabt")
    END IF

    CALL rs_pools_give_back_rs_vect(rs_pools, rs_v)

    IF (calculate_forces) THEN
      DEALLOCATE (atom_of_kind,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"atom_of_kind")
    END IF

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE integrate_v_rspace_full

! *****************************************************************************

    SUBROUTINE integrate_pgf_product_rspace(la_max,zeta,la_min,&
                                            lb_max,zetb,lb_min,&
                                            ra,rab,rab2,rsgrid,&
                                            cube_info,l_info,hab,pab,o1,o2,&
                                            eps_gvg_rspace,&
                                            calculate_forces,force_a,force_b,&
                                            ithread,mytimings,compute_tau)

    INTEGER, INTENT(IN)                      :: la_max
    REAL(KIND=dp), INTENT(IN)                :: zeta
    INTEGER, INTENT(IN)                      :: la_min, lb_max
    REAL(KIND=dp), INTENT(IN)                :: zetb
    INTEGER, INTENT(IN)                      :: lb_min
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: ra, rab
    REAL(KIND=dp), INTENT(IN)                :: rab2
    TYPE(realspace_grid_type), POINTER       :: rsgrid
    TYPE(cube_info_type), INTENT(IN)         :: cube_info
    TYPE(l_info_type), INTENT(IN)            :: l_info
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: hab
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: pab
    INTEGER, INTENT(IN)                      :: o1, o2
    REAL(KIND=dp), INTENT(IN)                :: eps_gvg_rspace
    LOGICAL, INTENT(IN)                      :: calculate_forces
    REAL(KIND=dp), DIMENSION(3), &
      INTENT(INOUT)                          :: force_a, force_b
    INTEGER, OPTIONAL                        :: ithread
    REAL(KIND=dp), DIMENSION(10), OPTIONAL   :: mytimings
    LOGICAL, INTENT(IN), OPTIONAL            :: compute_tau

    INTEGER :: ax, ay, az, bx, by, bz, cmax, coef_max, gridbounds(2,3), i, &
      ico, icoef, ig, ithread_l, jco, k, l, la, la_max_l, la_min_l, lb, &
      lb_cube_min, lb_max_l, lb_min_l, length, lx, lx_max, lxa, lxb, lxy, &
      lxy_max, lxyz, lxyz_max, lya, lyb, lza, lzb, offset, start, ub_cube_max
    INTEGER, DIMENSION(3)                    :: cubecenter, lb_cube, ng, &
                                                ub_cube
    INTEGER, DIMENSION(:), POINTER           :: ly_max, lz_max, sphere_bounds
    INTEGER, DIMENSION(:, :), POINTER        :: map
    INTEGER, POINTER                         :: ipzyx(:,:,:,:,:,:)
    LOGICAL                                  :: my_compute_tau
    REAL(KIND=dp) :: a, axpm0, b, binomial_k_lxa, binomial_l_lxb, cutoff, &
      exp_x0, exp_x1, exp_x2, f, ftza, ftzb, pabval, pg, prefactor, radius, &
      rpg, ya, yap, yb, ybp, za, zap, zb, zbp, zetp
    REAL(KIND=dp), DIMENSION(3)              :: dr, rap, rb, rbp, roffset, rp
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: grid
    REAL(KIND=dp), POINTER                   :: alpha(:,:), dpy(:,:), &
                                                dpz(:,:), polx(:,:), &
                                                poly(:,:), polz(:,:), pzyx(:)

!   ---------------------------------------------------------------------------

    IF (PRESENT(ithread)) THEN
       ithread_l=ithread
    ELSE
       ithread_l=0
    ENDIF

    ! my_compute_tau defaults to .FALSE.
    ! IF (.true.) it will compute 0.5 * (nabla x_a).(v(r) nabla x_b)
    IF (PRESENT(compute_tau)) THEN
       my_compute_tau=compute_tau
    ELSE
       my_compute_tau=.FALSE.
    ENDIF

    !IF (ithread_l.eq.0) t_a_1=m_walltime()

    IF (calculate_forces) THEN
      la_max_l=la_max+1  ! needed for the derivative of the gaussian, unimportant which one
      la_min_l=MAX(la_min-1,0) ! just in case the la_min,lb_min is not zero
      lb_min_l=MAX(lb_min-1,0)
      lb_max_l=lb_max
    ELSE
      la_max_l=la_max
      la_min_l=la_min
      lb_min_l=lb_min
      lb_max_l=lb_max
    END IF

    IF (my_compute_tau) THEN
      la_max_l=la_max_l+1
      lb_max_l=lb_max_l+1
      la_min_l=MAX(la_min_l-1,0)
      lb_min_l=MAX(lb_min_l-1,0)
    ENDIF

    coef_max=la_max_l+lb_max_l+1
    zetp = zeta + zetb
    f = zetb/zetp
    prefactor = EXP(-zeta*f*rab2)
!   *** position of the gaussian product
    rap(:) = f*rab(:)
    rbp(:) = rap(:) - rab(:)
    rp(:) = ra(:) + rap(:)  ! this is the gaussian center in real coordinates
    rb(:) = ra(:) + rab(:)

    cutoff=1.0_dp
    IF (PRESENT(pab)) THEN
    radius=exp_radius_very_extended(la_min,la_max,lb_min,lb_max,pab,o1,o2,ra,rb,rp,&
                                    zetp,eps_gvg_rspace,prefactor,cutoff)
    ELSE
    radius=exp_radius_very_extended(la_min,la_max,lb_min,lb_max,ra=ra,rb=rb,rp=rp,&
            zetp=zetp,eps=eps_gvg_rspace,prefactor=prefactor,cutoff=cutoff)
    ENDIF

    IF (radius == 0.0_dp) THEN
       !IF (ithread_l.EQ.0) t_a_2=m_walltime()
       !IF (ithread_l.EQ.0) mytimings(1)=mytimings(1)+(t_a_2-t_a_1)
       RETURN
    ENDIF
    CALL return_cube(cube_info,radius,lb_cube,ub_cube,sphere_bounds)
    CALL return_l_info(l_info,la_min_l,la_max_l,lb_min_l,lb_max_l,ithread_l,lx_max,lxy_max, &
                         lxyz_max,ly_max,lz_max,&
                         map,polx,poly,polz,dpy,dpz,alpha,pzyx,ipzyx,cmax)

!   *** properties of the grid ***
    !IF (ithread_l.eq.0) t_c_1=m_walltime()

    dr(:) = rsgrid%dr(:)
    ng(:) = rsgrid%npts(:)

    grid => rsgrid%r(:,:,:)
    gridbounds(1,1)=LBOUND(GRID,1)
    gridbounds(2,1)=UBOUND(GRID,1)
    gridbounds(1,2)=LBOUND(GRID,2)
    gridbounds(2,2)=UBOUND(GRID,2)
    gridbounds(1,3)=LBOUND(GRID,3)
    gridbounds(2,3)=UBOUND(GRID,3)
    !gridbounds(1,1)=LBOUND(rsgrid%r,1)
    !gridbounds(2,1)=UBOUND(rsgrid%r,1)
    !gridbounds(1,2)=LBOUND(rsgrid%r,2)
    !gridbounds(2,2)=UBOUND(rsgrid%r,2)
    !gridbounds(1,3)=LBOUND(rsgrid%r,3)
    !gridbounds(2,3)=UBOUND(rsgrid%r,3)

    cubecenter(:) = FLOOR(rp(:)/dr(:))
    roffset(:) = rp(:) - REAL(cubecenter(:),dp)*dr(:)
    lb_cube_min = MINVAL(lb_cube(:))
    ub_cube_max = MAXVAL(ub_cube(:))

!   *** a mapping so that the ig corresponds to the right grid point, also with pbc
    DO i=1,3
      IF ( rsgrid % perd ( i ) == 1 ) THEN
        start=lb_cube(i)
        DO
         offset=MODULO(cubecenter(i)+start,ng(i))+1-start
         length=MIN(ub_cube(i),ng(i)-offset)-start
         DO ig=start,start+length
            map(ig,i) = ig+offset
         END DO
         IF (start+length.GE.ub_cube(i)) EXIT
         start=start+length+1
        END DO
      ELSE
        ! this takes partial grid + border regions into account
        offset=MODULO(cubecenter(i),ng(i))+rsgrid%lb(i)
        offset=offset-rsgrid%lb_local(i)+1
        DO ig=lb_cube(i),ub_cube(i)
           map(ig,i) = ig+offset
        END DO
      END IF
    ENDDO

!   *** initialise the pol x,y,z terms
!  new trick to avoid to many exps (reuse the result from the previous gridpoint): 
!  exp( -a*(x+d)**2)=exp(-a*x**2)*(-2*a*x*d)*exp(-a*d**2)
!  exp(-2*a*(x+d)*d)=exp(-2*a*x*d)*exp(-2*a*d**2)

    rpg=REAL(lb_cube(3)-1,dp)*dr(3)-roffset(3)
    exp_x2=exp(-zetp*rpg**2)
    exp_x1=exp(-2*zetp*rpg*dr(3))
    exp_x0=exp(-zetp*dr(3)**2)

    DO ig=lb_cube(3),ub_cube(3)
      lxyz=0
      rpg = rpg + dr(3)
      ! zap = EXP(-zetp*rpg**2)
      exp_x1 = exp_x1 * exp_x0
      exp_x2 = exp_x2 * exp_x1
      exp_x1 = exp_x1 * exp_x0
      zap = exp_x2
      za  = rpg + rap(3)
      zb  = za  - rab(3)
      DO lza=0,la_max_l
       zbp=1.0_dp
       DO lzb=0,lb_max_l
          dpz(lzb,lza)=zap*zbp
          zbp=zbp*zb
       ENDDO
       zap=zap*za
      ENDDO
      DO lxa=0,la_max_l
      DO lxb=0,lb_max_l
       DO lya=0,la_max_l-lxa
       DO lyb=0,lb_max_l-lxb
          DO lza=MAX(la_min_l-lxa-lya,0),la_max_l-lxa-lya
          DO lzb=MAX(lb_min_l-lxb-lyb,0),lb_max_l-lxb-lyb
             lxyz=lxyz+1
             polz(lxyz,ig)=dpz(lzb,lza)
          ENDDO
          ENDDO
       ENDDO
       ENDDO
      ENDDO
      ENDDO
    ENDDO

    rpg=REAL(lb_cube(2)-1,dp)*dr(2)-roffset(2)
    exp_x2=exp(-zetp*rpg**2)
    exp_x1=exp(-2*zetp*rpg*dr(2))
    exp_x0=exp(-zetp*dr(2)**2)
    DO ig=lb_cube(2),ub_cube(2)
      rpg = rpg + dr(2)
      ! yap = EXP(-zetp*rpg**2)
      exp_x1 = exp_x1 * exp_x0
      exp_x2 = exp_x2 * exp_x1
      exp_x1 = exp_x1 * exp_x0
      yap = exp_x2
      ya  = rpg + rap(2)
      yb  = ya  - rab(2)
      DO lya=0,la_max_l
       ybp=1.0_dp
       DO lyb=0,lb_max_l
          dpy(lyb,lya)=yap*ybp
          ybp=ybp*yb
       ENDDO
       yap=yap*ya
      ENDDO

      lxy=0
      DO lxa=0,la_max_l
      DO lxb=0,lb_max_l
       DO lya=0,la_max_l-lxa
       DO lyb=0,lb_max_l-lxb
             lxy=lxy+1
             poly(lxy,ig)=dpy(lyb,lya)
       ENDDO
       ENDDO
      ENDDO
      ENDDO
    ENDDO


!   *** make the alpha matrix ***

    alpha(:,:)=0.0_dp
    lx=0
    DO lxa=0,la_max_l
    DO lxb=0,lb_max_l
       lx=lx+1
       binomial_k_lxa=1.0_dp
       a=1.0_dp
       DO k=0,lxa
        binomial_l_lxb=1.0_dp
        b=1.0_dp
        DO l=0,lxb
           alpha(lxa-l+lxb-k+1,lx)=alpha(lxa-l+lxb-k+1,lx)+ &
                             binomial_k_lxa*binomial_l_lxb*a*b
           binomial_l_lxb=binomial_l_lxb*REAL(lxb-l,dp)/REAL(l+1,dp)
           b=b*(rp(1)-(ra(1)+rab(1)))
        ENDDO
        binomial_k_lxa=binomial_k_lxa*REAL(lxa-k,dp)/REAL(k+1,dp)
        a=a*(rp(1)-ra(1))
       ENDDO
    ENDDO
    ENDDO

    rpg=REAL(lb_cube(1)-1,dp)*dr(1)-roffset(1)
    exp_x2=exp(-zetp*rpg**2)
    exp_x1=exp(-2*zetp*rpg*dr(1))
    exp_x0=exp(-zetp*dr(1)**2)
    DO ig=lb_cube(1),ub_cube(1)
      rpg = rpg + dr(1)
      ! pg  = EXP(-zetp*rpg**2)
      exp_x1 = exp_x1 * exp_x0
      exp_x2 = exp_x2 * exp_x1
      exp_x1 = exp_x1 * exp_x0
      pg = exp_x2
      DO icoef=1,coef_max
         polx(icoef,ig)=pg
         pg=pg*(rpg)
      ENDDO
    ENDDO

!   *** do the loop over the grid

    !IF (ithread_l.EQ.0)t_b_1=m_walltime()
    CALL integrate_core(polx(1,-cmax),poly(1,-cmax),polz(1,-cmax),grid(1,1,1),alpha(1,1), &
                        lx_max,lxy_max,lxyz_max,coef_max,cmax, &
                        ly_max(1),lz_max(1),gridbounds(1,1),map(-cmax,1),pzyx(1), &
                        la_min_l,la_max_l,lb_min_l,lb_max_l,sphere_bounds(1))
    !IF (ithread_l.EQ.0)t_b_2=m_walltime()

    DO i=1,lxyz_max
       pzyx(i)=pzyx(i)*prefactor
    ENDDO

!   *** pzyx contains all the information needed to find the elements of hab
!   *** and optionally of derivatives of these elements

    ftza = 2.0_dp*zeta
    ftzb = 2.0_dp*zetb

    DO la=la_min,la_max
      DO ax=0,la
        DO ay=0,la-ax
          az = la - ax - ay
          ico=coset(ax,ay,az)
          DO lb=lb_min,lb_max
            DO bx=0,lb
              DO by=0,lb-bx
                bz = lb - bx - by
                jco=coset(bx,by,bz)
                IF (my_compute_tau) THEN
                    axpm0 =  0.5_dp * ( ax * bx * pzyx(ipzyx(MAX(ax-1,0),ay,az,MAX(bx-1,0),by,bz)) +  &
                                        ay * by * pzyx(ipzyx(ax,MAX(ay-1,0),az,bx,MAX(by-1,0),bz)) +  &
                                        az * bz * pzyx(ipzyx(ax,ay,MAX(az-1,0),bx,by,MAX(bz-1,0)))  &
                                        - ftza * bx * pzyx(ipzyx(ax+1,ay,az,MAX(bx-1,0),by,bz))  &
                                        - ftza * by * pzyx(ipzyx(ax,ay+1,az,bx,MAX(by-1,0),bz))  &
                                        - ftza * bz * pzyx(ipzyx(ax,ay,az+1,bx,by,MAX(bz-1,0)))  &
                                        - ax * ftzb * pzyx(ipzyx(MAX(ax-1,0),ay,az,bx+1,by,bz))  &
                                        - ay * ftzb * pzyx(ipzyx(ax,MAX(ay-1,0),az,bx,by+1,bz))  &
                                        - az * ftzb * pzyx(ipzyx(ax,ay,MAX(az-1,0),bx,by,bz+1)) +  &
                                        ftza * ftzb * pzyx(ipzyx(ax+1,ay,az,bx+1,by,bz)) + &
                                        ftza * ftzb * pzyx(ipzyx(ax,ay+1,az,bx,by+1,bz)) + &
                                        ftza * ftzb * pzyx(ipzyx(ax,ay,az+1,bx,by,bz+1)) )
                ELSE
                    axpm0 = pzyx(ipzyx(ax,ay,az,bx,by,bz))
                ENDIF
                hab(o1+ico,o2+jco) = hab(o1+ico,o2+jco) + axpm0
                IF (calculate_forces) THEN
                  IF (my_compute_tau) THEN
                     pabval=pab(o1+ico,o2+jco)*0.5_dp*ax * bx
                     CALL force_update(force_a,force_b,rab,pabval,ftza,ftzb,MAX(ax-1,0),ay,az,MAX(bx-1,0),by,bz,pzyx,ipzyx)
                     pabval=pab(o1+ico,o2+jco)*0.5_dp*ay * by
                     CALL force_update(force_a,force_b,rab,pabval,ftza,ftzb,ax,MAX(ay-1,0),az,bx,MAX(by-1,0),bz,pzyx,ipzyx)
                     pabval=pab(o1+ico,o2+jco)*0.5_dp*az * bz
                     CALL force_update(force_a,force_b,rab,pabval,ftza,ftzb,ax,ay,MAX(az-1,0),bx,by,MAX(bz-1,0),pzyx,ipzyx)
                     pabval=pab(o1+ico,o2+jco)*0.5_dp*(- ftza * bx )
                     CALL force_update(force_a,force_b,rab,pabval,ftza,ftzb,ax+1,ay,az,MAX(bx-1,0),by,bz ,pzyx,ipzyx)
                     pabval=pab(o1+ico,o2+jco)*0.5_dp*(- ftza * by )
                     CALL force_update(force_a,force_b,rab,pabval,ftza,ftzb,ax,ay+1,az,bx,MAX(by-1,0),bz ,pzyx,ipzyx)
                     pabval=pab(o1+ico,o2+jco)*0.5_dp*(- ftza * bz  )
                     CALL force_update(force_a,force_b,rab,pabval,ftza,ftzb,ax,ay,az+1,bx,by,MAX(bz-1,0) ,pzyx,ipzyx)
                     pabval=pab(o1+ico,o2+jco)*0.5_dp*(- ax * ftzb )
                     CALL force_update(force_a,force_b,rab,pabval,ftza,ftzb,MAX(ax-1,0),ay,az,bx+1,by,bz ,pzyx,ipzyx)
                     pabval=pab(o1+ico,o2+jco)*0.5_dp*(- ay * ftzb )
                     CALL force_update(force_a,force_b,rab,pabval,ftza,ftzb,ax,MAX(ay-1,0),az,bx,by+1,bz ,pzyx,ipzyx)
                     pabval=pab(o1+ico,o2+jco)*0.5_dp*(- az * ftzb )
                     CALL force_update(force_a,force_b,rab,pabval,ftza,ftzb,ax,ay,MAX(az-1,0),bx,by,bz+1 ,pzyx,ipzyx)
                     pabval=pab(o1+ico,o2+jco)*0.5_dp*(ftza * ftzb )
                     CALL force_update(force_a,force_b,rab,pabval,ftza,ftzb,ax+1,ay,az,bx+1,by,bz ,pzyx,ipzyx)
                     pabval=pab(o1+ico,o2+jco)*0.5_dp*(ftza * ftzb )
                     CALL force_update(force_a,force_b,rab,pabval,ftza,ftzb,ax,ay+1,az,bx,by+1,bz ,pzyx,ipzyx)
                     pabval=pab(o1+ico,o2+jco)*0.5_dp*(ftza * ftzb )
                     CALL force_update(force_a,force_b,rab,pabval,ftza,ftzb,ax,ay,az+1,bx,by,bz+1 ,pzyx,ipzyx)
                  ELSE
                     pabval=pab(o1+ico,o2+jco)
                     CALL force_update(force_a,force_b,rab,pabval,ftza,ftzb,ax,ay,az,bx,by,bz,pzyx,ipzyx)
                  ENDIF
                END IF
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO
    !IF (ithread_l.eq.0) t_c_2=m_walltime()
    !IF (ithread_l.EQ.0) t_a_2=m_walltime()
    !IF (ithread_l.EQ.0) mytimings(1)=mytimings(1)+(t_a_2-t_a_1)
    !IF (ithread_l.EQ.0) mytimings(2)=mytimings(2)+(t_b_2-t_b_1)
    !IF (ithread_l.EQ.0) mytimings(3)=mytimings(3)+(t_c_2-t_c_1)

  END SUBROUTINE integrate_pgf_product_rspace

! *****************************************************************************
  SUBROUTINE force_update(force_a,force_b,rab,pab,ftza,ftzb,ax,ay,az,bx,by,bz,pzyx,ipzyx) 
    REAL(KIND=dp), DIMENSION(3), &
      INTENT(INOUT)                          :: force_a, force_b
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rab
    REAL(KIND=dp), INTENT(IN)                :: pab, ftza, ftzb
    INTEGER, INTENT(IN)                      :: ax, ay, az, bx, by, bz
    REAL(KIND=dp), POINTER                   :: pzyx(:)
    INTEGER, POINTER                         :: ipzyx(:,:,:,:,:,:)

    REAL(kind=dp)                            :: axm1, axp1, axpm0, aym1, &
                                                ayp1, azm1, azp1, bxm1, bym1, &
                                                bzm1

    axpm0 = pzyx(ipzyx(ax,ay,az,bx,by,bz))
    axp1=pzyx(ipzyx(ax+1,ay,az,bx,by,bz))
    axm1=pzyx(ipzyx(MAX(0,ax-1),ay,az,bx,by,bz))
    ayp1=pzyx(ipzyx(ax,ay+1,az,bx,by,bz))
    aym1=pzyx(ipzyx(ax,MAX(0,ay-1),az,bx,by,bz))
    azp1=pzyx(ipzyx(ax,ay,az+1,bx,by,bz))
    azm1=pzyx(ipzyx(ax,ay,MAX(0,az-1),bx,by,bz))
    bxm1=pzyx(ipzyx(ax,ay,az,MAX(0,bx-1),by,bz))
    bym1=pzyx(ipzyx(ax,ay,az,bx,MAX(0,by-1),bz))
    bzm1=pzyx(ipzyx(ax,ay,az,bx,by,MAX(0,bz-1)))
    force_a(1) = force_a(1) + pab*(ftza*axp1 - REAL(ax,dp)* axm1)
    force_a(2) = force_a(2) + pab*(ftza*ayp1 - REAL(ay,dp)* aym1)
    force_a(3) = force_a(3) + pab*(ftza*azp1 - REAL(az,dp)* azm1)
    force_b(1) = force_b(1) + pab*(ftzb*(axp1 - rab(1)*axpm0) - REAL(bx,dp)* bxm1)
    force_b(2) = force_b(2) + pab*(ftzb*(ayp1 - rab(2)*axpm0) - REAL(by,dp)* bym1)
    force_b(3) = force_b(3) + pab*(ftzb*(azp1 - rab(3)*axpm0) - REAL(bz,dp)* bzm1)

  END SUBROUTINE force_update

! *****************************************************************************

  SUBROUTINE redistribute_matrix ( dh, h, para_env)
    TYPE(real_matrix_type), POINTER          :: dh, h
    TYPE(cp_para_env_type), POINTER          :: para_env

    INTEGER                                  :: hmax, i, ii, il, isend, j, &
                                                jj, k, plength, stat
    INTEGER, DIMENSION(:, :), POINTER        :: dhlist, hlist, llist
    REAL(KIND=dp), DIMENSION(:), POINTER     :: pblock
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: dh_block, h_block

! local blocks

     CALL cp_sm_scale_and_add(matrix_a=h,matrix_b=dh)

     CALL cp_sm_get_block_list ( h, hlist )
     CALL cp_sm_get_block_list ( dh, dhlist )

     hmax = SIZE ( hlist, 2 )
     CALL mp_max ( hmax, para_env%group )
     ALLOCATE ( llist(2,hmax), STAT=stat )
     IF ( stat /= 0 ) CALL stop_memory ( "redistribute_matrix", "llist", 2*hmax )
     ALLOCATE ( pblock(10*hmax), STAT=stat )
     IF ( stat /= 0 ) CALL stop_memory ( "redistribute_matrix", "pblock", 10*hmax )
     DO isend = 1, para_env%num_pe - 1
        llist = 0
        llist(1:2,1:SIZE(hlist,2)) = hlist
        CALL mp_shift ( llist(1:2,1:hmax), para_env%group, isend )
        j = 0
        plength = 0
        DO il = 1, hmax
          IF ( llist(1,il) == 0 ) EXIT
          IF ( llist(1,il) <= llist(2,il) ) THEN
            ii = llist(1,il)
            jj = llist(2,il)
          ELSE
            ii = llist(2,il)
            jj = llist(1,il)
          END IF
          CALL get_block_node(matrix=dh,&
                              block_row=ii,&
                              block_col=jj,&
                              block=dh_block)
          IF ( ASSOCIATED ( dh_block ) ) THEN
            j = j + 1
            llist(:,j) = llist(:,il)
            DO i = 1, SIZE(dh_block,2)
              DO k = 1, SIZE(dh_block,1)
                plength = plength + 1
                IF ( plength > SIZE(pblock) ) CALL reallocate ( pblock, 1, 2*plength )
                pblock(plength) = dh_block(k,i)
              END DO
            END DO
          END IF
        END DO
        llist(:,j+1:) = 0
        CALL mp_max ( plength, para_env%group )
        IF ( plength > SIZE(pblock) ) CALL reallocate ( pblock, 1, 2*plength )
        CALL mp_shift ( llist(1:2,1:hmax), para_env%group, -isend )
        CALL mp_shift ( pblock(1:plength), para_env%group, -isend )
        plength = 0
        DO il = 1, hmax
          IF ( llist(1,il) == 0 ) EXIT
          IF ( llist(1,il) <= llist(2,il) ) THEN
            ii = llist(1,il)
            jj = llist(2,il)
          ELSE
            ii = llist(2,il)
            jj = llist(1,il)
          END IF
          CALL get_block_node(matrix=h,&
                              block_row=ii,&
                              block_col=jj,&
                              block=h_block)
          IF ( .NOT. ASSOCIATED ( h_block ) ) THEN
            CALL stop_program ( "redistribute_matrix", "h block missmatch" )
          END IF
          DO i = 1, SIZE(h_block,2)
            DO k = 1, SIZE(h_block,1)
              plength = plength + 1
              h_block(k,i) = h_block(k,i) + pblock(plength)
            END DO
          END DO
        END DO
     END DO

     DEALLOCATE ( hlist )
     DEALLOCATE ( dhlist )
     DEALLOCATE ( llist, STAT=stat )
     IF ( stat /= 0 ) CALL stop_memory ( "redistribute_matrix", "llist" )
     DEALLOCATE ( pblock, STAT=stat )
     IF ( stat /= 0 ) CALL stop_memory ( "redistribute_matrix", "pblock" )

  END SUBROUTINE redistribute_matrix

! *****************************************************************************

END MODULE qs_integrate_potential
