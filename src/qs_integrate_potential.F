!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Build up the plane wave density by collocating the primitive Gaussian
!>      functions (pgf).
!> \par History
!>      Joost VandeVondele (02.2002)
!>            1) rewrote collocate_pgf for increased accuracy and speed
!>            2) collocate_core hack for PGI compiler
!>            3) added multiple grid feature
!>            4) new way to go over the grid
!>      Joost VandeVondele (05.2002)
!>            1) prelim. introduction of the real space grid type
!>      JGH [30.08.02] multigrid arrays independent from potential
!>      JGH [17.07.03] distributed real space code
!>      JGH [23.11.03] refactoring and new loop ordering
!>      JGH [04.12.03] OpneMP parallelization of main loops
!>      Joost VandeVondele (12.2003)
!>           1) modified to compute tau
!>      Joost removed incremental build feature
!>      Joost introduced map consistent
!>      Rewrote grid integration/collocation routines, [Joost VandeVondele,03.2007]
!> \author Matthias Krack (03.04.2001)
! *****************************************************************************
MODULE qs_integrate_potential
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE cp_array_r_utils,                ONLY: cp_2d_r_p_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_rs_pool_types,                ONLY: cp_rs_pool_p_type,&
                                             cp_rs_pool_type,&
                                             rs_pool_create_rs,&
                                             rs_pool_give_back_rs,&
                                             rs_pools_create_rs_vect,&
                                             rs_pools_give_back_rs_vect
  USE cube_utils,                      ONLY: cube_info_type,&
                                             return_cube,&
                                             return_cube_nonortho
  USE d3_poly,                         ONLY: poly_d32cp2k
  USE distribution_2d_types,           ONLY: distribution_2d_release
  USE gauss_colloc,                    ONLY: integrateGaussFull
  USE gaussian_gridlevels,             ONLY: gaussian_gridlevel,&
                                             gridlevel_info_type
  USE input_constants,                 ONLY: pw_interp,&
                                             spline3_pbc_interp
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp,&
                                             int_8,&
                                             int_size
  USE mathconstants,                   ONLY: fac
  USE memory_utilities,                ONLY: reallocate
  USE message_passing,                 ONLY: mp_sync
  USE orbital_pointers,                ONLY: coset,&
                                             ncoset
  USE particle_types,                  ONLY: particle_type
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_methods,                      ONLY: pw_copy,&
                                             pw_transfer,&
                                             pw_zero
  USE pw_pool_types,                   ONLY: pw_pool_p_type,&
                                             pw_pools_create_pws,&
                                             pw_pools_give_back_pws
  USE pw_spline_utils,                 ONLY: pw_restrict_s3
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_p_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_interactions,                 ONLY: exp_radius_very_extended
  USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
  USE realspace_grid_types,            ONLY: pw2rs,&
                                             realspace_grid_p_type,&
                                             realspace_grid_type,&
                                             rs_pw_transfer
  USE realspace_task_selection,        ONLY: distribute_matrix,&
                                             int2pair
  USE sparse_matrix_types,             ONLY: &
       add_block_node, allocate_matrix, allocate_matrix_set, &
       deallocate_matrix, deallocate_matrix_set, get_block_node, &
       real_matrix_p_type, real_matrix_type, replicate_matrix
  USE task_list_types,                 ONLY: task_list_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE virial_types,                    ONLY: virial_type
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  INTEGER :: debug_count=0

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_integrate_potential'

! *** Public subroutines ***

  PUBLIC :: integrate_v_rspace,&
            integrate_v_core_rspace,&
            integrate_pgf_product_rspace,&
            potential_pw2rs

CONTAINS

! *****************************************************************************
!> \brief computes the forces/virial due to the ionic cores with a potential on
!>      grid
! *****************************************************************************
  SUBROUTINE integrate_v_core_rspace(v_rspace,qs_env,error)
    TYPE(pw_p_type), INTENT(INOUT)           :: v_rspace
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'integrate_v_core_rspace', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: atom_a, dir, handle, iatom, &
                                                ierr, ikind, j, &
                                                natom_of_kind, npme, tp(3)
    INTEGER, DIMENSION(:), POINTER           :: atom_list, cores
    LOGICAL                                  :: use_virial
    REAL(KIND=dp)                            :: alpha_core_charge, &
                                                ccore_charge, eps_rho_rspace
    REAL(KIND=dp), DIMENSION(3)              :: force_a, force_b, ra
    REAL(KIND=dp), DIMENSION(3, 3)           :: my_virial_a, my_virial_b
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
    TYPE(virial_type), POINTER               :: virial

!   --------------------------------------------------------------------

    CALL timeset(routineN,"I","",handle)

    NULLIFY(pw_env,auxbas_rs_pool,hab,pab,cores)
    ALLOCATE(hab(1,1),pab(1,1),STAT=ierr)
    IF(ierr/=0) CALL stop_memory("integrate_v_core_rspace","hab,pab",1)

    CALL get_qs_env(qs_env=qs_env,pw_env=pw_env,error=error)
    CALL pw_env_get(pw_env=pw_env,auxbas_rs_pool=auxbas_rs_pool,error=error)
    CALL rs_pool_create_rs(auxbas_rs_pool,rs_v, error=error)

    CALL rs_pw_transfer(rs_v,v_rspace%pw,pw2rs)

    CALL get_qs_env(qs_env=qs_env,&
         atomic_kind_set=atomic_kind_set,&
         cell=cell,&
         dft_control=dft_control,&
         particle_set=particle_set,&
         para_env=para_env,pw_env=pw_env,&
         force=force,virial=virial,error=error)

    use_virial = virial%pv_availability.AND.(.NOT.virial%pv_numer)

    eps_rho_rspace = dft_control%qs_control%eps_rho_rspace

    DO ikind=1,SIZE(atomic_kind_set)

       atomic_kind => atomic_kind_set(ikind)

       CALL get_atomic_kind(atomic_kind=atomic_kind,&
            natom=natom_of_kind,&
            atom_list=atom_list,&
            alpha_core_charge=alpha_core_charge,&
            ccore_charge=ccore_charge)

       pab(1,1) = -ccore_charge
      
       IF (alpha_core_charge == 0.0_dp .OR. pab(1,1)== 0.0_dp) CYCLE

       CALL reallocate ( cores, 1, natom_of_kind )
       npme = 0
       cores = 0

       DO iatom = 1, natom_of_kind
          atom_a = atom_list(iatom)
          ra(:) = pbc(particle_set(atom_a)%r,cell)
          DO dir =1,3
             tp(dir) = FLOOR(DOT_PRODUCT(cell%h_inv(dir,:),ra)*rs_v%npts(dir))
             tp(dir) = MODULO ( tp(dir), rs_v%npts(dir) )
             tp(dir) = tp(dir) + rs_v%lb(dir)
          END DO
          IF(rs_v%parallel .AND. .NOT. rs_v%distributed) THEN
              ! replicated realspace grid, split the atoms up between procs
              IF (MODULO(natom_of_kind,rs_v%group_size) == rs_v % my_pos ) THEN
                 npme = npme + 1
                 cores (npme) = iatom
              ENDIF
           ELSE
              npme = npme + 1
              cores (npme) = iatom
           ENDIF
       END DO
       
      DO j=1,npme

        iatom = cores(j)
        atom_a = atom_list(iatom)
        ra(:) = pbc(particle_set(atom_a)%r,cell)
        hab(1,1) = 0.0_dp
        force_a(:) = 0.0_dp
        force_b(:) = 0.0_dp
        IF (use_virial) THEN
          my_virial_a = 0.0_dp
          my_virial_b = 0.0_dp
        END IF

        CALL integrate_pgf_product_rspace(0,alpha_core_charge,0,&
             0,0.0_dp,0,ra,(/0.0_dp,0.0_dp,0.0_dp/),0.0_dp,&
             rs_v,cell,pw_env%cube_info(1),hab,pab=pab,o1=0,o2=0,&
             eps_gvg_rspace=eps_rho_rspace,&
             calculate_forces=.TRUE.,force_a=force_a,&
             force_b=force_b,use_virial=use_virial,my_virial_a=my_virial_a,&
             my_virial_b=my_virial_b,use_subpatch=.TRUE.,error=error)

        force(ikind)%rho_core(:,iatom) =&
          force(ikind)%rho_core(:,iatom) + force_a(:)

        IF (use_virial) THEN
          virial%pv_virial = virial%pv_virial + my_virial_a
        END IF
     END DO

    END DO

    CALL rs_pool_give_back_rs(auxbas_rs_pool,rs_v,error=error)

    DEALLOCATE(hab,pab,cores,STAT=ierr)
    IF(ierr /= 0) CALL stop_memory("integrate_v_core_rspace",&
                  "hab,pab,cores")

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE integrate_v_core_rspace

! *****************************************************************************
!> \brief computes matrix elements corresponding to a given potential
!> \note
!>     integrates a given potential (or other object on a real
!>     space grid) = v_rspace using a multi grid technique (mgrid_*)
!>     over the basis set producing a number for every element of h
!>     (should have the same sparsity structure of S)
!>     additional screening is available using the magnitude of the
!>     elements in p (? I'm not sure this is a very good idea)
!>     this argument is optional
!>     derivatives of these matrix elements with respect to the ionic
!>     coordinates can be computed as well
! *****************************************************************************
  SUBROUTINE integrate_v_rspace(v_rspace, p,h,qs_env,calculate_forces,compute_tau,gapw,&
       matrix_dv,error)

    TYPE(pw_p_type)                          :: v_rspace
    TYPE(real_matrix_p_type), INTENT(IN), &
      OPTIONAL                               :: p
    TYPE(real_matrix_p_type), INTENT(INOUT)  :: h
    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(IN)                      :: calculate_forces
    LOGICAL, INTENT(IN), OPTIONAL            :: compute_tau, gapw
    TYPE(real_matrix_p_type), DIMENSION(:), &
      OPTIONAL, TARGET                       :: matrix_dv
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'integrate_v_rspace', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: add_tasks = 1000, &
                                                max_tasks = 3000
    REAL(kind=dp), PARAMETER                 :: mult_tasks = 2.0_dp

    INTEGER :: atom_a, atom_b, bcol, brow, handle, i, iatom, iatom_new, &
      iatom_old, idir, igrid_level, ikind, ikind_old, ilevel, ipgf, ipgf_new, &
      iset, iset_new, iset_old, istat, itask, ithread, jatom, jatom_new, &
      jatom_old, jkind, jkind_old, jpgf, jpgf_new, jset, jset_new, jset_old, &
      maxco, maxpgf, maxset, maxsgf_set, n, na1, na2, natom, nb1, nb2, ncoa, &
      ncob, nkind, nseta, nsetb, ntasks, nthread, offs_dv, sgfa, sgfb, stat
    INTEGER(KIND=int_8), DIMENSION(:), &
      POINTER                                :: atom_pair_recv, atom_pair_send
    INTEGER(kind=int_8), DIMENSION(:, :), &
      POINTER                                :: tasks
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, lb_max, &
                                                lb_min, npgfa, npgfb, nsgfa, &
                                                nsgfb
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb
    LOGICAL :: atom_pair_changed, atom_pair_done, distributed_grids, failure, &
      h_duplicated, has_dv, map_consistent, my_compute_tau, my_gapw, &
      new_set_pair_coming, p_duplicated, pab_required, scatter, use_subpatch, &
      use_virial
    REAL(KIND=dp)                            :: dab, eps_gvg_rspace, &
                                                kind_radius_b, rab2, zetp
    REAL(KIND=dp), DIMENSION(3)              :: force_a, force_b, ra, rab, rb
    REAL(KIND=dp), DIMENSION(3, 3)           :: my_virial_a, my_virial_b
    REAL(KIND=dp), DIMENSION(:), POINTER     :: set_radius_a, set_radius_b
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: dist_ab, h_block, hab, &
                                                p_block, pab, rpgfa, rpgfb, &
                                                sphi_a, sphi_b, work, zeta, &
                                                zetb
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: habt, hadb, hdab, pabt, workt
    REAL(kind=dp), DIMENSION(:, :, :, :), &
      POINTER                                :: hadbt, hdabt
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_2d_r_p_type), DIMENSION(3)       :: dv_block
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_rs_pool_p_type), DIMENSION(:), &
      POINTER                                :: rs_pools
    TYPE(cube_info_type), DIMENSION(:), &
      POINTER                                :: cube_info
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(gridlevel_info_type), POINTER       :: gridlevel_info
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: ddv
    TYPE(real_matrix_type), POINTER          :: deltap, dh
    TYPE(realspace_grid_p_type), &
      DIMENSION(:), POINTER                  :: rs_v
    TYPE(section_vals_type), POINTER         :: input, interp_section
    TYPE(task_list_type), POINTER            :: task_list, task_list_soft
    TYPE(virial_type), POINTER               :: virial

!   ---------------------------------------------------------------------------

    failure=.FALSE.
    NULLIFY(pw_env, rs_pools, tasks, dist_ab)

    debug_count=debug_count+1

    offs_dv=0
    has_dv=PRESENT(matrix_dv)
    my_compute_tau = .FALSE.
    my_gapw = .FALSE.
    IF (PRESENT(compute_tau)) my_compute_tau = compute_tau
    IF (PRESENT(gapw)) my_gapw = gapw

    CALL get_qs_env(qs_env=qs_env,&
         atomic_kind_set=atomic_kind_set,&
         cell=cell,&
         dft_control=dft_control,&
         particle_set=particle_set,&
         sab_orb=sab_orb,&
         para_env=para_env,&
         input=input,&
         task_list=task_list,&
         task_list_soft=task_list_soft,&
         force=force,pw_env=pw_env,&
         virial=virial,error=error)

    CALL mp_sync(para_env%group)
    IF (my_compute_tau) THEN
       CALL timeset(routineN,"I","",handle)
    ELSE
       CALL timeset(routineN,"I","",handle)
    END IF

    ! get the task lists
    IF (my_gapw) task_list=>task_list_soft
    CPPrecondition(ASSOCIATED(task_list),cp_failure_level,routineP,error,failure)
    tasks  =>task_list%tasks
    dist_ab=>task_list%dist_ab
    atom_pair_send=>task_list%atom_pair_send
    atom_pair_recv=>task_list%atom_pair_recv
    ntasks=task_list%ntasks

    CPPrecondition(ASSOCIATED(pw_env),cp_failure_level,routineP,error,failure)
    CALL pw_env_get(pw_env, rs_pools=rs_pools,error=error)
    CALL rs_pools_create_rs_vect(rs_pools, rs_v, error=error)

    ! *** assign from pw_env
    gridlevel_info=>pw_env%gridlevel_info
    cube_info=>pw_env%cube_info

    interp_section => section_vals_get_subs_vals(input,"DFT%MGRID%INTERPOLATOR",&
         error=error)
    CALL potential_pw2rs(rs_v,v_rspace,pw_env,interp_section,error)

    !   *** having the potential on the rs_multigrids, just integrate ...
    nkind = SIZE(atomic_kind_set)
    natom = SIZE(particle_set)
    use_virial = virial%pv_availability.AND.(.NOT.virial%pv_numer)

    IF (calculate_forces) THEN
       ALLOCATE (atom_of_kind(natom),STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineP,"atom_of_kind",natom*int_size)
       CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
            atom_of_kind=atom_of_kind)
    END IF

    map_consistent=dft_control%qs_control%map_consistent
    IF (map_consistent) THEN
       eps_gvg_rspace = dft_control%qs_control%eps_rho_rspace ! needs to be consistent with rho_rspace
    ELSE
       eps_gvg_rspace = dft_control%qs_control%eps_gvg_rspace
    ENDIF

    pab_required = PRESENT(P) .AND. (calculate_forces .OR. .NOT. map_consistent)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
         maxco=maxco,&
         maxsgf_set=maxsgf_set)

    distributed_grids = .FALSE.
    DO igrid_level = 1, gridlevel_info%ngrid_levels
       IF ( rs_v(igrid_level)%rs_grid%distributed ) THEN
          distributed_grids = .TRUE.
       ENDIF
    ENDDO

    h_duplicated = .FALSE.
    dh => h%matrix
    IF (distributed_grids) THEN
       NULLIFY ( dh )
       h_duplicated = .TRUE.
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
            sparsity_id=-1, & ! unknown sparsity, as it will be modified in parallel
            matrix_symmetry=h%matrix%symmetry,error=error)
       IF (has_dv) THEN
          NULLIFY(ddv)
          CALL allocate_matrix_set(ddv,3,error=error)
          DO idir=1,3
             CALL allocate_matrix(matrix=ddv(idir)%matrix,&
                  nrow=h%matrix%nrow,&
                  ncol=h%matrix%ncol,&
                  nblock_row=h%matrix%nblock_row,&
                  nblock_col=h%matrix%nblock_col,&
                  first_row=h%matrix%first_row(:),&
                  last_row=h%matrix%last_row(:),&
                  first_col=h%matrix%first_col(:),&
                  last_col=h%matrix%last_col(:),&
                  matrix_name="LocalDV"//TRIM(ADJUSTL(cp_to_string(idir))),&
                  sparsity_id=-1, & ! unknown sparsity, as it will be modified in parallel
                  matrix_symmetry=h%matrix%symmetry,error=error)
          END DO
       END IF
    END IF

    p_duplicated = .FALSE.
    IF ( pab_required ) THEN
       deltap => p%matrix
       IF (distributed_grids) THEN
          p_duplicated = .TRUE.
          NULLIFY ( deltap )
          CALL replicate_matrix(p%matrix,deltap,target_name="LocalP",error=error)
          ! this matrix has no strict sparsity pattern in parallel
          deltap%sparsity_id=-1
          CALL distribution_2d_release(deltap%distribution_2d,error=error)
       END IF
    END IF

    !   *** Allocate work storage ***
    nthread = 1
    NULLIFY ( pabt, habt, workt )
    CALL reallocate(habt,1,maxco,1,maxco,0,nthread)
    CALL reallocate(workt,1,maxco,1,maxsgf_set,0,nthread)
    IF (pab_required) THEN
       CALL reallocate(pabt,1,maxco,1,maxco,0,nthread)
    ELSE
       IF (calculate_forces) CALL stop_program("integrate_v_rspace",&
            "need p for forces")
    ENDIF

    NULLIFY(hdabt,hadbt,hdab,hadb)

    IF (has_dv) THEN
       CALL reallocate(hdabt,1,3,1,maxco,1,1,0,nthread)
       CALL reallocate(hadbt,1,3,1,maxco,1,1,0,nthread)
    END IF

    !   get maximum numbers
    natom = SIZE( particle_set )
    maxset=0
    maxpgf=0
    DO ikind=1,nkind
       atomic_kind => atomic_kind_set(ikind)

       CALL get_atomic_kind(atomic_kind=atomic_kind,&
            softb = my_gapw, &
            orb_basis_set=orb_basis_set)

       IF (.NOT.ASSOCIATED(orb_basis_set)) CYCLE

       CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
            npgf=npgfa, nset=nseta )

       maxset=MAX(nseta,maxset)
       maxpgf=MAX(MAXVAL(npgfa),maxpgf)
    END DO

    IF (distributed_grids .AND. pab_required) THEN
        CALL distribute_matrix (rs_v(1)%rs_grid, deltap, atom_pair_send, atom_pair_recv, natom, scatter=.TRUE., error=error)
    ENDIF

    ithread = 0
    work => workt(:,:,ithread)
    hab => habt(:,:,ithread)
    IF (has_dv) THEN
       hdab => hdabt(:,:,:,ithread)
       hadb => hadbt(:,:,:,ithread)
    END IF
    IF (pab_required) THEN
       pab => pabt(:,:,ithread)
    ENDIF
    IF (has_dv) THEN
       hdab => hdabt(:,:,:,ithread)
       hadb => hadbt(:,:,:,ithread)
    END IF

    iatom_old = -1 ; jatom_old = -1 ; iset_old = -1 ; jset_old = -1 
    ikind_old = -1 ; jkind_old = -1

    loop_tasks: DO itask = 1, ntasks

       CALL int2pair(tasks(3,itask),ilevel,iatom,jatom,iset,jset,ipgf,jpgf,natom,maxset,maxpgf)

       IF (iatom .NE. iatom_old .OR. jatom .NE. jatom_old) THEN 

          ikind = particle_set(iatom)%atomic_kind%kind_number
          jkind = particle_set(jatom)%atomic_kind%kind_number

          IF (iatom .NE. iatom_old) ra(:) = pbc(particle_set(iatom)%r,cell)

          IF (iatom <= jatom) THEN
             brow = iatom
             bcol = jatom
          ELSE
             brow = jatom
             bcol = iatom
          END IF

          IF (ikind .NE. ikind_old ) THEN
             CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
                  softb = my_gapw, &
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
                  softb = my_gapw, &
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

          CALL get_block_node(matrix=dh,&
               block_row=brow,&
               block_col=bcol,&
               block=h_block)

          IF (.NOT.ASSOCIATED(h_block)) &
               CALL add_block_node ( dh, brow, bcol, h_block ,error=error)
          IF (has_dv) THEN
             DO idir=1,3
                CALL get_block_node(matrix=ddv(idir)%matrix,&
                     block_row=brow,&
                     block_col=bcol,&
                     BLOCK=dv_block(idir)%array)

                IF (.NOT.ASSOCIATED(dv_block(idir)%array)) &
                     CALL add_block_node ( ddv(idir)%matrix, brow, bcol, dv_block(idir)%array ,error=error)
             END DO
          END IF

          IF (pab_required) THEN
             CALL get_block_node(matrix=deltap,&
                  block_row=brow,&
                  block_col=bcol,&
                  BLOCK=p_block)
             IF (.NOT.ASSOCIATED(p_block)) &
                  CALL stop_program(routineP,"p_block not associated in deltap")
          END IF

          IF (calculate_forces) THEN
             atom_a = atom_of_kind(iatom)
             atom_b = atom_of_kind(jatom)
             force_a(:) = 0.0_dp
             force_b(:) = 0.0_dp
          ENDIF
          IF (use_virial) THEN
             my_virial_a = 0.0_dp
             my_virial_b = 0.0_dp
          ENDIF

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
          IF (pab_required) THEN
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

          IF (iatom<=jatom) THEN
             hab(1:ncoa,1:ncob) = 0._dp
          ELSE
             hab(1:ncob,1:ncoa) = 0._dp
          ENDIF

          ! expensive zero, presumably zeroing *far* too much in most cases.
          IF (has_dv) THEN
             hdab=0._dp
             hadb=0._dp
          END IF

          iset_old = iset
          jset_old = jset

       ENDIF

       rab(1) = dist_ab (1,itask)
       rab(2) = dist_ab (2,itask)
       rab(3) = dist_ab (3,itask)
       rab2  = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
       rb(1) = ra(1) + rab(1)
       rb(2) = ra(2) + rab(2)
       rb(3) = ra(3) + rab(3)
       zetp = zeta(ipgf,iset) + zetb(jpgf,jset)
       dab=SQRT(rab2)
       IF (dab.lt.0.1E0_dp .AND. dft_control%qs_control%map_paa) THEN
          igrid_level = 1
       ELSE
          igrid_level = gaussian_gridlevel(gridlevel_info,zetp)
       ENDIF

       na1 = (ipgf - 1)*ncoset(la_max(iset)) + 1
       na2 = ipgf*ncoset(la_max(iset))
       nb1 = (jpgf - 1)*ncoset(lb_max(jset)) + 1
       nb2 = jpgf*ncoset(lb_max(jset))

       ! check whether we need to use fawzi's generalised collocation scheme
       IF(rs_v(igrid_level)%rs_grid%distributed)THEN
          !tasks(4,:) is 0 for replicated, 1 for distributed 2 for exceptional distributed tasks
          IF (tasks(4,itask) .EQ. 2 ) THEN
             use_subpatch = .TRUE.
          ELSE
             use_subpatch = .FALSE.
          ENDIF 
       ELSE
          use_subpatch = .FALSE.
       ENDIF

       IF (pab_required) THEN
          IF (has_dv) THEN
             IF (iatom <= jatom) THEN
                CALL integrate_pgf_product_rspace(&
                     la_max(iset),zeta(ipgf,iset),la_min(iset),&
                     lb_max(jset),zetb(jpgf,jset),lb_min(jset),&
                     ra,rab,rab2,rs_v(igrid_level)%rs_grid,cell,&
                     cube_info(igrid_level),&
                     hab,pab=pab,o1=na1-1,o2=nb1-1, &
                     eps_gvg_rspace=eps_gvg_rspace,&
                     calculate_forces=calculate_forces,hdab=hdab,hadb=hadb,&
                     force_a=force_a,force_b=force_b,ithread=ithread,&
                     compute_tau=my_compute_tau,map_consistent=map_consistent,&
                     use_virial=use_virial,my_virial_a=my_virial_a,&
                     my_virial_b=my_virial_b, use_subpatch=use_subpatch,error=error)
             ELSE
                CALL integrate_pgf_product_rspace(&
                     lb_max(jset),zetb(jpgf,jset),lb_min(jset),&
                     la_max(iset),zeta(ipgf,iset),la_min(iset),&
                     rb,-rab,rab2,rs_v(igrid_level)%rs_grid,cell,&
                     cube_info(igrid_level),&
                     hab,pab=pab,o1=nb1-1,o2=na1-1, &
                     eps_gvg_rspace=eps_gvg_rspace,&
                     calculate_forces=calculate_forces,hdab=hadb,hadb=hdab,&
                     force_a=force_b,force_b=force_a,ithread=ithread,&
                     compute_tau=my_compute_tau,map_consistent=map_consistent,&
                     use_virial=use_virial,my_virial_a=my_virial_b,&
                     my_virial_b=my_virial_a,use_subpatch=use_subpatch,error=error)
             END IF
          ELSE

             IF (iatom <= jatom) THEN
                CALL integrate_pgf_product_rspace(&
                     la_max(iset),zeta(ipgf,iset),la_min(iset),&
                     lb_max(jset),zetb(jpgf,jset),lb_min(jset),&
                     ra,rab,rab2,rs_v(igrid_level)%rs_grid,cell,&
                     cube_info(igrid_level),&
                     hab,pab=pab,o1=na1-1,o2=nb1-1, &
                     eps_gvg_rspace=eps_gvg_rspace,&
                     calculate_forces=calculate_forces,&
                     force_a=force_a,force_b=force_b,ithread=ithread,&
                     compute_tau=my_compute_tau,map_consistent=map_consistent,&
                     use_virial=use_virial,my_virial_a=my_virial_a,&
                     my_virial_b=my_virial_b,use_subpatch=use_subpatch,error=error)
             ELSE
                CALL integrate_pgf_product_rspace(&
                     lb_max(jset),zetb(jpgf,jset),lb_min(jset),&
                     la_max(iset),zeta(ipgf,iset),la_min(iset),&
                     rb,-rab,rab2,rs_v(igrid_level)%rs_grid,cell,&
                     cube_info(igrid_level),&
                     hab,pab=pab,o1=nb1-1,o2=na1-1, &
                     eps_gvg_rspace=eps_gvg_rspace,&
                     calculate_forces=calculate_forces,&
                     force_a=force_b,force_b=force_a,ithread=ithread,&
                     compute_tau=my_compute_tau,map_consistent=map_consistent,&
                     use_virial=use_virial,my_virial_a=my_virial_b,&
                     my_virial_b=my_virial_a,use_subpatch=use_subpatch,error=error)
             END IF
          END IF
       ELSE
          IF (iatom <= jatom) THEN
             CALL integrate_pgf_product_rspace(&
                  la_max(iset),zeta(ipgf,iset),la_min(iset),&
                  lb_max(jset),zetb(jpgf,jset),lb_min(jset),&
                  ra,rab,rab2,rs_v(igrid_level)%rs_grid,cell,&
                  cube_info(igrid_level),&
                  hab,o1=na1-1,o2=nb1-1,&
                  eps_gvg_rspace=eps_gvg_rspace,&
                  calculate_forces=calculate_forces,&
                  force_a=force_a,force_b=force_b,ithread=ithread,&
                  compute_tau=my_compute_tau,&
                  map_consistent=map_consistent,use_subpatch=use_subpatch,error=error)
          ELSE
             CALL integrate_pgf_product_rspace(&
                  lb_max(jset),zetb(jpgf,jset),lb_min(jset),&
                  la_max(iset),zeta(ipgf,iset),la_min(iset),&
                  rb,-rab,rab2,rs_v(igrid_level)%rs_grid,cell,&
                  cube_info(igrid_level),&
                  hab,o1=nb1-1,o2=na1-1,&
                  eps_gvg_rspace=eps_gvg_rspace,&
                  calculate_forces=calculate_forces,&
                  force_a=force_b,force_b=force_a,ithread=ithread, &
                  compute_tau=my_compute_tau,&
                  map_consistent=map_consistent,use_subpatch=use_subpatch,error=error)
          END IF
       END IF

       new_set_pair_coming=.FALSE.
       atom_pair_done = .FALSE.
       IF (itask<ntasks)  THEN
          CALL int2pair(tasks(3,itask+1),ilevel,iatom_new,jatom_new,iset_new,jset_new,ipgf_new,jpgf_new,natom,maxset,maxpgf)
          IF (iatom_new .NE. iatom .OR. jatom_new .NE. jatom .OR. iset_new .NE. iset .OR. jset_new .NE. jset) THEN
             new_set_pair_coming=.TRUE.
          ENDIF
          IF(iatom_new .NE. iatom .OR. jatom_new .NE. jatom ) THEN
             atom_pair_done = .TRUE.
          ENDIF
       ELSE 
          ! do not forget the last block
          new_set_pair_coming=.TRUE.
          atom_pair_done = .TRUE.
       ENDIF

       ! contract the block into h if we're done with the current set pair
       IF (new_set_pair_coming) THEN
          IF (iatom <= jatom) THEN
             CALL dgemm("N","N",ncoa,nsgfb(jset),ncob,&
                  1.0_dp,hab(1,1),SIZE(hab,1),&
                  sphi_b(1,sgfb),SIZE(sphi_b,1),&
                  0.0_dp,work(1,1),SIZE(work,1))
             CALL dgemm("T","N",nsgfa(iset),nsgfb(jset),ncoa,&
                  1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                  work(1,1),SIZE(work,1),&
                  1.0_dp,h_block(sgfa,sgfb),SIZE(h_block,1))
             IF (has_dv) THEN
                DO idir=1,3
                   CALL dgemm("N","N",ncoa,nsgfb(jset),ncob,&
                        1.0_dp,hdab(idir,1,1),SIZE(hab,2),&
                        sphi_b(1,sgfb),SIZE(sphi_b,1),&
                        0.0_dp,work(1,1),SIZE(work,1))
                   CALL dgemm("T","N",nsgfa(iset),nsgfb(jset),ncoa,&
                        1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                        work(1,1),SIZE(work,1),&
                        1.0_dp,dv_block(idir)%array(sgfa,sgfb),SIZE(h_block,1))
                END DO
             END IF
          ELSE
             CALL dgemm("N","N",ncob,nsgfa(iset),ncoa,&
                  1.0_dp,hab(1,1),SIZE(hab,1),&
                  sphi_a(1,sgfa),SIZE(sphi_a,1),&
                  0.0_dp,work(1,1),SIZE(work,1))
             CALL dgemm("T","N",nsgfb(jset),nsgfa(iset),ncob,&
                  1.0_dp,sphi_b(1,sgfb),SIZE(sphi_b,1),&
                  work(1,1),SIZE(work,1),&
                  1.0_dp,h_block(sgfb,sgfa),SIZE(h_block,1))
             IF (has_dv) THEN
                DO idir=1,3
                   CALL dgemm("N","N",ncob,nsgfa(iset),ncoa,&
                        1.0_dp,hadb(idir,1,1),SIZE(hadb,2),&
                        sphi_a(1,sgfa),SIZE(sphi_a,1),&
                        0.0_dp,work(1,1),SIZE(work,1))
                   CALL dgemm("T","N",nsgfb(jset),nsgfa(iset),ncob,&
                        1.0_dp,sphi_b(1,sgfb),SIZE(sphi_b,1),&
                        work(1,1),SIZE(work,1),&
                        1.0_dp,dv_block(idir)%array(sgfb,sgfa),SIZE(h_block,1))
                END DO
             END IF
          END IF

       END IF

       IF (atom_pair_done) THEN
          IF (calculate_forces) THEN
             force(ikind)%rho_elec(:,atom_a) =&
                  force(ikind)%rho_elec(:,atom_a) + 2.0_dp*force_a(:)
             IF (iatom /= jatom ) THEN
                force(jkind)%rho_elec(:,atom_b) =&
                     force(jkind)%rho_elec(:,atom_b) + 2.0_dp*force_b(:)
             END IF
          ENDIF
          IF (use_virial) THEN
             IF (use_virial .AND. calculate_forces) THEN
                virial%pv_virial = virial%pv_virial + 2.0_dp*my_virial_a
                IF (iatom /= jatom) THEN
                   virial%pv_virial = virial%pv_virial + 2.0_dp*my_virial_b
                END IF
             END IF

          END IF

       ENDIF

    END DO loop_tasks

    IF ( h_duplicated ) THEN
       ! Reconstruct H matrix if using distributed RS grids
       ! note send and recv direction reversed WRT collocate
       scatter = .FALSE.
       CALL distribute_matrix (rs_v(1)%rs_grid, dh, atom_pair_recv, atom_pair_send,&
            natom, scatter, error, h%matrix)
       CALL deallocate_matrix ( dh ,error=error)

       IF (has_dv) THEN
          DO idir=1,3
             CALL distribute_matrix (rs_v(1)%rs_grid, ddv(idir)%matrix,&
                  atom_pair_recv, atom_pair_send, &
                  natom, scatter, error, matrix_dv(idir)%matrix)
          END DO
          CALL deallocate_matrix_set(ddv,error=error)
       END IF
    ELSE
       NULLIFY ( dh, ddv )
    END IF

    IF ( p_duplicated ) THEN
       CALL deallocate_matrix ( deltap ,error=error)
    ELSE
       NULLIFY ( deltap )
    END IF

    !   *** Release work storage ***

    DEALLOCATE (habt,workt,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"habt,workt")

    IF ( pab_required ) THEN
       DEALLOCATE (pabt,STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineP,"pabt")
    END IF

    CALL rs_pools_give_back_rs_vect(rs_pools, rs_v,error=error)

    IF (calculate_forces) THEN
       DEALLOCATE (atom_of_kind,STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineP,"atom_of_kind")
    END IF

    CALL mp_sync(para_env%group)

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE integrate_v_rspace

! *****************************************************************************
!> \brief transfers a potential from a pw_grid to a vector of 
!>      realspace multigrids
!> \param v_rspace INPUT : the potential on a planewave grid in Rspace
!> \param rs_v OUTPUT: the potential on the realspace multigrids
!> \note
!>      extracted from integrate_v_rspace
!>      should contain all parallel communication of integrate_v_rspace in the 
!>      case of replicated grids.
!> \par History
!>      09.2006 created [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE potential_pw2rs(rs_v,v_rspace,pw_env,interp_section,error)

    TYPE(realspace_grid_p_type), &
      DIMENSION(:), POINTER                  :: rs_v
    TYPE(pw_p_type), INTENT(IN)              :: v_rspace
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(section_vals_type), POINTER         :: interp_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'potential_pw2rs', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: auxbas_grid, handle, &
                                                igrid_level, interp_kind
    REAL(KIND=dp)                            :: scale
    TYPE(gridlevel_info_type), POINTER       :: gridlevel_info
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: mgrid_gspace, mgrid_rspace
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools

    CALL timeset(routineN,handle)

    ! *** set up of the potential on the multigrids
    CALL pw_env_get(pw_env, pw_pools=pw_pools, gridlevel_info=gridlevel_info, &
           auxbas_grid = auxbas_grid, error=error)

    CALL pw_pools_create_pws(pw_pools,mgrid_rspace,&
                use_data = REALDATA3D,&
                in_space = REALSPACE, error=error)

    ! use either realspace or fft techniques to get the potential on the rs multigrids
    CALL section_vals_val_get(interp_section,"KIND",i_val=interp_kind,error=error)
    SELECT CASE(interp_kind)
    CASE (pw_interp)
       CALL pw_pools_create_pws(pw_pools,mgrid_gspace,&
                                 use_data = COMPLEXDATA1D,&
                                 in_space = RECIPROCALSPACE, error=error)
       CALL pw_transfer(v_rspace%pw,mgrid_gspace(auxbas_grid)%pw)
       DO igrid_level=1,gridlevel_info%ngrid_levels
         IF ( igrid_level /= auxbas_grid ) THEN
              CALL pw_copy(mgrid_gspace(auxbas_grid)%pw,mgrid_gspace(igrid_level)%pw)
              CALL pw_transfer(mgrid_gspace(igrid_level)%pw,&
                                             mgrid_rspace(igrid_level)%pw)
         ELSE
              IF (mgrid_gspace(auxbas_grid)%pw%pw_grid%spherical) THEN
                  CALL pw_transfer(mgrid_gspace(auxbas_grid)%pw,&
                                             mgrid_rspace(auxbas_grid)%pw)
              ELSE ! fft forward + backward should be identical
                  CALL pw_copy(v_rspace%pw,mgrid_rspace(auxbas_grid)%pw)
              ENDIF
         ENDIF
         ! *** Multiply by the grid volume element ratio ***
         IF ( igrid_level /= auxbas_grid ) THEN
            scale = mgrid_rspace(igrid_level)%pw%pw_grid%dvol/&
                    mgrid_rspace(auxbas_grid)%pw%pw_grid%dvol
            mgrid_rspace(igrid_level)%pw%cr3d = &
                                      scale*mgrid_rspace(igrid_level)%pw%cr3d
         END IF
       END DO
       CALL pw_pools_give_back_pws(pw_pools,mgrid_gspace,error=error)
    CASE(spline3_pbc_interp)
       CALL pw_copy(v_rspace%pw,mgrid_rspace(1)%pw)
       DO igrid_level=1,gridlevel_info%ngrid_levels-1
          CALL pw_zero(mgrid_rspace(igrid_level+1)%pw)
          CALL pw_restrict_s3(mgrid_rspace(igrid_level)%pw,&
               mgrid_rspace(igrid_level+1)%pw,pw_pools(igrid_level+1)%pool,&
               interp_section,error=error)
          ! *** Multiply by the grid volume element ratio
          mgrid_rspace(igrid_level+1) % pw % cr3d = &
                 mgrid_rspace(igrid_level+1) % pw % cr3d * 8._dp
       END DO
    CASE default
       CALL cp_unimplemented_error(routineN,"interpolation not supported "//&
            cp_to_string(interp_kind),error=error)
    END SELECT

    DO igrid_level=1,gridlevel_info%ngrid_levels
       CALL rs_pw_transfer(rs_v(igrid_level)%rs_grid,&
                           mgrid_rspace(igrid_level)%pw,pw2rs)
    ENDDO
    ! *** give back the pw multi-grids
    CALL pw_pools_give_back_pws(pw_pools,mgrid_rspace,error=error)

    CALL timestop(handle)

  END SUBROUTINE potential_pw2rs

! *****************************************************************************
!> \brief low level function to compute matrix elements of primitive gaussian functions
! *****************************************************************************
    SUBROUTINE integrate_pgf_product_rspace(la_max,zeta,la_min,&
                                            lb_max,zetb,lb_min,&
                                            ra,rab,rab2,rsgrid,cell,&
                                            cube_info,hab,pab,o1,o2,&
                                            eps_gvg_rspace,&
                                            calculate_forces,hdab,hadb,force_a,force_b,&
                                            ithread,mytimings,compute_tau,map_consistent,&
                                            collocate_rho0,rpgf0_s,use_virial,my_virial_a,&
                                            my_virial_b,use_subpatch,error)

    INTEGER, INTENT(IN)                      :: la_max
    REAL(KIND=dp), INTENT(IN)                :: zeta
    INTEGER, INTENT(IN)                      :: la_min, lb_max
    REAL(KIND=dp), INTENT(IN)                :: zetb
    INTEGER, INTENT(IN)                      :: lb_min
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: ra, rab
    REAL(KIND=dp), INTENT(IN)                :: rab2
    TYPE(realspace_grid_type), POINTER       :: rsgrid
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cube_info_type), INTENT(IN)         :: cube_info
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: hab
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: pab
    INTEGER, INTENT(IN)                      :: o1, o2
    REAL(KIND=dp), INTENT(IN)                :: eps_gvg_rspace
    LOGICAL, INTENT(IN)                      :: calculate_forces
    REAL(KIND=dp), DIMENSION(:, :, :), &
      OPTIONAL, POINTER                      :: hdab, hadb
    REAL(KIND=dp), DIMENSION(3), &
      INTENT(INOUT), OPTIONAL                :: force_a, force_b
    INTEGER, OPTIONAL                        :: ithread
    REAL(KIND=dp), DIMENSION(10), OPTIONAL   :: mytimings
    LOGICAL, INTENT(IN), OPTIONAL            :: compute_tau, map_consistent, &
                                                collocate_rho0
    REAL(dp), INTENT(IN), OPTIONAL           :: rpgf0_s
    LOGICAL, INTENT(IN), OPTIONAL            :: use_virial
    REAL(KIND=dp), DIMENSION(3,3), OPTIONAL  :: my_virial_a, my_virial_b
    TYPE(cp_error_type), INTENT(inout)       :: error
    LOGICAL, OPTIONAL                        :: use_subpatch

    CHARACTER(len=*), PARAMETER :: routineN = 'integrate_pgf_product_rspace', &
      routineP = moduleN//':'//routineN

    INTEGER :: ax, ay, az, bx, by, bz, cmax, coef_max, gridbounds(2,3), i, &
      ico, icoef, ig, ithread_l, jco, k, l, la, la_max_local, la_min_local, lb, &
      lb_cube_min, lb_max_local, lb_min_local, length, lx, lx_max, lxa, lxb, lxy, &
      lxy_max, lxyz, lxyz_max, lya, lyb, lza, lzb, offset, start, ub_cube_max
    INTEGER, DIMENSION(3)                    :: cubecenter, lb_cube, ng, &
                                                ub_cube 
    INTEGER, DIMENSION(:), POINTER           :: ly_max, lz_max, sphere_bounds
    LOGICAL                                  :: my_collocate_rho0, &
                                                my_compute_tau, &
                                                my_map_consistent, &
                                                my_use_virial,&
                                                subpatch_integrate
    REAL(KIND=dp) :: a, axpm0, b, binomial_k_lxa, binomial_l_lxb, cutoff, &
      der_a(3), der_b(3), exp_x0, exp_x1, exp_x2, f, ftza, ftzb, pabval, pg, &
      prefactor, radius, rpg, ya, yap, yb, ybp, za, zap, zb, zbp, zetp
    REAL(KIND=dp), DIMENSION(3)              :: dr, rap, rb, rbp, roffset, rp
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: grid

    INTEGER :: lxp,lyp,lzp,lp,iaxis
    INTEGER,       ALLOCATABLE, DIMENSION(:,:) :: map
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:,:,:,:) :: alpha
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:) :: coef_xyz
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:) :: coef_xyt
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:) :: coef_xtt
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:,:) :: coef_ttz
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:,:,:,:) :: coef_tyz

    REAL(kind=dp), ALLOCATABLE, DIMENSION(:,:,:) :: pol_z
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:,:,:) :: pol_y
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:,:) :: pol_x
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:,:) :: vab
    REAL(KIND=dp) :: t_exp_1,t_exp_2,t_exp_min_1,t_exp_min_2,t_exp_plus_1,t_exp_plus_2

!   ---------------------------------------------------------------------------

    subpatch_integrate = .FALSE.

    IF(PRESENT(use_subpatch)) THEN 
       IF(use_subpatch)THEN 
         subpatch_integrate = .TRUE.
       ENDIF
    ENDIF  

    IF (PRESENT(ithread)) THEN
       ithread_l=ithread
    ELSE
       ithread_l=0
    ENDIF

    IF (PRESENT(use_virial)) THEN
       my_use_virial=use_virial
    ELSE
       my_use_virial=.FALSE.
    ENDIF

    ! my_compute_tau defaults to .FALSE.
    ! IF (.true.) it will compute 0.5 * (nabla x_a).(v(r) nabla x_b)
    IF (PRESENT(compute_tau)) THEN
       my_compute_tau=compute_tau
    ELSE
       my_compute_tau=.FALSE.
    ENDIF
    ! use identical radii for integrate and collocate ?
    IF (PRESENT(map_consistent)) THEN
       my_map_consistent=map_consistent
    ELSE
       my_map_consistent=.FALSE.
    ENDIF

    IF (PRESENT(collocate_rho0).AND.PRESENT(rpgf0_s)) THEN
       my_collocate_rho0=collocate_rho0
    ELSE
       my_collocate_rho0=.FALSE.
    END IF

    !IF (ithread_l.eq.0) t_a_1=m_walltime()

    IF (calculate_forces) THEN
      la_max_local=la_max+1  ! needed for the derivative of the gaussian, unimportant which one
      la_min_local=MAX(la_min-1,0) ! just in case the la_min,lb_min is not zero
      lb_min_local=MAX(lb_min-1,0)
      lb_max_local=lb_max
      IF (my_use_virial) THEN
         la_max_local=la_max_local+1
         lb_max_local=lb_max_local+1
      ENDIF
    ELSE
      la_max_local=la_max
      la_min_local=la_min
      lb_min_local=lb_min
      lb_max_local=lb_max
    END IF

    IF (my_compute_tau) THEN
      la_max_local=la_max_local+1
      lb_max_local=lb_max_local+1
      la_min_local=MAX(la_min_local-1,0)
      lb_min_local=MAX(lb_min_local-1,0)
    ENDIF

    coef_max=la_max_local+lb_max_local+1
    zetp = zeta + zetb
    f = zetb/zetp
    prefactor = EXP(-zeta*f*rab2)
!   *** position of the gaussian product
    rap(:) = f*rab(:)
    rbp(:) = rap(:) - rab(:)
    rp(:) = ra(:) + rap(:)  ! this is the gaussian center in real coordinates
    rb(:) = ra(:) + rab(:)

    IF (my_map_consistent) THEN ! still assumes that eps_gvg_rspace=eps_rho_rspace
       cutoff=1.0_dp
       radius=exp_radius_very_extended(la_min,la_max,lb_min,lb_max,ra=ra,rb=rb,rp=rp,&
               zetp=zetp,eps=eps_gvg_rspace,prefactor=prefactor,cutoff=cutoff)
    ELSE IF (my_collocate_rho0) THEN
       cutoff    = 0.0_dp
       prefactor = 1.0_dp
       radius = rpgf0_s
!       radius=exp_radius_very_extended(la_min,la_max,lb_min,lb_max,ra=ra,rb=rb,rp=rp,&
!                            zetp=zetp,eps=eps_gvg_rspace,prefactor=prefactor,cutoff=cutoff)
    ELSE
       cutoff=1.0_dp
       IF (PRESENT(pab)) THEN
          radius=exp_radius_very_extended(la_min,la_max,lb_min,lb_max,pab,o1,o2,ra,rb,rp,&
                                       zetp,eps_gvg_rspace,prefactor,cutoff)
       ELSE
          radius=exp_radius_very_extended(la_min,la_max,lb_min,lb_max,ra=ra,rb=rb,rp=rp,&
               zetp=zetp,eps=eps_gvg_rspace,prefactor=prefactor,cutoff=cutoff)
       ENDIF
    ENDIF

    IF (radius == 0.0_dp) THEN
       !IF (ithread_l.EQ.0) t_a_2=m_walltime()
       !IF (ithread_l.EQ.0) mytimings(1)=mytimings(1)+(t_a_2-t_a_1)
       RETURN
    ENDIF

    ng(:) = rsgrid%npts(:)
    grid => rsgrid%r(:,:,:)
    ALLOCATE(vab(ncoset(la_max_local),ncoset(lb_max_local)))
    vab=0.0_dp

    IF (subpatch_integrate) THEN
        CALL integrate_general_subpatch()
    ELSE
        IF (rsgrid%orthorhombic ) THEN
          CALL integrate_ortho()
        ELSE
          CALL integrate_general_wings()
          !CALL integrate_general_opt()
        END IF
    END IF

!   *** vab contains all the information needed to find the elements of hab
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
                    axpm0 =  0.5_dp * ( ax * bx * vab(coset(MAX(ax-1,0),ay,az),coset(MAX(bx-1,0),by,bz)) +  &
                                        ay * by * vab(coset(ax,MAX(ay-1,0),az),coset(bx,MAX(by-1,0),bz)) +  &
                                        az * bz * vab(coset(ax,ay,MAX(az-1,0)),coset(bx,by,MAX(bz-1,0)))  &
                                        - ftza * bx * vab(coset(ax+1,ay,az),coset(MAX(bx-1,0),by,bz))  &
                                        - ftza * by * vab(coset(ax,ay+1,az),coset(bx,MAX(by-1,0),bz))  &
                                        - ftza * bz * vab(coset(ax,ay,az+1),coset(bx,by,MAX(bz-1,0)))  &
                                        - ax * ftzb * vab(coset(MAX(ax-1,0),ay,az),coset(bx+1,by,bz))  &
                                        - ay * ftzb * vab(coset(ax,MAX(ay-1,0),az),coset(bx,by+1,bz))  &
                                        - az * ftzb * vab(coset(ax,ay,MAX(az-1,0)),coset(bx,by,bz+1)) +  &
                                        ftza * ftzb * vab(coset(ax+1,ay,az),coset(bx+1,by,bz)) + &
                                        ftza * ftzb * vab(coset(ax,ay+1,az),coset(bx,by+1,bz)) + &
                                        ftza * ftzb * vab(coset(ax,ay,az+1),coset(bx,by,bz+1)) )
                ELSE
                    axpm0 = vab(coset(ax,ay,az),coset(bx,by,bz))
                ENDIF
                hab(o1+ico,o2+jco) = hab(o1+ico,o2+jco) + axpm0
                IF (calculate_forces .AND. PRESENT(force_a)) THEN
                  IF (my_compute_tau) THEN
                     pabval=pab(o1+ico,o2+jco)*0.5_dp*ax * bx
                     CALL force_update(force_a,force_b,rab,pabval,ftza,ftzb,MAX(ax-1,0),ay,az,MAX(bx-1,0),by,bz,vab)
                     pabval=pab(o1+ico,o2+jco)*0.5_dp*ay * by
                     CALL force_update(force_a,force_b,rab,pabval,ftza,ftzb,ax,MAX(ay-1,0),az,bx,MAX(by-1,0),bz,vab)
                     pabval=pab(o1+ico,o2+jco)*0.5_dp*az * bz
                     CALL force_update(force_a,force_b,rab,pabval,ftza,ftzb,ax,ay,MAX(az-1,0),bx,by,MAX(bz-1,0),vab)
                     pabval=pab(o1+ico,o2+jco)*0.5_dp*(- ftza * bx )
                     CALL force_update(force_a,force_b,rab,pabval,ftza,ftzb,ax+1,ay,az,MAX(bx-1,0),by,bz ,vab)
                     pabval=pab(o1+ico,o2+jco)*0.5_dp*(- ftza * by )
                     CALL force_update(force_a,force_b,rab,pabval,ftza,ftzb,ax,ay+1,az,bx,MAX(by-1,0),bz ,vab)
                     pabval=pab(o1+ico,o2+jco)*0.5_dp*(- ftza * bz  )
                     CALL force_update(force_a,force_b,rab,pabval,ftza,ftzb,ax,ay,az+1,bx,by,MAX(bz-1,0) ,vab)
                     pabval=pab(o1+ico,o2+jco)*0.5_dp*(- ax * ftzb )
                     CALL force_update(force_a,force_b,rab,pabval,ftza,ftzb,MAX(ax-1,0),ay,az,bx+1,by,bz ,vab)
                     pabval=pab(o1+ico,o2+jco)*0.5_dp*(- ay * ftzb )
                     CALL force_update(force_a,force_b,rab,pabval,ftza,ftzb,ax,MAX(ay-1,0),az,bx,by+1,bz ,vab)
                     pabval=pab(o1+ico,o2+jco)*0.5_dp*(- az * ftzb )
                     CALL force_update(force_a,force_b,rab,pabval,ftza,ftzb,ax,ay,MAX(az-1,0),bx,by,bz+1 ,vab)
                     pabval=pab(o1+ico,o2+jco)*0.5_dp*(ftza * ftzb )
                     CALL force_update(force_a,force_b,rab,pabval,ftza,ftzb,ax+1,ay,az,bx+1,by,bz ,vab)
                     pabval=pab(o1+ico,o2+jco)*0.5_dp*(ftza * ftzb )
                     CALL force_update(force_a,force_b,rab,pabval,ftza,ftzb,ax,ay+1,az,bx,by+1,bz ,vab)
                     pabval=pab(o1+ico,o2+jco)*0.5_dp*(ftza * ftzb )
                     CALL force_update(force_a,force_b,rab,pabval,ftza,ftzb,ax,ay,az+1,bx,by,bz+1 ,vab)
                  ELSE
                     pabval=pab(o1+ico,o2+jco)
                     CALL force_update(force_a,force_b,rab,pabval,ftza,ftzb,ax,ay,az,bx,by,bz,vab)
                     IF (my_use_virial) THEN
                       CALL virial_update(my_virial_a,my_virial_b,rab,pabval,ftza,ftzb,ax,ay,az,bx,by,bz,vab)
                     ENDIF
                  ENDIF
                END IF
                IF (calculate_forces .AND. PRESENT(hdab)) THEN
                  der_a(1:3) = 0.0_dp
                  der_b(1:3) = 0.0_dp
                  CALL hab_derivatives(der_a,der_b,rab,ftza,ftzb,ax,ay,az,bx,by,bz,vab)
                  hdab(1:3,o1+ico,o2+jco) = der_a(1:3)
                  hadb(1:3,o1+ico,o2+jco) = der_b(1:3)
                END IF
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO
    DEALLOCATE(vab)

  CONTAINS

! *****************************************************************************
   SUBROUTINE integrate_ortho()

    CALL return_cube(cube_info,radius,lb_cube,ub_cube,sphere_bounds)    
    cmax=MAXVAL(ub_cube) 

    dr(1) = rsgrid%dh(1,1)
    dr(2) = rsgrid%dh(2,2)
    dr(3) = rsgrid%dh(3,3)

    gridbounds(1,1)=LBOUND(GRID,1)
    gridbounds(2,1)=UBOUND(GRID,1)
    gridbounds(1,2)=LBOUND(GRID,2)
    gridbounds(2,2)=UBOUND(GRID,2)
    gridbounds(1,3)=LBOUND(GRID,3)
    gridbounds(2,3)=UBOUND(GRID,3)

    cubecenter(:) = FLOOR(rp(:)/dr(:))
    roffset(:) = rp(:) - REAL(cubecenter(:),dp)*dr(:)
    lb_cube_min = MINVAL(lb_cube(:))
    ub_cube_max = MAXVAL(ub_cube(:))

!   *** a mapping so that the ig corresponds to the right grid point, also with pbc
    ALLOCATE(map(-cmax:cmax,3))
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

    lp=la_max_local+lb_max_local
    ALLOCATE(coef_xyz(((lp+1)*(lp+2)*(lp+3))/6))
    ALLOCATE(pol_z(1:2,0:lp,-cmax:0))
    ALLOCATE(pol_y(1:2,0:lp,-cmax:0))
    ALLOCATE(pol_x(0:lp,-cmax:cmax))
#include "prep.f90"

!   *** initialise the coefficient matrix, we transform the sum
!
!   sum_{lxa,lya,lza,lxb,lyb,lzb} P_{lxa,lya,lza,lxb,lyb,lzb} (x-a_x)**lxa (y-a_y)**lya (z-a_z)**lza (x-b_x)**lxb (y-a_y)**lya (z-a_z)**lza
!
!   into
!
!   sum_{lxp,lyp,lzp} P_{lxp,lyp,lzp} (x-p_x)**lxp (y-p_y)**lyp (z-p_z)**lzp
!
!   where p is center of the product gaussian, and lp = la_max + lb_max
!   (current implementation is l**7)
!

#include "call_integrate.f90"

    CALL xyz_to_vab()
    ! deallocation needed to pass around a pgi bug..
    DEALLOCATE(coef_xyz)
    DEALLOCATE(pol_z)
    DEALLOCATE(pol_y)
    DEALLOCATE(pol_x)
    DEALLOCATE(map)

    END SUBROUTINE integrate_ortho 

! *****************************************************************************
    SUBROUTINE xyz_to_vab

    coef_xyz=coef_xyz*prefactor

!   *** initialise the coefficient matrix, we transform the sum
!
!   sum_{lxa,lya,lza,lxb,lyb,lzb} P_{lxa,lya,lza,lxb,lyb,lzb} (x-a_x)**lxa (y-a_y)**lya (z-a_z)**lza (x-b_x)**lxb (y-a_y)**lya (z-a_z)**lza
!
!   into
!
!   sum_{lxp,lyp,lzp} P_{lxp,lyp,lzp} (x-p_x)**lxp (y-p_y)**lyp (z-p_z)**lzp
!
!   where p is center of the product gaussian, and lp = la_max + lb_max
!   (current implementation is l**7)
!
!
!   compute polynomial expansion coefs -> (x-a)**lxa (x-b)**lxb -> sum_{ls} alpha(ls,lxa,lxb,1)*(x-p)**ls
!
!   *** make the alpha matrix ***
    ALLOCATE(alpha(0:lp,0:la_max_local,0:lb_max_local,3))
    alpha(:,:,:,:)=0.0_dp
    DO iaxis=1,3
    DO lxa=0,la_max_local
    DO lxb=0,lb_max_local
       binomial_k_lxa=1.0_dp
       a=1.0_dp
       DO k=0,lxa
        binomial_l_lxb=1.0_dp
        b=1.0_dp
        DO l=0,lxb
           alpha(lxa-l+lxb-k,lxa,lxb,iaxis)=alpha(lxa-l+lxb-k,lxa,lxb,iaxis)+ &
                             binomial_k_lxa*binomial_l_lxb*a*b
           binomial_l_lxb=binomial_l_lxb*REAL(lxb-l,dp)/REAL(l+1,dp)
           b=b*(rp(iaxis)-(ra(iaxis)+rab(iaxis)))
        ENDDO
        binomial_k_lxa=binomial_k_lxa*REAL(lxa-k,dp)/REAL(k+1,dp)
        a=a*(-ra(iaxis)+rp(iaxis))
       ENDDO
    ENDDO
    ENDDO
    ENDDO

    !
    !   compute v_{lxa,lya,lza,lxb,lyb,lzb} given v_{lxp,lyp,lzp} and alpha(ls,lxa,lxb,1)
    !   use a three step procedure
    !

    ALLOCATE(coef_ttz(0:la_max_local,0:lb_max_local))
    ALLOCATE(coef_tyz(0:la_max_local,0:lb_max_local,0:la_max_local,0:lb_max_local))
    lxyz=0
    DO lzp=0,lp
       coef_tyz=0.0_dp
       DO lyp=0,lp-lzp
          coef_ttz=0.0_dp
          DO lxp=0,lp-lzp-lyp
             lxyz=lxyz+1
             DO lxb=0,lb_max_local
             DO lxa=0,la_max_local
                coef_ttz(lxa,lxb)=coef_ttz(lxa,lxb)+coef_xyz(lxyz)*alpha(lxp,lxa,lxb,1)
             ENDDO
             ENDDO

          ENDDO

          DO lyb=0,lb_max_local
          DO lya=0,la_max_local
             DO lxb=0,lb_max_local-lyb
             DO lxa=0,la_max_local-lya
                coef_tyz(lxa,lxb,lya,lyb)=coef_tyz(lxa,lxb,lya,lyb)+coef_ttz(lxa,lxb)*alpha(lyp,lya,lyb,2)
             ENDDO
             ENDDO
          ENDDO
          ENDDO

       ENDDO

       DO lzb=0,lb_max_local
       DO lza=0,la_max_local
          DO lyb=0,lb_max_local-lzb
          DO lya=0,la_max_local-lza
             DO lxb=MAX(lb_min_local-lzb-lyb,0),lb_max_local-lzb-lyb
             jco=coset(lxb,lyb,lzb)
             DO lxa=MAX(la_min_local-lza-lya,0),la_max_local-lza-lya
                ico=coset(lxa,lya,lza)
                vab(ico,jco)=vab(ico,jco)+coef_tyz(lxa,lxb,lya,lyb)*alpha(lzp,lza,lzb,3)
             ENDDO
             ENDDO
          ENDDO
          ENDDO
       ENDDO
       ENDDO

    ENDDO
    ! deallocation needed to pass around a pgi bug..
    DEALLOCATE(coef_tyz)
    DEALLOCATE(coef_ttz)
    DEALLOCATE(alpha)

    END SUBROUTINE xyz_to_vab

! *****************************************************************************
    SUBROUTINE integrate_general_opt()
    INTEGER :: dir, i, i_index, il, ilx, ily, ilz, index_max(3), &
      index_min(3), ismax, ismin, j, j_index, jl, jlx, jly, jlz, k, k_index, &
      kl, klx, kly, klz, lpx, lpy, lpz, lx, ly, lz, offset(3)
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: grid_map
    INTEGER, ALLOCATABLE, DIMENSION(:, :, :) :: coef_map
    REAL(KIND=dp) :: a, b, c, d, di, dip, dj, djp, dk, dkp, exp0i, exp1i, &
      exp2i, gp(3), gridval, hmatgrid(3,3), pointj(3), pointk(3), v(3)
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: coef_ijk
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: hmatgridp

! 
! transform P_{lxp,lyp,lzp} into a P_{lip,ljp,lkp} such that
! sum_{lxp,lyp,lzp} P_{lxp,lyp,lzp} (x-x_p)**lxp (y-y_p)**lyp (z-z_p)**lzp =
! sum_{lip,ljp,lkp} P_{lip,ljp,lkp} (i-i_p)**lip (j-j_p)**ljp (k-k_p)**lkp
! 

      lp=la_max_local+lb_max_local
      ALLOCATE(coef_xyz(((lp+1)*(lp+2)*(lp+3))/6))
      ALLOCATE(coef_ijk(((lp+1)*(lp+2)*(lp+3))/6))
      ALLOCATE(coef_xyt(((lp+1)*(lp+2))/2))
      ALLOCATE(coef_xtt(0:lp))

      ! aux mapping array to simplify life
      ALLOCATE(coef_map(0:lp,0:lp,0:lp))
      coef_map=HUGE(coef_map) 
      lxyz=0
      DO lzp=0,lp
      DO lyp=0,lp-lzp
      DO lxp=0,lp-lzp-lyp
          lxyz=lxyz+1
          coef_map(lxp,lyp,lzp)=lxyz
      ENDDO
      ENDDO
      ENDDO

      ! cell hmat in grid points
      hmatgrid(:,1)=cell%hmat(:,1)/ng(1)
      hmatgrid(:,2)=cell%hmat(:,2)/ng(2)
      hmatgrid(:,3)=cell%hmat(:,3)/ng(3)

      ! center in grid coords
      gp=MATMUL(cell%h_inv,rp)*ng

      ! added bt matt
      cubecenter(:) = FLOOR(gp)

      !t2=nanotime_ia32()
      !write(6,*) t2-t1
      !t1=nanotime_ia32()

      CALL return_cube_nonortho(cube_info,radius,index_min,index_max,rp)

      DO dir=1,3
         offset(dir)=MODULO(cubecenter(dir),ng(dir))+rsgrid%lb(dir)
         offset(dir)=offset(dir)-rsgrid%lb_local(dir)+1
      ENDDO

      ALLOCATE(grid_map(index_min(1):index_max(1)))
      DO i=index_min(1),index_max(1)
         grid_map(i)=MODULO(i,ng(1))+1
         IF (rsgrid % perd ( 1 )==1) THEN
            grid_map(i)=MODULO(i,ng(1))+1
         ELSE
            grid_map(i)=i-cubecenter(1)+offset(1)
         ENDIF
      ENDDO


      coef_ijk=0.0_dp

      ! go over the grid, but cycle if the point is not within the radius
      DO k=index_min(3),index_max(3)
        dk=k-gp(3)
        pointk=hmatgrid(:,3)*dk

        ! allow for generalised rs grids
        IF (rsgrid % perd ( 3 )==1) THEN
           k_index=MODULO(k,ng(3))+1
        ELSE
           k_index=k-cubecenter(3)+offset(3)
        ENDIF

        coef_xyt=0.0_dp

        DO j=index_min(2),index_max(2)
          dj=j-gp(2)
          pointj=pointk+hmatgrid(:,2)*dj
          IF (rsgrid % perd ( 2 )==1) THEN
             j_index=MODULO(j,ng(2))+1
          ELSE
             j_index=j-cubecenter(2)+offset(2)
          ENDIF

          coef_xtt=0.0_dp

          ! find bounds for the inner loop
          ! based on a quadratic equation in i
          ! a*i**2+b*i+c=radius**2
          v=pointj-gp(1)*hmatgrid(:,1)
          a=DOT_PRODUCT(hmatgrid(:,1),hmatgrid(:,1))
          b=2*DOT_PRODUCT(v,hmatgrid(:,1))
          c=DOT_PRODUCT(v,v)
          d=b*b-4*a*(c-radius**2)

          IF (d<0) THEN
              CYCLE
          ELSE
              d=SQRT(d)
              ismin=CEILING((-b-d)/(2*a))
              ismax=FLOOR((-b+d)/(2*a))
          ENDIF
          ! prepare for computing -zetp*rsq
          a=-zetp*a
          b=-zetp*b
          c=-zetp*c
          i=ismin-1
          exp2i=EXP((a*i+b)*i+c)
          exp1i=EXP(2*a*i+a+b)
          exp0i=EXP(2*a)

          coef_xtt=0.0_dp

          DO i=ismin,ismax
             di=i-gp(1)

             exp2i=exp2i*exp1i
             exp1i=exp1i*exp0i

             i_index=grid_map(i)
             gridval=grid(i_index,j_index,k_index)*exp2i

             dip=1.0_dp
             DO il=0,lp
                coef_xtt(il)=coef_xtt(il)+gridval*dip
                dip=dip*di
             ENDDO
          ENDDO

          lxy=0
          djp=1.0_dp
          DO jl=0,lp
            DO il=0,lp-jl
               lxy=lxy+1
               coef_xyt(lxy)=coef_xyt(lxy)+coef_xtt(il)*djp
            ENDDO
            djp=djp*dj
          ENDDO

        ENDDO

        lxyz = 0
        dkp=1.0_dp
        DO kl=0,lp
           lxy=0
           DO jl=0,lp-kl
              DO il=0,lp-kl-jl
                 lxyz=lxyz+1 ; lxy=lxy+1
                 coef_ijk(lxyz)=coef_ijk(lxyz)+dkp*coef_xyt(lxy)
              ENDDO
              lxy=lxy+kl
           ENDDO
           dkp=dkp*dk
        ENDDO

      ENDDO

      ! transform using multinomials
      ALLOCATE(hmatgridp(3,3,0:lp))
      hmatgridp(:,:,0)=1.0_dp
      DO k=1,lp
         hmatgridp(:,:,k)=hmatgridp(:,:,k-1)*hmatgrid(:,:)
      ENDDO

      coef_xyz=0.0_dp
      lpx=lp
      DO klx=0,lpx
      DO jlx=0,lpx-klx
      DO ilx=0,lpx-klx-jlx
         lx=ilx+jlx+klx
         lpy=lp-lx
         DO kly=0,lpy
         DO jly=0,lpy-kly
         DO ily=0,lpy-kly-jly
            ly=ily+jly+kly
            lpz=lp-lx-ly
            DO klz=0,lpz
            DO jlz=0,lpz-klz
            DO ilz=0,lpz-klz-jlz
               lz=ilz+jlz+klz

               il=ilx+ily+ilz
               jl=jlx+jly+jlz
               kl=klx+kly+klz
               coef_xyz(coef_map(lx,ly,lz))=coef_xyz(coef_map(lx,ly,lz))+ coef_ijk(coef_map(il,jl,kl))* &
                                            hmatgridp(1,1,ilx) * hmatgridp(1,2,jlx) * hmatgridp(1,3,klx) * &
                                            hmatgridp(2,1,ily) * hmatgridp(2,2,jly) * hmatgridp(2,3,kly) * &
                                            hmatgridp(3,1,ilz) * hmatgridp(3,2,jlz) * hmatgridp(3,3,klz) * &
                                            fac(lx)*fac(ly)*fac(lz)/ &
                        (fac(ilx)*fac(ily)*fac(ilz)*fac(jlx)*fac(jly)*fac(jlz)*fac(klx)*fac(kly)*fac(klz))
            ENDDO
            ENDDO
            ENDDO
         ENDDO
         ENDDO
         ENDDO
      ENDDO
      ENDDO
      ENDDO

      CALL xyz_to_vab()

      ! deallocation needed to pass around a pgi bug..
      DEALLOCATE(hmatgridp)
      DEALLOCATE(grid_map)
      DEALLOCATE(coef_map)
      DEALLOCATE(coef_xtt)
      DEALLOCATE(coef_xyt)
      DEALLOCATE(coef_ijk)
      DEALLOCATE(coef_xyz)
    END SUBROUTINE integrate_general_opt

! *****************************************************************************

    SUBROUTINE integrate_general_subpatch
    INTEGER                                  :: stat
    INTEGER, DIMENSION(2, 3)                 :: local_b
    INTEGER, DIMENSION(3)                    :: local_s, periodic
    REAL(dp), DIMENSION((&
      la_max_local+lb_max_local+1)*(&
      la_max_local+lb_max_local+2)*(&
      la_max_local+lb_max_local+3)/6)        :: poly_d3

        periodic=1 ! cell%perd         
        lp=la_max_local+lb_max_local
        local_b(1,:)=rsgrid%lb_real-rsgrid%lb
        local_b(2,:)=rsgrid%ub_real-rsgrid%lb
        local_s=rsgrid%lb_real-rsgrid%lb_local
        CALL integrateGaussFull(h=cell%hmat,h_inv=cell%h_inv,&
            grid=grid,poly=poly_d3,alphai=zetp,posi=rp,max_r2=radius*radius,&
            periodic=periodic,gdim=ng,local_bounds=local_b,local_shift=local_s,&
            error=error,scale=rsgrid%ngpts/ABS(cell%deth))
        ! defaults: local_shift=(/0,0,0/),poly_shift=(/0.0_dp,0.0_dp,0.0_dp/),scale=1.0_dp,
        ALLOCATE(coef_xyz(((lp+1)*(lp+2)*(lp+3))/6),stat=stat)
        CPPostconditionNoFail(stat==0,cp_failure_level,routineP,error)
        CALL poly_d32cp2k(coef_xyz,lp,poly_d3,error)
        CALL xyz_to_vab()
        DEALLOCATE(coef_xyz,stat=stat)
        CPPostconditionNoFail(stat==0,cp_failure_level,routineP,error)
    END SUBROUTINE

    SUBROUTINE integrate_general_wings()
    INTEGER                                  :: periodic(3), stat
    INTEGER, DIMENSION(2, 3)                 :: local_b
    REAL(dp), DIMENSION((&
      la_max_local+lb_max_local+1)*(&
      la_max_local+lb_max_local+2)*(&
      la_max_local+lb_max_local+3)/6)        :: poly_d3

        periodic=1 ! cell%perd
        lp=la_max_local+lb_max_local
        local_b(1,:)=rsgrid%lb_local-rsgrid%lb
        local_b(2,:)=MIN(local_b(1,:)+rsgrid%npts-1,rsgrid%ub_local-rsgrid%lb)
        CALL integrateGaussFull(h=cell%hmat,h_inv=cell%h_inv,&
            grid=grid,poly=poly_d3,alphai=zetp,posi=rp,max_r2=radius*radius,&
            periodic=periodic,gdim=ng,local_bounds=local_b,&
            error=error,scale=rsgrid%ngpts/ABS(cell%deth))
        ! defaults: local_shift=(/0,0,0/),poly_shift=(/0.0_dp,0.0_dp,0.0_dp/),scale=1.0_dp,
        ALLOCATE(coef_xyz(((lp+1)*(lp+2)*(lp+3))/6),stat=stat)
        CPPostconditionNoFail(stat==0,cp_failure_level,routineP,error)
        CALL poly_d32cp2k(coef_xyz,lp,poly_d3,error)
        CALL xyz_to_vab()
        DEALLOCATE(coef_xyz,stat=stat)
        CPPostconditionNoFail(stat==0,cp_failure_level,routineP,error)
    END SUBROUTINE

! *****************************************************************************
    SUBROUTINE integrate_general()
    INTEGER                                  :: i, index_max(3), &
                                                index_min(3), ipoint(3), j, k
    REAL(KIND=dp)                            :: gridval, point(3)

      CALL return_cube_nonortho(cube_info,radius,index_min,index_max,rp)

      ! go over the grid, but cycle if the point is not within the radius
      DO k=index_min(3),index_max(3)
      DO j=index_min(2),index_max(2)
      DO i=index_min(1),index_max(1)
         ! point in real space
         point=MATMUL(cell%hmat,REAL((/i,j,k/),KIND=dp)/ng)
         ! skip if outside of the sphere
         IF (SUM((point-rp)**2)>radius**2) CYCLE
         ! point on the grid (including pbc)
         ipoint=MODULO((/i,j,k/),ng)+1
         ! integrate on the grid
         gridval=grid(ipoint(1),ipoint(2),ipoint(3))
         CALL primitive_integrate(point,gridval) 
      ENDDO
      ENDDO
      ENDDO
    END SUBROUTINE integrate_general

! *****************************************************************************
    SUBROUTINE primitive_integrate(point,gridval)
    REAL(KIND=dp)                            :: point(3), gridval

    REAL(KIND=dp)                            :: dra(3), drap(3), drb(3), &
                                                drbp(3), myexp

       myexp=EXP(-zetp*SUM((point-rp)**2))*prefactor*gridval
        dra=point-ra
        drb=point-rb
        drap(1)=1.0_dp
        DO lxa=0,la_max_local
        drbp(1)=1.0_dp
        DO lxb=0,lb_max_local
           drap(2)=1.0_dp
           DO lya=0,la_max_local-lxa
           drbp(2)=1.0_dp
           DO lyb=0,lb_max_local-lxb
              drap(3)=1.0_dp
              DO lza=1,MAX(la_min_local-lxa-lya,0)
                 drap(3)=drap(3)*dra(3)
              ENDDO
              DO lza=MAX(la_min_local-lxa-lya,0),la_max_local-lxa-lya
              drbp(3)=1.0_dp
              DO lzb=1,MAX(lb_min_local-lxb-lyb,0)
                 drbp(3)=drbp(3)*drb(3)
              ENDDO
              DO lzb=MAX(lb_min_local-lxb-lyb,0),lb_max_local-lxb-lyb
                ico=coset(lxa,lya,lza)
                jco=coset(lxb,lyb,lzb)
                vab(ico,jco)=vab(ico,jco)+myexp*PRODUCT(drap)*PRODUCT(drbp)
                drbp(3)=drbp(3)*drb(3)
              ENDDO
              drap(3)=drap(3)*dra(3)
              ENDDO
           drbp(2)=drbp(2)*drb(2)
           ENDDO
           drap(2)=drap(2)*dra(2)
           ENDDO
        drbp(1)=drbp(1)*drb(1)
        ENDDO
        drap(1)=drap(1)*dra(1)
        ENDDO

    END SUBROUTINE

  END SUBROUTINE integrate_pgf_product_rspace

! *****************************************************************************
!> \brief given a set of matrix elements, perform the correct contraction to obtain the virial
! *****************************************************************************
  SUBROUTINE virial_update(my_virial_a,my_virial_b,rab,pab,&
                           ftza,ftzb,ax,ay,az,bx,by,bz,vab)
    REAL(KIND=dp), DIMENSION(3, 3), &
      INTENT(INOUT)                          :: my_virial_a, my_virial_b
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rab
    REAL(KIND=dp), INTENT(IN)                :: pab, ftza, ftzb
    INTEGER, INTENT(IN)                      :: ax, ay, az, bx, by, bz
    REAL(KIND=dp)                            :: vab(:,:)

    my_virial_a(1,1) = my_virial_a(1,1) &
        + pab*ftza*vab(coset(ax+2,ay,az),coset(bx,by,bz)) & 
        - pab*REAL(ax,dp)*vab(coset(MAX(0,ax-1)+1,ay,az),coset(bx,by,bz))
    my_virial_a(1,2) = my_virial_a(1,2) &
        + pab*ftza*vab(coset(ax+1,ay+1,az),coset(bx,by,bz)) & 
        - pab*REAL(ax,dp)*vab(coset(MAX(0,ax-1),ay+1,az),coset(bx,by,bz))
    my_virial_a(1,3) = my_virial_a(1,3) &
        + pab*ftza*vab(coset(ax+1,ay,az+1),coset(bx,by,bz)) & 
        - pab*REAL(ax,dp)*vab(coset(MAX(0,ax-1),ay,az+1),coset(bx,by,bz))
    my_virial_a(2,1) = my_virial_a(2,1) &
        + pab*ftza*vab(coset(ax+1,ay+1,az),coset(bx,by,bz)) & 
        - pab*REAL(ay,dp)*vab(coset(ax+1,MAX(0,ay-1),az),coset(bx,by,bz))
    my_virial_a(2,2) = my_virial_a(2,2) &
        + pab*ftza*vab(coset(ax,ay+2,az),coset(bx,by,bz)) & 
        - pab*REAL(ay,dp)*vab(coset(ax,MAX(0,ay-1)+1,az),coset(bx,by,bz))
    my_virial_a(2,3) = my_virial_a(2,3) &
        + pab*ftza*vab(coset(ax,ay+1,az+1),coset(bx,by,bz)) & 
        - pab*REAL(ay,dp)*vab(coset(ax,MAX(0,ay-1),az+1),coset(bx,by,bz))
    my_virial_a(3,1) = my_virial_a(3,1) &
        + pab*ftza*vab(coset(ax+1,ay,az+1),coset(bx,by,bz)) & 
        - pab*REAL(az,dp)*vab(coset(ax+1,ay,MAX(0,az-1)),coset(bx,by,bz))
    my_virial_a(3,2) = my_virial_a(3,2) &
        + pab*ftza*vab(coset(ax,ay+1,az+1),coset(bx,by,bz)) & 
        - pab*REAL(az,dp)*vab(coset(ax,ay+1,MAX(0,az-1)),coset(bx,by,bz))
    my_virial_a(3,3) = my_virial_a(3,3) &
        + pab*ftza*vab(coset(ax,ay,az+2),coset(bx,by,bz)) & 
        - pab*REAL(az,dp)*vab(coset(ax,ay,MAX(0,az-1)+1),coset(bx,by,bz))

    my_virial_b(1,1) = my_virial_b(1,1) + pab*ftzb* ( &
          vab(coset(ax+2,ay,az),coset(bx,by,bz)) &
        - vab(coset(ax+1,ay,az),coset(bx,by,bz))*rab(1) &  
        - vab(coset(ax+1,ay,az),coset(bx,by,bz))*rab(1) &  
        + vab(coset(ax,ay,az),coset(bx,by,bz))*rab(1)*rab(1) ) &
        - pab*REAL(bx,dp)*vab(coset(ax,ay,az),coset(MAX(0,bx-1)+1,by,bz))  
    my_virial_b(1,2) = my_virial_b(1,2) + pab*ftzb* ( &
          vab(coset(ax+1,ay+1,az),coset(bx,by,bz)) &
        - vab(coset(ax,ay+1,az),coset(bx,by,bz))*rab(1) &  
        - vab(coset(ax+1,ay,az),coset(bx,by,bz))*rab(2) &  
        + vab(coset(ax,ay,az),coset(bx,by,bz))*rab(1)*rab(2) ) &
        - pab*REAL(bx,dp)*vab(coset(ax,ay,az),coset(MAX(0,bx-1),by+1,bz))  
    my_virial_b(1,3) = my_virial_b(1,3) + pab*ftzb* ( &
          vab(coset(ax+1,ay,az+1),coset(bx,by,bz)) &
        - vab(coset(ax,ay,az+1),coset(bx,by,bz))*rab(1) &  
        - vab(coset(ax+1,ay,az),coset(bx,by,bz))*rab(3) &  
        + vab(coset(ax,ay,az),coset(bx,by,bz))*rab(1)*rab(3) ) &
        - pab*REAL(bx,dp)*vab(coset(ax,ay,az),coset(MAX(0,bx-1),by,bz+1))  
    my_virial_b(2,1) = my_virial_b(2,1) + pab*ftzb* ( &
          vab(coset(ax+1,ay+1,az),coset(bx,by,bz)) &
        - vab(coset(ax+1,ay,az),coset(bx,by,bz))*rab(2) &  
        - vab(coset(ax,ay+1,az),coset(bx,by,bz))*rab(1) &  
        + vab(coset(ax,ay,az),coset(bx,by,bz))*rab(2)*rab(1) ) &
        - pab*REAL(by,dp)*vab(coset(ax,ay,az),coset(bx+1,MAX(0,by-1),bz))  
    my_virial_b(2,2) = my_virial_b(2,2) + pab*ftzb* ( &
          vab(coset(ax,ay+2,az),coset(bx,by,bz)) &
        - vab(coset(ax,ay+1,az),coset(bx,by,bz))*rab(2) &  
        - vab(coset(ax,ay+1,az),coset(bx,by,bz))*rab(2) &  
        + vab(coset(ax,ay,az),coset(bx,by,bz))*rab(2)*rab(2) ) &
        - pab*REAL(by,dp)*vab(coset(ax,ay,az),coset(bx,MAX(0,by-1)+1,bz))  
    my_virial_b(2,3) = my_virial_b(2,3) + pab*ftzb* ( &
          vab(coset(ax,ay+1,az+1),coset(bx,by,bz)) &
        - vab(coset(ax,ay,az+1),coset(bx,by,bz))*rab(2) &  
        - vab(coset(ax,ay+1,az),coset(bx,by,bz))*rab(3) &  
        + vab(coset(ax,ay,az),coset(bx,by,bz))*rab(2)*rab(3) ) &
        - pab*REAL(by,dp)*vab(coset(ax,ay,az),coset(bx,MAX(0,by-1),bz+1))  
    my_virial_b(3,1) = my_virial_b(3,1) + pab*ftzb* ( &
          vab(coset(ax+1,ay,az+1),coset(bx,by,bz)) &
        - vab(coset(ax+1,ay,az),coset(bx,by,bz))*rab(3) &  
        - vab(coset(ax,ay,az+1),coset(bx,by,bz))*rab(1) &  
        + vab(coset(ax,ay,az),coset(bx,by,bz))*rab(3)*rab(1) ) &
        - pab*REAL(bz,dp)*vab(coset(ax,ay,az),coset(bx+1,by,MAX(0,bz-1)))  
    my_virial_b(3,2) = my_virial_b(3,2) + pab*ftzb* ( &
          vab(coset(ax,ay+1,az+1),coset(bx,by,bz)) &
        - vab(coset(ax,ay+1,az),coset(bx,by,bz))*rab(3) &  
        - vab(coset(ax,ay,az+1),coset(bx,by,bz))*rab(2) &  
        + vab(coset(ax,ay,az),coset(bx,by,bz))*rab(3)*rab(2) ) &
        - pab*REAL(bz,dp)*vab(coset(ax,ay,az),coset(bx,by+1,MAX(0,bz-1)))  
    my_virial_b(3,3) = my_virial_b(3,3) + pab*ftzb* ( &
          vab(coset(ax,ay,az+2),coset(bx,by,bz)) &
        - vab(coset(ax,ay,az+1),coset(bx,by,bz))*rab(3) &  
        - vab(coset(ax,ay,az+1),coset(bx,by,bz))*rab(3) &  
        + vab(coset(ax,ay,az),coset(bx,by,bz))*rab(3)*rab(3) ) &
        - pab*REAL(bz,dp)*vab(coset(ax,ay,az),coset(bx,by,MAX(0,bz-1)+1))  

  END SUBROUTINE virial_update

! *****************************************************************************
!> \brief given a bunch of matrix elements, performe the right contractions to obtain the forces
! *****************************************************************************
  SUBROUTINE force_update(force_a,force_b,rab,pab,ftza,ftzb,ax,ay,az,bx,by,bz,vab)
    REAL(KIND=dp), DIMENSION(3), &
      INTENT(INOUT)                          :: force_a, force_b
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rab
    REAL(KIND=dp), INTENT(IN)                :: pab, ftza, ftzb
    INTEGER, INTENT(IN)                      :: ax, ay, az, bx, by, bz
    REAL(KIND=dp)                            :: vab(:,:)

    REAL(KIND=dp)                            :: axm1, axp1, axpm0, aym1, &
                                                ayp1, azm1, azp1, bxm1, bym1, &
                                                bzm1

    axpm0 = vab(coset(ax,ay,az),coset(bx,by,bz))
    axp1=vab(coset(ax+1,ay,az),coset(bx,by,bz))
    axm1=vab(coset(MAX(0,ax-1),ay,az),coset(bx,by,bz))
    ayp1=vab(coset(ax,ay+1,az),coset(bx,by,bz))
    aym1=vab(coset(ax,MAX(0,ay-1),az),coset(bx,by,bz))
    azp1=vab(coset(ax,ay,az+1),coset(bx,by,bz))
    azm1=vab(coset(ax,ay,MAX(0,az-1)),coset(bx,by,bz))
    bxm1=vab(coset(ax,ay,az),coset(MAX(0,bx-1),by,bz))
    bym1=vab(coset(ax,ay,az),coset(bx,MAX(0,by-1),bz))
    bzm1=vab(coset(ax,ay,az),coset(bx,by,MAX(0,bz-1)))
    force_a(1) = force_a(1) + pab*(ftza*axp1 - REAL(ax,dp)* axm1)
    force_a(2) = force_a(2) + pab*(ftza*ayp1 - REAL(ay,dp)* aym1)
    force_a(3) = force_a(3) + pab*(ftza*azp1 - REAL(az,dp)* azm1)
    force_b(1) = force_b(1) + pab*(ftzb*(axp1 - rab(1)*axpm0) - REAL(bx,dp)* bxm1)
    force_b(2) = force_b(2) + pab*(ftzb*(ayp1 - rab(2)*axpm0) - REAL(by,dp)* bym1)
    force_b(3) = force_b(3) + pab*(ftzb*(azp1 - rab(3)*axpm0) - REAL(bz,dp)* bzm1)

  END SUBROUTINE force_update

! *****************************************************************************
!> \brief given a bunch of matrix elements perform the right contractions to obtain the
!>      derivatives of the hab matirx
! *****************************************************************************
  SUBROUTINE hab_derivatives(der_a,der_b,rab,ftza,ftzb,ax,ay,az,bx,by,bz,vab)
    REAL(KIND=dp), DIMENSION(3), &
      INTENT(INOUT)                          :: der_a, der_b
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rab
    REAL(KIND=dp), INTENT(IN)                :: ftza, ftzb
    INTEGER, INTENT(IN)                      :: ax, ay, az, bx, by, bz
    REAL(KIND=dp)                            :: vab(:,:)

    REAL(KIND=dp)                            :: axm1, axp1, axpm0, aym1, &
                                                ayp1, azm1, azp1, bxm1, bym1, &
                                                bzm1

    axpm0 = vab(coset(ax,ay,az),coset(bx,by,bz))
    axp1=vab(coset(ax+1,ay,az),coset(bx,by,bz))
    axm1=vab(coset(MAX(0,ax-1),ay,az),coset(bx,by,bz))
    ayp1=vab(coset(ax,ay+1,az),coset(bx,by,bz))
    aym1=vab(coset(ax,MAX(0,ay-1),az),coset(bx,by,bz))
    azp1=vab(coset(ax,ay,az+1),coset(bx,by,bz))
    azm1=vab(coset(ax,ay,MAX(0,az-1)),coset(bx,by,bz))
    bxm1=vab(coset(ax,ay,az),coset(MAX(0,bx-1),by,bz))
    bym1=vab(coset(ax,ay,az),coset(bx,MAX(0,by-1),bz))
    bzm1=vab(coset(ax,ay,az),coset(bx,by,MAX(0,bz-1)))
    der_a(1) =  (ftza*axp1 - REAL(ax,dp)* axm1)
    der_a(2) =  (ftza*ayp1 - REAL(ay,dp)* aym1)
    der_a(3) =  (ftza*azp1 - REAL(az,dp)* azm1)
    der_b(1) =  (ftzb*(axp1 - rab(1)*axpm0) - REAL(bx,dp)* bxm1)
    der_b(2) =  (ftzb*(ayp1 - rab(2)*axpm0) - REAL(by,dp)* bym1)
    der_b(3) =  (ftzb*(azp1 - rab(3)*axpm0) - REAL(bz,dp)* bzm1)

  END SUBROUTINE hab_derivatives

END MODULE qs_integrate_potential
