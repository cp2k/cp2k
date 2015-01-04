!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
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
MODULE qs_integrate_potential_product
  USE admm_types,                      ONLY: admm_type
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_add_block_node, cp_dbcsr_col_block_sizes, cp_dbcsr_copy, &
       cp_dbcsr_create, cp_dbcsr_deallocate_matrix, cp_dbcsr_distribution, &
       cp_dbcsr_finalize, cp_dbcsr_get_block_p, cp_dbcsr_get_data_size, &
       cp_dbcsr_get_matrix_type, cp_dbcsr_get_num_blocks, cp_dbcsr_init, &
       cp_dbcsr_p_type, cp_dbcsr_row_block_sizes, cp_dbcsr_type, &
       cp_dbcsr_work_create, dbcsr_distribution_has_threads, &
       dbcsr_distribution_obj
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cube_utils,                      ONLY: cube_info_type
  USE gaussian_gridlevels,             ONLY: gridlevel_info_type
  USE input_constants,                 ONLY: do_admm_exch_scaling_merlot,&
                                             use_aux_fit_basis_set,&
                                             use_orb_basis_set
  USE kinds,                           ONLY: dp,&
                                             int_8
  USE memory_utilities,                ONLY: reallocate
  USE orbital_pointers,                ONLY: ncoset
  USE particle_types,                  ONLY: particle_type
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_types,                        ONLY: pw_p_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_integrate_potential_low,      ONLY: integrate_pgf_product_rspace
  USE qs_kind_types,                   ONLY: get_qs_kind,&
                                             get_qs_kind_set,&
                                             qs_kind_type
  USE realspace_grid_types,            ONLY: realspace_grid_desc_p_type,&
                                             realspace_grid_p_type,&
                                             rs_grid_release,&
                                             rs_grid_retain
  USE rs_pw_interface,                 ONLY: potential_pw2rs
  USE task_list_methods,               ONLY: int2pair,&
                                             rs_distribute_matrix
  USE task_list_types,                 ONLY: task_list_type
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE virial_types,                    ONLY: virial_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  INTEGER :: debug_count=0

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.FALSE.

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_integrate_potential_product'

! *** Public subroutines ***
! *** Don't include this routines directly, use the interface to
! *** qs_integrate_potential

  PUBLIC :: integrate_v_rspace

CONTAINS

! *****************************************************************************
!> \brief computes matrix elements corresponding to a given potential
!> \param v_rspace ...
!> \param p ...
!> \param h ...
!> \param qs_env ...
!> \param calculate_forces ...
!> \param force_adm ... whether force of in aux. dens. matrix is calculated
!> \param ispin ...
!> \param compute_tau ...
!> \param gapw ...
!> \param basis_set_id ...
!> \param pw_env_external ...
!> \param task_list_external ...
!> \param error ...
!> \par History
!>      IAB (29-Apr-2010): Added OpenMP parallelisation to task loop
!>                         (c) The Numerical Algorithms Group (NAG) Ltd, 2010 on behalf of the HECToR project
!>      Some refactoring, get priorities for options correct (JGH, 04.2014)
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
  SUBROUTINE integrate_v_rspace(v_rspace, p, h, qs_env, calculate_forces, force_adm, ispin, &
       compute_tau, gapw, basis_set_id, pw_env_external, task_list_external, error)

    TYPE(pw_p_type)                          :: v_rspace
    TYPE(cp_dbcsr_p_type), INTENT(IN), &
      OPTIONAL                               :: p
    TYPE(cp_dbcsr_p_type), INTENT(INOUT)     :: h
    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(IN)                      :: calculate_forces
    LOGICAL, INTENT(IN), OPTIONAL            :: force_adm
    INTEGER, INTENT(IN), OPTIONAL            :: ispin
    LOGICAL, INTENT(IN), OPTIONAL            :: compute_tau, gapw
    INTEGER, INTENT(IN), OPTIONAL            :: basis_set_id
    TYPE(pw_env_type), OPTIONAL, POINTER     :: pw_env_external
    TYPE(task_list_type), OPTIONAL, POINTER  :: task_list_external
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'integrate_v_rspace', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: add_tasks = 1000, &
                                                max_tasks = 3000
    REAL(kind=dp), PARAMETER                 :: mult_tasks = 2.0_dp

    INTEGER :: atom_a, atom_b, bcol, brow, handle, i, iatom, igrid_level, &
      ikind, ikind_old, ilevel, ipair, ipgf, ipgf_new, iset, iset_new, &
      iset_old, itask, ithread, jatom, jkind, jkind_old, jpgf, jpgf_new, &
      jset, jset_new, jset_old, maxco, maxpgf, maxset, maxsgf_set, &
      my_basis_set_id, na1, na2, natom, nb1, nb2, ncoa, ncob, nkind, nseta, &
      nsetb, nthread, offs_dv, sgfa, sgfb, stat
    INTEGER(KIND=int_8), DIMENSION(:), &
      POINTER                                :: atom_pair_recv, atom_pair_send
    INTEGER(kind=int_8), DIMENSION(:, :), &
      POINTER                                :: tasks
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: block_touched
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, lb_max, &
                                                lb_min, npgfa, npgfb, nsgfa, &
                                                nsgfb
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb
    LOGICAL :: atom_pair_changed, atom_pair_done, distributed_grids, failure, &
      found, h_duplicated, map_consistent, my_compute_tau, my_force_adm, &
      my_gapw, new_set_pair_coming, p_duplicated, pab_required, scatter, &
      use_subpatch, use_virial
    REAL(KIND=dp)                            :: admm_scal_fac, dab, &
                                                eps_gvg_rspace, rab2, zetp
    REAL(KIND=dp), DIMENSION(3)              :: force_a, force_b, ra, rab, &
                                                rab_inv, rb
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
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: ddv
    TYPE(cp_dbcsr_type), POINTER             :: deltap, dh
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cube_info_type), DIMENSION(:), &
      POINTER                                :: cube_info
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(gridlevel_info_type), POINTER       :: gridlevel_info
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(realspace_grid_desc_p_type), &
      DIMENSION(:), POINTER                  :: rs_descs
    TYPE(realspace_grid_p_type), &
      DIMENSION(:), POINTER                  :: rs_v
    TYPE(task_list_type), POINTER            :: task_list, task_list_soft
    TYPE(virial_type), POINTER               :: virial

!$  TYPE(dbcsr_distribution_obj), POINTER    :: dist
!$  INTEGER                                  :: omp_get_max_threads, &
!$                                              omp_get_thread_num

    CALL timeset(routineN,handle)

    failure=.FALSE.
    NULLIFY(pw_env, rs_descs, tasks, dist_ab, admm_env)

    debug_count=debug_count+1

    offs_dv=0

    ! this routine works in two modes:
    ! normal mode : <a| V | b>
    ! tau mode    : < nabla a| V | nabla b>
    my_compute_tau = .FALSE.
    IF (PRESENT(compute_tau)) my_compute_tau = compute_tau

    my_force_adm = .FALSE. 
    IF (PRESENT(force_adm)) my_force_adm = force_adm


    ! this sets the basis set to be used. GAPW(==soft basis) overwrites basis_set_id
    ! default is "use_orb_basis_set"
    my_gapw = .FALSE.
    IF (PRESENT(gapw)) my_gapw = gapw
    IF (PRESENT(basis_set_id)) THEN
      my_basis_set_id = basis_set_id
    ELSE
     my_basis_set_id = use_orb_basis_set
    END IF

    ! get the task lists
    ! task lists have to be in sync with basis sets
    ! there is an option to provide the task list from outside (not through qs_env)
    ! outside option has highest priority
    SELECT CASE (my_basis_set_id)
    CASE (use_orb_basis_set)
      CALL get_qs_env(qs_env=qs_env,&
         task_list=task_list,&
         task_list_soft=task_list_soft,&
         error=error)
    CASE (use_aux_fit_basis_set)
      CALL get_qs_env(qs_env=qs_env,&
         task_list_aux_fit=task_list,&
         task_list_soft=task_list_soft,&
         error=error)
    END SELECT
    IF (my_gapw) task_list=>task_list_soft
    IF (PRESENT(task_list_external)) task_list=>task_list_external
    CPPrecondition(ASSOCIATED(task_list),cp_failure_level,routineP,error,failure)

    ! the information on the grids is provided through pw_env
    ! pw_env has to be the parent env for the potential grid (input)
    ! there is an option to provide an external grid
    CALL get_qs_env(qs_env=qs_env,pw_env=pw_env,error=error)
    IF (PRESENT(pw_env_external)) pw_env=>pw_env_external

    ! get all the general information on the system we are working on
    CALL get_qs_env(qs_env=qs_env,&
         atomic_kind_set=atomic_kind_set,&
         qs_kind_set=qs_kind_set,&
         cell=cell,&
         dft_control=dft_control,&
         particle_set=particle_set,&
         para_env=para_env,&
         force=force,&
         virial=virial,error=error)

    admm_scal_fac = 1.0_dp
    IF(my_force_adm) THEN
       CALL get_qs_env(qs_env=qs_env,admm_env=admm_env,error=error)
       ! Calculate bare scaling of force according to Merlot, 1. IF: ADMMP, 2. IF: ADMMS,
       IF( ( .NOT. admm_env%charge_constrain) .AND. &
            (admm_env%scaling_model == do_admm_exch_scaling_merlot ) ) THEN
            admm_scal_fac = admm_env%gsi(ispin)**2
       ELSE IF( admm_env%charge_constrain .AND. &
           (admm_env%scaling_model == do_admm_exch_scaling_merlot ) ) THEN
            admm_scal_fac = (admm_env%gsi(ispin))**(2.0_dp/3.0_dp)
       END IF
    END IF

    ! short cuts to task list variables
    tasks  =>task_list%tasks
    dist_ab=>task_list%dist_ab
    atom_pair_send=>task_list%atom_pair_send
    atom_pair_recv=>task_list%atom_pair_recv

    CPPrecondition(ASSOCIATED(pw_env),cp_failure_level,routineP,error,failure)
    CALL pw_env_get(pw_env, rs_descs=rs_descs, rs_grids=rs_v, error=error)
    DO i=1,SIZE(rs_v)
      CALL rs_grid_retain(rs_v(i)%rs_grid,error=error)
    END DO

    ! assign from pw_env
    gridlevel_info=>pw_env%gridlevel_info
    cube_info=>pw_env%cube_info

    ! transform the potential on the rs_multigrids
    CALL potential_pw2rs(rs_v,v_rspace,pw_env,error)

    nkind = SIZE(qs_kind_set)
    natom = SIZE(particle_set)
    use_virial = virial%pv_availability.AND.(.NOT.virial%pv_numer)

    IF (calculate_forces) THEN
       ALLOCATE (atom_of_kind(natom),STAT=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       CALL get_atomic_kind_set(atomic_kind_set, atom_of_kind=atom_of_kind)
    END IF

    map_consistent=dft_control%qs_control%map_consistent
    IF (map_consistent) THEN
       eps_gvg_rspace = dft_control%qs_control%eps_rho_rspace ! needs to be consistent with rho_rspace
    ELSE
       eps_gvg_rspace = dft_control%qs_control%eps_gvg_rspace
    ENDIF

    pab_required = PRESENT(p) .AND. (calculate_forces .OR. .NOT. map_consistent)

    CALL get_qs_kind_set(qs_kind_set=qs_kind_set,&
         maxco=maxco,&
         maxsgf_set=maxsgf_set,&
         basis_set_id=my_basis_set_id)

    distributed_grids = .FALSE.
    DO igrid_level = 1, gridlevel_info%ngrid_levels
       IF ( rs_v(igrid_level)%rs_grid%desc%distributed ) THEN
          distributed_grids = .TRUE.
       ENDIF
    ENDDO

    h_duplicated = .FALSE.
    dh => h%matrix
    IF (distributed_grids) THEN
       NULLIFY ( dh )
       h_duplicated = .TRUE.
       ALLOCATE(dh)
       CALL cp_dbcsr_init(dh, error=error)
       CALL cp_dbcsr_create(dh, 'LocalH', &
            cp_dbcsr_distribution (h%matrix),&
            cp_dbcsr_get_matrix_type (h%matrix), cp_dbcsr_row_block_sizes(h%matrix),&
            cp_dbcsr_col_block_sizes(h%matrix), cp_dbcsr_get_num_blocks(h%matrix),&
            cp_dbcsr_get_data_size(h%matrix),&
            error=error)

    END IF

    p_duplicated = .FALSE.
    IF ( pab_required ) THEN
       deltap => p%matrix
       IF (distributed_grids) THEN
          p_duplicated = .TRUE.
          NULLIFY ( deltap )
          ALLOCATE(deltap)
          CALL cp_dbcsr_init(deltap, error=error)
          CALL cp_dbcsr_copy(deltap,p%matrix,name="LocalP",error=error)
       END IF
    END IF

    nthread = 1
!$  nthread = omp_get_max_threads()

    !   *** Allocate work storage ***
    NULLIFY ( pabt, habt, workt )
    CALL reallocate(habt,1,maxco,1,maxco,0,nthread)
    CALL reallocate(workt,1,maxco,1,maxsgf_set,0,nthread)
    IF (pab_required) THEN
       CALL reallocate(pabt,1,maxco,1,maxco,0,nthread)
    ELSE
       IF (calculate_forces) THEN
          CALL stop_program(routineN,moduleN,__LINE__,&
                            "Need p for forces")
       END IF
    END IF

    NULLIFY(hdabt,hadbt,hdab,hadb)

    !   get maximum numbers
    natom = SIZE( particle_set )
    maxset=0
    maxpgf=0
    DO ikind=1,nkind
       SELECT CASE (my_basis_set_id)
       CASE (use_orb_basis_set)
         CALL get_qs_kind(qs_kind_set(ikind),&
              softb = my_gapw, &
              orb_basis_set=orb_basis_set)
       CASE (use_aux_fit_basis_set)
         CALL get_qs_kind(qs_kind_set(ikind),&
              softb = my_gapw, &
              aux_fit_basis_set=orb_basis_set,&
              basis_set_id=my_basis_set_id)
       END SELECT

       IF (.NOT.ASSOCIATED(orb_basis_set)) CYCLE

       CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
            npgf=npgfa, nset=nseta )

       maxset=MAX(nseta,maxset)
       maxpgf=MAX(MAXVAL(npgfa),maxpgf)
    END DO

    IF (distributed_grids .AND. pab_required) THEN
        CALL rs_distribute_matrix (rs_descs, deltap, atom_pair_send, atom_pair_recv, natom, scatter=.TRUE., error=error)
    ENDIF

    IF (debug_this_module) &
      ALLOCATE(block_touched(natom,natom))

!$omp parallel default(none), &
!$omp shared(workt,habt,hdabt,hadbt,pabt,tasks,particle_set,natom,maxset), &
!$omp shared(maxpgf,my_basis_set_id,my_gapw,dh,ddv,deltap,use_virial,admm_scal_fac), &
!$omp shared(pab_required,calculate_forces,ncoset,rs_v,cube_info,my_compute_tau), &
!$omp shared(map_consistent,eps_gvg_rspace,force,virial,cell,atom_of_kind,dist_ab), &
!$omp shared(gridlevel_info,task_list,failure,block_touched,nthread,qs_kind_set), &
!$omp private(ithread,work,hab,hdab,hadb,pab,iset_old,jset_old), &
!$omp private(ikind_old,jkind_old,iatom,jatom,iset,jset,ikind,jkind,ilevel,ipgf,jpgf), &
!$omp private(brow,bcol,orb_basis_set,first_sgfa,la_max,la_min,npgfa,nseta,nsgfa), &
!$omp private(rpgfa,set_radius_a,sphi_a,zeta,first_sgfb,lb_max,lb_min,npgfb), &
!$omp private(nsetb,nsgfb,rpgfb,set_radius_b,sphi_b,zetb,found,error,atom_a,atom_b), &
!$omp private(force_a,force_b,my_virial_a,my_virial_b,atom_pair_changed,h_block), &
!$omp private(p_block,ncoa,sgfa,ncob,sgfb,rab,rab2,ra,rb,zetp,dab,igrid_level), &
!$omp private(na1,na2,nb1,nb2,use_subpatch,rab_inv,new_set_pair_coming,atom_pair_done), &
!$omp private(iset_new,jset_new,ipgf_new,jpgf_new,dist), &
!$omp private(itask)

    ithread = 0
!$  ithread = omp_get_thread_num()
    work => workt(:,:,ithread)
    hab => habt(:,:,ithread)
    IF (pab_required) THEN
       pab => pabt(:,:,ithread)
    END IF

    iset_old = -1 ; jset_old = -1
    ikind_old = -1 ; jkind_old = -1


    ! Here we loop over gridlevels first, finalising the matrix after each grid level is
    ! completed.  On each grid level, we loop over atom pairs, which will only access
    ! a single block of each matrix, so with OpenMP, each matrix block is only touched
    ! by a single thread for each grid level
    loop_gridlevels: DO igrid_level = 1, gridlevel_info%ngrid_levels

    CALL cp_dbcsr_work_create(dh,work_mutable=.TRUE.,n=nthread,error=error)
!$  dist => dh%matrix%m%dist
!$  CALL cp_assert (dbcsr_distribution_has_threads(dist), cp_fatal_level,&
!$       cp_internal_error, routineN, "No thread distribution defined.",&
!$       error=error)
!$omp barrier

    IF (debug_this_module) THEN
!$omp single
      block_touched = -1
!$omp end single
!$omp flush
    END IF

!$omp do schedule (dynamic, MAX(1,task_list%npairs(igrid_level)/(nthread*50)))
    loop_pairs: DO ipair = 1, task_list%npairs(igrid_level)
    loop_tasks: DO itask = task_list%taskstart(ipair,igrid_level), task_list%taskstop(ipair,igrid_level)

       CALL int2pair(tasks(3,itask),ilevel,iatom,jatom,iset,jset,ipgf,jpgf,natom,maxset,maxpgf)

       ! At the start of a block of tasks, get atom data (and kind data, if needed)
       IF (itask .EQ. task_list%taskstart(ipair,igrid_level) ) THEN

          ikind = particle_set(iatom)%atomic_kind%kind_number
          jkind = particle_set(jatom)%atomic_kind%kind_number

          ra(:) = pbc(particle_set(iatom)%r,cell)

          IF (iatom <= jatom) THEN
             brow = iatom
             bcol = jatom
          ELSE
             brow = jatom
             bcol = iatom
          END IF

          IF (ikind .NE. ikind_old ) THEN
             SELECT CASE (my_basis_set_id)
             CASE (use_orb_basis_set)
               CALL get_qs_kind(qs_kind_set(ikind),&
                    softb = my_gapw, &
                    orb_basis_set=orb_basis_set)
             CASE (use_aux_fit_basis_set)
               CALL get_qs_kind(qs_kind_set(ikind),&
                    softb = my_gapw, &
                    aux_fit_basis_set=orb_basis_set,&
                    basis_set_id = my_basis_set_id)
             END SELECT

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
             SELECT CASE (my_basis_set_id)
             CASE (use_orb_basis_set)
               CALL get_qs_kind(qs_kind_set(jkind),&
                    softb = my_gapw, &
                    orb_basis_set=orb_basis_set)
             CASE (use_aux_fit_basis_set)
               CALL get_qs_kind(qs_kind_set(jkind),&
                    softb = my_gapw, &
                    aux_fit_basis_set=orb_basis_set,&
                    basis_set_id=my_basis_set_id)
             END SELECT
             CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                  first_sgf=first_sgfb,&
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

          IF (debug_this_module) THEN
!$omp critical (block_touched_critical)
            IF ((block_touched(brow,bcol).NE.ithread) .AND. (block_touched(brow,bcol).NE. -1) ) THEN
              CALL stop_program(routineN,moduleN,__LINE__,&
                                "Block has been modified by another thread")
            END IF
            block_touched(brow,bcol) = ithread
!$omp end critical (block_touched_critical)
          END IF

          NULLIFY(h_block)
          CALL cp_dbcsr_get_block_p(dh,brow,bcol,h_block,found)
          IF (.NOT.ASSOCIATED(h_block)) THEN
               CALL cp_dbcsr_add_block_node ( dh, brow, bcol, h_block ,error=error)
          END IF

          IF (pab_required) THEN
             CALL cp_dbcsr_get_block_p(matrix=deltap,&
                  row=brow,col=bcol,BLOCK=p_block,found=found)

             IF (.NOT.ASSOCIATED(p_block)) THEN
                CALL stop_program(routineN,moduleN,__LINE__,&
                                  "p_block not associated in deltap")
             END IF
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

       na1 = (ipgf - 1)*ncoset(la_max(iset)) + 1
       na2 = ipgf*ncoset(la_max(iset))
       nb1 = (jpgf - 1)*ncoset(lb_max(jset)) + 1
       nb2 = jpgf*ncoset(lb_max(jset))

       ! check whether we need to use fawzi's generalised collocation scheme
       IF(rs_v(igrid_level)%rs_grid%desc%distributed)THEN
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
          IF (iatom <= jatom) THEN
             CALL integrate_pgf_product_rspace(&
                  la_max(iset),zeta(ipgf,iset),la_min(iset),&
                  lb_max(jset),zetb(jpgf,jset),lb_min(jset),&
                  ra,rab,rab2,rs_v(igrid_level)%rs_grid,cell,&
                  cube_info(igrid_level),&
                  hab,pab=pab,o1=na1-1,o2=nb1-1, &
                  eps_gvg_rspace=eps_gvg_rspace,&
                  calculate_forces=calculate_forces,&
                  force_a=force_a,force_b=force_b,&
                  compute_tau=my_compute_tau,map_consistent=map_consistent,&
                  use_virial=use_virial,my_virial_a=my_virial_a,&
                  my_virial_b=my_virial_b,use_subpatch=use_subpatch,subpatch_pattern=tasks(6,itask),error=error)
          ELSE
             rab_inv=-rab
             CALL integrate_pgf_product_rspace(&
                  lb_max(jset),zetb(jpgf,jset),lb_min(jset),&
                  la_max(iset),zeta(ipgf,iset),la_min(iset),&
                  rb,rab_inv,rab2,rs_v(igrid_level)%rs_grid,cell,&
                  cube_info(igrid_level),&
                  hab,pab=pab,o1=nb1-1,o2=na1-1, &
                  eps_gvg_rspace=eps_gvg_rspace,&
                  calculate_forces=calculate_forces,&
                  force_a=force_b,force_b=force_a,&
                  compute_tau=my_compute_tau,map_consistent=map_consistent,&
                  use_virial=use_virial,my_virial_a=my_virial_b,&
                  my_virial_b=my_virial_a,use_subpatch=use_subpatch,subpatch_pattern=tasks(6,itask),error=error)
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
                  force_a=force_a,force_b=force_b,&
                  compute_tau=my_compute_tau,&
                  map_consistent=map_consistent,use_subpatch=use_subpatch,subpatch_pattern=tasks(6,itask),error=error)
          ELSE
             rab_inv=-rab
             CALL integrate_pgf_product_rspace(&
                  lb_max(jset),zetb(jpgf,jset),lb_min(jset),&
                  la_max(iset),zeta(ipgf,iset),la_min(iset),&
                  rb,rab_inv,rab2,rs_v(igrid_level)%rs_grid,cell,&
                  cube_info(igrid_level),&
                  hab,o1=nb1-1,o2=na1-1,&
                  eps_gvg_rspace=eps_gvg_rspace,&
                  calculate_forces=calculate_forces,&
                  force_a=force_b,force_b=force_a, &
                  compute_tau=my_compute_tau,&
                  map_consistent=map_consistent,use_subpatch=use_subpatch,subpatch_pattern=tasks(6,itask),error=error)
          END IF
       END IF

       new_set_pair_coming=.FALSE.
       atom_pair_done = .FALSE.
       IF (itask < task_list%taskstop(ipair,igrid_level))  THEN
          CALL int2pair(tasks(3,itask+1),ilevel,iatom,jatom,iset_new,jset_new,ipgf_new,jpgf_new,natom,maxset,maxpgf)
          IF (iset_new .NE. iset .OR. jset_new .NE. jset) THEN
             new_set_pair_coming=.TRUE.
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
       END IF

       IF (atom_pair_done) THEN
!$omp critical(force_critical)
          IF (calculate_forces) THEN
             force(ikind)%rho_elec(:,atom_a) =&
                  force(ikind)%rho_elec(:,atom_a) + 2.0_dp*admm_scal_fac*force_a(:)
             IF (iatom /= jatom ) THEN
                force(jkind)%rho_elec(:,atom_b) =&
                     force(jkind)%rho_elec(:,atom_b) + 2.0_dp*admm_scal_fac*force_b(:)
             END IF
          ENDIF
          IF (use_virial) THEN
             IF (use_virial .AND. calculate_forces) THEN
                virial%pv_virial = virial%pv_virial + 2.0_dp*admm_scal_fac*my_virial_a
                IF (iatom /= jatom) THEN
                   virial%pv_virial = virial%pv_virial + 2.0_dp*admm_scal_fac*my_virial_b
                END IF
             END IF
          END IF
!$omp end critical (force_critical)
       ENDIF
    END DO loop_tasks
    END DO loop_pairs
!$omp end do

    CALL cp_dbcsr_finalize(dh, error=error)

    END DO loop_gridlevels

!$omp end parallel

    IF (debug_this_module) &
      DEALLOCATE(block_touched)

    IF ( h_duplicated ) THEN
       ! Reconstruct H matrix if using distributed RS grids
       ! note send and recv direction reversed WRT collocate
       scatter = .FALSE.
       CALL rs_distribute_matrix (rs_descs, dh, atom_pair_recv, atom_pair_send,&
            natom, scatter, error, h%matrix)
       CALL cp_dbcsr_deallocate_matrix ( dh ,error=error)

    ELSE
       NULLIFY ( dh, ddv )
    END IF

    IF ( p_duplicated ) THEN
       CALL cp_dbcsr_deallocate_matrix ( deltap ,error=error)
    ELSE
       NULLIFY ( deltap )
    END IF

    !   *** Release work storage ***

    DEALLOCATE (habt,workt,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    IF ( pab_required ) THEN
       DEALLOCATE (pabt,STAT=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

    IF (ASSOCIATED(rs_v)) THEN
      DO i=1,SIZE(rs_v)
        CALL rs_grid_release(rs_v(i)%rs_grid, error=error)
      END DO
    END IF

    IF (calculate_forces) THEN
       DEALLOCATE (atom_of_kind,STAT=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

    CALL timestop(handle)

  END SUBROUTINE integrate_v_rspace

END MODULE qs_integrate_potential_product
