!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief generate the tasks lists used by collocate and integrate routines
!> \par History
!>      01.2008 [Joost VandeVondele] refactered out of qs_collocate / qs_integrate
!> \author Joost VandeVondele
! *****************************************************************************
MODULE task_list_methods
  USE ao_util,                         ONLY: exp_radius_very_extended
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_p_type,&
                                             gto_basis_set_type
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_interface,              ONLY: &
       convert_sizes_to_offsets, cp_dbcsr_col_block_sizes, cp_dbcsr_finalize, &
       cp_dbcsr_get_block_p, cp_dbcsr_get_info, cp_dbcsr_p_type, &
       cp_dbcsr_put_block, cp_dbcsr_row_block_sizes, cp_dbcsr_type, &
       cp_dbcsr_work_create
  USE cube_utils,                      ONLY: compute_cube_center,&
                                             cube_info_type,&
                                             return_cube,&
                                             return_cube_nonortho
  USE gaussian_gridlevels,             ONLY: gaussian_gridlevel,&
                                             gridlevel_info_type
  USE input_constants,                 ONLY: use_aux_fit_basis_set,&
                                             use_orb_basis_set
  USE kinds,                           ONLY: dp,&
                                             int_8
  USE kpoint_types,                    ONLY: get_kpoint_info,&
                                             kpoint_type
  USE memory_utilities,                ONLY: reallocate
  USE message_passing,                 ONLY: mp_allgather,&
                                             mp_alltoall,&
                                             mp_gather,&
                                             mp_scatter
  USE particle_types,                  ONLY: particle_type
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE qs_kind_types,                   ONLY: get_qs_kind,&
                                             qs_kind_type
  USE qs_ks_types,                     ONLY: get_ks_env,&
                                             qs_ks_env_type
  USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                             neighbor_list_iterate,&
                                             neighbor_list_iterator_create,&
                                             neighbor_list_iterator_p_type,&
                                             neighbor_list_iterator_release,&
                                             neighbor_list_set_p_type
  USE realspace_grid_types,            ONLY: realspace_grid_desc_p_type,&
                                             realspace_grid_desc_type,&
                                             realspace_grid_p_type,&
                                             rs_grid_create,&
                                             rs_grid_locate_rank,&
                                             rs_grid_release,&
                                             rs_grid_reorder_ranks
  USE task_list_types,                 ONLY: task_list_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: sort

!$ USE omp_lib
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'task_list_methods'

  PUBLIC :: generate_qs_task_list,&
            task_list_inner_loop
  PUBLIC :: distribute_tasks,&
            int2pair,&
            rs_distribute_matrix

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param ks_env ...
!> \param task_list ...
!> \param reorder_rs_grid_ranks Flag that indicates if this routine should
!>        or should not overwrite the rs descriptor (see comment below)
!> \param skip_load_balance_distributed ...
!> \param soft_valid ...
!> \param basis_set_id ...
!> \param pw_env_external ...
!> \param sab_orb_external ...
!> \param error ...
!> \par History
!>      01.2008 factored out of calculate_rho_elec [Joost VandeVondele]
!>      04.2010 divides tasks into grid levels and atom pairs for integrate/collocate [Iain Bethune]
!>              (c) The Numerical Algorithms Group (NAG) Ltd, 2010 on behalf of the HECToR project
!>      06.2015 adjusted to be used with multiple images (k-points) [JGH]
!> \note  If this routine is called several times with different task lists,
!>        the default behaviour is to re-optimize the grid ranks and overwrite
!>        the rs descriptor and grids. reorder_rs_grid_ranks = .FALSE. prevents the code
!>        of performing a new optimization by leaving the rank order in
!>        its current state.
! *****************************************************************************
  SUBROUTINE generate_qs_task_list(ks_env, task_list,&
             reorder_rs_grid_ranks, skip_load_balance_distributed, &
             soft_valid, basis_set_id, pw_env_external, sab_orb_external, error)

    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(task_list_type), POINTER            :: task_list
    LOGICAL, INTENT(IN)                      :: reorder_rs_grid_ranks, &
                                                skip_load_balance_distributed
    LOGICAL, INTENT(IN), OPTIONAL            :: soft_valid
    INTEGER, INTENT(IN), OPTIONAL            :: basis_set_id
    TYPE(pw_env_type), OPTIONAL, POINTER     :: pw_env_external
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), OPTIONAL, POINTER        :: sab_orb_external
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'generate_qs_task_list', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: max_tasks = 2000

    INTEGER :: cindex, curr_tasks, handle, i, iatom, iatom_old, igrid_level, &
      igrid_level_old, ikind, img, img_old, inode, ipair, ipgf, iset, jatom, &
      jatom_old, jkind, jpgf, jset, maxpgf, maxset, my_basis_set_id, natoms, &
      nimages, nkind, nseta, nsetb, stat
    INTEGER, DIMENSION(3)                    :: cellind
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, lb_max, &
                                                lb_min, npgfa, npgfb
    INTEGER, DIMENSION(:, :, :), POINTER     :: cell_to_index
    LOGICAL                                  :: dokp, failure, my_soft
    REAL(KIND=dp)                            :: kind_radius_a, kind_radius_b
    REAL(KIND=dp), DIMENSION(3)              :: ra, rab
    REAL(KIND=dp), DIMENSION(:), POINTER     :: set_radius_a, set_radius_b
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: rpgfa, rpgfb, zeta, zetb
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cube_info_type), DIMENSION(:), &
      POINTER                                :: cube_info
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(gridlevel_info_type), POINTER       :: gridlevel_info
    TYPE(gto_basis_set_p_type), &
      DIMENSION(:), POINTER                  :: basis_set_list
    TYPE(gto_basis_set_type), POINTER        :: basis_set_a, basis_set_b, &
                                                orb_basis_set
    TYPE(kpoint_type), POINTER               :: kpoints
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
      DIMENSION(:), POINTER                  :: rs_grids

    CALL timeset(routineN,handle)

    failure=.FALSE.

    ! by default, the full density is calculated
    my_soft=.FALSE.
    IF (PRESENT(soft_valid)) my_soft = soft_valid

    CALL get_ks_env(ks_env,&
         qs_kind_set=qs_kind_set,&
         cell=cell,&
         particle_set=particle_set,&
         dft_control=dft_control,&
         pw_env=pw_env,error=error)

    IF( PRESENT(basis_set_id)) THEN
      my_basis_set_id = basis_set_id
    ELSE
      my_basis_set_id = use_orb_basis_set
    END IF

    SELECT CASE (my_basis_set_id)
    CASE (use_orb_basis_set)
      CALL get_ks_env(ks_env,&
           sab_orb=sab_orb,error=error)
    CASE (use_aux_fit_basis_set)
      CALL get_ks_env(ks_env,&
           sab_aux_fit=sab_orb,error=error)
    END SELECT

    ! maybe just don't pas ks_env instead ??
    IF (PRESENT(pw_env_external)) pw_env=>pw_env_external
    IF (PRESENT(sab_orb_external)) sab_orb=>sab_orb_external

    CALL pw_env_get(pw_env, rs_descs=rs_descs, rs_grids=rs_grids, error=error)

    ! *** assign from pw_env
    gridlevel_info=>pw_env%gridlevel_info
    cube_info=>pw_env%cube_info

    ! find maximum numbers
    nkind = SIZE(qs_kind_set)
    natoms = SIZE( particle_set )
    maxset=0
    maxpgf=0
    DO ikind=1,nkind
       qs_kind => qs_kind_set(ikind)
       SELECT CASE (my_basis_set_id)
       CASE (use_orb_basis_set)
          CALL get_qs_kind(qs_kind=qs_kind,&
               softb = my_soft, &
               orb_basis_set=orb_basis_set)
       CASE (use_aux_fit_basis_set)
          CALL get_qs_kind(qs_kind=qs_kind,&
               softb = my_soft, &
               aux_fit_basis_set=orb_basis_set,&
               basis_set_id=my_basis_set_id)
       END SELECT

       IF (.NOT.ASSOCIATED(orb_basis_set)) CYCLE
       CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
            npgf=npgfa, nset=nseta )

       maxset=MAX(nseta,maxset)
       maxpgf=MAX(MAXVAL(npgfa),maxpgf)
    END DO

    ! kpoint related
    nimages = dft_control%nimages
    IF(nimages > 1) THEN
       dokp = .TRUE.
       NULLIFY(kpoints)
       CALL get_ks_env(ks_env=ks_env,kpoints=kpoints,error=error)
       CALL get_kpoint_info(kpoint=kpoints,cell_to_index=cell_to_index,error=error)
    ELSE
       dokp = .FALSE.
       NULLIFY(cell_to_index)
    END IF

    ! free the atom_pair lists if allocated
    IF (ASSOCIATED(task_list%atom_pair_send)) DEALLOCATE (task_list%atom_pair_send)
    IF (ASSOCIATED(task_list%atom_pair_recv)) DEALLOCATE (task_list%atom_pair_recv)

    ! construct a list of all tasks
    IF (.NOT.ASSOCIATED(task_list%tasks)) THEN
       CALL reallocate(task_list%tasks,1,7,1,max_tasks)
       CALL reallocate(task_list%dist_ab,1,3,1,max_tasks)
    ENDIF
    task_list%ntasks = 0
    curr_tasks = SIZE(task_list%tasks,2)

    ALLOCATE (basis_set_list(nkind))
    DO ikind=1,nkind
       qs_kind => qs_kind_set(ikind)
       SELECT CASE (my_basis_set_id)
       CASE (use_orb_basis_set)
          CALL get_qs_kind(qs_kind=qs_kind,softb=my_soft,orb_basis_set=basis_set_a)
       CASE (use_aux_fit_basis_set)
          CALL get_qs_kind(qs_kind=qs_kind,softb=my_soft,&
               aux_fit_basis_set=basis_set_a,basis_set_id=my_basis_set_id)
       END SELECT
       IF (ASSOCIATED(basis_set_a)) THEN
          basis_set_list(ikind)%gto_basis_set => basis_set_a
       ELSE
          NULLIFY(basis_set_list(ikind)%gto_basis_set)
       END IF
    END DO
    CALL neighbor_list_iterator_create(nl_iterator,sab_orb)
    DO WHILE (neighbor_list_iterate(nl_iterator)==0)
       CALL get_iterator_info(nl_iterator,ikind=ikind,jkind=jkind,inode=inode,&
                              iatom=iatom,jatom=jatom,r=rab,cell=cellind)
       basis_set_a => basis_set_list(ikind)%gto_basis_set
       IF (.NOT.ASSOCIATED(basis_set_a)) CYCLE
       basis_set_b => basis_set_list(jkind)%gto_basis_set
       IF (.NOT.ASSOCIATED(basis_set_b)) CYCLE
       ra(:) = pbc(particle_set(iatom)%r,cell)
       ! basis ikind
       la_max       =>  basis_set_a%lmax
       la_min       =>  basis_set_a%lmin
       npgfa        =>  basis_set_a%npgf
       nseta        =   basis_set_a%nset
       rpgfa        =>  basis_set_a%pgf_radius
       set_radius_a =>  basis_set_a%set_radius
       kind_radius_a=   basis_set_a%kind_radius
       zeta         =>  basis_set_a%zet
       ! basis jkind
       lb_max       =>  basis_set_b%lmax
       lb_min       =>  basis_set_b%lmin
       npgfb        =>  basis_set_b%npgf
       nsetb        =   basis_set_b%nset
       rpgfb        =>  basis_set_b%pgf_radius
       set_radius_b =>  basis_set_b%set_radius
       kind_radius_b=   basis_set_b%kind_radius
       zetb         =>  basis_set_b%zet

       IF(dokp) THEN
          cindex = cell_to_index(cellind(1),cellind(2),cellind(3))
       ELSE
          cindex = 1
       END IF

       CALL task_list_inner_loop(task_list%tasks,task_list%dist_ab,task_list%ntasks,curr_tasks,&
                rs_descs,dft_control,cube_info,gridlevel_info,cindex,&
                iatom,jatom,rpgfa,rpgfb,zeta,zetb,kind_radius_b,set_radius_a,set_radius_b,ra,rab,&
                la_max,la_min,lb_max,lb_min,npgfa,npgfb,maxpgf,maxset,natoms,nseta,nsetb,error)

    END DO
    CALL neighbor_list_iterator_release(nl_iterator)

    DEALLOCATE (basis_set_list)

    ! redistribute the task list so that all tasks map on the local rs grids
    CALL distribute_tasks (rs_descs, task_list%ntasks, natoms, maxset, maxpgf, nimages,&
         task_list%tasks, rval=task_list%dist_ab, atom_pair_send=task_list%atom_pair_send,&
         atom_pair_recv=task_list%atom_pair_recv, symmetric=.TRUE., &
         reorder_rs_grid_ranks=reorder_rs_grid_ranks, skip_load_balance_distributed=skip_load_balance_distributed,&
         error=error)

    ! If the rank order has changed, reallocate any of the distributed rs_grids

    IF (reorder_rs_grid_ranks) THEN
       DO i=1,gridlevel_info%ngrid_levels
          IF (rs_descs(i)%rs_desc%distributed) THEN
             CALL rs_grid_release(rs_grids(i)%rs_grid, error)
             NULLIFY(rs_grids(i)%rs_grid)
             CALL rs_grid_create(rs_grids(i)%rs_grid, rs_descs(i)%rs_desc, error)
          END IF
       END DO
    END IF

    ! Now we have the final list of tasks, setup the task_list with the
    ! data needed for the loops in integrate_v/calculate_rho

    IF (ASSOCIATED(task_list%taskstart)) THEN
       DEALLOCATE (task_list%taskstart,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    IF (ASSOCIATED(task_list%taskstop)) THEN
       DEALLOCATE (task_list%taskstop,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    IF (ASSOCIATED(task_list%npairs)) THEN
       DEALLOCATE (task_list%npairs,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

    ! First, count the number of unique atom pairs per grid level

    ALLOCATE (task_list%npairs(SIZE(rs_descs)),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    iatom_old = -1 ; jatom_old = -1 ; igrid_level_old = -1 ; img_old = -1
    ipair = 0
    task_list%npairs = 0

    DO i=1,task_list%ntasks
       CALL int2pair(task_list%tasks(3,i),igrid_level,iatom,jatom,iset,jset,ipgf,jpgf,&
                     natoms,maxset,maxpgf)
       img = INT(task_list%tasks(7,i))
       IF (igrid_level .NE. igrid_level_old) THEN
          IF (igrid_level_old .NE. -1) THEN
             task_list%npairs(igrid_level_old) = ipair
          END IF
          ipair = 1
          igrid_level_old = igrid_level
          iatom_old = iatom
          jatom_old = jatom
          img_old   = img
       ELSE IF (iatom .NE. iatom_old .OR. jatom .NE. jatom_old .OR. img .NE. img_old) THEN
          ipair = ipair + 1
          iatom_old = iatom
          jatom_old = jatom
          img_old   = img
       END IF
    END DO
    ! Take care of the last iteration
    IF (task_list%ntasks /= 0) THEN
       task_list%npairs(igrid_level) = ipair
    END IF

    ! Second, for each atom pair, find the indices in the task list
    ! of the first and last task

    ! Array sized for worst case
    ALLOCATE (task_list%taskstart(MAXVAL(task_list%npairs),SIZE(rs_descs)), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (task_list%taskstop(MAXVAL(task_list%npairs),SIZE(rs_descs)), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    iatom_old = -1 ; jatom_old = -1 ; igrid_level_old = -1 ; img_old = -1
    ipair = 0
    task_list%taskstart = 0
    task_list%taskstop = 0

    DO i=1,task_list%ntasks
       CALL int2pair(task_list%tasks(3,i),igrid_level,iatom,jatom,iset,jset,ipgf,jpgf,&
                     natoms,maxset,maxpgf)
       img = INT(task_list%tasks(7,i))
       IF (igrid_level .NE. igrid_level_old) THEN
          IF (igrid_level_old .NE. -1) THEN
             task_list%taskstop(ipair,igrid_level_old) = i - 1
          END IF
          ipair = 1
          task_list%taskstart(ipair,igrid_level) = i
          igrid_level_old = igrid_level
          iatom_old = iatom
          jatom_old = jatom
          img_old   = img
       ELSE IF (iatom .NE. iatom_old .OR. jatom .NE. jatom_old .OR. img .NE. img_old) THEN
          ipair = ipair + 1
          task_list%taskstart(ipair,igrid_level) = i
          task_list%taskstop(ipair-1,igrid_level) = i - 1
          iatom_old = iatom
          jatom_old = jatom
          img_old   = img
       END IF
    END DO
    ! Take care of the last iteration
    IF (task_list%ntasks /= 0) THEN
       task_list%taskstop(ipair,igrid_level) = task_list%ntasks
    END IF

    CALL timestop(handle)

  END SUBROUTINE generate_qs_task_list

! *****************************************************************************
!> \brief ...
!> \param tasks ...
!> \param dist_ab ...
!> \param ntasks ...
!> \param curr_tasks ...
!> \param rs_descs ...
!> \param dft_control ...
!> \param cube_info ...
!> \param gridlevel_info ...
!> \param cindex ...
!> \param iatom ...
!> \param jatom ...
!> \param rpgfa ...
!> \param rpgfb ...
!> \param zeta ...
!> \param zetb ...
!> \param kind_radius_b ...
!> \param set_radius_a ...
!> \param set_radius_b ...
!> \param ra ...
!> \param rab ...
!> \param la_max ...
!> \param la_min ...
!> \param lb_max ...
!> \param lb_min ...
!> \param npgfa ...
!> \param npgfb ...
!> \param maxpgf ...
!> \param maxset ...
!> \param natoms ...
!> \param nseta ...
!> \param nsetb ...
!> \param error ...
!> \par History
!>      Joost VandeVondele: 10.2008 refactored
! *****************************************************************************
  SUBROUTINE task_list_inner_loop(tasks, dist_ab, ntasks, curr_tasks, rs_descs,dft_control,&
             cube_info,gridlevel_info,cindex,&
             iatom,jatom,rpgfa,rpgfb,zeta,zetb,kind_radius_b,set_radius_a,set_radius_b,ra,rab,&
             la_max,la_min,lb_max,lb_min,npgfa,npgfb,maxpgf,maxset,natoms,nseta,nsetb,error)

    INTEGER(KIND=int_8), DIMENSION(:, :), &
      POINTER                                :: tasks
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: dist_ab
    INTEGER                                  :: ntasks, curr_tasks
    TYPE(realspace_grid_desc_p_type), &
      DIMENSION(:), POINTER                  :: rs_descs
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(cube_info_type), DIMENSION(:), &
      POINTER                                :: cube_info
    TYPE(gridlevel_info_type), POINTER       :: gridlevel_info
    INTEGER                                  :: cindex, iatom, jatom
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: rpgfa, rpgfb, zeta, zetb
    REAL(KIND=dp)                            :: kind_radius_b
    REAL(KIND=dp), DIMENSION(:), POINTER     :: set_radius_a, set_radius_b
    REAL(KIND=dp), DIMENSION(3)              :: ra, rab
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, lb_max, &
                                                lb_min, npgfa, npgfb
    INTEGER                                  :: maxpgf, maxset, natoms, &
                                                nseta, nsetb
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'task_list_inner_loop', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: cube_center(3), igrid_level, &
                                                ipgf, iset, jpgf, jset, &
                                                lb_cube(3), ub_cube(3)
    REAL(KIND=dp)                            :: dab, rab2, radius, zetp

    rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
    dab = SQRT(rab2)

    loop_iset: DO iset=1,nseta

       IF (set_radius_a(iset) + kind_radius_b < dab) CYCLE

       loop_jset: DO jset=1,nsetb

          IF (set_radius_a(iset) + set_radius_b(jset) < dab) CYCLE

          loop_ipgf: DO ipgf=1,npgfa(iset)

             IF (rpgfa(ipgf,iset) + set_radius_b(jset) < dab) CYCLE

             loop_jpgf: DO jpgf=1,npgfb(jset)

                IF (rpgfa(ipgf,iset) + rpgfb(jpgf,jset) < dab) CYCLE

                zetp = zeta(ipgf,iset) + zetb(jpgf,jset)
                igrid_level = gaussian_gridlevel(gridlevel_info,zetp)

                CALL compute_pgf_properties(cube_center,lb_cube,ub_cube,radius,&
                           rs_descs(igrid_level)%rs_desc,cube_info(igrid_level),&
                           la_max(iset),zeta(ipgf,iset),la_min(iset),&
                           lb_max(jset),zetb(jpgf,jset),lb_min(jset),&
                           ra,rab,rab2,dft_control%qs_control%eps_rho_rspace,error)

                CALL pgf_to_tasks(tasks,dist_ab,ntasks,curr_tasks,&
                          rab,cindex,iatom,jatom,iset,jset,ipgf,jpgf,natoms,maxset,maxpgf,&
                          la_max(iset),lb_max(jset),rs_descs(igrid_level)%rs_desc,&
                          igrid_level,gridlevel_info%ngrid_levels,cube_center,&
                          lb_cube,ub_cube,error)

             END DO loop_jpgf

          END DO loop_ipgf

       END DO loop_jset

    END DO loop_iset

  END SUBROUTINE task_list_inner_loop

! *****************************************************************************
!> \brief combines the calculation of several basic properties of a given pgf:
!>  its center, the bounding cube, the radius, the cost,
!>  tries to predict the time needed for processing this task
!>      in this way an improved load balance might be obtained
!> \param cube_center ...
!> \param lb_cube ...
!> \param ub_cube ...
!> \param radius ...
!> \param rs_desc ...
!> \param cube_info ...
!> \param la_max ...
!> \param zeta ...
!> \param la_min ...
!> \param lb_max ...
!> \param zetb ...
!> \param lb_min ...
!> \param ra ...
!> \param rab ...
!> \param rab2 ...
!> \param eps ...
!> \param error ...
!> \par History
!>      10.2008 refactored [Joost VandeVondele]
!> \note
!>      -) this requires the radius to be computed in the same way as
!>      collocate_pgf_product_rspace, we should factor that part into a subroutine
!>      -) we're assuming that integrate_pgf and collocate_pgf are the same cost for load balancing
!>         this is more or less true for map_consistent
!>      -) in principle, the computed radius could be recycled in integrate_pgf/collocate_pgf if it is certainly
!>         the same, this could lead to a small speedup
!>      -) the cost function is a fit through the median cost of mapping a pgf with a given l and a given radius (in grid points)
!>         fitting the measured data on an opteron/g95 using the expression
!>         a*(l+b)(r+c)**3+d which is based on the innerloop of the collocating routines
! *****************************************************************************
  SUBROUTINE compute_pgf_properties(cube_center,lb_cube,ub_cube,radius, &
               rs_desc,cube_info,la_max,zeta,la_min,lb_max,zetb,lb_min,ra,rab,rab2,eps,error)

    INTEGER, DIMENSION(3), INTENT(OUT)       :: cube_center, lb_cube, ub_cube
    REAL(KIND=dp), INTENT(OUT)               :: radius
    TYPE(realspace_grid_desc_type), POINTER  :: rs_desc
    TYPE(cube_info_type), INTENT(IN)         :: cube_info
    INTEGER, INTENT(IN)                      :: la_max
    REAL(KIND=dp), INTENT(IN)                :: zeta
    INTEGER, INTENT(IN)                      :: la_min, lb_max
    REAL(KIND=dp), INTENT(IN)                :: zetb
    INTEGER, INTENT(IN)                      :: lb_min
    REAL(KIND=dp), INTENT(IN)                :: ra(3), rab(3), rab2, eps
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_pgf_properties', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: extent(3)
    INTEGER, DIMENSION(:), POINTER           :: sphere_bounds
    REAL(KIND=dp)                            :: cutoff, f, prefactor, rb(3), &
                                                zetp
    REAL(KIND=dp), DIMENSION(3)              :: rp

! the radius for this task

    zetp = zeta + zetb
    rp(:) = ra(:) + zetb/zetp*rab(:)
    rb(:) = ra(:)+rab(:)
    cutoff    = 1.0_dp
    f         = zetb/zetp
    prefactor = EXP(-zeta*f*rab2)
    radius=exp_radius_very_extended(la_min,la_max,lb_min,lb_max,ra=ra,rb=rb,rp=rp,&
               zetp=zetp,eps=eps, prefactor=prefactor,cutoff=cutoff)

    CALL compute_cube_center(cube_center,rs_desc,zeta,zetb,ra,rab)
    ! compute cube_center, the center of the gaussian product to map (folded to within the unit cell)
    cube_center(:) = MODULO ( cube_center(:), rs_desc%npts(:) )
    cube_center(:) = cube_center(:) + rs_desc%lb(:)

    IF (rs_desc%orthorhombic ) THEN
       CALL return_cube(cube_info,radius,lb_cube,ub_cube,sphere_bounds)
    ELSE
       CALL return_cube_nonortho(cube_info,radius,lb_cube,ub_cube,rp)
       !? unclear if extent is computed correctly.
       extent(:)  = ub_cube(:) - lb_cube(:)
       lb_cube(:) = -extent(:) / 2 - 1
       ub_cube(:) =  extent(:) / 2
    ENDIF

  END SUBROUTINE compute_pgf_properties
! *****************************************************************************
!> \brief predicts the cost of a task in kcycles for a given task
!>        the model is based on a fit of actual data, and might need updating
!>        as collocate_pgf_product_rspace changes (or CPUs/compilers change)
!>        maybe some dynamic approach, improving the cost model on the fly could
!>        work as well
!>        the cost model does not yet take into account the fraction of space
!>        that is mapped locally for a given cube and rs_grid (generalised tasks)
!> \param lb_cube ...
!> \param ub_cube ...
!> \param fraction ...
!> \param lmax ...
!> \param is_ortho ...
!> \retval cost_model ...
! *****************************************************************************
  INTEGER FUNCTION cost_model(lb_cube,ub_cube,fraction,lmax,is_ortho)
    INTEGER, DIMENSION(3), INTENT(IN)       :: lb_cube, ub_cube
    REAL(KIND=dp), INTENT(IN)               :: fraction

    INTEGER :: cmax,lmax
    LOGICAL :: is_ortho
    REAL(KIND=dp) :: v1,v2,v3,v4,v5

    cmax=MAXVAL(((ub_cube-lb_cube)+1)/2)

    IF (is_ortho) THEN
       v1=1.504760E+00_dp
       v2=3.126770E+00_dp
       v3=5.074106E+00_dp
       v4=1.091568E+00_dp
       v5=1.070187E+00_dp
    ELSE
       v1=7.831105E+00_dp
       v2=2.675174E+00_dp
       v3=7.546553E+00_dp
       v4=6.122446E-01_dp
       v5=3.886382E+00_dp
    ENDIF
    cost_model = CEILING(((lmax+v1)*(cmax+v2)**3*v3*fraction+v4+v5*lmax**7)/1000.0_dp)

  END FUNCTION cost_model
! *****************************************************************************
!> \brief pgf_to_tasks converts a given pgf to one or more tasks, in particular
!>        this determines by which CPUs a given pgf gets collocated
!>        the format of the task array is as follows
!>        tasks(1,i) := destination
!>        tasks(2,i) := source
!>        tasks(3,i) := compressed type (iatom, jatom, ....)
!>        tasks(4,i) := type (0: replicated, 1: distributed local, 2: distributed generalised)
!>        tasks(5,i) := cost
!>        tasks(6,i) := alternate destination code (0 if none availabe)
!>        tasks(7,i) := cell index of matrix associated with this task (always 1 for no k-points)
!>
!> \param tasks ...
!> \param dist_ab ...
!> \param ntasks ...
!> \param curr_tasks ...
!> \param rab ...
!> \param cindex ...
!> \param iatom ...
!> \param jatom ...
!> \param iset ...
!> \param jset ...
!> \param ipgf ...
!> \param jpgf ...
!> \param natoms ...
!> \param maxset ...
!> \param maxpgf ...
!> \param la_max ...
!> \param lb_max ...
!> \param rs_desc ...
!> \param igrid_level ...
!> \param n_levels ...
!> \param cube_center ...
!> \param lb_cube ...
!> \param ub_cube ...
!> \param error ...
!> \par History
!>      10.2008 Refactored based on earlier routines by MattW [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE pgf_to_tasks(tasks,dist_ab,ntasks,curr_tasks,&
             rab,cindex,iatom,jatom,iset,jset,ipgf,jpgf,natoms,maxset,maxpgf,&
             la_max,lb_max,rs_desc,igrid_level,n_levels,cube_center,lb_cube,ub_cube,error)

    INTEGER(KIND=int_8), DIMENSION(:, :), &
      POINTER                                :: tasks
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: dist_ab
    INTEGER, INTENT(INOUT)                   :: ntasks, curr_tasks
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rab
    INTEGER, INTENT(IN)                      :: cindex, iatom, jatom, iset, &
                                                jset, ipgf, jpgf, natoms, &
                                                maxset, maxpgf, la_max, lb_max
    TYPE(realspace_grid_desc_type), POINTER  :: rs_desc
    INTEGER, INTENT(IN)                      :: igrid_level, n_levels
    INTEGER, DIMENSION(3), INTENT(IN)        :: cube_center, lb_cube, ub_cube
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER, PARAMETER                       :: add_tasks = 1000
    REAL(kind=dp), PARAMETER                 :: mult_tasks = 2.0_dp

    INTEGER                                  :: added_tasks, cost, j, lmax
    LOGICAL                                  :: is_ortho
    REAL(KIND=dp)                            :: tfraction

    ntasks = ntasks + 1
    IF ( ntasks > curr_tasks ) THEN
       curr_tasks = INT((curr_tasks+add_tasks)*mult_tasks)
       CALL reallocate(tasks,1,7,1,curr_tasks)
    END IF

    IF(rs_desc%distributed) THEN

       ! finds the node(s) that need to process this task
       ! on exit tasks(4,:) is 1 for distributed tasks and 2 for generalised tasks
       CALL rs_find_node(rs_desc,igrid_level,n_levels,cube_center, &
            ntasks=ntasks,tasks=tasks,lb_cube=lb_cube,ub_cube=ub_cube,added_tasks=added_tasks,error=error)

    ELSE
       tasks(1,ntasks) = encode_rank(rs_desc % my_pos, igrid_level, n_levels)
       tasks(4,ntasks) = 0
       tasks(6,ntasks) = 0
       added_tasks = 1
    ENDIF

    IF (SIZE(dist_ab,2).NE.SIZE(tasks,2)) THEN
       CALL reallocate(dist_ab,1,3,1,SIZE(tasks,2))
    ENDIF

    lmax=la_max+lb_max
    is_ortho=(tasks(4,ntasks)==0 .OR. tasks(4,ntasks)==1) .AND. rs_desc%orthorhombic
    ! we assume the load is shared equally between processes dealing with a generalised Gaussian.
    ! this could be refined in the future
    tfraction=1.0_dp/added_tasks

    cost=cost_model(lb_cube,ub_cube,tfraction,lmax,is_ortho)

    DO j= 1, added_tasks

       tasks (2,ntasks-added_tasks+j) = encode_rank(rs_desc % my_pos, igrid_level, n_levels)
       tasks (5,ntasks-added_tasks+j) = cost
       tasks (7,ntasks-added_tasks+j) = cindex

       !encode the atom pairs and basis info as a single long integer
       CALL pair2int(tasks(3,ntasks-added_tasks+j),igrid_level,&
            iatom,jatom,iset,jset,ipgf,jpgf,natoms,maxset,maxpgf)

       dist_ab (1,ntasks-added_tasks+j) = rab(1)
       dist_ab (2,ntasks-added_tasks+j) = rab(2)
       dist_ab (3,ntasks-added_tasks+j) = rab(3)

    ENDDO

  END SUBROUTINE pgf_to_tasks

! *****************************************************************************
!> \brief converts a pgf index pair (ipgf,iset,iatom),(jpgf,jset,jatom)
!>      to a unique integer.
!>      a list of integers can be sorted, and will result in a list of pgf pairs
!>      for which all atom pairs are grouped, and for each atom pair all set pairs are grouped
!>      and for each set pair, all pgfs are grouped
!> \param res ...
!> \param ilevel ...
!> \param iatom ...
!> \param jatom ...
!> \param iset ...
!> \param jset ...
!> \param ipgf ...
!> \param jpgf ...
!> \param natom ...
!> \param maxset ...
!> \param maxpgf ...
!> \par History
!>      11.2007 created [Joost]
!> \note
!>      will hopefully not overflow any more
! *****************************************************************************
  SUBROUTINE pair2int(res,ilevel,iatom,jatom,iset,jset,ipgf,jpgf,natom,maxset,maxpgf)
    INTEGER(KIND=int_8), INTENT(OUT)         :: res
    INTEGER, INTENT(IN)                      :: ilevel, iatom, jatom, iset, &
                                                jset, ipgf, jpgf, natom, &
                                                maxset, maxpgf

    INTEGER(KIND=int_8)                      :: maxpgf8, maxset8, natom8

    natom8=natom ; maxset8=maxset ; maxpgf8=maxpgf
    !
    ! this encoding yields the right order of the tasks for collocation after the sort
    ! in distribute_tasks. E.g. for a atom pair, all sets and pgfs are computed in one go.
    ! The exception is the gridlevel. Tasks are first ordered wrt to grid_level. This implies
    ! that a given density matrix block will be decontracted several times, but cache effects on the
    ! grid make up for this.
    !
    res=ilevel*natom8**2*maxset8**2*maxpgf8**2+&
        ((iatom-1)*natom8+jatom-1)*maxset8**2*maxpgf8**2 + &
        ((iset-1)*maxset8+jset-1)*maxpgf8**2+ &
         (ipgf-1)*maxpgf8+jpgf-1

  END SUBROUTINE pair2int

! *****************************************************************************
!> \brief ...
!> \param res ...
!> \param ilevel ...
!> \param iatom ...
!> \param jatom ...
!> \param iset ...
!> \param jset ...
!> \param ipgf ...
!> \param jpgf ...
!> \param natom ...
!> \param maxset ...
!> \param maxpgf ...
! *****************************************************************************
  SUBROUTINE int2pair(res,ilevel,iatom,jatom,iset,jset,ipgf,jpgf,natom,maxset,maxpgf)
    INTEGER(KIND=int_8), INTENT(IN)          :: res
    INTEGER, INTENT(OUT)                     :: ilevel, iatom, jatom, iset, &
                                                jset, ipgf, jpgf
    INTEGER, INTENT(IN)                      :: natom, maxset, maxpgf

    INTEGER(KIND=int_8)                      :: iatom8, ijatom, ijpgf, ijset, &
                                                ipgf8, iset8, jatom8, jpgf8, &
                                                jset8, maxpgf8, maxset8, &
                                                natom8, tmp

    natom8=natom ; maxset8=maxset ; maxpgf8=maxpgf

    ilevel=INT(res/(natom8**2*maxset8**2*maxpgf8**2))
    tmp=MOD(res,natom8**2*maxset8**2*maxpgf8**2)
    ijatom=tmp/(maxpgf8**2*maxset8**2)
    iatom8=ijatom/natom8+1
    jatom8=MOD(ijatom,natom8)+1
    tmp=MOD(tmp,maxpgf8**2*maxset8**2)
    ijset=tmp/maxpgf8**2
    iset8=ijset/maxset8+1
    jset8=MOD(ijset,maxset8)+1
    ijpgf=MOD(tmp,maxpgf8**2)
    ipgf8=ijpgf/maxpgf8+1
    jpgf8=MOD(ijpgf,maxpgf8)+1

    iatom=INT(iatom8) ; jatom=INT(jatom8); iset=INT(iset8) ; jset=INT(jset8)
    ipgf=INT(ipgf8) ; jpgf=INT(jpgf8)

  END SUBROUTINE int2pair

! *****************************************************************************
!> \brief performs load balancing of the tasks on the distributed grids
!> \param tasks ...
!> \param ntasks ...
!> \param rs_descs ...
!> \param grid_level ...
!> \param natoms ...
!> \param maxset ...
!> \param maxpgf ...
!> \param error ...
!> \par History
!>      created 2008-10-03 [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE load_balance_distributed (tasks, ntasks, rs_descs, grid_level, &
                                       natoms, maxset, maxpgf, error)

    INTEGER(KIND=int_8), DIMENSION(:, :), &
      POINTER                                :: tasks
    INTEGER                                  :: ntasks
    TYPE(realspace_grid_desc_p_type), &
      DIMENSION(:), POINTER                  :: rs_descs
    INTEGER                                  :: grid_level, natoms, maxset, &
                                                maxpgf
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'load_balance_distributed', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, stat
    INTEGER, DIMENSION(:, :, :), POINTER     :: list
    LOGICAL                                  :: failure

    CALL timeset(routineN,handle)

    failure=.FALSE.
    NULLIFY(list)
    ! here we create for each cpu (0:ncpu-1) a list of possible destinations.
    ! if a destination would not be in this list, it is a bug
    CALL create_destination_list(list,rs_descs,grid_level,error)

    ! now, walk over the tasks, filling in the loads of each destination
    CALL compute_load_list(list, rs_descs, grid_level, tasks, ntasks, natoms, maxset, maxpgf,&
                           create_list=.TRUE., error=error)

    ! optimize loads & fluxes
    CALL optimize_load_list(list,rs_descs(1)%rs_desc%group,rs_descs(1)%rs_desc%my_pos, error)

    ! now, walk over the tasks, using the list to set the destinations
    CALL compute_load_list(list, rs_descs,grid_level, tasks, ntasks, natoms, maxset, maxpgf,&
                           create_list=.FALSE., error=error)

    DEALLOCATE (list,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE load_balance_distributed


! *****************************************************************************
!> \brief this serial routine adjusts the fluxes in the global list
!>
!> \param list_global ...
!> \par History
!>      created 2008-10-06 [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE balance_global_list(list_global)
    INTEGER, DIMENSION(:, :, 0:)             :: list_global

    CHARACTER(LEN=*), PARAMETER :: routineN = 'balance_global_list', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: Max_Iter = 100
    REAL(KIND=dp), PARAMETER                 :: Tolerance_factor = 0.005_dp

    INTEGER                                  :: dest, handle, icpu, idest, &
                                                iflux, ilocal, k, maxdest, &
                                                Ncpu, Nflux
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: flux_connections
    LOGICAL                                  :: solution_optimal
    REAL(KIND=dp)                            :: average, load_shift, &
                                                max_load_shift, tolerance
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: load, optimized_flux, &
                                                optimized_load
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: flux_limits

    CALL timeset(routineN,handle)

     Ncpu=SIZE(list_global,3)
     maxdest=SIZE(list_global,2)
     ALLOCATE (load(0:Ncpu-1))
     load=0.0_dp
     ALLOCATE (optimized_load(0:Ncpu-1))

     ! figure out the number of fluxes
     ! we assume that the global_list is symmetric
     Nflux=0
     DO icpu=0,ncpu-1
        DO idest=1,maxdest
           dest=list_global(1,idest,icpu)
           IF (dest<ncpu .AND. dest>icpu) Nflux=Nflux+1
        ENDDO
     ENDDO
     ALLOCATE (optimized_flux(Nflux))
     ALLOCATE (flux_limits(2,Nflux))
     ALLOCATE (flux_connections(2,Nflux))

     ! reorder data
     flux_limits=0
     Nflux=0
     DO icpu=0,ncpu-1
        load(icpu)=SUM(list_global(2,:,icpu))
        DO idest=1,maxdest
           dest=list_global(1,idest,icpu)
           IF (dest<ncpu) THEN
              IF (dest.NE.icpu) THEN
                 IF (dest>icpu) THEN
                    Nflux=Nflux+1
                    flux_limits(2,Nflux)=list_global(2,idest,icpu)
                    flux_connections(1,Nflux)=icpu
                    flux_connections(2,Nflux)=dest
                 ELSE
                    DO iflux=1,Nflux
                       IF (flux_connections(1,iflux)==dest .AND. flux_connections(2,iflux)==icpu) THEN
                          flux_limits(1,iflux)=-list_global(2,idest,icpu)
                          EXIT
                       ENDIF
                    ENDDO
                 ENDIF
              ENDIF
           ENDIF
        ENDDO
     ENDDO

     solution_optimal=.FALSE.
     optimized_flux=0.0_dp

     ! an iterative solver, if iterated till convergence the maximum load is minimal
     ! we terminate before things are fully converged, since this does show up in the timings
     ! once the largest shift becomes less than a small fraction of the average load, we're done
     ! we're perfectly happy if the load balance is within 1 percent or so
     ! the maximum load normally converges even faster
     average=SUM(load)/SIZE(load)
     tolerance=Tolerance_factor*average

     optimized_load(:)=load
     DO k=1,Max_iter
       max_load_shift=0.0_dp
       DO iflux=1,Nflux
          load_shift=(optimized_load(flux_connections(1,iflux))-optimized_load(flux_connections(2,iflux)))/2
          load_shift=MAX(flux_limits(1,iflux)-optimized_flux(iflux),load_shift)
          load_shift=MIN(flux_limits(2,iflux)-optimized_flux(iflux),load_shift)
          max_load_shift=MAX(ABS(load_shift),max_load_shift)
          optimized_load(flux_connections(1,iflux))=optimized_load(flux_connections(1,iflux))-load_shift
          optimized_load(flux_connections(2,iflux))=optimized_load(flux_connections(2,iflux))+load_shift
          optimized_flux(iflux)=optimized_flux(iflux)+load_shift
       ENDDO
       IF (max_load_shift<tolerance) THEN
          solution_optimal=.TRUE.
          EXIT
       ENDIF
     ENDDO

     ! now adjust the load list to reflect the optimized fluxes
     ! reorder data
     Nflux=0
     DO icpu=0,ncpu-1
        DO idest=1,maxdest
           IF (list_global(1,idest,icpu)==icpu) ilocal=idest
        ENDDO
        DO idest=1,maxdest
           dest=list_global(1,idest,icpu)
           IF (dest<ncpu) THEN
              IF (dest.NE.icpu) THEN
                 IF (dest>icpu) THEN
                    Nflux=Nflux+1
                    IF (optimized_flux(Nflux)>0) THEN
                       list_global(2,ilocal,icpu)=list_global(2,ilocal,icpu)+&
                                                  list_global(2,idest,icpu)-NINT(optimized_flux(Nflux))
                       list_global(2,idest,icpu)=NINT(optimized_flux(Nflux))
                    ELSE
                       list_global(2,ilocal,icpu)=list_global(2,ilocal,icpu)+&
                                                  list_global(2,idest,icpu)
                       list_global(2,idest,icpu)=0
                    ENDIF
                 ELSE
                    DO iflux=1,Nflux
                       IF (flux_connections(1,iflux)==dest .AND. flux_connections(2,iflux)==icpu) THEN
                          IF (optimized_flux(iflux)>0) THEN
                             list_global(2,ilocal,icpu)=list_global(2,ilocal,icpu)+&
                                                        list_global(2,idest,icpu)
                             list_global(2,idest,icpu)=0
                          ELSE
                             list_global(2,ilocal,icpu)=list_global(2,ilocal,icpu)+&
                                                        list_global(2,idest,icpu)+NINT(optimized_flux(iflux))
                             list_global(2,idest,icpu)=-NINT(optimized_flux(iflux))
                          ENDIF
                          EXIT
                       ENDIF
                    ENDDO
                 ENDIF
              ENDIF
           ENDIF
        ENDDO
     ENDDO

    CALL timestop(handle)

  END SUBROUTINE balance_global_list

! *****************************************************************************
!> \brief this routine gets back optimized loads for all destinations
!>
!> \param list ...
!> \param group ...
!> \param my_pos ...
!> \param error ...
!> \par History
!>      created 2008-10-06 [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE optimize_load_list(list,group,my_pos,error)
    INTEGER, DIMENSION(:, :, 0:)             :: list
    INTEGER, INTENT(IN)                      :: group, my_pos
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'optimize_load_list', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, icpu, idest, itmp, &
                                                jcpu, maxdest, ncpu, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: load_all, load_per_source
    INTEGER, ALLOCATABLE, DIMENSION(:, :, :) :: list_global
    LOGICAL                                  :: failure

    CALL timeset(routineN,handle)

    failure=.FALSE.

    ncpu=SIZE(list,3)
    maxdest=SIZE(list,2)

    ALLOCATE (load_all(maxdest*ncpu),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    !
    ! XXXXX this allocation can cause an out-of-memory if the number of MPI tasks is sufficient
    !
    ALLOCATE (load_per_source(maxdest*ncpu*ncpu),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (list_global(2,maxdest,ncpu),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    load_all(:)=RESHAPE(list(2,:,:),(/maxdest*ncpu/))

    CALL mp_gather(load_all,load_per_source,0,group)

    IF (my_pos==0) THEN

       ! collect the load of all CPUs
       load_all=0
       DO icpu=1,ncpu
          load_all(:)=load_all+load_per_source(maxdest*ncpu*(icpu-1)+1:maxdest*ncpu*icpu)
       ENDDO

       list_global(1,:,:)=list(1,:,:)
       list_global(2,:,:)=RESHAPE(load_all,(/maxdest,ncpu/))

       ! optimize the fluxes
       CALL balance_global_list(list_global)

       ! figure out the amount of data a given CPU can send to a destination
       !
       i=0
       DO jcpu=1,ncpu
         DO  icpu=1,ncpu
           DO idest=1,maxdest
              i=i+1
              itmp=MIN(list_global(2,idest,icpu),load_per_source(i))
              load_per_source(i)=itmp
              list_global(2,idest,icpu)=list_global(2,idest,icpu)-itmp
           ENDDO
         ENDDO
       ENDDO

    ENDIF

    CALL mp_scatter(load_per_source,load_all,0,group)
    list(2,:,:)=RESHAPE(load_all,(/maxdest,ncpu/))

    DEALLOCATE (load_all,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (load_per_source,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (list_global,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE optimize_load_list
! *****************************************************************************
!> \brief fill the load list with values derived from the tasks array
!>        from the alternate locations, we select the alternate location that
!>        can be used without increasing the number of matrix blocks needed to
!>        distribute.
!>        Replicated tasks are not yet considered
!>
!> \param list ...
!> \param rs_descs ...
!> \param grid_level ...
!> \param tasks ...
!> \param ntasks ...
!> \param natoms ...
!> \param maxset ...
!> \param maxpgf ...
!> \param create_list ...
!> \param error ...
!> \par History
!>      created 2008-10-06 [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE compute_load_list(list,rs_descs,grid_level,tasks,&
                               ntasks,natoms,maxset,maxpgf,create_list,error)
    INTEGER, DIMENSION(:, :, 0:)             :: list
    TYPE(realspace_grid_desc_p_type), &
      DIMENSION(:), POINTER                  :: rs_descs
    INTEGER                                  :: grid_level
    INTEGER(KIND=int_8), DIMENSION(:, :), &
      POINTER                                :: tasks
    INTEGER                                  :: ntasks, natoms, maxset, maxpgf
    LOGICAL                                  :: create_list
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_load_list', &
      routineP = moduleN//':'//routineN

    INTEGER :: cost, dest, handle, i, iatom, ilevel, img, img_old, iopt, &
      ipgf, iset, itask, itask_start, itask_stop, jatom, jpgf, jset, li, &
      maxdest, ncpu, ndest_pair, nopt, nshort, rank
    INTEGER(KIND=int_8)                      :: bit_pattern, ipair, &
                                                ipair_old, natom8
    INTEGER(KIND=int_8), ALLOCATABLE, &
      DIMENSION(:)                           :: loads
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: all_dests, index
    INTEGER, DIMENSION(6)                    :: options
    LOGICAL                                  :: failure

    CALL timeset(routineN,handle)

    ALLOCATE (loads(0:rs_descs(grid_level)%rs_desc%group_size-1))
    CALL get_current_loads(loads, rs_descs, grid_level, ntasks, natoms, maxset, maxpgf, &
                           tasks, use_reordered_ranks=.FALSE. , error=error)

    failure=.FALSE.

    maxdest=SIZE(list,2)
    ncpu=SIZE(list,3)
    natom8 =natoms

    ! first find the tasks that deal with the same atom pair
    itask_stop=0
    ipair_old=HUGE(ipair_old)
    img_old = -1
    ALLOCATE (all_dests(0))
    ALLOCATE (INDEX(0))

    DO

       ! first find the range of tasks that deal with the same atom pair
       itask_start=itask_stop+1
       itask_stop =itask_start
       IF(itask_stop>ntasks) EXIT
       CALL int2pair(tasks(3,itask_stop),ilevel,iatom,jatom,iset,jset,ipgf,jpgf,natoms,maxset,maxpgf)
       img_old=INT(tasks(7,itask_stop))
       ipair_old=(iatom-1)*natom8+(jatom-1)
       DO
          IF(itask_stop+1>ntasks) EXIT
          CALL int2pair(tasks(3,itask_stop+1),ilevel,iatom,jatom,iset,jset,ipgf,jpgf,natoms,maxset,maxpgf)
          ipair=(iatom-1)*natom8+(jatom-1)
          img=INT(tasks(7,itask_stop+1))
          IF (ipair==ipair_old .AND. img==img_old) THEN
             itask_stop=itask_stop+1
          ELSE
             EXIT
          ENDIF
       ENDDO
       ipair=ipair_old
       nshort=itask_stop-itask_start+1

       ! find the unique list of destinations on this grid level only
       DEALLOCATE (all_dests)
       ALLOCATE (all_dests(nshort))
       DEALLOCATE (index)
       ALLOCATE (INDEX(nshort))
       DO i=1,nshort
         CALL int2pair(tasks(3,itask_start + i - 1),ilevel,iatom,jatom,iset,jset,ipgf,jpgf,natoms,maxset,maxpgf)
           IF (ilevel .EQ. grid_level) THEN
             all_dests(i)=decode_rank(tasks(1,itask_start + i - 1), SIZE(rs_descs))
           ELSE
             all_dests(i)=HUGE(all_dests(i))
           END IF
       ENDDO
       CALL sort(all_dests,nshort,index)
       ndest_pair=1
       DO i=2,nshort
          IF ((all_dests(ndest_pair).NE.all_dests(i) ) .AND. (all_dests(i) .NE. HUGE(all_dests(i)))) THEN
             ndest_pair=ndest_pair+1
             all_dests(ndest_pair)=all_dests(i)
          ENDIF
       ENDDO

       DO itask=itask_start,itask_stop

          dest=decode_rank(tasks(1,itask), SIZE(rs_descs)) ! notice that dest can be changed
          CALL int2pair(tasks(3,itask),ilevel,iatom,jatom,iset,jset,ipgf,jpgf,natoms,maxset,maxpgf)
          ! Only proceed with tasks which are on this grid level
          IF (ilevel .NE. grid_level) CYCLE
          ipair=(iatom-1)*natom8+(jatom-1)
          cost=INT(tasks(5,itask))

          SELECT CASE(tasks(4,itask))
          CASE(1)
            bit_pattern=tasks(6,itask)
            nopt=0
            IF (BTEST(bit_pattern,0)) THEN
               rank=rs_grid_locate_rank(rs_descs(ilevel)%rs_desc,dest,(/-1,0,0/))
               IF (ANY(all_dests(1:ndest_pair).EQ.rank)) THEN
                  nopt=nopt+1
                  options(nopt)=rank
               ENDIF
            ENDIF
            IF (BTEST(bit_pattern,1)) THEN
               rank=rs_grid_locate_rank(rs_descs(ilevel)%rs_desc,dest,(/+1,0,0/))
               IF (ANY(all_dests(1:ndest_pair).EQ.rank)) THEN
                  nopt=nopt+1
                  options(nopt)=rank
               ENDIF
            ENDIF
            IF (BTEST(bit_pattern,2)) THEN
               rank=rs_grid_locate_rank(rs_descs(ilevel)%rs_desc,dest,(/0,-1,0/))
               IF (ANY(all_dests(1:ndest_pair).EQ.rank)) THEN
                  nopt=nopt+1
                  options(nopt)=rank
               ENDIF
            ENDIF
            IF (BTEST(bit_pattern,3)) THEN
               rank=rs_grid_locate_rank(rs_descs(ilevel)%rs_desc,dest,(/0,+1,0/))
               IF (ANY(all_dests(1:ndest_pair).EQ.rank)) THEN
                  nopt=nopt+1
                  options(nopt)=rank
               ENDIF
            ENDIF
            IF (BTEST(bit_pattern,4)) THEN
               rank=rs_grid_locate_rank(rs_descs(ilevel)%rs_desc,dest,(/0,0,-1/))
               IF (ANY(all_dests(1:ndest_pair).EQ.rank)) THEN
                  nopt=nopt+1
                  options(nopt)=rank
               ENDIF
            ENDIF
            IF (BTEST(bit_pattern,5)) THEN
               rank=rs_grid_locate_rank(rs_descs(ilevel)%rs_desc,dest,(/0,0,+1/))
               IF (ANY(all_dests(1:ndest_pair).EQ.rank)) THEN
                  nopt=nopt+1
                  options(nopt)=rank
               ENDIF
            ENDIF
            IF (nopt>0) THEN
              ! set it to the rank with the lowest load
              rank=options(1)
              DO iopt=2,nopt
                 IF (loads(rank)>loads(options(iopt))) rank=options(iopt)
              ENDDO
            ELSE
              rank=dest
            ENDIF
            li=list_index(list,rank,dest)
            IF (create_list) THEN
               list(2,li,dest)=list(2,li,dest)+cost
            ELSE
               IF (list(1,li,dest)==dest) THEN
                   tasks(1,itask)=encode_rank(dest, ilevel, SIZE(rs_descs))
               ELSE
                 IF (list(2,li,dest)>=cost) THEN
                     list(2,li,dest)=list(2,li,dest)-cost
                     tasks(1,itask)=encode_rank(list(1,li,dest), ilevel, SIZE(rs_descs))
                 ELSE
                     tasks(1,itask)=encode_rank(dest, ilevel, SIZE(rs_descs))
                 ENDIF
               ENDIF
            ENDIF
          CASE(2) ! generalised
            li=list_index(list,dest,dest)
            IF (create_list) THEN
               list(2,li,dest)=list(2,li,dest)+cost
            ELSE
               IF (list(1,li,dest)==dest) THEN
                   tasks(1,itask)=encode_rank(dest, ilevel, SIZE(rs_descs))
               ELSE
                 IF (list(2,li,dest)>=cost) THEN
                     list(2,li,dest)=list(2,li,dest)-cost
                     tasks(1,itask)=encode_rank(list(1,li,dest), ilevel, SIZE(rs_descs))
                 ELSE
                     tasks(1,itask)=encode_rank(dest, ilevel, SIZE(rs_descs))
                 ENDIF
               ENDIF
            ENDIF
          CASE DEFAULT
            CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
          END SELECT


       ENDDO

    ENDDO

    CALL timestop(handle)

  END SUBROUTINE compute_load_list
! *****************************************************************************
!> \brief small helper function to return the proper index in the list array
!>
!> \param list ...
!> \param rank ...
!> \param dest ...
!> \retval list_index ...
!> \par History
!>      created 2008-10-06 [Joost VandeVondele]
! *****************************************************************************
  INTEGER FUNCTION list_index(list,rank,dest)
    INTEGER, DIMENSION(:, :,0:), INTENT(IN)  :: list
    INTEGER, INTENT(IN)                      :: rank,dest
    list_index=1
    DO
      IF (list(1,list_index,dest)==rank) EXIT
      list_index=list_index+1
    ENDDO
  END FUNCTION list_index
! *****************************************************************************
!> \brief create a list with possible destinations (i.e. the central cpu and neighbors) for each cpu
!>        note that we allocate it with an additional field to store the load of this destination
!>
!> \param list ...
!> \param rs_descs ...
!> \param grid_level ...
!> \param error ...
!> \par History
!>      created 2008-10-06 [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE create_destination_list(list,rs_descs,grid_level,error)
    INTEGER, DIMENSION(:, :, :), POINTER     :: list
    TYPE(realspace_grid_desc_p_type), &
      DIMENSION(:), POINTER                  :: rs_descs
    INTEGER, INTENT(IN)                      :: grid_level
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'create_destination_list', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, icpu, j, maxcount, &
                                                ncpu, stat, ultimate_max
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: index, sublist
    LOGICAL                                  :: failure

    failure=.FALSE.

    CALL timeset(routineN,handle)

    CPPrecondition(.NOT.ASSOCIATED(list),cp_failure_level,routineP,error,failure)
    ncpu=rs_descs(grid_level)%rs_desc%group_size
    ultimate_max=7

    ALLOCATE (list(2,ultimate_max,0:ncpu-1),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE (INDEX(ultimate_max),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (sublist(ultimate_max),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    sublist=HUGE(sublist)

    maxcount=1
    DO icpu=0,ncpu-1
       sublist(1)=icpu
       sublist(2)=rs_grid_locate_rank(rs_descs(grid_level)%rs_desc,icpu,(/-1,0,0/))
       sublist(3)=rs_grid_locate_rank(rs_descs(grid_level)%rs_desc,icpu,(/+1,0,0/))
       sublist(4)=rs_grid_locate_rank(rs_descs(grid_level)%rs_desc,icpu,(/0,-1,0/))
       sublist(5)=rs_grid_locate_rank(rs_descs(grid_level)%rs_desc,icpu,(/0,+1,0/))
       sublist(6)=rs_grid_locate_rank(rs_descs(grid_level)%rs_desc,icpu,(/0,0,-1/))
       sublist(7)=rs_grid_locate_rank(rs_descs(grid_level)%rs_desc,icpu,(/0,0,+1/))
       ! only retain unique values of the destination
       CALL sort(sublist,ultimate_max,index)
       j=1
       DO i=2,7
          IF(sublist(i).NE.sublist(j)) THEN
             j=j+1
             sublist(j)=sublist(i)
          ENDIF
       ENDDO
       maxcount=MAX(maxcount,j)
       sublist(j+1:ultimate_max)=HUGE(sublist)
       list(1,:,icpu)=sublist
       list(2,:,icpu)=0
    ENDDO

    CALL reallocate(list,1,2,1,maxcount,0,ncpu-1)

    CALL timestop(handle)

  END SUBROUTINE create_destination_list

! *****************************************************************************
!> \brief given a task list, compute the load of each process everywhere
!>        giving this function the ability to loop over a (sub)set of rs_grids, 
!>        and do all the communication in one shot, would speed it up
!> \param loads ...
!> \param rs_descs ...
!> \param grid_level ...
!> \param ntasks ...
!> \param natom ...
!> \param maxset ...
!> \param maxpgf ...
!> \param tasks ...
!> \param use_reordered_ranks ...
!> \param error ...
!> \par History
!>      none
!> \author MattW 21/11/2007
! *****************************************************************************
  SUBROUTINE get_current_loads(loads, rs_descs, grid_level, ntasks, natom, maxset, &
                               maxpgf, tasks, use_reordered_ranks, error)
    INTEGER(KIND=int_8), DIMENSION(:)        :: loads
    TYPE(realspace_grid_desc_p_type), &
      DIMENSION(:), POINTER                  :: rs_descs
    INTEGER                                  :: grid_level, ntasks, natom, &
                                                maxset, maxpgf
    INTEGER(KIND=int_8), DIMENSION(:, :), &
      POINTER                                :: tasks
    LOGICAL, INTENT(IN)                      :: use_reordered_ranks
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_current_loads', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, iatom, ilevel, &
                                                ipgf, iset, jatom, jpgf, &
                                                jset, stat
    INTEGER(KIND=int_8)                      :: total_cost_local
    INTEGER(KIND=int_8), ALLOCATABLE, &
      DIMENSION(:)                           :: recv_buf_i, send_buf_i
    LOGICAL                                  :: failure
    TYPE(realspace_grid_desc_type), POINTER  :: desc

    CALL timeset(routineN,handle)

    failure = .FALSE.

    desc => rs_descs(grid_level)%rs_desc

    ! allocate local arrays
    ALLOCATE (send_buf_i(desc%group_size),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (recv_buf_i(desc%group_size),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    ! communication step 1 : compute the total local cost of the tasks
    !                        each proc needs to know the amount of work he will receive

    ! send buffer now contains for each target the cost of the tasks it will receive
    send_buf_i=0
    DO i = 1, ntasks
       CALL int2pair(tasks(3,i),ilevel,iatom,jatom,iset,jset,ipgf,jpgf,natom,maxset,maxpgf)
       IF (ilevel .NE. grid_level) CYCLE
       IF (use_reordered_ranks) THEN
         send_buf_i( rs_descs(ilevel)%rs_desc%virtual2real(decode_rank(tasks(1,i), SIZE(rs_descs))) +1 ) = &
         send_buf_i( rs_descs(ilevel)%rs_desc%virtual2real(decode_rank(tasks(1,i), SIZE(rs_descs))) +1 ) &
                   + tasks(5,i)
       ELSE
         send_buf_i( decode_rank(tasks(1,i), SIZE(rs_descs)) +1 ) = &
         send_buf_i( decode_rank(tasks(1,i), SIZE(rs_descs)) +1 ) &
                   + tasks(5,i)
       END IF
    END DO
    CALL mp_alltoall(send_buf_i, recv_buf_i, 1, desc % group)

    ! communication step 2 : compute the global cost of the tasks
    total_cost_local = SUM(recv_buf_i)

    ! after this step, the recv buffer contains the local cost for each CPU
    CALL mp_allgather(total_cost_local, loads, desc % group)

    CALL timestop(handle)

  END SUBROUTINE get_current_loads
! *****************************************************************************
!> \brief performs load balancing shifting tasks on the replicated grids
!>        this modifies the destination of some of the tasks on replicated
!>        grids, and in this way balances the load
!> \param rs_descs ...
!> \param ntasks ...
!> \param tasks ...
!> \param natoms ...
!> \param maxset ...
!> \param maxpgf ...
!> \param error ...
!> \par History
!>      none
!> \author MattW 21/11/2007
! *****************************************************************************
  SUBROUTINE load_balance_replicated ( rs_descs, ntasks, tasks, natoms, maxset, maxpgf, error)

    TYPE(realspace_grid_desc_p_type), &
      DIMENSION(:), POINTER                  :: rs_descs
    INTEGER                                  :: ntasks
    INTEGER(KIND=int_8), DIMENSION(:, :), &
      POINTER                                :: tasks
    INTEGER, INTENT(IN)                      :: natoms, maxset, maxpgf
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'load_balance_replicated', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, i, iatom, ilevel, ipgf, iset, j, jatom, jpgf, jset, &
      no_overloaded, no_underloaded, proc_receiving, stat
    INTEGER(KIND=int_8)                      :: average_cost, cost_task_rep, &
                                                count, offset, &
                                                total_cost_global
    INTEGER(KIND=int_8), ALLOCATABLE, &
      DIMENSION(:)                           :: load_imbalance, loads, &
                                                recv_buf_i
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: index
    LOGICAL                                  :: failure
    TYPE(realspace_grid_desc_type), POINTER  :: desc

    CALL timeset(routineN,handle)

    failure = .FALSE.

    desc=>rs_descs(1)%rs_desc

    ! allocate local arrays
    ALLOCATE (recv_buf_i(desc%group_size),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (loads(desc%group_size),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    recv_buf_i = 0
    DO i = 1, SIZE(rs_descs)
      CALL get_current_loads ( loads, rs_descs, i, ntasks, natoms, maxset, maxpgf, &
                               tasks, use_reordered_ranks=.TRUE. ,error=error )
      recv_buf_i(:) = recv_buf_i + loads
    END DO

    total_cost_global = SUM(recv_buf_i)
    average_cost = total_cost_global / desc%group_size

    !
    ! compute how to redistribute the replicated tasks so that the average cost is reached
    !

    ! load imbalance measures the load of a given CPU relative 
    ! to the optimal load distribution (load=average)
    ALLOCATE (load_imbalance(desc%group_size),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (INDEX(desc%group_size),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    load_imbalance(:) = recv_buf_i - average_cost
    no_overloaded = 0
    no_underloaded = 0

    DO i = 1, desc%group_size
       IF (load_imbalance(i) .GT. 0 ) no_overloaded = no_overloaded + 1
       IF (load_imbalance(i) .LT. 0 ) no_underloaded = no_underloaded + 1
    ENDDO

    ! sort the recv_buffer on number of tasks, gives us index which provides a
    ! mapping between processor ranks and how overloaded the processor
    CALL sort(recv_buf_i,SIZE(recv_buf_i),index)

    ! find out the number of replicated tasks each proc has
    ! but only those tasks which have not yet been assigned
    cost_task_rep = 0
    DO i = 1, ntasks
       IF(tasks(4,i) .EQ. 0 &
          .AND. decode_rank(tasks(1,i),SIZE(rs_descs))==decode_rank(tasks(2,i), SIZE(rs_descs)) ) THEN
          cost_task_rep = cost_task_rep + tasks(5,i)
       ENDIF
    ENDDO

    ! now, correct the load imbalance for the overloaded CPUs
    ! they will send away not more than the total load of replicated tasks
    CALL mp_allgather(cost_task_rep,recv_buf_i,desc % group)

    DO i = 1, desc%group_size
       ! At the moment we can only offload replicated tasks
       IF(load_imbalance(i) .GT. 0)&
            load_imbalance(i) = MIN( load_imbalance(i),recv_buf_i(i))
    ENDDO

    ! simplest algorithm I can think of of is that the processor with the most
    ! excess tasks fills up the process needing most, then moves on to next most.
    ! At the moment if we've got less replicated tasks than we're overloaded then
    ! task balancing will be incomplete

    ! only need to do anything if I've excess tasks
    IF (load_imbalance( desc%my_pos + 1 ) .GT. 0 ) THEN

       count = 0 ! weighted amount of tasks offloaded
       offset = 0 ! no of underloaded processes already filled by other more overloaded procs

       ! calculate offset
       DO i = desc%group_size, desc%group_size-no_overloaded+1, -1
          IF ( INDEX(i) .EQ. desc%my_pos+1 ) THEN
             EXIT
          ELSE
             offset = offset + load_imbalance ( INDEX (i) )
          ENDIF
       ENDDO

       ! find my starting processor to send to
       proc_receiving = HUGE(proc_receiving)
       DO i = 1, no_underloaded
          offset = offset + load_imbalance ( INDEX (i) )
          IF ( offset .LE. 0 ) THEN
             proc_receiving = i
             EXIT
          ENDIF
       ENDDO

       ! offset now contains minus the number of tasks proc_receiving requires
       ! we fill this up by adjusting the destination of tasks on the replicated grid, 
       ! then move to next most underloaded proc
       DO j = 1, ntasks
          IF (tasks(4,j) .EQ. 0 &
              .AND. decode_rank(tasks(1,j),SIZE(rs_descs))==decode_rank(tasks(2,j),SIZE(rs_descs)) ) THEN
             ! just avoid sending to non existing procs due to integer truncation 
             ! in the computation of the average
             IF(proc_receiving .GT. no_underloaded) EXIT
             ! set new destination
             CALL int2pair(tasks(3,j),ilevel,iatom,jatom,iset,jset,ipgf,jpgf,natoms,maxset,maxpgf)
             tasks (1,j) = encode_rank( INDEX ( proc_receiving ) - 1, ilevel, SIZE(rs_descs))
             offset = offset + tasks(5,j)
             count = count + tasks(5,j)
             IF(count .GE. load_imbalance (desc%my_pos + 1 ) ) EXIT
             IF(offset .GT. 0 ) THEN
                proc_receiving = proc_receiving + 1
                ! just avoid sending to non existing procs due to integer truncation 
                ! in the computation of the average
                IF(proc_receiving .GT. no_underloaded) EXIT
                offset = load_imbalance ( INDEX ( proc_receiving ) )
             ENDIF
          ENDIF
       ENDDO
    ENDIF

    DEALLOCATE (index,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (load_imbalance,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE load_balance_replicated

! *****************************************************************************
!> \brief given an input task list, redistribute so that all tasks can be processed locally,
!>        i.e. dest equals rank
!> \param rs_descs ...
!> \param ntasks ...
!> \param tasks ...
!> \param rval ...
!> \param ntasks_recv ...
!> \param tasks_recv ...
!> \param rval_recv ...
!> \param error ...
!> \par History
!>      none
!> \author MattW 21/11/2007
! *****************************************************************************
  SUBROUTINE create_local_tasks ( rs_descs, ntasks, tasks, rval, &
                                  ntasks_recv, tasks_recv, rval_recv, error)

    TYPE(realspace_grid_desc_p_type), &
      DIMENSION(:), POINTER                  :: rs_descs
    INTEGER                                  :: ntasks
    INTEGER(KIND=int_8), DIMENSION(:, :), &
      POINTER                                :: tasks
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: rval
    INTEGER                                  :: ntasks_recv
    INTEGER(KIND=int_8), DIMENSION(:, :), &
      POINTER                                :: tasks_recv
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: rval_recv
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'create_local_tasks', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, j, k, l, rank, &
                                                stat, task_dim
    INTEGER(KIND=int_8), ALLOCATABLE, &
      DIMENSION(:)                           :: recv_buf_i, send_buf_i
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: recv_disps, recv_sizes, &
                                                send_disps, send_sizes
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: recv_buf_r, send_buf_r
    TYPE(realspace_grid_desc_type), POINTER  :: desc

    CALL timeset(routineN,handle)

    failure = .FALSE.

    desc=>rs_descs(1)%rs_desc

    ! allocate local arrays
    ALLOCATE (send_sizes(desc%group_size),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (recv_sizes(desc%group_size),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (send_disps(desc%group_size),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (recv_disps(desc%group_size),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (send_buf_i(desc%group_size),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (recv_buf_i(desc%group_size),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    ! fill send buffer, now counting how many tasks will be send (stored in an int8 array for convenience only).
    send_buf_i=0
    DO i = 1, ntasks
       rank = rs_descs(decode_level(tasks(1,i), SIZE(rs_descs)))%rs_desc%virtual2real(decode_rank(tasks(1,i), SIZE(rs_descs)))
       send_buf_i( rank+1 ) = send_buf_i( rank+1 ) + 1
    END DO

    CALL mp_alltoall(send_buf_i, recv_buf_i, 1, desc % group)

    ! pack the tasks, and send them around

    task_dim = SIZE(tasks,1)

    send_sizes = 0
    send_disps = 0
    recv_sizes = 0
    recv_disps = 0

    send_sizes(1) = INT(send_buf_i(1) * task_dim)
    recv_sizes(1) = INT(recv_buf_i(1) * task_dim)
    DO i = 2,desc%group_size
       send_sizes(i) = INT(send_buf_i(i) * task_dim)
       recv_sizes(i) = INT(recv_buf_i(i) * task_dim)
       send_disps(i)=send_disps(i-1)+send_sizes(i-1)
       recv_disps(i)=recv_disps(i-1)+recv_sizes(i-1)
    ENDDO

    ! deallocate old send/recv buffers
    DEALLOCATE (send_buf_i,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (recv_buf_i,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    ! allocate them with new sizes
    ALLOCATE (send_buf_i(SUM(send_sizes)),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (recv_buf_i(SUM(recv_sizes)),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    ! do packing
    send_buf_i = 0
    send_sizes = 0
    DO j = 1,ntasks
       i = rs_descs(decode_level(tasks(1,j), SIZE(rs_descs)))%rs_desc%virtual2real(decode_rank(tasks(1,j), SIZE(rs_descs))) +1
       DO k=1,task_dim
          send_buf_i(send_disps(i)+send_sizes(i)+k)=tasks(k,j)
       ENDDO
       send_sizes(i)=send_sizes(i) + task_dim
    ENDDO

    ! do communication
    CALL mp_alltoall(send_buf_i, send_sizes, send_disps, recv_buf_i, recv_sizes, recv_disps, desc % group)

    DEALLOCATE (send_buf_i,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    ntasks_recv=SUM(recv_sizes)/task_dim
    ALLOCATE (tasks_recv(task_dim,ntasks_recv),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    ! do unpacking
    l = 0
    DO i = 1, desc % group_size
       DO j = 0,recv_sizes(i)/task_dim-1
          l = l + 1
          DO k=1,task_dim
             tasks_recv(k,l)=recv_buf_i(recv_disps(i)+j*task_dim+k)
          ENDDO
       ENDDO
    ENDDO

    DEALLOCATE (recv_buf_i,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    ! send rvals (to be removed :-)

    ! take care of the new task_dim in the send_sizes
    send_sizes(:) = (send_sizes(:)/task_dim)*SIZE(rval,1)
    recv_sizes(:) = (recv_sizes(:)/task_dim)*SIZE(rval,1)
    send_disps(:) = (send_disps(:)/task_dim)*SIZE(rval,1)
    recv_disps(:) = (recv_disps(:)/task_dim)*SIZE(rval,1)
    task_dim=SIZE(rval,1)

    ALLOCATE (send_buf_r(SUM(send_sizes)),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (recv_buf_r(SUM(recv_sizes)),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    !do packing
    send_sizes = 0
    DO j = 1,ntasks
       i = rs_descs(decode_level(tasks(1,j), SIZE(rs_descs)))%rs_desc%virtual2real(decode_rank(tasks(1,j), SIZE(rs_descs))) +1
       DO k=1,task_dim
          send_buf_r(send_disps(i)+send_sizes(i)+k)=rval(k,j)
       ENDDO
       send_sizes(i)=send_sizes(i) + task_dim
    ENDDO

    ! do communication
    CALL mp_alltoall(send_buf_r, send_sizes, send_disps,&
         recv_buf_r, recv_sizes, recv_disps, desc % group)

    DEALLOCATE (send_buf_r,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (rval_recv(task_dim,SUM(recv_sizes)/task_dim),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    ! do unpacking
    l = 0
    DO i = 1,desc % group_size
       DO j = 0,recv_sizes(i)/task_dim-1
          l = l + 1
          DO k=1,task_dim
             rval_recv(k,l)=recv_buf_r(recv_disps(i)+j*task_dim+k)
          ENDDO
       ENDDO
    ENDDO

    DEALLOCATE (recv_buf_r,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (send_sizes,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (recv_sizes,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (send_disps,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (recv_disps,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE create_local_tasks

! *****************************************************************************
!> \brief Assembles tasks to be performed on local grid
!> \param rs_descs the grids
!> \param ntasks Number of tasks for local processing
!> \param natoms ...
!> \param maxset ...
!> \param maxpgf ...
!> \param nimages ...
!> \param tasks the task set generated on this processor
!> \param rval ...
!> \param atom_pair_send ...
!> \param atom_pair_recv ...
!> \param symmetric ...
!> \param reorder_rs_grid_ranks ...
!> \param skip_load_balance_distributed ...
!> \param error ...
!> \par History
!>      none
!> \author MattW 21/11/2007
! *****************************************************************************
  SUBROUTINE distribute_tasks ( rs_descs, ntasks, natoms, maxset, maxpgf, nimages,&
             tasks, rval, atom_pair_send, atom_pair_recv,&
             symmetric, reorder_rs_grid_ranks, skip_load_balance_distributed, error)

    TYPE(realspace_grid_desc_p_type), &
      DIMENSION(:), POINTER                  :: rs_descs
    INTEGER                                  :: ntasks, natoms, maxset, &
                                                maxpgf, nimages
    INTEGER(KIND=int_8), DIMENSION(:, :), &
      POINTER                                :: tasks
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: rval
    INTEGER(KIND=int_8), DIMENSION(:), &
      POINTER                                :: atom_pair_send, atom_pair_recv
    LOGICAL, INTENT(IN)                      :: symmetric, &
                                                reorder_rs_grid_ranks, &
                                                skip_load_balance_distributed
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'distribute_tasks', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, igrid_level, irank, &
                                                k, ntasks_recv, stat
    INTEGER(KIND=int_8)                      :: load_gap, max_load, &
                                                replicated_load
    INTEGER(KIND=int_8), ALLOCATABLE, &
      DIMENSION(:)                           :: taskid, total_loads, &
                                                total_loads_tmp, trial_loads
    INTEGER(KIND=int_8), DIMENSION(:, :), &
      POINTER                                :: loads, tasks_recv
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: index, real2virtual, &
                                                total_index
    LOGICAL                                  :: distributed_grids, failure, &
                                                fixed_first_grid
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: rval_recv
    TYPE(realspace_grid_desc_type), POINTER  :: desc

    CALL timeset(routineN,handle)

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(tasks),cp_failure_level,routineP,error,failure)

    ! *** figure out if we have distributed grids
    distributed_grids=.FALSE.
    DO igrid_level=1,SIZE(rs_descs)
       IF ( rs_descs(igrid_level)%rs_desc%distributed ) THEN
          distributed_grids=.TRUE.
       ENDIF
    END DO
    desc=>rs_descs(1)%rs_desc

    IF (distributed_grids) THEN

       ALLOCATE (loads(0:desc%group_size-1,SIZE(rs_descs)),STAT=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE (total_loads(0:desc%group_size-1),STAT=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

       total_loads = 0

       ! First round of balancing on the distributed grids
       ! we just balance each of the distributed grids independently
       ! XXXX this call to load_balance_distributed can be a memory hog XXXX
       DO igrid_level=1,SIZE(rs_descs)
         IF ( rs_descs(igrid_level)%rs_desc%distributed ) THEN

           IF (.NOT.skip_load_balance_distributed) &
               CALL load_balance_distributed(tasks, ntasks, rs_descs, igrid_level, natoms, maxset,maxpgf, error)

           CALL get_current_loads( loads(:,igrid_level), rs_descs, igrid_level, ntasks, natoms, maxset,maxpgf, &
                                   tasks, use_reordered_ranks=.FALSE. , error=error)

           total_loads(:)=total_loads+loads(:,igrid_level)

         END IF
       END DO

       ! calculate the total load of replicated tasks, so we can decide if it is worth
       ! reordering the distributed grid levels

       replicated_load=0
       DO igrid_level=1,SIZE(rs_descs)
        IF (.NOT. rs_descs(igrid_level)%rs_desc%distributed) THEN
          CALL get_current_loads( loads(:,igrid_level), rs_descs, igrid_level, ntasks, natoms, maxset,maxpgf, &
                                     tasks, use_reordered_ranks=.FALSE., error=error)
          replicated_load = replicated_load + SUM(loads(:,igrid_level))
        END IF
       END DO

       !IF (desc%my_pos==0) THEN
       ! WRITE(*,*) "Total replicated load is ",replicated_load
       !END IF

       ! Now we adjust the rank ordering based on the current loads
       ! we leave the first distributed level and all the replicated levels in the default order
       IF (reorder_rs_grid_ranks) THEN
         fixed_first_grid=.FALSE.
         DO igrid_level=1,SIZE(rs_descs)
           IF (rs_descs(igrid_level)%rs_desc%distributed ) THEN
             IF (fixed_first_grid .EQV. .FALSE.) THEN
               total_loads(:) = loads(:,igrid_level)
               fixed_first_grid = .TRUE.
             ELSE
               ALLOCATE (trial_loads(0:desc%group_size-1),STAT=stat)
               CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

               trial_loads(:) = total_loads + loads(:,igrid_level)
               max_load = MAXVAL(trial_loads)
               load_gap = 0
               DO irank=0, desc%group_size-1
                 load_gap = load_gap + max_load - trial_loads(irank)
               END DO

               ! If there is not enough replicated load to load balance well enough
               ! then we will reorder this grid level
               IF (load_gap > replicated_load * 1.05_dp) THEN

                 ALLOCATE (INDEX(0:desc%group_size-1),STAT=stat)
                 CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
                 ALLOCATE (total_index(0:desc%group_size-1),STAT=stat)
                 CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
                 ALLOCATE (total_loads_tmp(0:desc%group_size-1),STAT=stat)
                 CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
                 ALLOCATE (real2virtual(0:desc%group_size-1),STAT=stat)
                 CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

                 total_loads_tmp(:) = total_loads
                 CALL sort(total_loads_tmp, desc%group_size, total_index)
                 CALL sort(loads(:,igrid_level), desc%group_size, index)

                 ! Reorder so that the rank with smallest load on this grid level is paired with
                 ! the highest load in total
                 DO irank=0, desc%group_size-1
                   total_loads(total_index(irank)-1) = total_loads(total_index(irank)-1) + &
                               loads(desc%group_size - irank - 1, igrid_level)
                   real2virtual(total_index(irank)-1)=INDEX(desc%group_size - irank - 1) - 1
                 END DO

                 CALL rs_grid_reorder_ranks(rs_descs(igrid_level)%rs_desc, real2virtual, error)

                 DEALLOCATE (index,STAT=stat)
                 CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
                 DEALLOCATE (total_index,STAT=stat)
                 CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
                 DEALLOCATE (total_loads_tmp,STAT=stat)
                 CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
                 DEALLOCATE (real2virtual,STAT=stat)
                 CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
               ELSE
                 total_loads(:) = trial_loads
               ENDIF

               DEALLOCATE (trial_loads,STAT=stat)
               CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

             ENDIF
           ENDIF
         END DO
       END IF

       ! Now we use the replicated tasks to balance out the rest of the load
       CALL load_balance_replicated(rs_descs, ntasks, tasks, natoms, maxset, maxpgf, error)

       !total_loads = 0
       !DO igrid_level=1,SIZE(rs_descs)
       !  CALL get_current_loads(loads(:,igrid_level), rs_descs, igrid_level, ntasks, natoms, maxset,maxpgf, &
       !                         tasks, use_reordered_ranks=.TRUE. , error=error)
       !  total_loads = total_loads + loads(:, igrid_level)
       !END DO

       !IF (desc%my_pos==0) THEN
       !  WRITE(*,*) ""
       !  WRITE(*,*) "At the end of the load balancing procedure"
       !  WRITE(*,*) "Maximum  load:",MAXVAL(total_loads)
       !  WRITE(*,*) "Average  load:",SUM(total_loads)/SIZE(total_loads)
       !  WRITE(*,*) "Minimum  load:",MINVAL(total_loads)
       !ENDIF

       ! given a list of tasks, this will do the needed reshuffle so that all tasks will be local
       CALL create_local_tasks ( rs_descs, ntasks, tasks, rval, ntasks_recv, tasks_recv, rval_recv, error)

       !
       ! tasks list are complete, we can compute the list of atomic blocks (atom pairs)
       ! we will be sending. These lists are needed for redistribute_matrix.
       !
       ALLOCATE (atom_pair_send(ntasks),STAT=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       CALL get_atom_pair ( atom_pair_send, tasks, send=.TRUE.,&
            symmetric=symmetric, natoms=natoms, maxset=maxset, maxpgf=maxpgf, &
            nimages=nimages, rs_descs=rs_descs, error=error )

       ! natom_send=SIZE(atom_pair_send)
       ! CALL mp_sum(natom_send,desc%group)
       ! IF (desc%my_pos==0) THEN
       !     WRITE(*,*) ""
       !     WRITE(*,*) "Total number of atomic blocks to be send:",natom_send
       ! ENDIF

       ALLOCATE (atom_pair_recv(ntasks_recv),STAT=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       CALL get_atom_pair ( atom_pair_recv, tasks_recv, send=.FALSE.,&
            symmetric=symmetric, natoms=natoms, maxset=maxset, maxpgf=maxpgf, &
            nimages=nimages, rs_descs=rs_descs, error=error )

       ! cleanup, at this point we  don't need the original tasks & rvals anymore
       DEALLOCATE (tasks,STAT=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE (rval,STAT=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE (loads,STAT=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE (total_loads,STAT=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    ELSE

       tasks_recv =>tasks
       ntasks_recv=ntasks
       rval_recv=>rval

    ENDIF

    !
    ! here we sort the task list we will process locally.
    ! the sort is determined by pair2int
    !
    rval => rval_recv

    ALLOCATE (taskid(ntasks_recv),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (INDEX(ntasks_recv),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    taskid(:)=tasks_recv(3,1:ntasks_recv)
    CALL sort(taskid,SIZE(taskid),index)

    DO k=1,SIZE(tasks_recv,1)
       tasks_recv(k,1:ntasks_recv)=tasks_recv(k,index)
    ENDDO

    DO k=1,SIZE(rval,1)
       rval(k,1:ntasks_recv)=rval(k,index)
    ENDDO

    DEALLOCATE (taskid,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (index,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    !
    ! final lists are ready
    !

    tasks=>tasks_recv
    ntasks=ntasks_recv

    CALL timestop(handle)

  END SUBROUTINE distribute_tasks

! *****************************************************************************
!> \brief ...
!> \param atom_pair ...
!> \param my_tasks ...
!> \param send ...
!> \param symmetric ...
!> \param natoms ...
!> \param maxset ...
!> \param maxpgf ...
!> \param nimages ...
!> \param rs_descs ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE get_atom_pair ( atom_pair, my_tasks, send, symmetric, natoms, maxset, &
                             maxpgf, nimages, rs_descs, error)

    INTEGER(KIND=int_8), DIMENSION(:), &
      POINTER                                :: atom_pair
    INTEGER(KIND=int_8), DIMENSION(:, :), &
      POINTER                                :: my_tasks
    LOGICAL                                  :: send, symmetric
    INTEGER                                  :: natoms, maxset, maxpgf, &
                                                nimages
    TYPE(realspace_grid_desc_p_type), &
      DIMENSION(:), POINTER                  :: rs_descs
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_atom_pair', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: acol, arow, i, iatom, ilevel, &
                                                img, ipgf, iset, j, jatom, &
                                                jpgf, jset, stat
    INTEGER(KIND=int_8)                      :: natom8, nim8
    INTEGER, DIMENSION(:), POINTER           :: index
    LOGICAL                                  :: failure

! calculate list of atom pairs
! fill pair list taking into acount symmetry

  failure = .FALSE.

  natom8 = natoms
  nim8 = nimages

  atom_pair = 0
  DO i = 1,SIZE(atom_pair)
     CALL int2pair(my_tasks(3,i),ilevel,iatom,jatom,iset,jset,ipgf,jpgf,natoms,maxset,maxpgf)
     img = INT(my_tasks(7,i))
     IF( symmetric ) THEN
        IF(iatom.LE.jatom) THEN
           arow = iatom
           acol = jatom
        ELSE
           arow = jatom
           acol = iatom
        ENDIF
     ELSE
        arow = iatom
        acol = jatom
     ENDIF
     IF ( send ) THEN
        ! If sending, we need to use the 'real rank' as the pair has to be sent to the process which
        ! actually has the correct part of the rs_grid to do the mapping
        atom_pair(i) = rs_descs(ilevel)%rs_desc%virtual2real(decode_rank(my_tasks(1,i),SIZE(rs_descs))) &
                        *natom8*natom8*nim8+(arow-1)*natom8*nim8 + (acol-1)*nim8 + (img-1)
     ELSE
        ! If we are recieving, then no conversion is needed as the rank is that of the process with the
        ! required matrix block, and the ordering of the rs grid is irrelevant
        atom_pair(i) = decode_rank(my_tasks(2,i),SIZE(rs_descs)) &
                        *natom8*natom8*nim8+(arow-1)*natom8*nim8 + (acol-1)*nim8 + (img-1)
     ENDIF
  ENDDO

  ALLOCATE (INDEX(SIZE(atom_pair)),STAT=stat)
  CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

  ! find unique atom pairs that I'm sending/receiving
  CALL sort(atom_pair,SIZE(atom_pair),index)

  DEALLOCATE (index,STAT=stat)
  CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

  IF (SIZE(atom_pair)>1) THEN
     j=1
     ! first atom pair must be allowed
     DO i = 2,SIZE(atom_pair)
        IF( atom_pair(i) .GT. atom_pair(i-1) ) THEN
           j = j + 1
           atom_pair(j) = atom_pair(i)
        ENDIF
     ENDDO
     ! reallocate the atom pair list
     CALL  reallocate(atom_pair,1,j)
  ENDIF

END SUBROUTINE get_atom_pair

! *****************************************************************************
!> \brief redistributes the matrix so that it can be used in realspace operations
!>        i.e. according to the task lists for collocate and integrate.
!>        This routine can become a bottleneck in large calculations.
!> \param rs_descs ...
!> \param pmats ...
!> \param atom_pair_send ...
!> \param atom_pair_recv ...
!> \param natoms ...
!> \param nimages ...
!> \param scatter ...
!> \param hmats ...
!> \param error ...
! *****************************************************************************
SUBROUTINE rs_distribute_matrix( rs_descs, pmats, atom_pair_send, atom_pair_recv, &
     natoms, nimages, scatter, hmats, error )

    TYPE(realspace_grid_desc_p_type), &
      DIMENSION(:), POINTER                  :: rs_descs
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: pmats
    INTEGER(KIND=int_8), DIMENSION(:), &
      POINTER                                :: atom_pair_send, atom_pair_recv
    INTEGER                                  :: natoms, nimages
    LOGICAL                                  :: scatter
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: hmats
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'rs_distribute_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER :: acol, arow, handle, i, img, j, k, l, me, nblkcols_total, &
      nblkrows_total, ncol, nrow, nthread, nthread_left, stat
    INTEGER(KIND=int_8)                      :: natom8, nim8, pair
    INTEGER, ALLOCATABLE, DIMENSION(:) :: first_col, first_row, last_col, &
      last_row, recv_disps, recv_pair_count, recv_pair_disps, recv_sizes, &
      send_disps, send_pair_count, send_pair_disps, send_sizes
    LOGICAL                                  :: failure, found
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: recv_buf_r, send_buf_r
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: h_block, p_block
    TYPE(cp_dbcsr_type), POINTER             :: hmat, pmat
    TYPE(realspace_grid_desc_type), POINTER  :: desc

!$  INTEGER(kind=omp_lock_kind), ALLOCATABLE, DIMENSION(:) :: locks

  CALL timeset(routineN,handle)

  failure = .FALSE.
  IF (.not.scatter) THEN
    CPPrecondition(PRESENT(hmats),cp_failure_level,routineP,error,failure)
  END IF

  desc => rs_descs(1)%rs_desc
  me = desc%my_pos+1

  ! allocate local arrays
  ALLOCATE (send_sizes(desc%group_size),STAT=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  ALLOCATE (recv_sizes(desc%group_size),STAT=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  ALLOCATE (send_disps(desc%group_size),STAT=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  ALLOCATE (recv_disps(desc%group_size),STAT=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  ALLOCATE (send_pair_count(desc%group_size),STAT=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  ALLOCATE (recv_pair_count(desc%group_size),STAT=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  ALLOCATE (send_pair_disps(desc%group_size),STAT=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  ALLOCATE (recv_pair_disps(desc%group_size),STAT=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

  pmat => pmats(1)%matrix
  CALL cp_dbcsr_get_info (pmat, &
       nblkrows_total=nblkrows_total,&
       nblkcols_total=nblkcols_total)
  ALLOCATE (first_row(nblkrows_total), last_row(nblkrows_total),&
            first_col(nblkcols_total), last_col(nblkcols_total),STAT=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  CALL convert_sizes_to_offsets (cp_dbcsr_row_block_sizes(pmat), first_row, last_row)
  CALL convert_sizes_to_offsets (cp_dbcsr_col_block_sizes(pmat), first_col, last_col)
  ! set up send buffer sizes
  natom8 = natoms
  nim8 = nimages

  send_sizes = 0
  send_pair_count = 0
  DO i = 1, SIZE(atom_pair_send)

     ! proc we're sending this block to
     k = INT(atom_pair_send(i) / (nim8*natom8**2)) + 1

     pair = MOD(atom_pair_send(i),nim8*natom8**2)

     arow = INT(pair / (nim8*natom8)) + 1
     acol = INT(MOD(pair, nim8*natom8) / nim8) + 1

     nrow = last_row(arow) - first_row(arow) + 1
     ncol = last_col(acol) - first_col(acol) + 1

     send_sizes(k) = send_sizes(k) + nrow * ncol
     send_pair_count(k) = send_pair_count(k) + 1

  ENDDO

  send_disps = 0
  send_pair_disps = 0
  DO i = 2, desc % group_size
     send_disps(i) = send_disps(i-1) + send_sizes(i-1)
     send_pair_disps(i) = send_pair_disps(i-1) + send_pair_count(i-1)
  ENDDO

  ALLOCATE (send_buf_r(SUM(send_sizes)),STAT=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

  ! set up recv buffer

  recv_sizes = 0
  recv_pair_count = 0
  DO i = 1, SIZE(atom_pair_recv)

     ! proc we're receiving this data from
     k = INT(atom_pair_recv(i) / (nim8*natom8**2)) + 1

     pair = MOD(atom_pair_recv(i),nim8*natom8**2)

     arow = INT(pair / (nim8*natom8)) + 1
     acol = INT(MOD(pair,nim8*natom8) / nim8) + 1

     nrow = last_row(arow) - first_row(arow) + 1
     ncol = last_col(acol) - first_col(acol) + 1

     recv_sizes(k) = recv_sizes(k) + nrow * ncol
     recv_pair_count(k) = recv_pair_count(k) + 1

  ENDDO

  recv_disps = 0
  recv_pair_disps = 0
  DO i = 2, desc % group_size
     recv_disps(i) = recv_disps(i-1) + recv_sizes(i-1)
     recv_pair_disps(i) = recv_pair_disps(i-1) + recv_pair_count(i-1)
  ENDDO
  ALLOCATE (recv_buf_r(SUM(recv_sizes)),STAT=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

!$omp parallel default(none), &
!$omp          shared(desc,send_pair_count,send_pair_disps,natom8,nim8,nimages),&
!$omp          shared(last_row,first_row,last_col,first_col),&
!$omp          shared(pmats,send_buf_r,send_disps,send_sizes),&
!$omp          shared(atom_pair_send,me,hmats,nblkrows_total),&
!$omp          shared(atom_pair_recv,recv_buf_r,scatter,recv_pair_disps), &
!$omp          shared(recv_sizes,recv_disps,recv_pair_count,locks), &
!$omp          private(i,img,pair,arow,acol,nrow,ncol,p_block,found,j,k,l),&
!$omp          private(nthread,h_block,error,nthread_left,hmat,pmat,failure),&
!$omp          private(stat)

  nthread = 1
!$ nthread = omp_get_num_threads()
  nthread_left = 1
!$ nthread_left = MAX(1, nthread-1)

  ! do packing
!$omp do schedule(guided)
  DO l =1, desc%group_size
   IF (l .EQ. me) CYCLE
    send_sizes(l)=0
    DO i = 1, send_pair_count(l)
     pair = MOD(atom_pair_send(send_pair_disps(l)+i),nim8*natom8**2)

     arow = INT(pair / (nim8*natom8)) + 1
     acol = INT(MOD(pair, nim8*natom8) / nim8) + 1
     img  = INT(MOD(pair, nim8)) + 1

     nrow = last_row(arow)-first_row(arow) + 1
     ncol = last_col(acol)-first_col(acol) + 1

     pmat => pmats(img)%matrix
     CALL cp_dbcsr_get_block_p(matrix=pmat,row=arow,col=acol,&
          block=p_block,found=found)
     CPPostcondition(found,cp_failure_level,routineP,error,failure)

     DO k = 1, ncol
        DO j = 1, nrow
           send_buf_r(send_disps(l)+send_sizes(l)+j+(k-1)*nrow) = p_block(j,k)
        ENDDO
     ENDDO
     send_sizes(l)=send_sizes(l)+nrow*ncol
    ENDDO
  ENDDO
!$omp end do

  IF (.NOT. scatter) THEN
  ! We need locks to protect concurrent summation into H
!$omp single
!$  ALLOCATE (locks(nthread*10))
!$omp end single

!$omp do
!$  do i=1,nthread*10
!$    call omp_init_lock(locks(i))
!$  end do
!$omp end do 
  END IF

!$omp master
  ! do communication
  CALL mp_alltoall(send_buf_r, send_sizes, send_disps,&
       recv_buf_r, recv_sizes, recv_disps, desc % group)
!$omp end master

  ! If this is a scatter, then no need to copy local blocks,
  ! If not, we sum them directly into H (bypassing the alltoall)
  IF (.NOT. scatter) THEN

  ! Distribute work over remaining threads assuming one is still in the alltoall
!$omp do schedule(dynamic,MAX(1,send_pair_count(me)/nthread_left))
    DO i = 1, send_pair_count(me)
     pair = MOD(atom_pair_send(send_pair_disps(me)+i),natom8**2)

     arow = INT(pair / (nim8*natom8)) + 1
     acol = INT(MOD(pair, nim8*natom8) / nim8) + 1
     img  = INT(MOD(pair, nim8)) + 1

     nrow = last_row(arow) - first_row(arow) + 1
     ncol = last_col(acol) - first_col(acol) + 1

     hmat => hmats(img)%matrix
     pmat => pmats(img)%matrix
     CALL cp_dbcsr_get_block_p(matrix=hmat,row=arow,col=acol,&
          BLOCK=h_block,found=found)
     CPPostcondition(found,cp_failure_level,routineP,error,failure)
     CALL cp_dbcsr_get_block_p(matrix=pmat,row=arow,col=acol,&
          BLOCK=p_block,found=found)
     CPPostcondition(found,cp_failure_level,routineP,error,failure)

!$  call omp_set_lock(locks((arow-1)*nthread*10/nblkrows_total+1))
        DO k = 1, ncol
           DO j = 1, nrow
              h_block(j,k) = h_block(j,k) + p_block(j,k)
           ENDDO
        ENDDO
!$  call omp_unset_lock(locks((arow-1)*nthread*10/nblkrows_total+1))
    ENDDO
!$omp end do
  ELSE
    ! We will insert new blocks into P, so create mutable work matrices
    DO img=1,nimages
       pmat => pmats(img)%matrix
       CALL cp_dbcsr_work_create(pmat,work_mutable=.TRUE.,&
         nblks_guess=SIZE(atom_pair_recv)/nthread,sizedata_guess=SIZE(recv_buf_r)/nthread,&
         n=nthread,error=error)
    END DO
  END IF

! wait for comm and setup to finish
!$omp barrier

  !do unpacking
!$omp do schedule(guided)
  DO l = 1, desc%group_size
    IF (l .EQ. me) CYCLE
    recv_sizes(l) = 0
    DO i = 1, recv_pair_count(l)
     pair = MOD(atom_pair_recv(recv_pair_disps(l)+i),natom8**2)

     arow = INT(pair / (nim8*natom8)) + 1
     acol = INT(MOD(pair, nim8*natom8) / nim8) + 1
     img  = INT(MOD(pair, nim8)) + 1

     nrow = last_row(arow)-first_row(arow) + 1
     ncol = last_col(acol)-first_col(acol) + 1

     pmat => pmats(img)%matrix
     NULLIFY(p_block)
     CALL cp_dbcsr_get_block_p(matrix=pmat,row=arow,col=acol,&
          BLOCK=p_block,found=found)

     IF ( PRESENT ( hmats ) ) THEN
        hmat => hmats(img)%matrix
        CALL cp_dbcsr_get_block_p(matrix=hmat,row=arow,col=acol,&
             BLOCK=h_block,found=found)
        CPPostcondition(found,cp_failure_level,routineP,error,failure)
     ENDIF

     IF(scatter .AND. .NOT. ASSOCIATED(p_block)) THEN
        CALL cp_dbcsr_put_block(pmat, arow, acol,&
           block=recv_buf_r(recv_disps(l)+recv_sizes(l)+1:recv_disps(l)+recv_sizes(l)+nrow*ncol))
     END IF
     IF (.not.scatter) THEN
!$  call omp_set_lock(locks((arow-1)*nthread*10/nblkrows_total+1))
        DO k = 1, ncol
           DO j = 1, nrow
              h_block(j,k) = h_block(j,k) + recv_buf_r( recv_disps(l) + recv_sizes(l)+ j + (k-1)*nrow )
           ENDDO
        ENDDO
!$  call omp_unset_lock(locks((arow-1)*nthread*10/nblkrows_total+1))
      ENDIF
     recv_sizes(l)=recv_sizes(l)+nrow*ncol
    ENDDO
  ENDDO
!$omp end do

!$ IF (.not.scatter) THEN
!$omp do
!$  do i=1,nthread*10
!$    call omp_destroy_lock(locks(i))
!$  end do
!$omp end do
!$ END IF

!$omp single
!$ IF (.not.scatter) THEN
!$  DEALLOCATE (locks,STAT=stat)
!$  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
!$ END IF
!$omp end single nowait

  IF (scatter) THEN
    ! Blocks were added to P
    DO img=1,nimages
       pmat => pmats(img)%matrix
       CALL cp_dbcsr_finalize(pmat, error=error)
    END DO
  END IF
!$omp end parallel

  DEALLOCATE (send_buf_r,STAT=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  DEALLOCATE (recv_buf_r,STAT=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

  DEALLOCATE (send_sizes,STAT=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  DEALLOCATE (recv_sizes,STAT=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  DEALLOCATE (send_disps,STAT=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  DEALLOCATE (recv_disps,STAT=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  DEALLOCATE (send_pair_count,STAT=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  DEALLOCATE (recv_pair_count,STAT=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  DEALLOCATE (send_pair_disps,STAT=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  DEALLOCATE (recv_pair_disps,STAT=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

  DEALLOCATE (first_row, last_row, first_col, last_col, STAT=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

  CALL timestop(handle)

END SUBROUTINE rs_distribute_matrix

! *****************************************************************************
!> \brief determines the rank of the processor who's real rs grid contains point
!> \param rs_desc ...
!> \param igrid_level ...
!> \param n_levels ...
!> \param cube_center ...
!> \param ntasks ...
!> \param tasks ...
!> \param lb_cube ...
!> \param ub_cube ...
!> \param added_tasks ...
!> \param error ...
!> \par History
!>      11.2007 created [MattW]
!>      10.2008 rewritten [Joost VandeVondele]
!> \author MattW
! *****************************************************************************
  SUBROUTINE rs_find_node(rs_desc,igrid_level,n_levels,cube_center,ntasks,tasks,&
                          lb_cube,ub_cube,added_tasks,error)

    TYPE(realspace_grid_desc_type), POINTER  :: rs_desc
    INTEGER, INTENT(IN)                      :: igrid_level, n_levels
    INTEGER, DIMENSION(3), INTENT(IN)        :: cube_center
    INTEGER, INTENT(INOUT)                   :: ntasks
    INTEGER(KIND=int_8), DIMENSION(:, :), &
      POINTER                                :: tasks
    INTEGER, DIMENSION(3), INTENT(IN)        :: lb_cube, ub_cube
    INTEGER, INTENT(OUT)                     :: added_tasks
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'rs_find_node', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: add_tasks = 1000
    REAL(kind=dp), PARAMETER                 :: mult_tasks = 2.0_dp

    INTEGER :: bit_index, coord(3), curr_tasks, dest, i, icoord(3), idest, &
      itask, ix, iy, iz, lb_coord(3), lb_domain(3), lbc(3), ub_coord(3), &
      ub_domain(3), ubc(3)
    INTEGER(KIND=int_8)                      :: bit_pattern
    LOGICAL                                  :: dir_periodic(3), failure

    coord(1)=rs_desc % x2coord(cube_center(1))
    coord(2)=rs_desc % y2coord(cube_center(2))
    coord(3)=rs_desc % z2coord(cube_center(3))
    dest=rs_desc%coord2rank(coord(1),coord(2),coord(3))

    failure=.FALSE.

    ! the real cube coordinates
    lbc=lb_cube+cube_center
    ubc=ub_cube+cube_center

    IF ( ALL ( (rs_desc % lb_global(:,dest)-rs_desc%border) .LE. lbc )  .AND. &
         ALL ( (rs_desc % ub_global(:,dest)+rs_desc%border) .GE. ubc )) THEN
       !standard distributed collocation/integration
       tasks(1,ntasks)=encode_rank(dest, igrid_level, n_levels)
       tasks(4,ntasks)=1
       tasks(6,ntasks)=0
       added_tasks = 1

       ! here we figure out if there is an alternate location for this task
       ! i.e. even though the cube_center is not in the real local domain,
       ! it might fully fit in the halo of the neighbor
       ! if its radius is smaller than the maximum radius
       ! the list of possible neighbors is stored here as a bit pattern
       ! to reconstruct what the rank of the neigbor is,
       ! we can use (note this requires the correct rs_grid)
       !  IF (BTEST(bit_pattern,0)) rank=rs_grid_locate_rank(rs_desc,tasks(1,ntasks),(/-1,0,0/))
       !  IF (BTEST(bit_pattern,1)) rank=rs_grid_locate_rank(rs_desc,tasks(1,ntasks),(/+1,0,0/))
       !  IF (BTEST(bit_pattern,2)) rank=rs_grid_locate_rank(rs_desc,tasks(1,ntasks),(/0,-1,0/))
       !  IF (BTEST(bit_pattern,3)) rank=rs_grid_locate_rank(rs_desc,tasks(1,ntasks),(/0,+1,0/))
       !  IF (BTEST(bit_pattern,4)) rank=rs_grid_locate_rank(rs_desc,tasks(1,ntasks),(/0,0,-1/))
       !  IF (BTEST(bit_pattern,5)) rank=rs_grid_locate_rank(rs_desc,tasks(1,ntasks),(/0,0,+1/))
       bit_index=0
       bit_pattern=0
       DO i=1,3
          IF (rs_desc % perd ( i ) ==  1) THEN
             bit_pattern=IBCLR(bit_pattern,bit_index)
             bit_index=bit_index+1
             bit_pattern=IBCLR(bit_pattern,bit_index)
             bit_index=bit_index+1
          ELSE
             ! fits the left neighbor ?
             IF (ubc(i)<=rs_desc%lb_global(i,dest)-1+rs_desc%border) THEN
                bit_pattern=IBSET(bit_pattern,bit_index)
                bit_index=bit_index+1
             ELSE
                bit_pattern=IBCLR(bit_pattern,bit_index)
                bit_index=bit_index+1
             ENDIF
             ! fits the right neighbor ?
             IF (lbc(i)>=rs_desc%ub_global(i,dest)+1-rs_desc%border) THEN
                bit_pattern=IBSET(bit_pattern,bit_index)
                bit_index=bit_index+1
             ELSE
                bit_pattern=IBCLR(bit_pattern,bit_index)
                bit_index=bit_index+1
             ENDIF
          ENDIF
       ENDDO
       tasks(6,ntasks)=bit_pattern

    ELSE
       ! generalised collocation/integration needed
       ! first we figure out how many neighbors we have to add to include the lbc/ubc
       ! in the available domains (inclusive of halo points)
       ! first we 'ignore' periodic boundary conditions
       ! i.e. ub_coord-lb_coord+1 might be larger than group_dim
       lb_coord=coord
       ub_coord=coord
       lb_domain=rs_desc % lb_global(:,dest) - rs_desc % border
       ub_domain=rs_desc % ub_global(:,dest) + rs_desc % border
       DO i=1,3
          ! only if the grid is not periodic in this direction we need to take care of adding neighbors
          IF (rs_desc%perd(i)==0) THEN
             ! if the domain lower bound is greater than the lbc we need to add the size of the neighbor domain
             DO
               IF (lb_domain(i)>lbc(i)) THEN
                  lb_coord(i)=lb_coord(i)-1
                  icoord=MODULO(lb_coord,rs_desc%group_dim)
                  idest=rs_desc%coord2rank(icoord(1),icoord(2),icoord(3))
                  lb_domain(i)=lb_domain(i)-(rs_desc % ub_global(i,idest)-rs_desc % lb_global(i,idest)+1)
               ELSE
                  EXIT
               ENDIF
             ENDDO
             ! same for the upper bound
             DO
               IF (ub_domain(i)<ubc(i)) THEN
                  ub_coord(i)=ub_coord(i)+1
                  icoord=MODULO(ub_coord,rs_desc%group_dim)
                  idest=rs_desc%coord2rank(icoord(1),icoord(2),icoord(3))
                  ub_domain(i)=ub_domain(i)+(rs_desc % ub_global(i,idest)-rs_desc % lb_global(i,idest)+1)
               ELSE
                  EXIT
               ENDIF
             ENDDO
          ENDIF
       ENDDO

       ! some care is needed for the periodic boundaries
       DO i=1,3
          IF (ub_domain(i)-lb_domain(i)+1>=rs_desc%npts(i)) THEN
             dir_periodic(i)=.TRUE.
             lb_coord(i)=0
             ub_coord(i)=rs_desc%group_dim(i)-1
          ELSE
             dir_periodic(i)=.FALSE.
          ENDIF
       ENDDO

       added_tasks=PRODUCT(ub_coord-lb_coord+1)
       itask=ntasks
       ntasks=ntasks+added_tasks-1
       IF ( ntasks > SIZE(tasks,2)) THEN
          curr_tasks = INT((SIZE(tasks,2)+add_tasks)*mult_tasks)
          CALL reallocate(tasks,1,7,1,curr_tasks)
       END IF
       DO iz=lb_coord(3),ub_coord(3)
       DO iy=lb_coord(2),ub_coord(2)
       DO ix=lb_coord(1),ub_coord(1)
          icoord=MODULO((/ix,iy,iz/),rs_desc%group_dim)
          idest=rs_desc%coord2rank(icoord(1),icoord(2),icoord(3))
          tasks(1,itask)=encode_rank(idest, igrid_level, n_levels)
          tasks(4,itask)=2
          tasks(6,itask)=0
          ! encode the domain size for this task
          ! if the bit is set, we need to add the border in that direction
          IF (ix==lb_coord(1).AND..NOT.dir_periodic(1)) tasks(6,itask)=IBSET(tasks(6,itask),0)
          IF (ix==ub_coord(1).AND..NOT.dir_periodic(1)) tasks(6,itask)=IBSET(tasks(6,itask),1)
          IF (iy==lb_coord(2).AND..NOT.dir_periodic(2)) tasks(6,itask)=IBSET(tasks(6,itask),2)
          IF (iy==ub_coord(2).AND..NOT.dir_periodic(2)) tasks(6,itask)=IBSET(tasks(6,itask),3)
          IF (iz==lb_coord(3).AND..NOT.dir_periodic(3)) tasks(6,itask)=IBSET(tasks(6,itask),4)
          IF (iz==ub_coord(3).AND..NOT.dir_periodic(3)) tasks(6,itask)=IBSET(tasks(6,itask),5)
          itask=itask+1
       ENDDO
       ENDDO
       ENDDO
    ENDIF

  END SUBROUTINE rs_find_node

! *****************************************************************************
!> \brief utility functions for encoding the grid level with a rank, allowing
!>        different grid levels to maintain a different rank ordering without
!>        losing information.  These encoded_ints are stored in the tasks(1:2,:) array
!> \param rank ...
!> \param grid_level ...
!> \param n_levels ...
!> \retval encoded_int ...
!> \par History
!>      4.2009 created [Iain Bethune]
!>        (c) The Numerical Algorithms Group (NAG) Ltd, 2009 on behalf of the HECToR project
! *****************************************************************************
  FUNCTION encode_rank (rank, grid_level, n_levels) RESULT(encoded_int)

    INTEGER, INTENT(IN)                      :: rank, grid_level, n_levels
    INTEGER(KIND=int_8)                      :: encoded_int

! ordered so can still sort by rank

    encoded_int = rank*n_levels + grid_level - 1

  END FUNCTION

! *****************************************************************************
!> \brief ...
!> \param encoded_int ...
!> \param n_levels ...
!> \retval rank ...
! *****************************************************************************
  FUNCTION decode_rank (encoded_int, n_levels) RESULT(rank)

    INTEGER(KIND=int_8), INTENT(IN)          :: encoded_int
    INTEGER, INTENT(IN)                      :: n_levels
    INTEGER                                  :: rank

    rank = INT(encoded_int / n_levels)

  END FUNCTION

! *****************************************************************************
!> \brief ...
!> \param encoded_int ...
!> \param n_levels ...
!> \retval grid_level ...
! *****************************************************************************
  FUNCTION decode_level (encoded_int, n_levels) RESULT(grid_level)

    INTEGER(KIND=int_8), INTENT(IN)          :: encoded_int
    INTEGER, INTENT(IN)                      :: n_levels
    INTEGER                                  :: grid_level

    INTEGER(KIND=int_8)                      :: n_levels8

    n_levels8 = n_levels

    grid_level = INT(MODULO (encoded_int, n_levels8)) + 1

  END FUNCTION decode_level

END MODULE task_list_methods
