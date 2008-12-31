!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief generate the tasks lists used by collocate and integrate routines
!> \par History
!>      01.2008 [Joost VandeVondele] refactered out of qs_collocate / qs_integrate
!> \author Joost VandeVondele 
! *****************************************************************************
MODULE task_list_methods
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_rs_pool_types,                ONLY: cp_rs_pool_p_type,&
                                             rs_pools_create_rs_vect,&
                                             rs_pools_give_back_rs_vect
  USE cube_utils,                      ONLY: compute_cube_center,&
                                             cube_info_type,&
                                             return_cube,&
                                             return_cube_nonortho
  USE f77_blas
  USE gaussian_gridlevels,             ONLY: gaussian_gridlevel,&
                                             gridlevel_info_type
  USE input_constants,                 ONLY: use_aux_fit_basis_set,&
                                             use_orb_basis_set
  USE kinds,                           ONLY: dp,&
                                             int_8
  USE memory_utilities,                ONLY: reallocate
  USE message_passing,                 ONLY: mp_alltoall,&
                                             mp_gather,&
                                             mp_scatter
  USE particle_types,                  ONLY: particle_type
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_interactions,                 ONLY: exp_radius_very_extended
  USE qs_neighbor_list_types,          ONLY: &
       first_list, first_node, get_neighbor_list, get_neighbor_list_set, &
       get_neighbor_node, neighbor_list_set_p_type, neighbor_list_type, &
       neighbor_node_type, next
  USE realspace_grid_types,            ONLY: realspace_grid_p_type,&
                                             realspace_grid_type,&
                                             rs_grid_locate_rank
  USE sparse_matrix_types,             ONLY: add_block_node,&
                                             get_block_node,&
                                             real_matrix_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: sort
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'task_list_methods'

  PUBLIC :: generate_qs_task_list,&
            task_list_inner_loop
  PUBLIC :: distribute_tasks,&
            pair2int,&
            int2pair,&
            distribute_matrix

CONTAINS

! *****************************************************************************
!> \par History
!>      01.2008 factored out of calculate_rho_elec [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE generate_qs_task_list(qs_env, tasks, dist_ab, ntasks, atom_pair_send, atom_pair_recv,soft_valid, &
                                   basis_set_id, error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    INTEGER(kind=int_8), DIMENSION(:, :), &
      POINTER                                :: tasks
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: dist_ab
    INTEGER                                  :: ntasks
    INTEGER(kind=int_8), DIMENSION(:), &
      POINTER                                :: atom_pair_send, atom_pair_recv
    LOGICAL, INTENT(IN), OPTIONAL            :: soft_valid
    INTEGER, INTENT(IN), OPTIONAL            :: basis_set_id
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'generate_qs_task_list', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: add_tasks = 1000, &
                                                max_tasks = 2000
    REAL(kind=dp), PARAMETER                 :: mult_tasks = 2.0_dp

    INTEGER :: ab, curr_tasks, handle, iatom, ikind, ilist, inode, jatom, &
      jkind, maxpgf, maxset, my_basis_set_id, natoms, nkind, nlist, nnode, &
      nseta, nsetb
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, lb_max, &
                                                lb_min, npgfa, npgfb
    LOGICAL                                  :: failure, my_soft
    REAL(KIND=dp)                            :: kind_radius_a, kind_radius_b
    REAL(KIND=dp), DIMENSION(3)              :: ra, rab
    REAL(KIND=dp), DIMENSION(:), POINTER     :: set_radius_a, set_radius_b
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: rpgfa, rpgfb, zeta, zetb
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_rs_pool_p_type), DIMENSION(:), &
      POINTER                                :: rs_pools
    TYPE(cube_info_type), DIMENSION(:), &
      POINTER                                :: cube_info
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(gridlevel_info_type), POINTER       :: gridlevel_info
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb
    TYPE(neighbor_list_type), POINTER        :: sab_orb_neighbor_list
    TYPE(neighbor_node_type), POINTER        :: sab_orb_neighbor_node
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(realspace_grid_p_type), &
      DIMENSION(:), POINTER                  :: rs_rho

    CALL timeset(routineN,handle)

    failure=.FALSE.

    ! by default, the full density is calculated
    my_soft=.FALSE.
    IF (PRESENT(soft_valid)) my_soft = soft_valid

    IF( PRESENT(basis_set_id)) THEN
      my_basis_set_id = basis_set_id
    ELSE
      my_basis_set_id = use_orb_basis_set
    END IF

    SELECT CASE (my_basis_set_id)
    CASE (use_orb_basis_set)
      CALL get_qs_env(qs_env=qs_env,&
           atomic_kind_set=atomic_kind_set,&
           cell=cell,&
           particle_set=particle_set,&
           dft_control=dft_control,&
           sab_orb=sab_orb,&
           pw_env=pw_env,error=error)
    CASE (use_aux_fit_basis_set)
      CALL get_qs_env(qs_env=qs_env,&
           atomic_kind_set=atomic_kind_set,&
           cell=cell,&
           particle_set=particle_set,&
           dft_control=dft_control,&
           sab_aux_fit=sab_orb,&
           pw_env=pw_env,error=error)
    END SELECT
    CALL pw_env_get(pw_env, rs_pools=rs_pools, error=error)

    ! *** assign from pw_env
    gridlevel_info=>pw_env%gridlevel_info
    cube_info=>pw_env%cube_info

    CALL rs_pools_create_rs_vect(rs_pools, rs_rho, error=error)

    ! find maximum numbers
    nkind = SIZE(atomic_kind_set)
    natoms = SIZE( particle_set )
    maxset=0
    maxpgf=0
    DO ikind=1,nkind
       atomic_kind => atomic_kind_set(ikind)
       SELECT CASE (my_basis_set_id)
       CASE (use_orb_basis_set)
         CALL get_atomic_kind(atomic_kind=atomic_kind,&
              softb = my_soft, &
              orb_basis_set=orb_basis_set)
       CASE (use_aux_fit_basis_set)
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
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

    ! free the atom_pair lists if allocated
    IF (ASSOCIATED(atom_pair_send)) DEALLOCATE(atom_pair_send)
    IF (ASSOCIATED(atom_pair_recv)) DEALLOCATE(atom_pair_recv)

    ! construct a list of all tasks
    IF (.NOT.ASSOCIATED(tasks)) THEN
       CALL reallocate(tasks,1,6,1,max_tasks)
       CALL reallocate(dist_ab,1,3,1,max_tasks)
    ENDIF
    ntasks = 0
    curr_tasks = SIZE(tasks,2)

    loop_ikind: DO ikind=1,nkind

       atomic_kind => atomic_kind_set(ikind)
       SELECT CASE (my_basis_set_id)
       CASE (use_orb_basis_set)
         CALL get_atomic_kind(atomic_kind=atomic_kind,&
              softb = my_soft, orb_basis_set=orb_basis_set)
       CASE (use_aux_fit_basis_set)
         CALL get_atomic_kind(atomic_kind=atomic_kind,&
              softb = my_soft, aux_fit_basis_set=orb_basis_set,&
              basis_set_id=my_basis_set_id)
       END SELECT

       IF (.NOT.ASSOCIATED(orb_basis_set)) CYCLE

       CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
            kind_radius=kind_radius_a,&
            npgf=npgfa,&
            nset=nseta,&
            lmax=la_max,&
            lmin=la_min,&
            pgf_radius=rpgfa,&
            set_radius=set_radius_a,&
            zet=zeta)

       loop_jkind: DO jkind=1,nkind

          atomic_kind => atomic_kind_set(jkind)
          SELECT CASE (my_basis_set_id)
          CASE (use_orb_basis_set)
            CALL get_atomic_kind(atomic_kind=atomic_kind,&
                 softb = my_soft, orb_basis_set=orb_basis_set)
          CASE (use_aux_fit_basis_set)
            CALL get_atomic_kind(atomic_kind=atomic_kind,&
                 softb = my_soft, aux_fit_basis_set=orb_basis_set,&
                 basis_set_id=my_basis_set_id)
          END SELECT

          IF (.NOT.ASSOCIATED(orb_basis_set)) CYCLE

          CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
               kind_radius=kind_radius_b,&
               npgf=npgfb,&
               nset=nsetb,&
               lmax=lb_max,&
               lmin=lb_min,&
               pgf_radius=rpgfb,&
               set_radius=set_radius_b,&
               zet=zetb)
          ab = ikind + nkind*(jkind - 1)

          IF (ASSOCIATED(sab_orb(ab)%neighbor_list_set)) THEN
             CALL get_neighbor_list_set(neighbor_list_set=&
                  sab_orb(ab)%neighbor_list_set, nlist=nlist)
             sab_orb_neighbor_list => first_list(sab_orb(ab)%neighbor_list_set)
          ELSE
             nlist=0
          END IF

          loop_ilist: DO ilist = 1, nlist

             CALL get_neighbor_list(neighbor_list=sab_orb_neighbor_list,&
                  atom=iatom,nnode=nnode)

             ra(:) = pbc(particle_set(iatom)%r,cell)

             sab_orb_neighbor_node => first_node(sab_orb_neighbor_list)

             loop_inode: DO inode = 1, nnode

                CALL get_neighbor_node(neighbor_node=sab_orb_neighbor_node,&
                     neighbor=jatom, r=rab(:))

                CALL task_list_inner_loop(tasks, dist_ab, ntasks, curr_tasks, rs_rho,dft_control,cube_info,gridlevel_info,cell,&
                                  iatom,jatom,rpgfa,rpgfb,zeta,zetb,kind_radius_b,set_radius_a,set_radius_b,ra,rab,&
                                  la_max,la_min,lb_max,lb_min,npgfa,npgfb,maxpgf,maxset,natoms,nseta,nsetb,error)

                sab_orb_neighbor_node => next(sab_orb_neighbor_node)

             END DO loop_inode

             sab_orb_neighbor_list => next(sab_orb_neighbor_list)

          END DO loop_ilist

       END DO loop_jkind

    END DO loop_ikind

    ! redistribute the task list so that all tasks map on the local rs grids
    CALL distribute_tasks (rs_rho,ntasks, natoms,&
            maxset,maxpgf, tasks, rval=dist_ab, atom_pair_send=atom_pair_send,&
            atom_pair_recv=atom_pair_recv, symmetric=.TRUE., error=error)

    CALL rs_pools_give_back_rs_vect(rs_pools, rs_rho, error=error)

    CALL timestop(handle)

  END SUBROUTINE generate_qs_task_list

! *****************************************************************************
!> \par History
!>      Joost VandeVondele: 10.2008 refactored
! *****************************************************************************
  SUBROUTINE task_list_inner_loop(tasks, dist_ab, ntasks, curr_tasks, rs_rho,dft_control,cube_info,gridlevel_info,cell,&
                                  iatom,jatom,rpgfa,rpgfb,zeta,zetb,kind_radius_b,set_radius_a,set_radius_b,ra,rab,&
                                  la_max,la_min,lb_max,lb_min,npgfa,npgfb,maxpgf,maxset,natoms,nseta,nsetb,error)

    INTEGER(kind=int_8), DIMENSION(:, :), &
      POINTER                                :: tasks
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: dist_ab
    INTEGER                                  :: ntasks, curr_tasks
    TYPE(realspace_grid_p_type), &
      DIMENSION(:), POINTER                  :: rs_rho
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(cube_info_type), DIMENSION(:), &
      POINTER                                :: cube_info
    TYPE(gridlevel_info_type), POINTER       :: gridlevel_info
    TYPE(cell_type), POINTER                 :: cell
    INTEGER                                  :: iatom, jatom
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
                IF (dab.lt.0.1E0_dp .AND. dft_control%qs_control%map_paa) THEN
                   igrid_level = 1
                ELSE
                   igrid_level = gaussian_gridlevel(gridlevel_info,zetp)
                ENDIF

                CALL compute_pgf_properties(cube_center,lb_cube,ub_cube,radius,&
                           cell,rs_rho(igrid_level)%rs_grid,cube_info(igrid_level),&
                           la_max(iset),zeta(ipgf,iset),la_min(iset),&
                           lb_max(jset),zetb(jpgf,jset),lb_min(jset),&
                           ra,rab,rab2,dft_control%qs_control%eps_rho_rspace)
   
                CALL pgf_to_tasks(tasks,dist_ab,ntasks,curr_tasks,&
                          rab,iatom,jatom,iset,jset,ipgf,jpgf,natoms,maxset,maxpgf,la_max(iset),lb_max(jset),&
                          rs_rho(igrid_level)%rs_grid,igrid_level,cube_center,lb_cube,ub_cube,error)

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
!> \par History
!>      10.2008 refactored [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE compute_pgf_properties(cube_center,lb_cube,ub_cube,radius, &
               cell,rs_grid,cube_info,la_max,zeta,la_min,lb_max,zetb,lb_min,ra,rab,rab2,eps)

    INTEGER, DIMENSION(3), INTENT(OUT)       :: cube_center, lb_cube, ub_cube
    REAL(KIND=dp), INTENT(OUT)               :: radius
    TYPE(cell_type), POINTER                 :: cell
    TYPE(realspace_grid_type), POINTER       :: rs_grid
    TYPE(cube_info_type), INTENT(IN)         :: cube_info
    INTEGER, INTENT(IN)                      :: la_max
    REAL(KIND=dp), INTENT(IN)                :: zeta
    INTEGER, INTENT(IN)                      :: la_min, lb_max
    REAL(KIND=dp), INTENT(IN)                :: zetb
    INTEGER, INTENT(IN)                      :: lb_min
    REAL(KIND=dp), INTENT(IN)                :: ra(3), rab(3), rab2, eps

    INTEGER, DIMENSION(:), POINTER           :: sphere_bounds
    REAL(KIND=dp)                            :: cutoff, extent(3), f, &
                                                prefactor, rb(3), zetp
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

    CALL compute_cube_center(cube_center,rs_grid,zeta,zetb,ra,rab)
    ! compute cube_center, the center of the gaussian product to map (folded to within the unit cell)
    cube_center(:) = MODULO ( cube_center(:), rs_grid%npts(:) )
    cube_center(:) = cube_center(:) + rs_grid%lb(:)

    IF (rs_grid%orthorhombic ) THEN
       CALL return_cube(cube_info,radius,lb_cube,ub_cube,sphere_bounds)
    ELSE
       ! untested so far
       CALL return_cube_nonortho(cube_info,radius,lb_cube,ub_cube,rp)
       extent(:)  = ub_cube(:) - lb_cube(:)          
       lb_cube(:) = -extent(:) / 2 - 1          
       ub_cube(:) = extent(:) / 2
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
!>
!> \par History
!>      10.2008 Refactored based on earlier routines by MattW [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE pgf_to_tasks(tasks,dist_ab,ntasks,curr_tasks,&
                          rab,iatom,jatom,iset,jset,ipgf,jpgf,natoms,maxset,maxpgf,la_max,lb_max,&
                          rs_grid,igrid_level,cube_center,lb_cube,ub_cube,error)

    INTEGER(kind=int_8), DIMENSION(:, :), &
      POINTER                                :: tasks
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: dist_ab
    INTEGER, INTENT(INOUT)                   :: ntasks, curr_tasks
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rab
    INTEGER, INTENT(IN)                      :: iatom, jatom, iset, jset, &
                                                ipgf, jpgf, natoms, maxset, &
                                                maxpgf, la_max, lb_max
    TYPE(realspace_grid_type), POINTER       :: rs_grid
    INTEGER, INTENT(IN)                      :: igrid_level
    INTEGER, DIMENSION(3), INTENT(IN)        :: cube_center, lb_cube, ub_cube
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER, PARAMETER                       :: add_tasks = 1000
    REAL(kind=dp), PARAMETER                 :: mult_tasks = 2.0_dp

    INTEGER                                  :: added_tasks, cost, j, lmax
    LOGICAL                                  :: is_ortho
    REAL(KIND=dp)                            :: fraction

    ntasks = ntasks + 1
    IF ( ntasks > curr_tasks ) THEN
       curr_tasks = (curr_tasks+add_tasks)*mult_tasks
       CALL reallocate(tasks,1,6,1,curr_tasks)
    END IF

    IF(rs_grid%distributed) THEN

       ! finds the node(s) that need to process this task
       ! on exit tasks(4,:) is 1 for distributed tasks and 2 for generalised tasks
       CALL rs_find_node(rs_grid,cube_center, &
            ntasks=ntasks,tasks=tasks,lb_cube=lb_cube,ub_cube=ub_cube,added_tasks=added_tasks,error=error)

    ELSE
       tasks(1,ntasks) = rs_grid%my_pos
       tasks(4,ntasks) = 0
       tasks(6,ntasks) = 0
       added_tasks = 1
    ENDIF

    IF (SIZE(dist_ab,2).NE.SIZE(tasks,2)) THEN
       CALL reallocate(dist_ab,1,3,1,SIZE(tasks,2))
    ENDIF

    lmax=la_max+lb_max
    is_ortho=(tasks(4,ntasks)==0 .OR. tasks(4,ntasks)==1) .AND. rs_grid%orthorhombic
    ! we assume the load is shared equally between processes dealing with a generalised Gaussian.
    ! this could be refined in the future
    fraction=1.0_dp/added_tasks

    cost=cost_model(lb_cube,ub_cube,fraction,lmax,is_ortho)

    DO j= 1, added_tasks

       tasks (2,ntasks-added_tasks+j) = rs_grid % my_pos
       tasks (5,ntasks-added_tasks+j) = cost

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
!> \note
!>      will hopefully not overflow any more
!> \par History
!>      11.2007 created [Joost]
! *****************************************************************************
  SUBROUTINE pair2int(res,ilevel,iatom,jatom,iset,jset,ipgf,jpgf,natom,maxset,maxpgf)
    INTEGER(kind=int_8)                      :: res
    INTEGER                                  :: ilevel, iatom, jatom, iset, &
                                                jset, ipgf, jpgf, natom, &
                                                maxset, maxpgf

    INTEGER(kind=int_8)                      :: maxpgf8, maxset8, natom8

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
  SUBROUTINE int2pair(res,ilevel,iatom,jatom,iset,jset,ipgf,jpgf,natom,maxset,maxpgf)
    INTEGER(kind=int_8)                      :: res
    INTEGER                                  :: ilevel, iatom, jatom, iset, &
                                                jset, ipgf, jpgf, natom, &
                                                maxset, maxpgf

    INTEGER(kind=int_8)                      :: iatom8, ijatom, ijpgf, ijset, &
                                                ipgf8, iset8, jatom8, jpgf8, &
                                                jset8, maxpgf8, maxset8, &
                                                natom8, tmp

    natom8=natom ; maxset8=maxset ; maxpgf8=maxpgf
    
    ilevel=res/(natom8**2*maxset8**2*maxpgf8**2)
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

    iatom=iatom8 ; jatom=jatom8; iset=iset8 ; jset=jset8 
    ipgf=ipgf8 ; jpgf=jpgf8

  END SUBROUTINE int2pair

! *****************************************************************************
!> \brief performs load balancing of the tasks on the distributed grids
!> \par History
!>      created 2008-10-03 [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE load_balance_distributed (tasks, ntasks, rs_rho, natoms, maxset,maxpgf, error)

    INTEGER(kind=int_8), DIMENSION(:, :), &
      POINTER                                :: tasks
    INTEGER                                  :: ntasks
    TYPE(realspace_grid_p_type), &
      DIMENSION(:), POINTER                  :: rs_rho
    INTEGER                                  :: natoms, maxset, maxpgf
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'load_balance_distributed', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, stat
    INTEGER, DIMENSION(:, :, :), POINTER     :: list
    LOGICAL                                  :: failure

    CALL timeset(routineN,handle)

    failure=.FALSE.
    NULLIFY(list)
    ! here we create for each cpu (0:ncpu-1) a list of possible destinations.
    ! if a destination would not be in this list, it is a bug
    CALL create_destination_list(list,rs_rho,error)

    ! now, walk over the tasks, filling in the loads of each destination
    CALL compute_load_list(list,rs_rho,tasks,ntasks,natoms, maxset,maxpgf, create_list=.TRUE., error=error)

    ! optimize loads & fluxes
    CALL optimize_load_list(list,rs_rho(1)%rs_grid%group,rs_rho(1)%rs_grid%my_pos, error)

    ! now, walk over the tasks, using the list to set the destinations
    CALL compute_load_list(list,rs_rho,tasks,ntasks,natoms, maxset,maxpgf, create_list=.FALSE., error=error)

    DEALLOCATE(list,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE load_balance_distributed


! *****************************************************************************
!> \brief this serial routine adjusts the fluxes in the global list
!>
!> \par History
!>      created 2008-10-06 [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE balance_global_list(list_global)
    INTEGER, DIMENSION(:, :, 0:)             :: list_global

    INTEGER, PARAMETER                       :: Max_Iter = 100
    REAL(KIND=dp), PARAMETER                 :: Tolerance_factor = 0.005_dp

    INTEGER                                  :: dest, icpu, idest, iflux, &
                                                ilocal, k, maxdest, Ncpu, &
                                                Nflux
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: flux_connections
    LOGICAL                                  :: solution_ok, solution_optimal
    REAL(KIND=dp)                            :: average, load_shift, &
                                                max_load_shift, tolerance
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: load, optimized_flux, &
                                                optimized_load
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: flux_limits

     Ncpu=SIZE(list_global,3)
     maxdest=SIZE(list_global,2)
     ALLOCATE(load(0:Ncpu-1))
     load=0.0_dp
     ALLOCATE(optimized_load(0:Ncpu-1))


     ! figure out the number of fluxes
     ! we assume that the global_list is symmetric
     Nflux=0
     DO icpu=0,ncpu-1
        DO idest=1,maxdest
           dest=list_global(1,idest,icpu)
           IF (dest<ncpu .AND. dest>icpu) Nflux=Nflux+1
        ENDDO
     ENDDO
     ALLOCATE(optimized_flux(Nflux))
     ALLOCATE(flux_limits(2,Nflux))
     ALLOCATE(flux_connections(2,Nflux))

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

     solution_ok=.TRUE.
     solution_optimal=.FALSE.
     optimized_flux=0.0_dp

     ! an iterative solver, if iterated till convergence the maximum load is minimal
     ! we terminate before things are fully converged, since this does show up in the timings
     ! once the largest shift becomes less than a small fraction of the average load, we're done 
     ! we're perfectly happy if the load balance is within 1 percent or so
     ! the maximum load normally converges even faster
     average=SUM(load)/SIZE(load)
     tolerance=Tolerance_factor*average

     optimized_load=load
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
                                                  list_global(2,idest,icpu)-optimized_flux(Nflux)
                       list_global(2,idest,icpu)=optimized_flux(Nflux)
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
                                                        list_global(2,idest,icpu)+optimized_flux(iflux)
                             list_global(2,idest,icpu)=-optimized_flux(iflux)
                          ENDIF
                          EXIT
                       ENDIF
                    ENDDO
                 ENDIF
              ENDIF
           ENDIF
        ENDDO
     ENDDO

  END SUBROUTINE balance_global_list

! *****************************************************************************
!> \brief this routine gets back optimized loads for all destinations
!>
!> \par History
!>      created 2008-10-06 [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE optimize_load_list(list,group,my_pos,error)
    INTEGER, DIMENSION(:, :, 0:)             :: list
    INTEGER, INTENT(IN)                      :: group, my_pos
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'optimize_load_list', &
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

    ALLOCATE(load_all(maxdest*ncpu),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(load_per_source(maxdest*ncpu*ncpu),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(list_global(2,maxdest,ncpu),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    load_all=RESHAPE(list(2,:,:),(/maxdest*ncpu/))
    CALL mp_gather(load_all,load_per_source,0,group)

    IF (my_pos==0) THEN

       ! collect the load of all CPUs
       load_all=0
       DO icpu=1,ncpu
          load_all=load_all+load_per_source(maxdest*ncpu*(icpu-1)+1:maxdest*ncpu*icpu)
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

    DEALLOCATE(load_all,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(load_per_source,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(list_global,STAT=stat)
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
!> \par History
!>      created 2008-10-06 [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE compute_load_list(list,rs_rho,tasks,ntasks,natoms, maxset,maxpgf, create_list, error)
    INTEGER, DIMENSION(:, :, 0:)             :: list
    TYPE(realspace_grid_p_type), &
      DIMENSION(:)                           :: rs_rho
    INTEGER(kind=int_8), DIMENSION(:, :), &
      POINTER                                :: tasks
    INTEGER                                  :: ntasks, natoms, maxset, maxpgf
    LOGICAL                                  :: create_list
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'compute_load_list', &
      routineP = moduleN//':'//routineN

    INTEGER :: dest, handle, i, iatom, ilevel, iopt, ipair_old, ipgf, iset, &
      itask, itask_start, itask_stop, jatom, jpgf, jset, li, maxdest, ncpu, &
      ndest_pair, nopt, nshort, rank
    INTEGER(KIND=int_8)                      :: bit_pattern, ipair, natom8
    INTEGER(KIND=int_8), ALLOCATABLE, &
      DIMENSION(:)                           :: loads
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: all_dests, index
    INTEGER, DIMENSION(6)                    :: options
    LOGICAL                                  :: failure

    CALL timeset(routineN,handle)

    ALLOCATE(loads(0:rs_rho(1)%rs_grid%group_size-1))
    CALL get_current_loads(loads, rs_rho, ntasks, tasks, skip_distributed=.TRUE., error=error)

    failure=.FALSE.

    maxdest=SIZE(list,2)
    ncpu=SIZE(list,3)
    natom8 =natoms

    ! first find the tasks that deal with the same atom pair
    itask_stop=0
    ipair_old=HUGE(ipair_old)
    ALLOCATE(all_dests(0))
    ALLOCATE(INDEX(0))

    DO 

       ! first find the range of tasks that deal with the same atom pair
       itask_start=itask_stop+1
       itask_stop =itask_start
       IF(itask_stop>ntasks) EXIT
       CALL int2pair(tasks(3,itask_stop),ilevel,iatom,jatom,iset,jset,ipgf,jpgf,natoms,maxset,maxpgf)
       ipair_old=(iatom-1)*natom8+(jatom-1)
       DO
          IF(itask_stop+1>ntasks) EXIT
          CALL int2pair(tasks(3,itask_stop+1),ilevel,iatom,jatom,iset,jset,ipgf,jpgf,natoms,maxset,maxpgf)
          ipair=(iatom-1)*natom8+(jatom-1)
          IF (ipair==ipair_old) THEN
             itask_stop=itask_stop+1
          ELSE
             EXIT
          ENDIF
       ENDDO
       ipair=ipair_old
       nshort=itask_stop-itask_start+1

       ! find the unique list of destinations
       DEALLOCATE(all_dests)
       ALLOCATE(all_dests(nshort))
       DEALLOCATE(index)
       ALLOCATE(INDEX(nshort))
       all_dests=tasks(1,itask_start:itask_stop)
       CALL sort(all_dests,nshort,index)
       ndest_pair=1
       DO i=2,nshort
          IF (all_dests(ndest_pair).NE.all_dests(i)) THEN
             ndest_pair=ndest_pair+1
             all_dests(ndest_pair)=all_dests(i)
          ENDIF
       ENDDO

       DO itask=itask_start,itask_stop

          dest=tasks(1,itask) ! notice that dest can be changed
          CALL int2pair(tasks(3,itask),ilevel,iatom,jatom,iset,jset,ipgf,jpgf,natoms,maxset,maxpgf)
          ipair=(iatom-1)*natom8+(jatom-1)

          SELECT CASE(tasks(4,itask))
          CASE(0)
            IF (.FALSE.) THEN
               ! set the dest to the one with lowest distributed load
               dest=all_dests(1)
               DO iopt=2,ndest_pair
                  IF (loads(dest)>loads(all_dests(iopt))) dest=all_dests(iopt)
               ENDDO

               ! set rank to the one with the lowest distributed load, but not dest if possible
               ! pick a feasible
               rank=dest
               DO iopt=1,maxdest
                  IF (list(1,iopt,dest)==dest .OR. list(1,iopt,dest)>ncpu-1) CYCLE
                  IF (ANY(all_dests(1:ndest_pair).EQ.list(1,iopt,dest))) THEN
                     rank=list(1,iopt,dest)
                  ENDIF
               ENDDO
               ! an optimal one
               DO iopt=1,maxdest
                  IF (list(1,iopt,dest)==dest .OR. list(1,iopt,dest)>ncpu-1) CYCLE
                  IF (ANY(all_dests(1:ndest_pair).EQ.list(1,iopt,dest))) THEN
                     IF (loads(rank)>loads(list(1,iopt,dest))) rank=list(1,iopt,dest)
                  ENDIF
               ENDDO
               li=list_index(list,rank,dest)
               IF (create_list) THEN
                  list(2,li,dest)=list(2,li,dest)+tasks(5,itask)
               ELSE
                  IF (list(1,li,dest)==dest) THEN
                      tasks(1,itask)=dest
                  ELSE
                    IF (list(2,li,dest)>=tasks(5,itask)) THEN
                        list(2,li,dest)=list(2,li,dest)-tasks(5,itask)
                        tasks(1,itask)=list(1,li,dest)
                    ELSE
                        tasks(1,itask)=dest
                    ENDIF
                  ENDIF
               ENDIF
            ELSE
               ! just ignore these replicated tasks   
            ENDIF
          CASE(1)
            bit_pattern=tasks(6,itask)
            nopt=0
            IF (BTEST(bit_pattern,0)) THEN
               rank=rs_grid_locate_rank(rs_rho(ilevel)%rs_grid,dest,(/-1,0,0/))
               IF (ANY(all_dests(1:ndest_pair).EQ.rank)) THEN
                  nopt=nopt+1
                  options(nopt)=rank
               ENDIF
            ENDIF
            IF (BTEST(bit_pattern,1)) THEN
               rank=rs_grid_locate_rank(rs_rho(ilevel)%rs_grid,dest,(/+1,0,0/))
               IF (ANY(all_dests(1:ndest_pair).EQ.rank)) THEN
                  nopt=nopt+1
                  options(nopt)=rank
               ENDIF
            ENDIF
            IF (BTEST(bit_pattern,2)) THEN
               rank=rs_grid_locate_rank(rs_rho(ilevel)%rs_grid,dest,(/0,-1,0/))
               IF (ANY(all_dests(1:ndest_pair).EQ.rank)) THEN
                  nopt=nopt+1
                  options(nopt)=rank
               ENDIF
            ENDIF
            IF (BTEST(bit_pattern,3)) THEN
               rank=rs_grid_locate_rank(rs_rho(ilevel)%rs_grid,dest,(/0,+1,0/))
               IF (ANY(all_dests(1:ndest_pair).EQ.rank)) THEN
                  nopt=nopt+1
                  options(nopt)=rank
               ENDIF
            ENDIF
            IF (BTEST(bit_pattern,4)) THEN
               rank=rs_grid_locate_rank(rs_rho(ilevel)%rs_grid,dest,(/0,0,-1/))
               IF (ANY(all_dests(1:ndest_pair).EQ.rank)) THEN
                  nopt=nopt+1
                  options(nopt)=rank
               ENDIF
            ENDIF
            IF (BTEST(bit_pattern,5)) THEN
               rank=rs_grid_locate_rank(rs_rho(ilevel)%rs_grid,dest,(/0,0,+1/))
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
               list(2,li,dest)=list(2,li,dest)+tasks(5,itask)
            ELSE
               IF (list(1,li,dest)==dest) THEN
                   tasks(1,itask)=dest
               ELSE
                 IF (list(2,li,dest)>=tasks(5,itask)) THEN
                     list(2,li,dest)=list(2,li,dest)-tasks(5,itask)
                     tasks(1,itask)=list(1,li,dest)
                 ELSE
                     tasks(1,itask)=dest
                 ENDIF
               ENDIF
            ENDIF
          CASE(2) ! generalised
            li=list_index(list,dest,dest)
            IF (create_list) THEN
               list(2,li,dest)=list(2,li,dest)+tasks(5,itask)
            ELSE
               IF (list(1,li,dest)==dest) THEN
                   tasks(1,itask)=dest
               ELSE
                 IF (list(2,li,dest)>=tasks(5,itask)) THEN
                     list(2,li,dest)=list(2,li,dest)-tasks(5,itask)
                     tasks(1,itask)=list(1,li,dest)
                 ELSE
                     tasks(1,itask)=dest
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
!> \par History
!>      created 2008-10-06 [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE create_destination_list(list,rs_rho,error)
    INTEGER, DIMENSION(:, :, :), POINTER     :: list
    TYPE(realspace_grid_p_type), &
      DIMENSION(:), POINTER                  :: rs_rho
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_destination_list', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, icount, icpu, &
                                                igrid_level, j, maxcount, &
                                                ncpu, ndistributed, stat, &
                                                ultimate_max
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: index, sublist
    INTEGER, DIMENSION(3)                    :: coord
    LOGICAL                                  :: failure

    failure=.FALSE.

    CALL timeset(routineN,handle)

    CPPrecondition(.NOT.ASSOCIATED(list),cp_failure_level,routineP,error,failure)
    ncpu=rs_rho(1)%rs_grid%group_size
    ndistributed=0
    DO igrid_level=1,SIZE(rs_rho)
       IF (rs_rho(igrid_level)%rs_grid%distributed) ndistributed=ndistributed+1
    ENDDO
    ultimate_max=7*ndistributed

    ! each grid is distributed and yields different neighbors
    ALLOCATE(list(2,ultimate_max,0:ncpu-1),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE(INDEX(ultimate_max),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(sublist(ultimate_max),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    sublist=HUGE(sublist)

    maxcount=1
    DO icpu=0,ncpu-1
       icount=1
       sublist(1)=icpu
       DO igrid_level=1,SIZE(rs_rho)
          IF (rs_rho(igrid_level)%rs_grid%distributed) THEN
             coord=rs_rho(igrid_level)%rs_grid%rank2coord(:,icpu)  
             sublist(icount+1)=rs_grid_locate_rank(rs_rho(igrid_level)%rs_grid,icpu,(/-1,0,0/))    
             sublist(icount+2)=rs_grid_locate_rank(rs_rho(igrid_level)%rs_grid,icpu,(/+1,0,0/))    
             sublist(icount+3)=rs_grid_locate_rank(rs_rho(igrid_level)%rs_grid,icpu,(/0,-1,0/))    
             sublist(icount+4)=rs_grid_locate_rank(rs_rho(igrid_level)%rs_grid,icpu,(/0,+1,0/))    
             sublist(icount+5)=rs_grid_locate_rank(rs_rho(igrid_level)%rs_grid,icpu,(/0,0,-1/))    
             sublist(icount+6)=rs_grid_locate_rank(rs_rho(igrid_level)%rs_grid,icpu,(/0,0,+1/))    
             icount=icount+6
          ENDIF
       END DO
       ! only retain unique values of the destination
       CALL sort(sublist,icount,index)
       j=1
       DO i=2,icount
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
!>        XXXXXXXXXX should be implemented with scatter / gather 
!> \par History
!>      none
!> \author MattW 21/11/2007
! *****************************************************************************
  SUBROUTINE get_current_loads(loads, rs_rho, ntasks, tasks, skip_distributed, error)
    INTEGER(kind=int_8), DIMENSION(:)        :: loads
    TYPE(realspace_grid_p_type), &
      DIMENSION(:)                           :: rs_rho
    INTEGER                                  :: ntasks
    INTEGER(kind=int_8), DIMENSION(:, :), &
      POINTER                                :: tasks
    LOGICAL                                  :: skip_distributed
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_current_loads', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, stat
    INTEGER(KIND=int_8)                      :: total_cost_local
    INTEGER(kind=int_8), ALLOCATABLE, &
      DIMENSION(:)                           :: recv_buf_i, send_buf_i
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: recv_disps, recv_sizes, &
                                                send_disps, send_sizes
    TYPE(realspace_grid_type), POINTER       :: rs

    CALL timeset(routineN,handle)

    rs=>rs_rho(1)%rs_grid

    ! allocate local arrays
    ALLOCATE(send_sizes(rs%group_size),STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","send_sizes",rs%group_size)
    ALLOCATE(recv_sizes(rs%group_size),STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","recv_sizes",rs%group_size)
    ALLOCATE(send_disps(rs%group_size),STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","send_disps",rs%group_size)
    ALLOCATE(recv_disps(rs%group_size),STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","recv_disps",rs%group_size)
    ALLOCATE(send_buf_i(rs%group_size),STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","send_buf",rs%group_size)
    ALLOCATE(recv_buf_i(rs%group_size),STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","recv_buf",rs%group_size)

!
! part 1, load balance on the replicated tasks, i.e. after this step replicated tasks can have
!         a destination that is different from their source
!

    send_sizes = 1
    recv_sizes = 1
    send_disps(1)=0
    recv_disps(1)=0
    DO i = 2, rs % group_size
       send_disps(i) = send_disps(i-1) + send_sizes(i-1)
       recv_disps(i) = recv_disps(i-1) + recv_sizes(i-1)
    ENDDO

! communication step 1 : compute the total local cost of the tasks
!                        each proc needs to know the amount of work he will receive


    ! send buffer now contains for each target the cost of the tasks it will receive
    send_buf_i=0
    DO i = 1, ntasks
       IF (skip_distributed .AND. tasks(4,i)==0) CYCLE
       send_buf_i( tasks(1,i)+1 ) = send_buf_i( tasks(1,i)+1 ) + tasks(5,i)
    END DO
    CALL mp_alltoall(send_buf_i, send_sizes, send_disps, recv_buf_i, recv_sizes, recv_disps, rs % group)
    total_cost_local = SUM(recv_buf_i)

! communication step 2 : compute the global cost of the tasks

    ! after this step, the recv buffer contains the local cost for each CPU
    send_buf_i = total_cost_local
    CALL mp_alltoall(send_buf_i, send_sizes, send_disps, recv_buf_i, recv_sizes, recv_disps, rs % group)

    loads=recv_buf_i

    DEALLOCATE(send_sizes,STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","send_sizes",rs%group_size)
    DEALLOCATE(recv_sizes,STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","recv_sizes",rs%group_size)
    DEALLOCATE(send_disps,STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","send_disps",rs%group_size)
    DEALLOCATE(recv_disps,STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","recv_disps",rs%group_size)

    CALL timestop(handle)

  END SUBROUTINE get_current_loads
! *****************************************************************************
!> \brief performs load balancing shifting tasks on the replicated grids
!>        this modifies the destination of some of the tasks on replicated
!>        grids, and in this way balances the load
!>        XXXXXXXXXXXXXXX should be implemented with scatter gather instead of all to all
!> \par History
!>      none
!> \author MattW 21/11/2007
! *****************************************************************************
  SUBROUTINE load_balance_replicated ( rs_rho, ntasks, tasks, error)

    TYPE(realspace_grid_p_type), &
      DIMENSION(:), POINTER                  :: rs_rho
    INTEGER                                  :: ntasks
    INTEGER(kind=int_8), DIMENSION(:, :), &
      POINTER                                :: tasks
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'load_balance_replicated', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, j, no_overloaded, &
                                                no_underloaded, &
                                                proc_receiving, stat
    INTEGER(KIND=int_8)                      :: average_cost, cost_task_rep, &
                                                count, offset, &
                                                total_cost_global
    INTEGER(kind=int_8), ALLOCATABLE, &
      DIMENSION(:)                           :: load_imbalance, recv_buf_i, &
                                                send_buf_i
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: index, recv_disps, &
                                                recv_sizes, send_disps, &
                                                send_sizes
    TYPE(realspace_grid_type), POINTER       :: rs

    CALL timeset(routineN,handle)

    rs=>rs_rho(1)%rs_grid

    ! allocate local arrays
    ALLOCATE(send_sizes(rs%group_size),STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","send_sizes",rs%group_size)
    ALLOCATE(recv_sizes(rs%group_size),STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","recv_sizes",rs%group_size)
    ALLOCATE(send_disps(rs%group_size),STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","send_disps",rs%group_size)
    ALLOCATE(recv_disps(rs%group_size),STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","recv_disps",rs%group_size)
    ALLOCATE(send_buf_i(rs%group_size),STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","send_buf",rs%group_size)
    ALLOCATE(recv_buf_i(rs%group_size),STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","recv_buf",rs%group_size)
    send_sizes = 1
    recv_sizes = 1
    send_disps(1)=0
    recv_disps(1)=0
    DO i = 2, rs % group_size
       send_disps(i) = send_disps(i-1) + send_sizes(i-1)
       recv_disps(i) = recv_disps(i-1) + recv_sizes(i-1)
    ENDDO

    CALL get_current_loads ( recv_buf_i, rs_rho, ntasks, tasks, skip_distributed=.FALSE., error=error )

    total_cost_global = SUM(recv_buf_i)
    average_cost = total_cost_global / rs%group_size

!
! compute how to redistribute the replicated tasks so that the average cost is reached
!

    ! load imbalance measures the load of a given CPU relative to the optimal load distribution (load=average)
    ALLOCATE(load_imbalance(rs%group_size),STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","index",rs%group_size)
    ALLOCATE(INDEX(rs%group_size),STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","load_imbalance",rs%group_size)

    load_imbalance = recv_buf_i - average_cost
    no_overloaded = 0
    no_underloaded = 0

    DO i = 1, rs%group_size
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
       IF(tasks(4,i) .EQ. 0 .AND. tasks(1,i)==tasks(2,i) ) THEN
          cost_task_rep = cost_task_rep + tasks(5,i)
       ENDIF
    ENDDO

    ! now, correct the load imbalance for the overloaded CPUs
    ! they will send away not more than the total load of replicated tasks
    send_buf_i = cost_task_rep
    CALL mp_alltoall(send_buf_i, send_sizes, send_disps, recv_buf_i, recv_sizes, recv_disps, rs % group)
    DO i = 1, rs%group_size
       ! At the moment we can only offload replicated tasks
       IF(load_imbalance(i) .GT. 0)&
            load_imbalance(i) = MIN( load_imbalance(i),recv_buf_i(i))
    ENDDO

    ! simplest algorithm I can think of of is that the processor with the most
    ! excess tasks fills up the process needing most, then moves on to next most.
    ! At the moment if we've got less replicated tasks than we're overloaded then 
    ! task balancing will be incomplete

    ! only need to do anything if I've excess tasks
    IF (load_imbalance( rs%my_pos + 1 ) .GT. 0 ) THEN

       count = 0 ! weighted amount of tasks offloaded
       offset = 0 ! no of underloaded processes already filled by other more overloaded procs

       ! calculate offset
       DO i = rs%group_size, rs%group_size-no_overloaded+1, -1
          IF ( INDEX(i) .EQ. rs%my_pos+1 ) THEN
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
       ! we fill this up by adjusting the destination of tasks on the replicated grid, then move to next most underloaded proc
       DO j = 1, ntasks
          IF (tasks(4,j) .EQ. 0 .AND. tasks(1,j)==tasks(2,j) ) THEN
             !just avoid sending to non existing procs due to integer truncation in the computation of the average
             IF(proc_receiving .GT. no_underloaded) EXIT
             ! set new destination
             tasks (1,j) = INDEX ( proc_receiving ) - 1
             offset = offset + tasks(5,j)
             count = count + tasks(5,j)
             IF(count .GE. load_imbalance (rs%my_pos + 1 ) ) EXIT
             IF(offset .GT. 0 ) THEN
                proc_receiving = proc_receiving + 1
                !just avoid sending to non existing procs due to integer truncation in the computation of the average
                IF(proc_receiving .GT. no_underloaded) EXIT
                offset = load_imbalance ( INDEX ( proc_receiving ) )
             ENDIF
          ENDIF
       ENDDO
    ENDIF

    DEALLOCATE(index,STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","index",rs%group_size)
    DEALLOCATE(load_imbalance,STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","load_imbalance",rs%group_size)

    DEALLOCATE(send_sizes,STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","send_sizes",rs%group_size)
    DEALLOCATE(recv_sizes,STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","recv_sizes",rs%group_size)
    DEALLOCATE(send_disps,STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","send_disps",rs%group_size)
    DEALLOCATE(recv_disps,STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","recv_disps",rs%group_size)

    CALL timestop(handle)

  END SUBROUTINE load_balance_replicated

! *****************************************************************************
!> \brief given an input task list, redistribute so that all tasks can be processed locally,
!>        i.e. dest equals rank
!> \par History
!>      none
!> \author MattW 21/11/2007
! *****************************************************************************
  SUBROUTINE create_local_tasks ( rs_rho, ntasks, tasks, rval, ntasks_recv, tasks_recv, rval_recv, error)

    TYPE(realspace_grid_p_type), &
      DIMENSION(:), POINTER                  :: rs_rho
    INTEGER                                  :: ntasks
    INTEGER(kind=int_8), DIMENSION(:, :), &
      POINTER                                :: tasks
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: rval
    INTEGER                                  :: ntasks_recv
    INTEGER(kind=int_8), DIMENSION(:, :), &
      POINTER                                :: tasks_recv
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: rval_recv
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_local_tasks', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, j, k, l, stat, &
                                                task_dim
    INTEGER(kind=int_8), ALLOCATABLE, &
      DIMENSION(:)                           :: recv_buf_i, send_buf_i
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: recv_disps, recv_sizes, &
                                                send_disps, send_sizes
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: recv_buf_r, send_buf_r
    TYPE(realspace_grid_type), POINTER       :: rs

    CALL timeset(routineN,handle)

    rs=>rs_rho(1)%rs_grid

    ! allocate local arrays
    ALLOCATE(send_sizes(rs%group_size),STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","send_sizes",rs%group_size)
    ALLOCATE(recv_sizes(rs%group_size),STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","recv_sizes",rs%group_size)
    ALLOCATE(send_disps(rs%group_size),STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","send_disps",rs%group_size)
    ALLOCATE(recv_disps(rs%group_size),STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","recv_disps",rs%group_size)

    send_sizes = 1
    recv_sizes = 1
    send_disps(1)=0
    recv_disps(1)=0
    DO i = 2, rs % group_size
       send_disps(i) = send_disps(i-1) + send_sizes(i-1)
       recv_disps(i) = recv_disps(i-1) + recv_sizes(i-1)
    ENDDO
    ALLOCATE(send_buf_i(rs%group_size),STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","send_buf",rs%group_size)
    ALLOCATE(recv_buf_i(rs%group_size),STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","recv_buf",rs%group_size)

    ! fill send buffer, now counting how many tasks will be send
    send_buf_i=0
    DO i = 1, ntasks
       send_buf_i( tasks(1,i)+1 ) = send_buf_i( tasks(1,i)+1 ) + 1 
    END DO

    CALL mp_alltoall(send_buf_i, send_sizes, send_disps, recv_buf_i, recv_sizes, recv_disps, rs % group)
    
    ! communication step 3 : pack the tasks, and send them around

    task_dim = SIZE(tasks,1)

    send_sizes = 0
    send_disps = 0
    recv_sizes = 0
    recv_disps = 0

    send_sizes(1) = send_buf_i(1) * task_dim
    recv_sizes(1) = recv_buf_i(1) * task_dim
    DO i = 2,rs%group_size
       send_sizes(i) = send_buf_i(i) * task_dim
       recv_sizes(i) = recv_buf_i(i) * task_dim
       send_disps(i)=send_disps(i-1)+send_sizes(i-1)           
       recv_disps(i)=recv_disps(i-1)+recv_sizes(i-1)
    ENDDO

    ! deallocate old send/recv buffers
    DEALLOCATE(send_buf_i,STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","send_buf_i",rs%group_size)
    DEALLOCATE(recv_buf_i,STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","recv_buf_i",rs%group_size)

    ! allocate them with new sizes
    ALLOCATE(send_buf_i(SUM(send_sizes)),STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","send_buf",rs%group_size)
    ALLOCATE(recv_buf_i(SUM(recv_sizes)),STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","send_buf",rs%group_size)

    ! do packing
    send_buf_i = 0
    send_sizes = 0
    DO j = 1,ntasks
       i = tasks(1,j)+1
       DO k=1,task_dim
          send_buf_i(send_disps(i)+send_sizes(i)+k)=tasks(k,j)
       ENDDO
       send_sizes(i)=send_sizes(i) + task_dim
    ENDDO

    ! do communication
    CALL mp_alltoall(send_buf_i, send_sizes, send_disps, recv_buf_i, recv_sizes, recv_disps, rs % group)

    DEALLOCATE(send_buf_i,STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","send_buf",rs%group_size)

    ntasks_recv=SUM(recv_sizes)/task_dim
    ALLOCATE(tasks_recv(task_dim,ntasks_recv),STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","tasks_recv",rs%group_size)

    ! do unpacking
    l = 0
    DO i = 1,rs % group_size
       DO j = 0,recv_sizes(i)/task_dim-1
          l = l + 1
          DO k=1,task_dim 
             tasks_recv(k,l)=recv_buf_i(recv_disps(i)+j*task_dim+k)
          ENDDO
       ENDDO
    ENDDO

    DEALLOCATE(recv_buf_i,STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","recv_buf",rs%group_size)

    ! send rvals (to be removed :-)

    ! take care of the new task_dim in the send_sizes
    send_sizes=(send_sizes/task_dim)*SIZE(rval,1)
    recv_sizes=(recv_sizes/task_dim)*SIZE(rval,1)
    send_disps=(send_disps/task_dim)*SIZE(rval,1)
    recv_disps=(recv_disps/task_dim)*SIZE(rval,1)
    task_dim=SIZE(rval,1)

    ALLOCATE(send_buf_r(SUM(send_sizes)),STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","send_buf_r",rs%group_size)
    ALLOCATE(recv_buf_r(SUM(recv_sizes)),STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","recv_buf_r",rs%group_size)

    !do packing
    send_sizes = 0
    DO j = 1,ntasks
       i = tasks(1,j)+1
       DO k=1,task_dim
          send_buf_r(send_disps(i)+send_sizes(i)+k)=rval(k,j)
       ENDDO
       send_sizes(i)=send_sizes(i) + task_dim
    ENDDO

    ! do communication
    CALL mp_alltoall(send_buf_r, send_sizes, send_disps,&
         recv_buf_r, recv_sizes, recv_disps, rs % group)

    DEALLOCATE(send_buf_r,STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","send_buf_r",rs%group_size)
    ALLOCATE(rval_recv(task_dim,SUM(recv_sizes)/task_dim),STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","rval_recv",rs%group_size)

    ! do unpacking
    l = 0
    DO i = 1,rs % group_size
       DO j = 0,recv_sizes(i)/task_dim-1
          l = l + 1
          DO k=1,task_dim
             rval_recv(k,l)=recv_buf_r(recv_disps(i)+j*task_dim+k)
          ENDDO
       ENDDO
    ENDDO

    DEALLOCATE(recv_buf_r,STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","recv_buf_r",rs%group_size)

    DEALLOCATE(send_sizes,STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","send_sizes",rs%group_size)
    DEALLOCATE(recv_sizes,STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","recv_sizes",rs%group_size)
    DEALLOCATE(send_disps,STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","send_disps",rs%group_size)
    DEALLOCATE(recv_disps,STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","recv_disps",rs%group_size)

    CALL timestop(handle)

  END SUBROUTINE create_local_tasks

! *****************************************************************************
!> \brief Assembles tasks to be performed on local grid
!> \param rs the grid
!> \param tasks the task set generated on this processor
!> \param npme Number of tasks for local processing
!> \par History
!>      none
!> \author MattW 21/11/2007
! *****************************************************************************
  SUBROUTINE distribute_tasks ( rs_rho, ntasks, natoms,&
       maxset,maxpgf,  tasks, rval, atom_pair_send, atom_pair_recv,&
       symmetric, error)

    TYPE(realspace_grid_p_type), &
      DIMENSION(:), POINTER                  :: rs_rho
    INTEGER                                  :: ntasks, natoms, maxset, maxpgf
    INTEGER(kind=int_8), DIMENSION(:, :), &
      POINTER                                :: tasks
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: rval
    INTEGER(kind=int_8), DIMENSION(:), &
      POINTER                                :: atom_pair_send, atom_pair_recv
    LOGICAL, INTENT(IN)                      :: symmetric
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'distribute_tasks', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, igrid_level, k, &
                                                ntasks_recv, stat
    INTEGER(kind=int_8), ALLOCATABLE, &
      DIMENSION(:)                           :: taskid
    INTEGER(kind=int_8), DIMENSION(:, :), &
      POINTER                                :: tasks_recv
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: index
    LOGICAL                                  :: distributed_grids
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: rval_recv
    TYPE(realspace_grid_type), POINTER       :: rs

    CALL timeset(routineN,handle)

    IF ( .NOT. ASSOCIATED ( tasks ) ) &
         CALL stop_program ( "get_my_tasks", "tasks not associated" )

    ! *** figure out if we have distributed grids
    distributed_grids=.FALSE.
    DO igrid_level=1,SIZE(rs_rho)
       IF ( rs_rho(igrid_level)%rs_grid%distributed ) THEN
          distributed_grids=.TRUE.
       ENDIF
    END DO
    rs=>rs_rho(1)%rs_grid

    IF (distributed_grids) THEN

       ! first round of balancing on the distributed grids
       CALL load_balance_distributed(tasks, ntasks, rs_rho, natoms, maxset,maxpgf, error)

       ! adjusting the replicated task to improve the load balance
       CALL load_balance_replicated(rs_rho, ntasks, tasks, error)

       ! given a list of tasks, this will do the needed reshuffle so that all tasks will be local
       CALL create_local_tasks ( rs_rho, ntasks, tasks, rval, ntasks_recv, tasks_recv, rval_recv, error)

       ! CALL get_current_loads( loads, rs_rho, ntasks, tasks, skip_distributed=.FALSE., error=error)
       ! IF (rs%my_pos==0) THEN
       !     WRITE(6,*) ""
       !     WRITE(6,*) "At the end of the load balancing procedure"
       !     WRITE(6,*) "Maximum load:",MAXVAL(loads)
       !     WRITE(6,*) "Average load:",SUM(loads)/SIZE(loads)
       ! ENDIF

       !
       ! tasks list are complete, we can compute the list of atomic blocks (atom pairs)
       ! we will be sending. These lists are needed for redistribute_matrix.
       !
       ALLOCATE(atom_pair_send(ntasks),STAT=stat)
       IF (stat/=0) CALL stop_memory("get_my_tasks","atom_pair_send",rs%group_size)
       CALL get_atom_pair ( atom_pair_send, tasks, send=.TRUE.,&
            symmetric=symmetric, natoms=natoms, maxset=maxset, maxpgf=maxpgf )
 
       ! natom_send=SIZE(atom_pair_send)
       ! CALL mp_sum(natom_send,rs%group)
       ! IF (rs%my_pos==0) THEN
       !     WRITE(6,*) ""
       !     WRITE(6,*) "Total number of atomic blocks to be send:",natom_send
       ! ENDIF

       ALLOCATE(atom_pair_recv(ntasks_recv),STAT=stat)
       IF (stat/=0) CALL stop_memory("get_my_tasks","atom_pair_recv",rs%group_size)
       CALL get_atom_pair ( atom_pair_recv, tasks_recv, send=.FALSE.,&
            symmetric=symmetric, natoms=natoms, maxset=maxset, maxpgf=maxpgf )

       ! cleanup, at this point we  don't need the original tasks & rvals anymore
       DEALLOCATE(tasks,STAT=stat)
       IF (stat/=0) CALL stop_memory("get_my_tasks","tasks",rs%group_size)
       DEALLOCATE(rval,STAT=stat)
       IF (stat/=0) CALL stop_memory("get_my_tasks","rval",rs%group_size)

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

    ALLOCATE(taskid(ntasks_recv),STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","taskid",rs%group_size)
    ALLOCATE(INDEX(ntasks_recv),STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","index",rs%group_size)

    taskid=tasks_recv(3,1:ntasks_recv)
    CALL sort(taskid,SIZE(taskid),index)

    DO k=1,SIZE(tasks_recv,1)
       tasks_recv(k,1:ntasks_recv)=tasks_recv(k,index)
    ENDDO

    DO k=1,SIZE(rval,1)
       rval(k,1:ntasks_recv)=rval(k,index)
    ENDDO

    DEALLOCATE(taskid,STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","task_id",rs%group_size)
    DEALLOCATE(index,STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","index",rs%group_size)

    !
    ! final lists are ready 
    !

    tasks=>tasks_recv
    ntasks=ntasks_recv

    CALL timestop(handle)

  END SUBROUTINE distribute_tasks

! *****************************************************************************
SUBROUTINE get_atom_pair ( atom_pair, my_tasks, send, symmetric, natoms, maxset, maxpgf)

    INTEGER(kind=int_8), DIMENSION(:), &
      POINTER                                :: atom_pair
    INTEGER(kind=int_8), DIMENSION(:, :), &
      POINTER                                :: my_tasks
    LOGICAL                                  :: send, symmetric
    INTEGER                                  :: natoms, maxset, maxpgf

    INTEGER                                  :: acol, arow, i, iatom, ilevel, &
                                                ipgf, iset, j, jatom, jpgf, &
                                                jset, stat
    INTEGER(kind=int_8)                      :: natom8
    INTEGER, DIMENSION(:), POINTER           :: index

! calculate list of atom pairs
! fill pair list taking into acount symmetry

  natom8=natoms

  atom_pair = 0
  DO i = 1,SIZE(atom_pair)
     CALL int2pair(my_tasks(3,i),ilevel,iatom,jatom,iset,jset,ipgf,jpgf,natoms,maxset,maxpgf)
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
        atom_pair(i) = my_tasks(1,i)*natom8*natom8+(arow-1)*natom8 + (acol-1)
     ELSE
        atom_pair(i) = my_tasks(2,i)*natom8*natom8+(arow-1)*natom8 + (acol-1)
     ENDIF
  ENDDO

  ALLOCATE(INDEX(SIZE(atom_pair)),STAT=stat)
  IF (stat/=0) CALL stop_memory("get_atom_pair","index",SIZE(atom_pair))  

  ! find unique atom pairs that I'm sending/receiving
  CALL sort(atom_pair,SIZE(atom_pair),index)

  DEALLOCATE(index,STAT=stat)
  IF (stat/=0) CALL stop_memory("get_atom_pair","index",SIZE(atom_pair))  

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
SUBROUTINE distribute_matrix( rs, pmat, atom_pair_send, atom_pair_recv, natoms, scatter, error, hmat )

    TYPE(realspace_grid_type), POINTER       :: rs
    TYPE(real_matrix_type), POINTER          :: pmat
    INTEGER(KIND=int_8), DIMENSION(:), &
      POINTER                                :: atom_pair_send, atom_pair_recv
    INTEGER                                  :: natoms
    LOGICAL                                  :: scatter
    TYPE(cp_error_type), INTENT(inout)       :: error
    TYPE(real_matrix_type), OPTIONAL, &
      POINTER                                :: hmat

    CHARACTER(len=*), PARAMETER :: routineN = 'distribute_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: acol, arow, handle, i, j, k, &
                                                l, ncol, nrow, stat
    INTEGER(KIND=int_8)                      :: natom8, pair
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: recv_disps, recv_sizes, &
                                                send_disps, send_sizes
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: recv_buf_r, send_buf_r
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: h_block, p_block

  CALL timeset(routineN,handle)

  ! allocate local arrays
  ALLOCATE(send_sizes(rs%group_size),STAT=stat)
  IF (stat/=0) CALL stop_memory("get_my_tasks","send_sizes",rs%group_size)
  ALLOCATE(recv_sizes(rs%group_size),STAT=stat)
  IF (stat/=0) CALL stop_memory("get_my_tasks","recv_sizes",rs%group_size)
  ALLOCATE(send_disps(rs%group_size),STAT=stat)
  IF (stat/=0) CALL stop_memory("get_my_tasks","send_disps",rs%group_size)
  ALLOCATE(recv_disps(rs%group_size),STAT=stat)
  IF (stat/=0) CALL stop_memory("get_my_tasks","recv_disps",rs%group_size)

  ! set up send buffer sizes
  natom8=natoms

  send_sizes = 0
  DO i = 1, SIZE(atom_pair_send)

     ! proc we're sending this block to
     k = atom_pair_send(i) / natom8**2 + 1

     pair = MOD(atom_pair_send(i),natom8**2)

     arow = pair / natom8 + 1
     acol = MOD(pair, natom8) + 1

     nrow = pmat%last_row(arow) - pmat%first_row(arow) + 1
     ncol = pmat%last_col(acol) - pmat%first_col(acol) + 1       

     send_sizes(k) = send_sizes(k) + nrow * ncol        

  ENDDO

  send_disps = 0
  DO i = 2, rs % group_size
     send_disps(i) = send_disps(i-1) + send_sizes(i-1)
  ENDDO

  ALLOCATE(send_buf_r(SUM(send_sizes)),STAT=stat)
  IF (stat/=0) CALL stop_memory("get_my_tasks","send_buf_r",SUM(send_sizes))

  send_buf_r = 0

  ! do packing
  send_sizes=0
  DO i = 1, SIZE(atom_pair_send)

     l = atom_pair_send(i) / natom8**2 + 1

     pair = MOD(atom_pair_send(i),natom8**2)

     arow = pair / natom8 + 1
     acol = MOD(pair, natom8) + 1

     nrow = pmat%last_row(arow)-pmat%first_row(arow) + 1
     ncol = pmat%last_col(acol)-pmat%first_col(acol) + 1

     CALL get_block_node(matrix=pmat,&
          block_row=arow,&
          block_col=acol,&
          BLOCK=p_block)
     IF ( .NOT. ASSOCIATED ( p_block ) ) THEN 
        CALL stop_program ( "pack_matrix almost there", "Matrix block not found" )
     ENDIF

     DO k = 1, ncol
        DO j = 1, nrow
           send_buf_r(send_disps(l)+send_sizes(l)+j+(k-1)*nrow) = p_block(j,k)
        ENDDO
     ENDDO
     send_sizes(l)=send_sizes(l)+nrow*ncol
  ENDDO

  ! set up recv buffer

  recv_sizes = 0
  DO i = 1, SIZE(atom_pair_recv)

     ! proc we're receiving this data from
     k = atom_pair_recv(i) / natom8**2 + 1

     pair = MOD(atom_pair_recv(i),natom8**2)

     arow = pair / natom8 + 1
     acol = MOD(pair,natom8) + 1

     nrow = pmat%last_row(arow) - pmat%first_row(arow) + 1
     ncol = pmat%last_col(acol) - pmat%first_col(acol) + 1       

     recv_sizes(k) = recv_sizes(k) + nrow * ncol        

  ENDDO

  recv_disps = 0     
  DO i = 2, rs % group_size
     recv_disps(i) = recv_disps(i-1) + recv_sizes(i-1)
  ENDDO

  ALLOCATE(recv_buf_r(SUM(recv_sizes)),STAT=stat)
  IF (stat/=0) CALL stop_memory("get_my_tasks","recv_buf_r",SUM(recv_sizes))

  recv_buf_r = 0

  ! do communication
  CALL mp_alltoall(send_buf_r, send_sizes, send_disps,&
       recv_buf_r, recv_sizes, recv_disps, rs % group)

  !do unpacking
  recv_sizes=0
  DO i = 1, SIZE(atom_pair_recv)
     l = atom_pair_recv(i) / natom8**2 + 1
     pair = MOD(atom_pair_recv(i),natom8**2)

     arow = pair / natom8 + 1
     acol = MOD(pair, natom8) + 1

     nrow = pmat%last_row(arow)-pmat%first_row(arow) + 1
     ncol = pmat%last_col(acol)-pmat%first_col(acol) + 1

     CALL get_block_node(matrix=pmat,&
          block_row=arow,&
          block_col=acol,&
          BLOCK=p_block)

     IF ( PRESENT ( hmat ) ) THEN
     CALL get_block_node(matrix=hmat,&
          block_row=arow,&
          block_col=acol,&
          BLOCK=h_block)
     ENDIF

     IF ( ASSOCIATED ( p_block ) .AND. l .NE. rs%my_pos+1) THEN
        !IF( scatter ) THEN
        !   WRITE(6,*) arow,acol,l 
        !   CALL stop_program ( "unpack_matrix here and new", "Matrix block already present" )
        !ENDIF
     ELSE
        CALL add_block_node ( pmat, arow, acol, p_block ,error=error)
     ENDIF

     DO k = 1, ncol
        DO j = 1, nrow
           IF ( scatter ) THEN
              p_block(j,k) = recv_buf_r( recv_disps(l) + recv_sizes(l)+ j + (k-1)*nrow )
           ELSE
              h_block(j,k) = h_block(j,k) + recv_buf_r( recv_disps(l) + recv_sizes(l)+ j + (k-1)*nrow )
           ENDIF
        ENDDO
     ENDDO
     recv_sizes(l)=recv_sizes(l)+nrow*ncol
  ENDDO

  DEALLOCATE(send_buf_r,STAT=stat)
  IF (stat/=0) CALL stop_memory("get_my_tasks","send_buf_r",rs%group_size)
  DEALLOCATE(recv_buf_r,STAT=stat)
  IF (stat/=0) CALL stop_memory("get_my_tasks","recv_buf_r",rs%group_size)

  DEALLOCATE(send_sizes,STAT=stat)
  IF (stat/=0) CALL stop_memory("get_my_tasks","send_sizes",rs%group_size)
  DEALLOCATE(recv_sizes,STAT=stat)
  IF (stat/=0) CALL stop_memory("get_my_tasks","recv_sizes",rs%group_size)
  DEALLOCATE(send_disps,STAT=stat)
  IF (stat/=0) CALL stop_memory("get_my_tasks","send_disps",rs%group_size)
  DEALLOCATE(recv_disps,STAT=stat)
  IF (stat/=0) CALL stop_memory("get_my_tasks","recv_disps",rs%group_size)

  CALL timestop(handle)
  
END SUBROUTINE distribute_matrix

! *****************************************************************************
!> \brief determines the rank of the processor who's real rs grid contains point
!> \par History
!>      11.2007 created [MattW]
!>      10.2008 rewritten [Joost VandeVondele]
!> \author MattW
! *****************************************************************************
SUBROUTINE rs_find_node(rs_grid,cube_center,ntasks,tasks,lb_cube,ub_cube,added_tasks,error)

    TYPE(realspace_grid_type), POINTER       :: rs_grid
    INTEGER, DIMENSION(3), INTENT(IN)        :: cube_center
    INTEGER, INTENT(INOUT)                   :: ntasks
    INTEGER(kind=int_8), DIMENSION(:, :), &
      POINTER                                :: tasks
    INTEGER, DIMENSION(3), INTENT(IN)        :: lb_cube, ub_cube
    INTEGER, INTENT(OUT)                     :: added_tasks
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rs_find_node', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: add_tasks = 1000
    REAL(kind=dp), PARAMETER                 :: mult_tasks = 2.0_dp

    INTEGER :: bit_index, coord(3), curr_tasks, dest, i, icoord(3), idest, &
      itask, ix, iy, iz, lb_coord(3), lb_domain(3), lbc(3), ub_coord(3), &
      ub_domain(3), ubc(3)
    INTEGER(KIND=int_8)                      :: bit_pattern
    LOGICAL                                  :: dir_periodic(3), failure

    coord(1)=rs_grid % x2coord(cube_center(1))
    coord(2)=rs_grid % y2coord(cube_center(2))
    coord(3)=rs_grid % z2coord(cube_center(3))
    dest=rs_grid%coord2rank(coord(1),coord(2),coord(3))

    failure=.FALSE.

    ! the real cube coordinates
    lbc=lb_cube+cube_center
    ubc=ub_cube+cube_center

    IF ( ALL ( (rs_grid % lb_global(:,dest)-rs_grid%border) .LE. lbc )  .AND. &
         ALL ( (rs_grid % ub_global(:,dest)+rs_grid%border) .GE. ubc )) THEN
       !standard distributed collocation/integration
       tasks(1,ntasks)=dest
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
       !  IF (BTEST(bit_pattern,0)) rank=rs_grid_locate_rank(rs_grid,tasks(1,ntasks),(/-1,0,0/))
       !  IF (BTEST(bit_pattern,1)) rank=rs_grid_locate_rank(rs_grid,tasks(1,ntasks),(/+1,0,0/))
       !  IF (BTEST(bit_pattern,2)) rank=rs_grid_locate_rank(rs_grid,tasks(1,ntasks),(/0,-1,0/))
       !  IF (BTEST(bit_pattern,3)) rank=rs_grid_locate_rank(rs_grid,tasks(1,ntasks),(/0,+1,0/))
       !  IF (BTEST(bit_pattern,4)) rank=rs_grid_locate_rank(rs_grid,tasks(1,ntasks),(/0,0,-1/))
       !  IF (BTEST(bit_pattern,5)) rank=rs_grid_locate_rank(rs_grid,tasks(1,ntasks),(/0,0,+1/))
       bit_index=0
       bit_pattern=0
       DO i=1,3
          IF (rs_grid % perd ( i ) ==  1) THEN
             bit_pattern=IBCLR(bit_pattern,bit_index)
             bit_index=bit_index+1 
             bit_pattern=IBCLR(bit_pattern,bit_index)
             bit_index=bit_index+1 
          ELSE
             ! fits the left neighbor ?
             IF (ubc(i)<=rs_grid%lb_global(i,dest)-1+rs_grid%border) THEN
                bit_pattern=IBSET(bit_pattern,bit_index)
                bit_index=bit_index+1
             ELSE
                bit_pattern=IBCLR(bit_pattern,bit_index)
                bit_index=bit_index+1
             ENDIF
             ! fits the right neighbor ? 
             IF (lbc(i)>=rs_grid%ub_global(i,dest)+1-rs_grid%border) THEN
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
       !generalised collocation/integration needed
       ! first we figure out how many neighbors we have to add to include the lbc/ubc
       ! in the available domains (inclusive of halo points)
       ! first we 'ignore' periodic boundary conditions
       ! i.e. ub_coord-lb_coord+1 might be larger than group_dim
       lb_coord=coord
       ub_coord=coord
       lb_domain=rs_grid % lb_global(:,dest) - rs_grid % border
       ub_domain=rs_grid % ub_global(:,dest) + rs_grid % border
       DO i=1,3
          ! only if the grid is not periodic in this direction we need to take care of adding neighbors
          IF (rs_grid%perd(i)==0) THEN
             ! if the domain lower bound is greater than the lbc we need to add the size of the neighbor domain
             DO
               IF (lb_domain(i)>lbc(i)) THEN
                  lb_coord(i)=lb_coord(i)-1
                  icoord=MODULO(lb_coord,rs_grid%group_dim)
                  idest=rs_grid%coord2rank(icoord(1),icoord(2),icoord(3))
                  lb_domain(i)=lb_domain(i)-(rs_grid % ub_global(i,idest)-rs_grid % lb_global(i,idest)+1)
               ELSE
                  EXIT
               ENDIF
             ENDDO
             ! same for the upper bound
             DO
               IF (ub_domain(i)<ubc(i)) THEN
                  ub_coord(i)=ub_coord(i)+1
                  icoord=MODULO(ub_coord,rs_grid%group_dim)
                  idest=rs_grid%coord2rank(icoord(1),icoord(2),icoord(3))
                  ub_domain(i)=ub_domain(i)+(rs_grid % ub_global(i,idest)-rs_grid % lb_global(i,idest)+1)
               ELSE
                  EXIT
               ENDIF
             ENDDO
          ENDIF
       ENDDO

       ! some care is needed for the periodic boundaries
       DO i=1,3
          IF (ub_domain(i)-lb_domain(i)+1>=rs_grid%npts(i)) THEN
             dir_periodic(i)=.TRUE. 
             lb_coord(i)=0 
             ub_coord(i)=rs_grid%group_dim(i)-1
          ELSE
             dir_periodic(i)=.FALSE.
          ENDIF
       ENDDO

       added_tasks=PRODUCT(ub_coord-lb_coord+1)
       itask=ntasks
       ntasks=ntasks+added_tasks-1
       IF ( ntasks > SIZE(tasks,2)) THEN
          curr_tasks = (SIZE(tasks,2)+add_tasks)*mult_tasks
          CALL reallocate(tasks,1,6,1,curr_tasks)
       END IF
       DO iz=lb_coord(3),ub_coord(3)
       DO iy=lb_coord(2),ub_coord(2)
       DO ix=lb_coord(1),ub_coord(1)
          icoord=MODULO((/ix,iy,iz/),rs_grid%group_dim)
          idest =rs_grid%coord2rank(icoord(1),icoord(2),icoord(3))
          tasks(1,itask)=idest
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

END MODULE task_list_methods
