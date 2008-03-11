!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
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
  USE cube_utils,                      ONLY: cube_info_type,&
                                             return_cube,&
                                             return_cube_nonortho
  USE gaussian_gridlevels,             ONLY: gaussian_gridlevel,&
                                             gridlevel_info_type
  USE kinds,                           ONLY: dp,&
                                             int_8
  USE memory_utilities,                ONLY: reallocate
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
                                             rs_find_node
  USE realspace_task_selection,        ONLY: pair2int,&
                                             rs_get_my_tasks
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'task_list_methods'

  PUBLIC :: generate_qs_task_list

CONTAINS

! *****************************************************************************
!> \par History
!>      01.2008 factored out of calculate_rho_elec [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE generate_qs_task_list(qs_env, soft_valid, tasks, dist_ab, ntasks, atom_pair_send, atom_pair_recv, error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(IN), OPTIONAL            :: soft_valid
    INTEGER(kind=int_8), DIMENSION(:, :), &
      POINTER                                :: tasks
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: dist_ab
    INTEGER                                  :: ntasks
    INTEGER(kind=int_8), DIMENSION(:), &
      POINTER                                :: atom_pair_send, atom_pair_recv
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'generate_qs_task_list', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: add_tasks = 1000, &
                                                max_tasks = 2000
    REAL(kind=dp), PARAMETER                 :: mult_tasks = 2.0_dp

    INTEGER :: ab, bcol, brow, cost, curr_tasks, dest, dir, handle, iatom, &
      igrid_level, ikind, ilist, inode, ipgf, iset, jatom, jkind, jpgf, jset, &
      maxpgf, maxset, natoms, nkind, nlist, nnode, nseta, nsetb, tp(3)
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, lb_max, &
                                                lb_min, npgfa, npgfb
    LOGICAL                                  :: distributed_rs_grids, &
                                                failure, my_soft
    REAL(KIND=dp)                            :: dab, kind_radius_a, &
                                                kind_radius_b, rab2, zetp
    REAL(KIND=dp), DIMENSION(3)              :: ra, rab, rp
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

    CALL get_qs_env(qs_env=qs_env,&
         atomic_kind_set=atomic_kind_set,&
         cell=cell,&
         particle_set=particle_set,&
         dft_control=dft_control,&
         sab_orb=sab_orb,&
         pw_env=pw_env,error=error)
    CALL pw_env_get(pw_env, rs_pools=rs_pools, error=error)

    ! *** assign from pw_env
    gridlevel_info=>pw_env%gridlevel_info
    cube_info=>pw_env%cube_info

    ! *** figure out if we have distributed grids
    distributed_rs_grids=.FALSE.
    CALL rs_pools_create_rs_vect(rs_pools, rs_rho, error=error)
    DO igrid_level=1,gridlevel_info%ngrid_levels
       IF ( rs_rho(igrid_level)%rs_grid%distributed ) THEN
          distributed_rs_grids=.TRUE.
       ENDIF
    END DO

    ! find maximum numbers
    nkind = SIZE(atomic_kind_set)
    natoms = SIZE( particle_set )
    maxset=0
    maxpgf=0
    DO ikind=1,nkind
       atomic_kind => atomic_kind_set(ikind)

       CALL get_atomic_kind(atomic_kind=atomic_kind,&
            softb = my_soft, &
            orb_basis_set=orb_basis_set)

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
       CALL reallocate(tasks,1,5,1,max_tasks)
       CALL reallocate(dist_ab,1,3,1,max_tasks)
    ENDIF
    ntasks = 0
    curr_tasks = SIZE(tasks,2)

    loop_ikind: DO ikind=1,nkind

       atomic_kind => atomic_kind_set(ikind)

       CALL get_atomic_kind(atomic_kind=atomic_kind,&
            softb = my_soft, orb_basis_set=orb_basis_set)

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

          CALL get_atomic_kind(atomic_kind=atomic_kind,&
               softb = my_soft, orb_basis_set=orb_basis_set)

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

                IF (iatom <= jatom) THEN
                   brow = iatom
                   bcol = jatom
                ELSE
                   brow = jatom
                   bcol = iatom
                END IF

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

                            ntasks = ntasks + 1
                            IF ( ntasks > curr_tasks ) THEN
                               curr_tasks = (curr_tasks+add_tasks)*mult_tasks
                               CALL reallocate(tasks,1,5,1,curr_tasks)
                               CALL reallocate(dist_ab,1,3,1,curr_tasks)
                            END IF

                            DO dir = 1,3
                               rp(:) = ra(:) + zetb(jpgf,jset)/zetp*rab(:)
                               rp(:) = pbc(rp,cell)
                               tp(dir) = FLOOR(DOT_PRODUCT(cell%h_inv(dir,:),rp)*rs_rho(igrid_level)%rs_grid%npts(dir))
                               tp(dir) = MODULO ( tp(dir), rs_rho(igrid_level)%rs_grid%npts(dir) )
                               tp(dir) = tp(dir) + rs_rho(igrid_level)%rs_grid%lb(dir)
                            END DO

                            IF(rs_rho(igrid_level)%rs_grid%distributed) THEN
                               CALL rs_find_node(rs_rho(igrid_level)%rs_grid,tp,dest)
                               tasks(4,ntasks) = 1
                            ELSE
                               dest =  rs_rho(igrid_level)%rs_grid%my_pos
                               tasks(4,ntasks) = 0 
                            ENDIF

                            ! cost_factor for this task
                            CALL estimate_task_cost(cost,rs_rho(igrid_level)%rs_grid,cube_info(igrid_level),&
                                                    la_max(iset),zeta(ipgf,iset),la_min(iset),&
                                                    lb_max(jset),zetb(jpgf,jset),lb_min(jset),&
                                                    ra,rab,rab2,dft_control%qs_control%eps_rho_rspace)

                            tasks (1,ntasks) = dest
                            tasks (2,ntasks) = rs_rho(igrid_level)% rs_grid % my_pos
                            tasks (5,ntasks) = cost

                            !encode the atom pairs and basis info as a single long integer
                            CALL pair2int(tasks(3,ntasks),igrid_level,iatom,jatom,iset,jset,ipgf,jpgf,natoms,maxset,maxpgf)

                            dist_ab (1,ntasks) = rab(1)
                            dist_ab (2,ntasks) = rab(2)
                            dist_ab (3,ntasks) = rab(3)

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

    ! redistribute the task list so that all tasks map on the local rs grids
    CALL rs_get_my_tasks (rs_rho(1)%rs_grid,distributed_rs_grids,ntasks, natoms,&
            maxset,maxpgf, tasks, rval=dist_ab, atom_pair_send=atom_pair_send,&
            atom_pair_recv=atom_pair_recv, error=error)

    CALL rs_pools_give_back_rs_vect(rs_pools, rs_rho, error=error)

    CALL timestop(handle)

  END SUBROUTINE generate_qs_task_list

! *****************************************************************************
!> \brief tries to predict the time needed for processing this task
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
!>      01.2008 created [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE estimate_task_cost(cost,rsgrid,cube_info,la_max,zeta,la_min,lb_max,zetb,lb_min,ra,rab,rab2,eps)

    INTEGER, INTENT(OUT)                     :: cost
    TYPE(realspace_grid_type), POINTER       :: rsgrid
    TYPE(cube_info_type), INTENT(IN)         :: cube_info
    INTEGER, INTENT(IN)                      :: la_max
    REAL(KIND=dp), INTENT(IN)                :: zeta
    INTEGER, INTENT(IN)                      :: la_min, lb_max
    REAL(KIND=dp), INTENT(IN)                :: zetb
    INTEGER, INTENT(IN)                      :: lb_min
    REAL(KIND=dp), INTENT(IN)                :: ra(3), rab(3), rab2, eps

    INTEGER                                  :: cmax, lb_cube(3), ub_cube(3)
    INTEGER, DIMENSION(:), POINTER           :: sphere_bounds
    REAL(KIND=dp)                            :: cutoff, f, prefactor, radius, &
                                                rap(3), rb(3), rbp(3), rp(3), &
                                                zetp

    zetp = zeta + zetb
    f    = zetb/zetp
    rap(:)    = f*rab(:)
    rbp(:)    = rap(:) - rab(:)
    rp(:)     = ra(:) + rap(:)
    rb(:)     = ra(:)+rab(:)
    cutoff    = 1.0_dp
    prefactor = EXP(-zeta*f*rab2)
    radius=exp_radius_very_extended(la_min,la_max,lb_min,lb_max,ra=ra,rb=rb,rp=rp,&
               zetp=zetp,eps=eps,&
               prefactor=prefactor,cutoff=cutoff)

    IF (radius .LE. 0.0_dp) THEN
       ! nevertheless a small cost (this is rare, but might need tuning)
       cost=20
    ELSE
       ! a bit overkill, we only need cmax
       IF (rsgrid%orthorhombic ) THEN
          CALL return_cube(cube_info,radius,lb_cube,ub_cube,sphere_bounds)
       ELSE
          ! untested so far
          CALL return_cube_nonortho(cube_info,radius,lb_cube,ub_cube,rp)
       ENDIF
       cmax=MAXVAL(((ub_cube-lb_cube)+1)/2)

       cost=CEILING((REAL(la_max+lb_max+1+2,KIND=dp)*(cmax+2)**3*7.0_dp+30000.0_dp)/1000.0_dp)
    ENDIF

  END SUBROUTINE estimate_task_cost

END MODULE task_list_methods
