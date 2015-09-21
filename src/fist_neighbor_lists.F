!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Generate the atomic neighbor lists for FIST.
!> \par History
!>      - build and update merged (11.02.2005,MK)
!>      - bug fix for PERIODIC NONE (24.02.06,MK)
!>      - Major rewriting (light memory neighbor lists): teo and joost (05.2006)
!>      - Completely new algorithm for the neighbor lists
!>        (faster and memory lighter) (Teo 08.2006)
!> \author MK (19.11.2002,24.07.2003)
!>      Teodoro Laino (08.2006) - MAJOR REWRITING
! *****************************************************************************
MODULE fist_neighbor_lists

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE cell_types,                      ONLY: cell_type,&
                                             get_cell,&
                                             pbc,&
                                             plane_distance,&
                                             scaled_to_real
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_units,                        ONLY: cp_unit_from_cp2k
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE exclusion_types,                 ONLY: exclusion_type
  USE fist_neighbor_list_types,        ONLY: fist_neighbor_add,&
                                             fist_neighbor_deallocate,&
                                             fist_neighbor_init,&
                                             fist_neighbor_type,&
                                             neighbor_kind_pairs_type
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE mathlib,                         ONLY: matvec_3x3
  USE memory_utilities,                ONLY: reallocate
  USE particle_types,                  ONLY: particle_type
  USE qmmm_ff_fist,                    ONLY: qmmm_ff_precond_only_qm
  USE string_utilities,                ONLY: compress
  USE subcell_types,                   ONLY: allocate_subcell,&
                                             deallocate_subcell,&
                                             give_ijk_subcell,&
                                             reorder_atoms_subcell,&
                                             subcell_type
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: sort
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  ! Global parameters (in this module)
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'fist_neighbor_lists'

  TYPE local_atoms_type
    INTEGER, DIMENSION(:), POINTER                   :: list,&
                                                        list_local_a_index
  END TYPE local_atoms_type

  ! Public subroutines
  PUBLIC :: build_fist_neighbor_lists

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param atomic_kind_set ...
!> \param particle_set ...
!> \param local_particles ...
!> \param cell ...
!> \param r_max ...
!> \param r_minsq ...
!> \param ei_scale14 ...
!> \param vdw_scale14 ...
!> \param nonbonded ...
!> \param para_env ...
!> \param build_from_scratch ...
!> \param geo_check ...
!> \param mm_section ...
!> \param full_nl ...
!> \param exclusions ...
!> \par History
!>      08.2006 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE build_fist_neighbor_lists(atomic_kind_set, particle_set, &
        local_particles, cell, r_max, r_minsq, ei_scale14, vdw_scale14, &
        nonbonded, para_env, build_from_scratch, geo_check, mm_section, &
        full_nl, exclusions)

    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(distribution_1d_type), OPTIONAL, &
      POINTER                                :: local_particles
    TYPE(cell_type), POINTER                 :: cell
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: r_max, r_minsq
    REAL(KIND=DP), INTENT(IN)                :: ei_scale14, vdw_scale14
    TYPE(fist_neighbor_type), POINTER        :: nonbonded
    TYPE(cp_para_env_type), POINTER          :: para_env
    LOGICAL, INTENT(IN)                      :: build_from_scratch, geo_check
    TYPE(section_vals_type), POINTER         :: mm_section
    LOGICAL, DIMENSION(:, :), OPTIONAL, &
      POINTER                                :: full_nl
    TYPE(exclusion_type), DIMENSION(:), &
      OPTIONAL                               :: exclusions

    CHARACTER(LEN=*), PARAMETER :: routineN = 'build_fist_neighbor_lists', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: kind_name, print_key_path, &
                                                unit_str
    INTEGER                                  :: atom_a, handle, iatom_local, &
                                                ikind, iw, maxatom, &
                                                natom_local_a, nkind, &
                                                output_unit
    LOGICAL                                  :: failure, &
                                                present_local_particles, &
                                                print_subcell_grid
    LOGICAL, DIMENSION(:), POINTER           :: skip_kind
    LOGICAL, DIMENSION(:, :), POINTER        :: my_full_nl
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(local_atoms_type), ALLOCATABLE, &
      DIMENSION(:)                           :: atom

    failure = .FALSE.
    CALL timeset(routineN,handle)
    NULLIFY(logger)
    logger => cp_get_default_logger()

    print_subcell_grid = .FALSE.
    output_unit = cp_print_key_unit_nr(logger,mm_section,"PRINT%SUBCELL",&
         extension=".Log")
    IF (output_unit > 0)  print_subcell_grid = .TRUE.

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
         maxatom=maxatom)

    present_local_particles = PRESENT(local_particles)

    ! if exclusions matters local particles are present. Seems like only the exclusions
    ! for the local particles are needed, which would imply a huge memory savings for fist
    IF (PRESENT(exclusions)) THEN
      CPASSERT(present_local_particles)
    ENDIF

    ! Allocate work storage
    nkind = SIZE(atomic_kind_set)
    ALLOCATE (atom(nkind))
    ALLOCATE (skip_kind(nkind))
    ! full_nl
    IF (PRESENT(full_nl)) THEN
       my_full_nl => full_nl
    ELSE
       ALLOCATE(my_full_nl(nkind,nkind))
       my_full_nl = .FALSE.
    END IF
    ! Initialize the local data structures
    DO ikind=1,nkind
       atomic_kind => atomic_kind_set(ikind)
       NULLIFY (atom(ikind)%list)
       NULLIFY (atom(ikind)%list_local_a_index)

       CALL get_atomic_kind(atomic_kind=atomic_kind,&
            atom_list=atom(ikind)%list,name=kind_name)
       skip_kind(ikind)=qmmm_ff_precond_only_qm(kind_name)
       IF(present_local_particles) THEN
          natom_local_a = local_particles%n_el(ikind)
       ELSE
          natom_local_a = SIZE(atom(ikind)%list)
       END IF
       IF (natom_local_a > 0) THEN
          ALLOCATE (atom(ikind)%list_local_a_index(natom_local_a))
          ! Build index vector for mapping
          DO iatom_local=1,natom_local_a
             IF (present_local_particles) THEN
                atom_a = local_particles%list(ikind)%array(iatom_local)
             ELSE
                atom_a = atom(ikind)%list(iatom_local)
             END IF
             atom(ikind)%list_local_a_index(iatom_local) = atom_a
          END DO

       END IF
    END DO

    IF (build_from_scratch) THEN
       IF (ASSOCIATED(nonbonded)) THEN
          CALL fist_neighbor_deallocate(nonbonded)
       END IF
    END IF

    ! Build the nonbonded neighbor lists
    CALL build_neighbor_lists(nonbonded, particle_set, atom, cell, &
         print_subcell_grid, output_unit, r_max, r_minsq, &
         ei_scale14, vdw_scale14, geo_check,"NONBONDED", skip_kind, &
         my_full_nl, exclusions)

    ! Sort the list according kinds for each cell
    CALL sort_neighbor_lists(nonbonded, nkind)

    print_key_path = "PRINT%NEIGHBOR_LISTS"

    IF (BTEST(cp_print_key_should_output(logger%iter_info,mm_section,print_key_path),&
         cp_p_file)) THEN
       iw = cp_print_key_unit_nr(logger=logger,&
            basis_section=mm_section,&
            print_key_path=print_key_path,&
            extension=".out",&
            middle_name="nonbonded_nl",&
            local=.TRUE.,&
            log_filename=.FALSE.,&
            file_position="REWIND")
       CALL section_vals_val_get(mm_section,TRIM(print_key_path)//"%UNIT",c_val=unit_str)
       CALL write_neighbor_lists(nonbonded,particle_set,cell,para_env,iw,&
            "NONBONDED NEIGHBOR LISTS",unit_str)
       CALL cp_print_key_finished_output(unit_nr=iw,&
            logger=logger,&
            basis_section=mm_section,&
            print_key_path=print_key_path,&
            local=.TRUE.)
    END IF

    ! Release work storage
    DO ikind=1,nkind
       NULLIFY (atom(ikind)%list)
       IF (ASSOCIATED(atom(ikind)%list_local_a_index)) THEN
          DEALLOCATE (atom(ikind)%list_local_a_index)
       END IF
    END DO
    IF (PRESENT(full_nl)) THEN
       NULLIFY(my_full_nl)
    ELSE
       DEALLOCATE(my_full_nl)
    END IF
    DEALLOCATE (atom)

    DEALLOCATE (skip_kind)

    CALL cp_print_key_finished_output(unit_nr=output_unit,&
         logger=logger,&
         basis_section=mm_section,&
         print_key_path="PRINT%SUBCELL")
    CALL timestop(handle)

  END SUBROUTINE build_fist_neighbor_lists

! *****************************************************************************
!> \brief ...
!> \param nonbonded ...
!> \param particle_set ...
!> \param atom ...
!> \param cell ...
!> \param print_subcell_grid ...
!> \param output_unit ...
!> \param r_max ...
!> \param r_minsq ...
!> \param ei_scale14 ...
!> \param vdw_scale14 ...
!> \param geo_check ...
!> \param name ...
!> \param skip_kind ...
!> \param full_nl ...
!> \param exclusions ...
!> \par History
!>      08.2006 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE build_neighbor_lists(nonbonded, particle_set, atom,cell, &
        print_subcell_grid, output_unit, r_max, r_minsq, &
        ei_scale14, vdw_scale14, geo_check, name, skip_kind, full_nl, exclusions)

    TYPE(fist_neighbor_type), POINTER        :: nonbonded
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(local_atoms_type), DIMENSION(:), &
      INTENT(IN)                             :: atom
    TYPE(cell_type), POINTER                 :: cell
    LOGICAL, INTENT(IN)                      :: print_subcell_grid
    INTEGER, INTENT(IN)                      :: output_unit
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: r_max, r_minsq
    REAL(KIND=dp), INTENT(IN)                :: ei_scale14, vdw_scale14
    LOGICAL, INTENT(IN)                      :: geo_check
    CHARACTER(LEN=*), INTENT(IN)             :: name
    LOGICAL, DIMENSION(:), POINTER           :: skip_kind
    LOGICAL, DIMENSION(:, :), POINTER        :: full_nl
    TYPE(exclusion_type), DIMENSION(:), &
      OPTIONAL                               :: exclusions

    CHARACTER(LEN=*), PARAMETER :: routineN = 'build_neighbor_lists', &
      routineP = moduleN//':'//routineN

    INTEGER :: a_i, a_j, a_k, atom_a, atom_b, b_i, b_j, b_k, b_pi, b_pj, &
      b_pk, bg_i, bg_j, bg_k, handle, i, i1, iatom_local, icell, icellmap, &
      id_kind, ii, ii_start, ij, ij_start, ik, ik_start, ikind, imap, &
      imax_cell, invcellmap, iw, ix, j, j1, jatom_local, jcell, jkind, jx, k, &
      kcell, kx, natom_local_a, ncellmax, nkind, nkind00, tmpdim, xdim, ydim, &
      zdim
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: kind_of, work
    INTEGER, ALLOCATABLE, DIMENSION(:, :, :) :: cellmap
    INTEGER, DIMENSION(3)                    :: isubcell, ncell, nsubcell, &
                                                periodic
    LOGICAL                                  :: any_full, atom_order, &
                                                check_spline, failure, &
                                                is_full, subcell000
    LOGICAL, ALLOCATABLE, DIMENSION(:, :, :) :: sphcub
    REAL(dp)                                 :: rab2, rab2_max, rab2_min, &
                                                rab_max
    REAL(dp), DIMENSION(3)                   :: abc, cell_v, cv_b, rab, rb, &
                                                sab_max
    REAL(KIND=dp)                            :: ic(3), icx(3), radius, vv
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: coord
    TYPE(neighbor_kind_pairs_type), POINTER  :: inv_neighbor_kind_pair, &
                                                neighbor_kind_pair
    TYPE(subcell_type), DIMENSION(:, :, :), &
      POINTER                                :: subcell_a, subcell_b

    CALL timeset(routineN,handle)
    failure = .FALSE.
    nkind = SIZE(atom)
    nsubcell = 1
    isubcell = 0
    ncell    = 0
    any_full = ANY(full_nl)
    CALL get_cell(cell=cell,&
                  abc=abc,&
                  periodic=periodic)
    ! Determines the number of subcells
    DO ikind=1,nkind
       DO jkind=ikind,nkind
          ! Calculate the square of the maximum interaction distance
          rab_max    = r_max(ikind,jkind)
          IF (skip_kind(ikind).AND.skip_kind(jkind)) CYCLE
          nsubcell(1) = MAX(nsubcell(1),CEILING(plane_distance(1,0,0,cell)/rab_max))
          nsubcell(2) = MAX(nsubcell(2),CEILING(plane_distance(0,1,0,cell)/rab_max))
          nsubcell(3) = MAX(nsubcell(3),CEILING(plane_distance(0,0,1,cell)/rab_max))
       END DO
    END DO
    ! Determines the number of periodic images and the number of interacting subcells
    DO ikind=1,nkind
       DO jkind=ikind,nkind
          IF (skip_kind(ikind).AND.skip_kind(jkind)) CYCLE
          ! Calculate the square of the maximum interaction distance
          rab_max    = r_max(ikind,jkind)
          sab_max(1) = rab_max/plane_distance(1,0,0,cell)
          sab_max(2) = rab_max/plane_distance(0,1,0,cell)
          sab_max(3) = rab_max/plane_distance(0,0,1,cell)
          ncell= MAX(ncell(:),CEILING(sab_max(:)*periodic(:)))
          isubcell=MAX(isubcell(:),CEILING(sab_max(:)*REAL(nsubcell(:),KIND=dp)))
       END DO
    END DO
    CALL fist_neighbor_init(nonbonded, ncell)
    ! Print headline
    IF (print_subcell_grid) THEN
       WRITE (UNIT=output_unit,FMT="(/,/,T2,A,/)")&
            "SUBCELL GRID  INFO FOR THE "//TRIM(name)//" NEIGHBOR LISTS"
       WRITE (UNIT=output_unit,FMT="(T4,A,10X,3I10)")" NUMBER OF SUBCELLS             ::",nsubcell
       WRITE (UNIT=output_unit,FMT="(T4,A,10X,3I10)")" NUMBER OF PERIODIC      IMAGES ::",ncell
       WRITE (UNIT=output_unit,FMT="(T4,A,10X,3I10)")" NUMBER OF INTERACTING SUBCELLS ::",isubcell
    END IF
    ! Allocate subcells
    CALL allocate_subcell(subcell_a,nsubcell,cell=cell)
    CALL allocate_subcell(subcell_b,nsubcell,cell=cell)
    ! Let's map the sequence of the periodic images
    ncellmax = MAXVAL(ncell)
    ALLOCATE(cellmap(-ncellmax:ncellmax,-ncellmax:ncellmax,-ncellmax:ncellmax))
    cellmap=-1
    imap=0
    nkind00 = nkind*(nkind+1)/2
    DO imax_cell=0,ncellmax
       DO kcell=-imax_cell,imax_cell
          DO jcell=-imax_cell,imax_cell
             DO icell=-imax_cell,imax_cell
                IF (cellmap(icell,jcell,kcell)==-1) THEN
                   imap=imap+1
                   cellmap(icell,jcell,kcell)=imap
                   CPASSERT(imap<=nonbonded%nlists)
                   neighbor_kind_pair => nonbonded%neighbor_kind_pairs(imap)

                   neighbor_kind_pair%cell_vector(1) = icell
                   neighbor_kind_pair%cell_vector(2) = jcell
                   neighbor_kind_pair%cell_vector(3) = kcell
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    ! Mapping the spherical interaction between subcells
    ALLOCATE(sphcub(-isubcell(1):isubcell(1),&
                    -isubcell(2):isubcell(2),&
                    -isubcell(3):isubcell(3)))
    sphcub = .FALSE.
    IF (ALL(isubcell/=0)) THEN
       radius = REAL(isubcell(1),KIND=dp)**2+ REAL(isubcell(2),KIND=dp)**2+&
                REAL(isubcell(3),KIND=dp)**2
       loop1: DO k = -isubcell(3),isubcell(3)
          loop2: DO j = -isubcell(2),isubcell(2)
             loop3: DO i = -isubcell(1),isubcell(1)
                ic = REAL((/i,j,k/),KIND=dp)
                ! subcell cube vertex
                DO kx = -1,1
                   icx(3) = ic(3) + SIGN(0.5_dp,REAL(kx,KIND=dp))
                   DO jx = -1,1
                      icx(2) = ic(2) + SIGN(0.5_dp,REAL(jx,KIND=dp))
                      DO ix = -1,1
                         icx(1) = ic(1) + SIGN(0.5_dp,REAL(ix,KIND=dp))
                         vv=icx(1)*icx(1)+icx(2)*icx(2)+icx(3)*icx(3)
                         vv=vv/radius
                         IF (vv <= 1.0_dp) THEN
                            sphcub(i,j,k) = .TRUE.
                            CYCLE loop3
                         END IF
                      END DO
                   END DO
                END DO
             END DO loop3
          END DO loop2
       END DO loop1
    END IF
    ! Mapping locally all atoms in the zeroth cell
    ALLOCATE(coord(3,SIZE(particle_set)))
    DO atom_a = 1, SIZE(particle_set)
       coord(:,atom_a) = pbc(particle_set(atom_a)%r,cell)
    END DO
    ! Associate particles to subcells (local particles)
    DO ikind=1,nkind
       IF (.NOT.ASSOCIATED(atom(ikind)%list_local_a_index)) CYCLE
       natom_local_a = SIZE(atom(ikind)%list_local_a_index)
       DO iatom_local=1,natom_local_a
          atom_a = atom(ikind)%list_local_a_index(iatom_local)
          CALL give_ijk_subcell(coord(:,atom_a),i,j,k,cell,nsubcell)
          subcell_a(i,j,k)%natom   = subcell_a(i,j,k)%natom + 1
       END DO
    END DO
    DO k = 1,nsubcell(3)
       DO j = 1,nsubcell(2)
          DO i = 1,nsubcell(1)
             ALLOCATE(subcell_a(i,j,k)%atom_list(subcell_a(i,j,k)%natom))
             subcell_a(i,j,k)%natom   = 0
          END DO
       END DO
    END DO
    DO ikind=1,nkind
       IF (.NOT.ASSOCIATED(atom(ikind)%list_local_a_index)) CYCLE
       natom_local_a = SIZE(atom(ikind)%list_local_a_index)
       DO iatom_local=1,natom_local_a
          atom_a = atom(ikind)%list_local_a_index(iatom_local)
          CALL give_ijk_subcell(coord(:,atom_a),i,j,k,cell,nsubcell)
          subcell_a(i,j,k)%natom   = subcell_a(i,j,k)%natom + 1
          subcell_a(i,j,k)%atom_list(subcell_a(i,j,k)%natom) = atom_a
       END DO
    END DO
    ! Associate particles to subcells (distributed particles)
    DO atom_b = 1, SIZE(particle_set)
       CALL give_ijk_subcell(coord(:,atom_b),i,j,k,cell,nsubcell)
       subcell_b(i,j,k)%natom   = subcell_b(i,j,k)%natom + 1
    END DO
    DO k = 1,nsubcell(3)
       DO j = 1,nsubcell(2)
          DO i = 1,nsubcell(1)
             ALLOCATE(subcell_b(i,j,k)%atom_list(subcell_b(i,j,k)%natom))
             subcell_b(i,j,k)%natom   = 0
          END DO
       END DO
    END DO
    DO atom_b = 1, SIZE(particle_set)
       CALL give_ijk_subcell(coord(:,atom_b),i,j,k,cell,nsubcell)
       subcell_b(i,j,k)%natom   = subcell_b(i,j,k)%natom + 1
       subcell_b(i,j,k)%atom_list(subcell_b(i,j,k)%natom) = atom_b
    END DO
    ! Reorder atoms associated to subcells
    tmpdim = MAXVAL(subcell_a(:,:,:)%natom)
    tmpdim = MAX(tmpdim,MAXVAL(subcell_b(:,:,:)%natom))
    ALLOCATE(work(3*tmpdim))
    ALLOCATE(kind_of(SIZE(particle_set)))
    DO i = 1, SIZE(particle_set)
       kind_of(i) = particle_set(i)%atomic_kind%kind_number
    END DO
    DO k = 1,nsubcell(3)
       DO j = 1,nsubcell(2)
          DO i = 1,nsubcell(1)
             CALL reorder_atoms_subcell(subcell_a(i,j,k)%atom_list,kind_of,work)
             CALL reorder_atoms_subcell(subcell_b(i,j,k)%atom_list,kind_of,work)
          END DO
       END DO
    END DO
    DEALLOCATE(work,kind_of)
    zdim  = nsubcell(3)
    ydim  = nsubcell(2)
    xdim  = nsubcell(1)
    is_full = .FALSE.
    ! We can skip until ik>=0.. this prescreens the order of the subcells
    ik_start = -isubcell(3)
    IF (.NOT.any_full) ik_start = 0
    ! Loop over first subcell
    loop_a_k: DO a_k=1,nsubcell(3)
    loop_a_j: DO a_j=1,nsubcell(2)
    loop_a_i: DO a_i=1,nsubcell(1)
       IF (subcell_a(a_i,a_j,a_k)%natom == 0) CYCLE
       ! Loop over second subcell
       loop_b_k: DO ik=ik_start,isubcell(3)
          bg_k = a_k+ik
          b_k  = MOD(bg_k,zdim)
          b_pk = bg_k/zdim
          IF (b_k<=0) THEN
             b_k  = zdim+b_k
             b_pk = b_pk - 1
          END IF
          IF ((periodic(3)==0).AND.(ABS(b_pk)>0)) CYCLE
          ! Setup the starting point.. this prescreens the order of the subcells
          ij_start = -isubcell(2)
          IF ((ik==0).AND.(ik_start==0)) ij_start = 0
       loop_b_j: DO ij=ij_start,isubcell(2)
          bg_j = a_j+ij
          b_j  = MOD(bg_j,ydim)
          b_pj = bg_j/ydim
          IF (b_j<=0) THEN
             b_j  = ydim+b_j
             b_pj = b_pj - 1
          END IF
          IF ((periodic(2)==0).AND.(ABS(b_pj)>0)) CYCLE
          ! Setup the starting point.. this prescreens the order of the subcells
          ii_start = -isubcell(1)
          IF ((ij==0).AND.(ij_start==0)) ii_start = 0
       loop_b_i: DO ii=ii_start,isubcell(1)
          ! Ellipsoidal screening of subcells
          IF (.NOT.sphcub(ii,ij,ik)) CYCLE
          bg_i  = a_i+ii
          b_i  = MOD(bg_i,xdim)
          b_pi = bg_i/xdim
          IF (b_i<=0) THEN
             b_i  = xdim+b_i
             b_pi = b_pi - 1
          END IF
          IF ((periodic(1)==0).AND.(ABS(b_pi)>0)) CYCLE
          IF ( subcell_b(b_i,b_j,b_k)%natom == 0) CYCLE
          ! Find the proper neighbor kind pair
          icellmap = cellmap(b_pi,b_pj,b_pk)
          neighbor_kind_pair => nonbonded%neighbor_kind_pairs(icellmap)
          ! Find the replica vector
          cell_v = 0.0_dp
          IF ((b_pi/=0).OR.(b_pj/=0).OR.(b_pk/=0)) THEN
             cv_b(1)=b_pi ; cv_b(2)=b_pj ; cv_b(3)=b_pk
             CALL scaled_to_real(cell_v,cv_b,cell)
          END IF
          subcell000 = (a_k==bg_k).AND.(a_j==bg_j).AND.(a_i==bg_i)
          ! Loop over particles inside subcell_a and subcell_b
          DO jatom_local = 1, subcell_b(b_i,b_j,b_k)%natom
             atom_b = subcell_b(b_i,b_j,b_k)%atom_list(jatom_local)
             jkind  = particle_set(atom_b)%atomic_kind%kind_number
             rb(1) = coord(1,atom_b)+cell_v(1)
             rb(2) = coord(2,atom_b)+cell_v(2)
             rb(3) = coord(3,atom_b)+cell_v(3)
             DO iatom_local = 1, subcell_a(a_i,a_j,a_k)%natom
                atom_a = subcell_a(a_i,a_j,a_k)%atom_list(iatom_local)
                ikind  = particle_set(atom_a)%atomic_kind%kind_number
                ! Screen interaction to avoid double counting
                atom_order = (atom_a<=atom_b)
                ! Special case for kind combination requiring the full NL
                IF (any_full) THEN
                   is_full = full_nl(ikind,jkind)
                   IF (is_full) THEN
                      atom_order = (atom_a==atom_b)
                   ELSE
                      IF (ik<0)           CYCLE
                      IF (ik==0.AND.ij<0) CYCLE
                      IF (ij==0.AND.ii<0) CYCLE
                   END IF
                END IF
                IF (subcell000.AND.atom_order) CYCLE
                rab(1) = rb(1) - coord(1,atom_a)
                rab(2) = rb(2) - coord(2,atom_a)
                rab(3) = rb(3) - coord(3,atom_a)
                rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
                rab_max  = r_max(ikind,jkind)
                rab2_max = rab_max*rab_max
                IF (rab2 < rab2_max) THEN
                   ! Diagonal storage
                   j1  = MIN(ikind,jkind)
                   i1  = MAX(ikind,jkind)-j1+1
                   j1  = nkind - j1 + 1
                   id_kind = nkind00-(j1*(j1+1)/2)+i1
                   ! Store the pair
                   CALL fist_neighbor_add(neighbor_kind_pair, atom_a, atom_b, &
                        rab=rab, &
                        check_spline=check_spline, id_kind=id_kind, &
                        skip=(skip_kind(ikind).AND.skip_kind(jkind)), &
                        cell=cell, ei_scale14=ei_scale14, &
                        vdw_scale14=vdw_scale14, exclusions=exclusions)
                   ! This is to handle properly when interaction radius is larger than cell size
                   IF ((atom_a==atom_b).AND.(ik_start==0)) THEN
                      invcellmap = cellmap(-b_pi,-b_pj,-b_pk)
                      inv_neighbor_kind_pair => nonbonded%neighbor_kind_pairs(invcellmap)
                      rab = rab - 2.0_dp * cell_v
                      CALL fist_neighbor_add(inv_neighbor_kind_pair, atom_a, atom_b, &
                           rab=rab, &
                           check_spline=check_spline, id_kind=id_kind, &
                           skip=(skip_kind(ikind).AND.skip_kind(jkind)), &
                           cell=cell, ei_scale14=ei_scale14, &
                           vdw_scale14=vdw_scale14, exclusions=exclusions)
                   END IF
                   ! Check for too close hits
                   IF (check_spline) THEN
                      rab2_min = r_minsq(ikind,jkind)
                      IF (rab2 < rab2_min) THEN
                         iw=cp_logger_get_default_unit_nr()
                         WRITE(iw,'(T2,A,2I7,2(A,F15.8),A)')"WARNING| Particles: ",&
                              atom_a,atom_b,&
                              " at distance [au]:",SQRT(rab2)," less than: ",&
                              SQRT(rab2_min),&
                              "; increase EMAX_SPLINE."
                         IF (rab2 < rab2_min/(1.06_dp)**2) THEN
                            IF (geo_check) THEN
                               CALL stop_program(routineN,moduleN,__LINE__,&
                                                 "GEOMETRY wrong or EMAX_SPLINE too small!")
                            END IF
                         END IF
                      END IF
                   END IF
                END IF
             END DO
          END DO
       END DO loop_b_i
       END DO loop_b_j
       END DO loop_b_k
    END DO loop_a_i
    END DO loop_a_j
    END DO loop_a_k
    DEALLOCATE(coord)
    DEALLOCATE(cellmap)
    DEALLOCATE(sphcub)
    CALL deallocate_subcell(subcell_a)
    CALL deallocate_subcell(subcell_b)

    CALL timestop(handle)
  END SUBROUTINE build_neighbor_lists

! *****************************************************************************
!> \brief Write a set of neighbor lists to the output unit.
!> \param nonbonded ...
!> \param particle_set ...
!> \param cell ...
!> \param para_env ...
!> \param output_unit ...
!> \param name ...
!> \param unit_str ...
!> \par History
!>      08.2006 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE write_neighbor_lists(nonbonded,particle_set,cell,para_env,output_unit,&
    name,unit_str)

    TYPE(fist_neighbor_type), POINTER        :: nonbonded
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER, INTENT(IN)                      :: output_unit
    CHARACTER(LEN=*), INTENT(IN)             :: name, unit_str

    CHARACTER(LEN=default_string_length)     :: string
    INTEGER                                  :: atom_a, atom_b, iab, ilist, &
                                                mype, nneighbor
    LOGICAL                                  :: print_headline
    REAL(dp)                                 :: conv, dab
    REAL(dp), DIMENSION(3)                   :: cell_v, ra, rab, rb
    TYPE(neighbor_kind_pairs_type), POINTER  :: neighbor_kind_pair

    mype = para_env%mepos
    ! Print headline
    string = ""
    WRITE (UNIT=string,FMT="(A,I5,A)")&
      TRIM(name)//" IN "//TRIM(unit_str)//" (PROCESS",mype,")"
    CALL compress(string)
    IF (output_unit > 0) WRITE (UNIT=output_unit,FMT="(/,/,T2,A)") TRIM(string)

    print_headline = .TRUE.
    nneighbor = 0
    conv = cp_unit_from_cp2k(1.0_dp,TRIM(unit_str))
    DO iab=1,SIZE(nonbonded%neighbor_kind_pairs)
       neighbor_kind_pair => nonbonded%neighbor_kind_pairs(iab)
       CALL matvec_3x3(cell_v, cell%hmat,REAL(neighbor_kind_pair%cell_vector,KIND=dp))
       DO ilist=1,neighbor_kind_pair%npairs
          nneighbor = nneighbor + 1
          IF (output_unit > 0) THEN
             ! Print second part of the headline
             atom_a = neighbor_kind_pair%list(1,ilist)
             atom_b = neighbor_kind_pair%list(2,ilist)
             IF (print_headline) THEN
                WRITE (UNIT=output_unit,FMT="(T3,2(A6,3(5X,A,5X)),1X,A11,10X,A8,A5,A10,A9)")&
                     "Atom-A","X","Y","Z","Atom-B","X","Y","Z","Cell(i,j,k)",&
                     "Distance","ONFO","VDW-scale","EI-scale"
                print_headline = .FALSE.
             END IF

             ra(:) = pbc(particle_set(atom_a)%r,cell)
             rb(:) = pbc(particle_set(atom_b)%r,cell)
             rab   = rb(:)-ra(:)+cell_v
             dab   = SQRT(DOT_PRODUCT(rab,rab))
             IF (ilist <= neighbor_kind_pair%nscale) THEN
                WRITE (UNIT=output_unit,FMT="(T3,2(I6,3(1X,F10.6)),3(1X,I3),10X,F8.4,L4,F11.5,F9.5)")&
                     atom_a,ra(1:3)*conv,&
                     atom_b,rb(1:3)*conv,&
                     neighbor_kind_pair%cell_vector,&
                     dab*conv,&
                     neighbor_kind_pair%is_onfo(ilist),&
                     neighbor_kind_pair%vdw_scale(ilist),&
                     neighbor_kind_pair%ei_scale(ilist)
             ELSE
                WRITE (UNIT=output_unit,FMT="(T3,2(I6,3(1X,F10.6)),3(1X,I3),10X,F8.4)")&
                     atom_a,ra(1:3)*conv,&
                     atom_b,rb(1:3)*conv,&
                     neighbor_kind_pair%cell_vector,&
                     dab*conv
             END IF
          END IF
       END DO ! ilist
    END DO ! iab

    string = ""
    WRITE (UNIT=string,FMT="(A,I12,A,I12)")&
      "Total number of neighbor interactions for process",mype,":",&
      nneighbor
    CALL compress(string)
    IF (output_unit > 0) WRITE (UNIT=output_unit,FMT="(/,T2,A)") TRIM(string)

  END SUBROUTINE write_neighbor_lists

! *****************************************************************************
!> \brief Sort the generated neighbor list according the kind
!> \param nonbonded ...
!> \param nkinds ...
!> \par History
!>      09.2007 created [tlaino] University of Zurich - Reducing memory usage
!>              for the FIST neighbor lists
!> \author Teodoro Laino - University of Zurich
! *****************************************************************************
  SUBROUTINE sort_neighbor_lists(nonbonded, nkinds)

    TYPE(fist_neighbor_type), POINTER        :: nonbonded
    INTEGER, INTENT(IN)                      :: nkinds

    CHARACTER(LEN=*), PARAMETER :: routineN = 'sort_neighbor_lists', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iab, id_kind, ikind, &
                                                ipair, jkind, max_alloc_size, &
                                                npairs, nscale, tmp
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: indj
    INTEGER, DIMENSION(:), POINTER           :: work
    INTEGER, DIMENSION(:, :), POINTER        :: list_copy
    LOGICAL                                  :: failure
    TYPE(neighbor_kind_pairs_type), POINTER  :: neighbor_kind_pair

    failure = .FALSE.
    NULLIFY(neighbor_kind_pair)
    CALL timeset(routineN,handle)
    ! define a lookup table to get jkind for a given id_kind
    ALLOCATE(indj(nkinds*(nkinds+1)/2))
    id_kind = 0
    DO jkind = 1, nkinds
       DO ikind = jkind, nkinds
          id_kind = id_kind + 1
          indj(id_kind)=jkind
       END DO
    END DO
    ! loop over all nlists and sort the pairs within each list.
    DO iab=1, nonbonded%nlists
       neighbor_kind_pair => nonbonded%neighbor_kind_pairs(iab)
       npairs = neighbor_kind_pair%npairs
       nscale = neighbor_kind_pair%nscale
       IF (npairs/=0) THEN
          IF (npairs>nscale) THEN
             ! 1) Sort the atom pairs by id_kind. Pairs whose interactions are
             ! scaled (possibly to zero for exclusion) are not touched. They
             ! stay packed in the beginning. Sorting is skipped altogether when
             ! all pairs have scaled interactions.
             ALLOCATE(work(1:npairs-nscale))
             ALLOCATE(list_copy(2,1:npairs-nscale))
             ! Copy of the pair list is required to perform the permutation below
             ! correctly.
             list_copy = neighbor_kind_pair%list(:,nscale+1:npairs)
             CALL sort(neighbor_kind_pair%id_kind(nscale+1:npairs),npairs-nscale,work)
             ! Reorder atoms using the same permutation that was used to sort
             ! the array id_kind.
             DO ipair = nscale+1, npairs
                tmp = work(ipair-nscale)
                neighbor_kind_pair%list(1,ipair) = list_copy(1,tmp)
                neighbor_kind_pair%list(2,ipair) = list_copy(2,tmp)
             END DO
             DEALLOCATE(work)
             DEALLOCATE(list_copy)
          END IF
          ! 2) determine the intervals (groups) in the pair list that correspond
          !    to a certain id_kind. also store the corresponding ikind and
          !    jkind. Note that this part does not assume ikind to be sorted,
          !    but it only makes sense when contiguous blobs of the same ikind
          !    are present.
          ! Allocate sufficient memory in case all pairs of atom kinds are
          ! present, and also provide storage for the pairs with exclusion
          ! flags, which are unsorted.
          max_alloc_size = nkinds*(nkinds+1)/2 + nscale
          IF (ASSOCIATED(neighbor_kind_pair%grp_kind_start)) THEN
             DEALLOCATE(neighbor_kind_pair%grp_kind_start)
          END IF
          ALLOCATE(neighbor_kind_pair%grp_kind_start(max_alloc_size))
          IF (ASSOCIATED(neighbor_kind_pair%grp_kind_end)) THEN
             DEALLOCATE(neighbor_kind_pair%grp_kind_end)
          END IF
          ALLOCATE(neighbor_kind_pair%grp_kind_end(max_alloc_size))
          IF (ASSOCIATED(neighbor_kind_pair%ij_kind)) THEN
             DEALLOCATE(neighbor_kind_pair%ij_kind)
          END IF
          ALLOCATE(neighbor_kind_pair%ij_kind(2,max_alloc_size))
          ! Start the first interval.
          ipair = 1
          neighbor_kind_pair%ngrp_kind = 1
          neighbor_kind_pair%grp_kind_start(neighbor_kind_pair%ngrp_kind)=ipair
          ! Get ikind and jkind corresponding to id_kind.
          id_kind = neighbor_kind_pair%id_kind(ipair)
          jkind = indj(id_kind)
          tmp   = nkinds - jkind
          ikind = nkinds + id_kind - nkinds*(nkinds+1)/2 + (tmp*(tmp+1)/2)
          neighbor_kind_pair%ij_kind(1,neighbor_kind_pair%ngrp_kind) = ikind
          neighbor_kind_pair%ij_kind(2,neighbor_kind_pair%ngrp_kind) = jkind
          ! Define the remaining intervals.
          DO ipair = 2, npairs
             IF (neighbor_kind_pair%id_kind(ipair)/=neighbor_kind_pair%id_kind(ipair-1)) THEN
                neighbor_kind_pair%grp_kind_end(neighbor_kind_pair%ngrp_kind) = ipair-1
                neighbor_kind_pair%ngrp_kind = neighbor_kind_pair%ngrp_kind + 1
                neighbor_kind_pair%grp_kind_start(neighbor_kind_pair%ngrp_kind)=ipair
                ! Get ikind and jkind corresponding to id_kind.
                id_kind = neighbor_kind_pair%id_kind(ipair)
                jkind = indj(id_kind)
                tmp   = nkinds - jkind
                ikind = nkinds + id_kind - nkinds*(nkinds+1)/2 + (tmp*(tmp+1)/2)
                neighbor_kind_pair%ij_kind(1,neighbor_kind_pair%ngrp_kind) = ikind
                neighbor_kind_pair%ij_kind(2,neighbor_kind_pair%ngrp_kind) = jkind
             END IF
          END DO
          ! Finish the last interval.
          neighbor_kind_pair%grp_kind_end(neighbor_kind_pair%ngrp_kind) = npairs
          ! Reduce the grp arrays to the actual size because not all pairs of
          ! atom types have to be present in this pair list.
          CALL reallocate(neighbor_kind_pair%grp_kind_start,1,neighbor_kind_pair%ngrp_kind)
          CALL reallocate(neighbor_kind_pair%grp_kind_end,  1,neighbor_kind_pair%ngrp_kind)
          CALL reallocate(neighbor_kind_pair%ij_kind,   1,2,1,neighbor_kind_pair%ngrp_kind)
       END IF
       ! Clean the memory..
       DEALLOCATE(neighbor_kind_pair%id_kind)
    END DO
    DEALLOCATE(indj)
    CALL timestop(handle)
  END SUBROUTINE sort_neighbor_lists

END MODULE fist_neighbor_lists
