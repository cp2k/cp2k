!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!****** cp2k/fist_neighbor_lists [1.0] *
!!
!!   NAME
!!     fist_neighbor_lists
!!
!!   FUNCTION
!!     Generate the atomic neighbor lists for FIST.
!!
!!   AUTHOR
!!     MK (19.11.2002,24.07.2003)
!!     Teodoro Laino (08.2006) - MAJOR REWRITING
!!
!!   MODIFICATION HISTORY
!!     - build and update merged (11.02.2005,MK)
!!     - bug fix for PERIODIC NONE (24.02.06,MK)
!!     - Major rewriting (light memory neighbor lists): teo and joost (05.2006)
!!     - Completely new algorithm for the neighbor lists 
!!       (faster and memory lighter) (Teo 08.2006)
!!
!!   SOURCE
!******************************************************************************

MODULE fist_neighbor_lists

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE cell_types,                      ONLY: cell_type,&
                                             get_cell,&
                                             pbc,&
                                             plane_distance,&
                                             real_to_scaled,&
                                             scaled_to_real
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE fist_neighbor_list_types,        ONLY: fist_neighbor_add,&
                                             fist_neighbor_deallocate,&
                                             fist_neighbor_init,&
                                             fist_neighbor_type,&
                                             neighbor_kind_pairs_type
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE mathlib,                         ONLY: matvec_3x3
  USE memory_utilities,                ONLY: reallocate
  USE particle_types,                  ONLY: particle_type
  USE qmmm_ff_fist,                    ONLY: qmmm_ff_precond_only_qm
  USE qs_neighbor_lists,               ONLY: allocate_subcell,&
                                             deallocate_subcell,&
                                             local_atoms_type,&
                                             subcell_type
  USE string_utilities,                ONLY: compress
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: sort
  USE f77_blas
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  ! Global parameters (in this module)
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'fist_neighbor_lists'

  ! Public subroutines
  PUBLIC :: build_fist_neighbor_lists,&
            build_neighbor_lists,&
            write_neighbor_lists

CONTAINS

!!****f* fist_neighbor_lists/build_fist_neighbor_lists [1.0] *
!!
!!   NAME
!!     build_fist_neighbor_lists
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     08.2006 created [tlaino]
!!
!!*** **********************************************************************
  SUBROUTINE build_fist_neighbor_lists(atomic_kind_set,particle_set,&
       local_particles,cell,r_max, r_minsq, nonbonded,para_env,&
       build_from_scratch,geo_check,force_env_section,pot_type,error)

    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(distribution_1d_type), OPTIONAL, &
      POINTER                                :: local_particles
    TYPE(cell_type), POINTER                 :: cell
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: r_max, r_minsq
    TYPE(fist_neighbor_type), POINTER        :: nonbonded
    TYPE(cp_para_env_type), POINTER          :: para_env
    LOGICAL, INTENT(IN)                      :: build_from_scratch, geo_check
    TYPE(section_vals_type), POINTER         :: force_env_section
    LOGICAL, DIMENSION(:, :), OPTIONAL, &
      POINTER                                :: pot_type
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'build_fist_neighbor_lists', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: kind_name, print_key_path
    INTEGER :: atom_a, handle, iatom_local, ikind, iw, maxatom, &
      maxatom_local, natom_local_a, nkind, output_unit, stat
    LOGICAL                                  :: failure, &
                                                present_local_particles, &
                                                print_subcell_grid
    LOGICAL, DIMENSION(:), POINTER           :: skip_kind
    LOGICAL, DIMENSION(:, :), POINTER        :: my_pot_type
    REAL(dp), DIMENSION(3)                   :: r_pbc
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(local_atoms_type), ALLOCATABLE, &
      DIMENSION(:)                           :: atom

    failure = .FALSE.
    CALL timeset(routineN,"I","",handle)
    NULLIFY(logger)
    logger => cp_error_get_logger(error)

    print_subcell_grid = .FALSE.
    output_unit = cp_print_key_unit_nr(logger,force_env_section,"MM%PRINT%SUBCELL",&
         extension=".Log",error=error)
    IF (output_unit > 0)  print_subcell_grid = .TRUE.

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
         maxatom=maxatom)

    ! Allocate work storage
    present_local_particles = PRESENT(local_particles)
    nkind = SIZE(atomic_kind_set)
    ALLOCATE (atom(nkind),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (skip_kind(nkind),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ! Pot_type
    IF (PRESENT(pot_type)) THEN
       my_pot_type => pot_type
    ELSE
       ALLOCATE(my_pot_type(nkind,nkind),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       my_pot_type = .FALSE.
    END IF
    ! Initialize the local data structures
    maxatom_local = 0
    DO ikind=1,nkind
       atomic_kind => atomic_kind_set(ikind)
       NULLIFY (atom(ikind)%list)
       NULLIFY (atom(ikind)%list_local_a_index)
       NULLIFY (atom(ikind)%r_pbc)
       NULLIFY (atom(ikind)%s_pbc)

       CALL get_atomic_kind(atomic_kind=atomic_kind,&
            atom_list=atom(ikind)%list,name=kind_name)
       skip_kind(ikind)=qmmm_ff_precond_only_qm(kind_name)
       IF(present_local_particles) THEN
          natom_local_a = local_particles%n_el(ikind)
       ELSE
          natom_local_a = SIZE(atom(ikind)%list)
       END IF
       maxatom_local = MAX(maxatom_local,natom_local_a)
       IF (natom_local_a > 0) THEN
          ALLOCATE (atom(ikind)%list_local_a_index(natom_local_a),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
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
          CALL fist_neighbor_deallocate(nonbonded,error)
       END IF
    END IF

    ! Build the nonbonded neighbor lists
    CALL build_neighbor_lists(nonbonded,particle_set,atom,cell,maxatom_local,&
         print_subcell_grid,output_unit,r_max,r_minsq, geo_check,"NONBONDED",&
         skip_kind,my_pot_type,error)

    ! Sort the list according kinds for each cell
    CALL sort_neighbor_lists(nonbonded, nkind, error)

    print_key_path = "MM%PRINT%NEIGHBOR_LISTS"

    IF (BTEST(cp_print_key_should_output(logger%iter_info,force_env_section,&
         print_key_path,error=error),&
         cp_p_file)) THEN
       iw = cp_print_key_unit_nr(logger=logger,&
            basis_section=force_env_section,&
            print_key_path=print_key_path,&
            extension=".out",&
            middle_name="nonbonded_nl",&
            local=.TRUE.,&
            log_filename=.FALSE.,&
            file_position="REWIND",&
            error=error)
       CALL write_neighbor_lists(nonbonded,particle_set,cell,para_env,iw,&
            "NONBONDED NEIGHBOR LISTS",error)
       CALL cp_print_key_finished_output(unit_nr=iw,&
            logger=logger,&
            basis_section=force_env_section,&
            print_key_path=print_key_path,&
            local=.TRUE.,&
            error=error)
    END IF

    ! Release work storage
    DO ikind=1,nkind
       NULLIFY (atom(ikind)%list)
       IF (ASSOCIATED(atom(ikind)%list_local_a_index)) THEN
          DEALLOCATE (atom(ikind)%list_local_a_index,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END DO
    IF (PRESENT(pot_type)) THEN
       NULLIFY(my_pot_type)
    ELSE
       DEALLOCATE(my_pot_type,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    DEALLOCATE (atom,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    DEALLOCATE (skip_kind,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL cp_print_key_finished_output(unit_nr=output_unit,&
         logger=logger,&
         basis_section=force_env_section,&
         print_key_path="MM%PRINT%SUBCELL",&
         error=error)
    CALL timestop(0.0_dp,handle)

  END SUBROUTINE build_fist_neighbor_lists

!!****f* fist_neighbor_lists/build_neighbor_lists [1.0] *
!!
!!   NAME
!!     build_neighbor_lists
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     08.2006 created [tlaino]
!!
!!*** **********************************************************************
  SUBROUTINE build_neighbor_lists(nonbonded,particle_set,atom,cell,maxatom_local,&
       print_subcell_grid,output_unit,r_max,r_minsq,geo_check,name,skip_kind,&
       full_nl,error)

    TYPE(fist_neighbor_type), POINTER        :: nonbonded
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(local_atoms_type), DIMENSION(:), &
      INTENT(IN)                             :: atom
    TYPE(cell_type), POINTER                 :: cell
    INTEGER, INTENT(IN)                      :: maxatom_local
    LOGICAL, INTENT(IN)                      :: print_subcell_grid
    INTEGER, INTENT(IN)                      :: output_unit
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: r_max, r_minsq
    LOGICAL, INTENT(IN)                      :: geo_check
    CHARACTER(LEN=*), INTENT(IN)             :: name
    LOGICAL, DIMENSION(:), POINTER           :: skip_kind
    LOGICAL, DIMENSION(:, :), POINTER        :: full_nl
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'build_neighbor_lists', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=8)                         :: unit_of_length_name
    INTEGER :: a_i, a_j, a_k, atom_a, atom_b, b_i, b_j, b_k, b_pi, b_pj, &
      b_pk, bg_i, bg_j, bg_k, handle, i, i1, iatom_local, icell, icellmap, &
      id_kind, ii, ii_start, ij, ij_start, ik, ik_start, ikind, imap, &
      imax_cell, invcellmap, iw, ix, j, j1, jatom_local, jcell, jkind, jx, k, &
      kcell, kx, natom_local_a, ncellmax, nkind, nkind00, stat, tmpdim, xdim, &
      ydim, zdim
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: work
    INTEGER, ALLOCATABLE, DIMENSION(:, :, :) :: cellmap
    INTEGER, DIMENSION(3)                    :: isubcell, ncell, nsubcell, &
                                                periodic
    LOGICAL                                  :: added, any_full, atom_order, &
                                                failure, is_full, subcell000
    LOGICAL, ALLOCATABLE, DIMENSION(:, :, :) :: sphcub
    REAL(dp)                                 :: rab2, rab2_max, rab2_min, &
                                                rab_max, unit_of_length
    REAL(dp), DIMENSION(3)                   :: abc, cell_v, cv_b, r, rab, &
                                                rb, sab_max
    REAL(KIND=dp)                            :: ic(3), icx(3), radius, vv
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: coord
    TYPE(neighbor_kind_pairs_type), POINTER  :: inv_neighbor_kind_pair, &
                                                neighbor_kind_pair
    TYPE(subcell_type), DIMENSION(:, :, :), &
      POINTER                                :: subcell_a, subcell_b

    CALL timeset(routineN,"I","",handle)
    failure = .FALSE.
    nkind = SIZE(atom)
    nsubcell = 1
    isubcell = 0
    ncell    = 0
    any_full = ANY(full_nl)
    CALL get_cell(cell=cell,&
                  abc=abc,&
                  periodic=periodic,&
                  unit_of_length=unit_of_length,&
                  unit_of_length_name=unit_of_length_name)
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
    CALL fist_neighbor_init(nonbonded, ncell, error)
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
    ALLOCATE(cellmap(-ncellmax:ncellmax,-ncellmax:ncellmax,-ncellmax:ncellmax),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
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
                   CPPostcondition(imap<=nonbonded%nlists,cp_failure_level,routineP,error,failure)
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
                    -isubcell(3):isubcell(3)),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
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
    ALLOCATE(coord(3,SIZE(particle_set)),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
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
             ALLOCATE(subcell_a(i,j,k)%atom_list(subcell_a(i,j,k)%natom),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
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
             ALLOCATE(subcell_b(i,j,k)%atom_list(subcell_b(i,j,k)%natom),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
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
    ALLOCATE(work(3*tmpdim),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DO k = 1,nsubcell(3)
       DO j = 1,nsubcell(2)
          DO i = 1,nsubcell(1)
             CALL reorder_atoms_subcell(subcell_a(i,j,k)%atom_list,particle_set,work)
             CALL reorder_atoms_subcell(subcell_b(i,j,k)%atom_list,particle_set,work)
          END DO
       END DO
    END DO
    DEALLOCATE(work,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
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
                IF (skip_kind(ikind).AND.skip_kind(jkind)) CYCLE
                ! Screen interaction to avoid double counting
                atom_order = (atom_a<=atom_b)
                ! Special case for kind combination requiring the full NL
                IF (any_full) THEN
                   is_full = full_nl(ikind,jkind)
                   IF (is_full) THEN
                      atom_order = (atom_a==atom_b)
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
                   CALL fist_neighbor_add(neighbor_kind_pair, atom_a, atom_b,&
                        exclusion_list = particle_set(atom_a)%list_exclude_vdw,&
                        particle_set=particle_set, rab=rab, added=added, id_kind=id_kind,&
                        is_full=is_full,cell=cell,error=error)
                   ! This is to handle properly when interaction radius is larger than cell size
                   IF ((atom_a==atom_b).AND.(ik_start==0)) THEN
                      invcellmap = cellmap(-b_pi,-b_pj,-b_pk)
                      inv_neighbor_kind_pair => nonbonded%neighbor_kind_pairs(invcellmap)
                      rab = rab - 2.0_dp * cell_v
                      CALL fist_neighbor_add(inv_neighbor_kind_pair, atom_a, atom_b,&
                           exclusion_list = particle_set(atom_a)%list_exclude_vdw,&
                           particle_set=particle_set, rab=rab, added=added, id_kind=id_kind,&
                           is_full=is_full,cell=cell,error=error) 
                   END IF
                   ! Check for too close hits
                   IF (added) THEN
                      rab2_min = r_minsq(ikind,jkind)
                      IF (rab2 < rab2_min) THEN
                         iw=cp_logger_get_default_unit_nr()
                         WRITE(iw,'(T2,A,2I7,2(A,F15.8),A)')"WARNING| Particles: ",&
                              atom_a,atom_b,&
                              " at distance [au]:",SQRT(rab2)," less than: ",&
                              SQRT(rab2_min),&
                              "; increase EMAX_SPLINE."
                         IF (rab2 < rab2_min/(1.06_dp)**2) THEN
                            IF (geo_check) CALL stop_program(routineN, &
                                 "GEOMETRY wrong or EMAX_SPLINE too small!")
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
    DEALLOCATE(coord,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(cellmap,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(sphcub,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL deallocate_subcell(subcell_a)
    CALL deallocate_subcell(subcell_b)

    CALL timestop(0.0_dp,handle)
  END SUBROUTINE build_neighbor_lists

!!****f* fist_neighbor_lists/give_ijk_subcell [1.0] *
!!
!!   NAME
!!     give_ijk_subcell
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     08.2006 created [tlaino]
!!
!!*** **********************************************************************
  SUBROUTINE give_ijk_subcell(r,i,j,k,cell,nsubcell)
    REAL(KIND=dp)                            :: r(3)
    INTEGER, INTENT(OUT)                     :: i, j, k
    TYPE(cell_type), POINTER                 :: cell
    INTEGER, DIMENSION(3), INTENT(IN)        :: nsubcell

    REAL(KIND=dp)                            :: r_pbc(3), s(3), s_pbc(3)

    r_pbc = r
    CALL real_to_scaled(s_pbc, r_pbc, cell)
    s(:) = s_pbc + 0.5_dp
    i = INT(s(1)*REAL(nsubcell(1),KIND=dp)) + 1
    j = INT(s(2)*REAL(nsubcell(2),KIND=dp)) + 1
    k = INT(s(3)*REAL(nsubcell(3),KIND=dp)) + 1
    i = MIN(MAX(i,1),nsubcell(1))
    j = MIN(MAX(j,1),nsubcell(2))
    k = MIN(MAX(k,1),nsubcell(3))
    
  END SUBROUTINE give_ijk_subcell

!!****f* fist_neighbor_lists/reorder_atoms_subcell [1.0] *
!!
!!   NAME
!!     reorder_atoms_subcell
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     08.2006 created [tlaino]
!!
!!*** **********************************************************************
  SUBROUTINE reorder_atoms_subcell(atom_list,particle_set,work)
    ! work needs to be dimensioned 3xSIZE(atom_list)
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    INTEGER, DIMENSION(:)                    :: work

    INTEGER                                  :: i, i0, i1, i2, j0, j1, j2

    i0 = 1
    j0 = SIZE(atom_list)
    i1 = j0+1
    j1 = 2*j0
    i2 = j1+1
    j2 = 3*j0
    ! Sort kind
    DO i = 1, SIZE(atom_list)
       work(i0+i-1) = particle_set(atom_list(i))%atomic_kind%kind_number
    END DO
    CALL sort(work(i0:j0),SIZE(atom_list),work(i1:j1))
    work(i2:j2) = atom_list
    DO i = 1, SIZE(atom_list)
       atom_list(i) = work(i2+work(i1+i-1)-1)
    END DO
  END SUBROUTINE reorder_atoms_subcell

!!****f* fist_neighbor_lists/write_neighbor_lists [1.0] *
!!
!!   NAME
!!     write_neighbor_lists
!!
!!   FUNCTION
!!     Write a set of neighbor lists to the output unit.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     08.2006 created [tlaino]
!!
!!*** **********************************************************************
  SUBROUTINE write_neighbor_lists(nonbonded,particle_set,cell,para_env,output_unit,&
                                  name,error)

    TYPE(fist_neighbor_type), POINTER        :: nonbonded
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER, INTENT(IN)                      :: output_unit
    CHARACTER(LEN=*), INTENT(IN)             :: name
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=8)                         :: unit_of_length_name
    CHARACTER(LEN=default_string_length)     :: string
    INTEGER                                  :: atom_a, atom_b, iab, ilist, &
                                                mype, nneighbor
    LOGICAL                                  :: print_headline
    REAL(dp)                                 :: dab, unit_of_length
    REAL(dp), DIMENSION(3)                   :: cell_v, ra, rab, rb
    TYPE(neighbor_kind_pairs_type), POINTER  :: neighbor_kind_pair

    mype = para_env%mepos
    CALL get_cell(cell=cell,&
                  unit_of_length=unit_of_length,&
                  unit_of_length_name=unit_of_length_name)

    ! Print headline
    string = ""
    WRITE (UNIT=string,FMT="(A,I5,A)")&
      TRIM(name)//" IN "//TRIM(unit_of_length_name)//" (PROCESS",mype,")"
    CALL compress(string)
    IF (output_unit > 0) WRITE (UNIT=output_unit,FMT="(/,/,T2,A)") TRIM(string)

    print_headline = .TRUE.
    nneighbor = 0
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
                WRITE (UNIT=output_unit,FMT="(/,T3,2(A6,3(5X,A,5X)),1X,A11,10X,A8)")&
                     "Atom-A","X","Y","Z","Atom-B","X","Y","Z","Cell(i,j,k)","Distance"
                print_headline = .FALSE.
             END IF
             
             ra(:) = pbc(particle_set(atom_a)%r,cell)
             rb(:) = pbc(particle_set(atom_b)%r,cell)
             rab   = rb(:)-ra(:)+cell_v
             dab   = SQRT(DOT_PRODUCT(rab,rab))
             WRITE (UNIT=output_unit,FMT="(/,T3,2(I6,3(1X,F10.6)),3(1X,I3),10X,F8.4)")&
                  atom_a,ra(1:3)/unit_of_length,&
                  atom_b,rb(1:3)/unit_of_length,&
                  neighbor_kind_pair%cell_vector,&
                  dab                              
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

!!****f* fist_neighbor_lists/sort_neighbor_lists [1.0] *
!!
!!   NAME
!!     sort_neighbor_lists
!!
!!   FUNCTION
!!     Sort the generated neighbor list according the kind
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Teodoro Laino - University of Zurich
!!
!!   MODIFICATION HISTORY
!!     09.2007 created [tlaino] University of Zurich - Reducing memory usage 
!!             for the FIST neighbor lists
!!
!!*** **********************************************************************
  SUBROUTINE sort_neighbor_lists(nonbonded, nkinds, error)

    TYPE(fist_neighbor_type), POINTER        :: nonbonded
    INTEGER, INTENT(IN)                      :: nkinds
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'sort_neighbor_lists', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iab, id_kind, ikind, &
                                                ipair, jkind, max_alloc_size, &
                                                npairs, stat, tmp
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: indj
    INTEGER, DIMENSION(:), POINTER           :: work
    INTEGER, DIMENSION(:, :), POINTER        :: ab_work
    LOGICAL                                  :: failure
    TYPE(neighbor_kind_pairs_type), POINTER  :: neighbor_kind_pair

    failure = .FALSE.
    NULLIFY(neighbor_kind_pair)
    CALL timeset(routineN,"I","",handle)
    max_alloc_size = nkinds*(nkinds+1)/2
    ALLOCATE(indj(max_alloc_size),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    id_kind = 0
    Ext_Loop: DO jkind = 1, nkinds
       DO ikind = jkind, nkinds
          id_kind = id_kind + 1
          indj(id_kind)=jkind
       END DO
    END DO Ext_Loop
    DO iab=1, nonbonded%nlists
       neighbor_kind_pair => nonbonded%neighbor_kind_pairs(iab)
       npairs = neighbor_kind_pair%npairs
       IF (npairs/=0) THEN 
          ALLOCATE(work(1:npairs),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(ab_work(2,1:npairs),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DO ipair = 1, npairs
             ab_work(1,ipair) = neighbor_kind_pair%list(1,ipair)
             ab_work(2,ipair) = neighbor_kind_pair%list(2,ipair)
          END DO
          CALL sort(neighbor_kind_pair%id_kind,npairs,work)
          ! Remap atoms
          DO ipair = 1, npairs
             tmp = work(ipair)
             neighbor_kind_pair%list(1,ipair) = ab_work(1,tmp)
             neighbor_kind_pair%list(2,ipair) = ab_work(2,tmp)
          END DO
          DEALLOCATE(work,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(ab_work,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          IF (ASSOCIATED(neighbor_kind_pair%grp_kind_start)) THEN
             DEALLOCATE(neighbor_kind_pair%grp_kind_start,stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF
          ALLOCATE(neighbor_kind_pair%grp_kind_start(max_alloc_size),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          IF (ASSOCIATED(neighbor_kind_pair%grp_kind_end)) THEN
             DEALLOCATE(neighbor_kind_pair%grp_kind_end,stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF
          ALLOCATE(neighbor_kind_pair%grp_kind_end(max_alloc_size),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          IF (ASSOCIATED(neighbor_kind_pair%ij_kind)) THEN
             DEALLOCATE(neighbor_kind_pair%ij_kind,stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF
          ALLOCATE(neighbor_kind_pair%ij_kind(2,max_alloc_size),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ipair = 1 
          neighbor_kind_pair%ngrp_kind = 1
          neighbor_kind_pair%grp_kind_start(neighbor_kind_pair%ngrp_kind)=ipair
          jkind = indj(neighbor_kind_pair%id_kind(ipair))
          tmp   = nkinds - jkind
          ikind = nkinds + neighbor_kind_pair%id_kind(ipair) - max_alloc_size +&
               (tmp*(tmp+1)/2)
          neighbor_kind_pair%ij_kind(1,neighbor_kind_pair%ngrp_kind) = ikind
          neighbor_kind_pair%ij_kind(2,neighbor_kind_pair%ngrp_kind) = jkind
          ! Define groups
          DO ipair = 2, npairs
             IF (neighbor_kind_pair%id_kind(ipair)/=neighbor_kind_pair%id_kind(ipair-1)) THEN
                neighbor_kind_pair%grp_kind_end(neighbor_kind_pair%ngrp_kind) = ipair-1
                neighbor_kind_pair%ngrp_kind = neighbor_kind_pair%ngrp_kind + 1
                neighbor_kind_pair%grp_kind_start(neighbor_kind_pair%ngrp_kind)=ipair
                ! 
                jkind = indj(neighbor_kind_pair%id_kind(ipair))
                tmp   = nkinds - jkind
                ikind = nkinds + neighbor_kind_pair%id_kind(ipair) - max_alloc_size +&
                     (tmp*(tmp+1)/2)
                neighbor_kind_pair%ij_kind(1,neighbor_kind_pair%ngrp_kind) = ikind
                neighbor_kind_pair%ij_kind(2,neighbor_kind_pair%ngrp_kind) = jkind
             END IF
          END DO
          neighbor_kind_pair%grp_kind_end(neighbor_kind_pair%ngrp_kind) = npairs
          ! Get the real size..
          CALL reallocate(neighbor_kind_pair%grp_kind_start,1,neighbor_kind_pair%ngrp_kind)
          CALL reallocate(neighbor_kind_pair%grp_kind_end,  1,neighbor_kind_pair%ngrp_kind)
          CALL reallocate(neighbor_kind_pair%ij_kind,   1,2,1,neighbor_kind_pair%ngrp_kind)
       END IF
       ! Clean the memory..
       DEALLOCATE(neighbor_kind_pair%id_kind,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END DO
    DEALLOCATE(indj,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE sort_neighbor_lists

END MODULE fist_neighbor_lists
